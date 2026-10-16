"""Shared machinery: a cursor over lines, tolerant number parsing, and the parser base class.

Design rule: a parser is a list of *section handlers*. The main loop walks the file once;
whenever a handler's trigger matches the current line, the handler consumes what it needs.
Each handler runs inside ``guard``: an exception inside one section becomes a
``ParseWarning`` and the walk continues at the next line, so a surprising format change in one
table can never take the whole file down. That is the failure mode of regex parsers that this
package is built to avoid.
"""

from __future__ import annotations

import re
import traceback
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from typing import Any

from ..schema import Calculation, Program, ProgramInfo, Termination

_FLOAT_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][-+]?\d+)?")
_FORTRAN_FLOAT_RE = re.compile(r"[-+]?\d*\.\d+(?:[dD][-+]?\d+)?|[-+]?\d+(?:\.\d*)?(?:[dD][-+]?\d+)?")

PERIODIC_TABLE = (
    "X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr "
    "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb "
    "Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr "
    "Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
).split()
SYMBOL_TO_Z = {s: i for i, s in enumerate(PERIODIC_TABLE)}
SYMBOL_TO_Z_LOWER = {s.lower(): i for i, s in enumerate(PERIODIC_TABLE)}


def symbol_to_z(symbol: str) -> int | None:
    """'C' -> 6, 'c' -> 6, 'Cl' -> 17, 'C1'/'C(iso=13)' -> 6, 'Cl3' -> 17; None when unknown."""
    s = re.match(r"[A-Za-z]{1,2}", symbol.strip())
    if not s:
        return None
    token = s.group(0)
    z = SYMBOL_TO_Z.get(token.capitalize())
    if z is None and len(token) == 2:
        z = SYMBOL_TO_Z.get(token[0].upper())
    return z or None  # 'X' (dummy atom, index 0) is not an element


def z_to_symbol(z: int) -> str:
    return PERIODIC_TABLE[z] if 0 <= z < len(PERIODIC_TABLE) else "X"


def to_float(token: str) -> float:
    """Parse a number as any of these programs print it: Fortran D exponents, '********' overflow."""
    t = token.strip().replace("D", "E").replace("d", "e")
    if set(t) <= {"*"}:
        return float("nan")
    return float(t)


def floats_in(line: str) -> list[float]:
    return [to_float(m) for m in _FLOAT_RE.findall(line)]


class Lines:
    """Cursor over the lines of a file with cheap lookahead."""

    __slots__ = ("lines", "pos")

    def __init__(self, lines: list[str]):
        self.lines = lines
        self.pos = 0

    def __len__(self) -> int:
        return len(self.lines)

    @property
    def done(self) -> bool:
        return self.pos >= len(self.lines)

    @property
    def current(self) -> str:
        return self.lines[self.pos] if self.pos < len(self.lines) else ""

    def peek(self, offset: int = 1) -> str:
        i = self.pos + offset
        return self.lines[i] if 0 <= i < len(self.lines) else ""

    def next(self) -> str:
        line = self.current
        self.pos += 1
        return line

    def skip(self, n: int = 1) -> None:
        self.pos = min(self.pos + n, len(self.lines))

    def skip_blank(self) -> None:
        while not self.done and not self.current.strip():
            self.pos += 1

    def find(self, predicate: Callable[[str], bool] | str, limit: int | None = None) -> bool:
        """Advance to the first line satisfying ``predicate`` (or containing ``predicate`` if str).

        Returns False and leaves the cursor at the limit/EOF when not found."""
        if isinstance(predicate, str):
            needle = predicate
            predicate = lambda s: needle in s  # noqa: E731
        end = len(self.lines) if limit is None else min(len(self.lines), self.pos + limit)
        while self.pos < end:
            if predicate(self.lines[self.pos]):
                return True
            self.pos += 1
        return False

    def take_while(self, predicate: Callable[[str], bool]) -> list[str]:
        out: list[str] = []
        while not self.done and predicate(self.current):
            out.append(self.next())
        return out

    def take_until_blank(self) -> list[str]:
        return self.take_while(lambda s: bool(s.strip()))

    def __iter__(self) -> Iterator[str]:
        while not self.done:
            yield self.next()


Handler = Callable[["Parser", str], None]


class Parser:
    """Base class. Subclasses register handlers with ``@section`` and may override ``finalize``."""

    program: Program = Program.UNKNOWN
    _handlers: list[tuple[Callable[[str], bool], Handler, str]]

    def __init_subclass__(cls, **kw: Any) -> None:
        super().__init_subclass__(**kw)
        cls._handlers = []
        for base in reversed(cls.__mro__):
            for name, attr in vars(base).items():
                trig = getattr(attr, "_trigger", None)
                if trig is not None:
                    cls._handlers.append((trig, attr, name))

    def __init__(self, lines: Lines):
        self.lines = lines
        self.calc = Calculation(program=ProgramInfo(name=self.program))
        self.state: dict[str, Any] = {}

    # -- warnings ----------------------------------------------------------------
    def warn(self, section: str, message: str) -> None:
        self.calc.provenance.warn(section, message, line=self.lines.pos + 1)

    @contextmanager
    def guard(self, section: str):
        start = self.lines.pos
        try:
            yield
        except Exception as exc:  # noqa: BLE001 - every failure must be recorded, never raised
            tb = traceback.extract_tb(exc.__traceback__)[-1]
            self.calc.provenance.warn(
                section,
                f"{type(exc).__name__}: {exc} (parser {tb.name}:{tb.lineno})",
                line=start + 1,
            )
            # Make sure we make progress even if the handler failed before consuming anything.
            if self.lines.pos == start:
                self.lines.pos = start + 1

    # -- main loop -----------------------------------------------------------------
    def run(self) -> Calculation:
        lines = self.lines
        handlers = self._handlers
        while not lines.done:
            line = lines.current
            matched = False
            for trigger, handler, name in handlers:
                if trigger(line):
                    matched = True
                    with self.guard(name):
                        handler(self, line)
                    break
            if not matched:
                lines.pos += 1
        with self.guard("finalize"):
            self.finalize()
        return self.calc

    def finalize(self) -> None:
        """Derive final fields; default marks termination INCOMPLETE if nothing decided."""
        if self.calc.termination is Termination.UNKNOWN:
            self.calc.termination = Termination.INCOMPLETE


def section(trigger: str | Callable[[str], bool]):
    """Decorator registering a method as a handler. ``trigger`` is a substring or predicate."""

    if isinstance(trigger, str):
        needle = trigger
        pred: Callable[[str], bool] = lambda s: needle in s  # noqa: E731
    else:
        pred = trigger

    def deco(fn: Handler) -> Handler:
        fn._trigger = pred  # type: ignore[attr-defined]
        return fn

    return deco


def regex_section(pattern: str):
    rx = re.compile(pattern)
    return section(lambda s: rx.search(s) is not None)
