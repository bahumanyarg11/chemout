"""Program detection and dispatch."""

from __future__ import annotations

import hashlib
import os
from collections.abc import Iterable, Iterator
from pathlib import Path

from .parsers.base import Lines
from .schema import Calculation, Program

# How many leading lines to look at when sniffing the program.
_SNIFF_LINES = 400


def read_text(path: str | os.PathLike[str]) -> str:
    """Read an output file tolerantly: UTF-8 first, then Latin-1, never raising on bytes."""
    data = Path(path).read_bytes()
    try:
        return data.decode("utf-8")
    except UnicodeDecodeError:
        return data.decode("latin-1")


def sniff_program(text_or_lines: str | list[str]) -> Program:
    lines = text_or_lines.splitlines() if isinstance(text_or_lines, str) else text_or_lines
    head = lines[:_SNIFF_LINES]
    joined = "\n".join(head)
    if "* O   R   C   A *" in joined or "Program Version" in joined and "ORCA" in joined:
        return Program.ORCA
    if "Entering Gaussian System" in joined or "Gaussian, Inc." in joined or " Gaussian 16" in joined:
        return Program.GAUSSIAN
    if "x T B" in joined or "|                           x T B" in joined or "xtb version" in joined:
        return Program.XTB
    if "Psi4: An Open-Source Ab Initio Electronic Structure Package" in joined:
        return Program.PSI4
    # Fall back to a scan of the whole file for the banner; some wrappers prepend noise.
    joined = "\n".join(lines)
    if "* O   R   C   A *" in joined:
        return Program.ORCA
    if "Entering Gaussian System" in joined or "Gaussian, Inc." in joined:
        return Program.GAUSSIAN
    if "x T B" in joined and "Grimme" in joined:
        return Program.XTB
    if "Psi4: An Open-Source" in joined:
        return Program.PSI4
    return Program.UNKNOWN


def _parser_for(program: Program):
    if program is Program.ORCA:
        from .parsers.orca import OrcaParser

        return OrcaParser
    if program is Program.GAUSSIAN:
        from .parsers.gaussian import GaussianParser

        return GaussianParser
    if program is Program.XTB:
        from .parsers.xtb import XtbParser

        return XtbParser
    raise ValueError(f"no parser for program {program.value!r}")


def parse_text(text: str, program: Program | str | None = None, path: str | None = None) -> Calculation:
    """Parse the contents of an output file already in memory."""
    lines = text.splitlines()
    if program is None:
        program = sniff_program(lines)
    program = Program(program)
    if program is Program.UNKNOWN:
        calc = Calculation()
        calc.provenance.path = path
        calc.provenance.warn("sniff", "could not identify the program that wrote this file")
        return calc
    parser = _parser_for(program)(Lines(lines))
    calc = parser.run()
    calc.provenance.path = path
    return calc


def parse(path: str | os.PathLike[str], program: Program | str | None = None) -> Calculation:
    """Parse one output file. Never raises on malformed content; see ``calc.provenance.warnings``."""
    p = Path(path)
    data = p.read_bytes()
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        text = data.decode("latin-1")
    calc = parse_text(text, program=program, path=str(p))
    calc.provenance.sha256 = hashlib.sha256(data).hexdigest()
    calc.provenance.size_bytes = len(data)
    return calc


def parse_many(paths: Iterable[str | os.PathLike[str]], program: Program | str | None = None) -> Iterator[Calculation]:
    for p in paths:
        yield parse(p, program=program)
