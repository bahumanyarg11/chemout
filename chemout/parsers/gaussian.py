"""Gaussian 98/03/09/16 output parser.

Gaussian prints geometries in "Input orientation" and (unless NoSymm) "Standard orientation".
We keep one geometry per optimisation point, preferring the standard orientation when it
exists, because forces, normal modes and multipoles are printed in that frame.
"""

from __future__ import annotations

import re

from ..schema import (
    HARTREE_TO_KCAL,
    ExcitedState,
    Geometry,
    JobType,
    Multipoles,
    Optimization,
    OptimizationStep,
    Orbitals,
    PopulationAnalysis,
    Program,
    SCFCycle,
    Termination,
    Thermochemistry,
    Timing,
    Transition,
    Vibrations,
)
from .base import Parser, floats_in, section, to_float, z_to_symbol

_VERSION_RE = re.compile(r"Gaussian\s+(\d+):\s+(\S+)\s+(\S+)")
_SCF_DONE_RE = re.compile(r"SCF Done:\s+E\((\S+)\)\s*=\s*(\S+)\s+A\.U\.\s+after\s+(\d+)\s+cycles")
_EIGEN_RE = re.compile(r"-?\d+\.\d+")
_EXC_RE = re.compile(
    r"Excited State\s+(\d+):\s+(\S+)\s+(-?\d+\.\d+)\s+eV\s+(-?\d+\.\d+)\s+nm\s+f=\s*(-?\d+\.\d+)"
    r"(?:\s+<S\*\*2>=\s*(-?\d+\.\d+))?"
)
_TRANS_RE = re.compile(r"^\s*(\d+)([AB]?)\s*(->|<-)\s*(\d+)([AB]?)\s+(-?\d+\.\d+)")
_TIME_RE = re.compile(r"(\d+)\s+days\s+(\d+)\s+hours\s+(\d+)\s+minutes\s+([\d.]+)\s+seconds")
_ROUTE_JOBS = {
    "opt": JobType.OPTIMIZATION,
    "freq": JobType.FREQUENCY,
    "td": JobType.EXCITED_STATES,
    "tda": JobType.EXCITED_STATES,
    "cis": JobType.EXCITED_STATES,
    "eomccsd": JobType.EXCITED_STATES,
    "nmr": JobType.NMR,
    "scan": JobType.SCAN,
    "irc": JobType.IRC,
    "force": JobType.GRADIENT,
    "sp": JobType.ENERGY,
}
_ROUTE_OPTIONS = {
    "symmetry",
    "symm",
    "nosymm",
    "nosymmetry",
    "population",
    "pop",
    "geom",
    "guess",
    "scf",
    "int",
    "integral",
    "iop",
    "scrf",
    "gfinput",
    "gfprint",
    "density",
    "output",
    "punch",
    "units",
    "test",
    "name",
    "maxdisk",
    "cphf",
    "pseudo",
    "temperature",
    "pressure",
    "chkbasis",
    "gen",
    "genecp",
    "extrabasis",
    "counterpoise",
    "field",
    "polar",
    "prop",
    "volume",
    "nmr",
    "stable",
    "charge",
    "sparse",
    "transformation",
    "trans",
    "window",
    "fmm",
    "nofmm",
    "extralinks",
    "nodensityfit",
    "densityfit",
    "denfit",
    "nocaf",
    "connectivity",
    "readisotopes",
    "externaliteration",
    "complex",
    "restart",
    "archive",
    "noraff",
    "fchk",
    "formcheck",
    "emprical",
    "dispersion",
    "noexternal",
    "scfcyc",
    "massage",
    "oniom",
    "admp",
    "bomd",
    "irc",
    "scan",
    "opt",
    "freq",
    "force",
    "sp",
    "td",
    "tda",
}
_METHOD_WITH_OPTIONS = {
    "eomccsd",
    "ccsd",
    "cis",
    "td",
    "tda",
    "casscf",
    "mp2",
    "mp4",
    "qcisd",
    "oniom",
    "sac-ci",
    "bd",
    "cid",
    "cisd",
    "w1u",
    "w1bd",
    "cbs-qb3",
    "g3",
    "g4",
    "g2",
    "g1",
}
_GAUSSIAN_BASIS_RE = re.compile(
    r"^(sto-|[3-6]-[0-9]+\+*g|d95|shc|cep-|lanl2|sdd|sddall|cc-pv|aug-cc|jun-cc|may-cc|apr-cc|def2|epr-|ugbs|mtsmall|dgdzvp|dgtzvp|cbsb7|gen$|chkbasis$|genecp$|svp$|tzvp$|qzvp$|tzv$|sv$)",
    re.I,
)
_DISPERSION_RE = re.compile(r"empiricaldispersion\s*=\s*\(?(\w+)", re.I)
_SCRF_RE = re.compile(r"scrf\s*(?:=\s*\(?([^)\s]+(?:,[^)\s]+)*)\)?)?", re.I)


def _same_structure(a: list[list[float]], b: list[list[float]], tol: float = 1e-4) -> bool:
    """True when two geometries have the same interatomic distances (orientation-independent)."""
    if not a or not b or len(a) != len(b):
        return False
    import numpy as np

    xa, xb = np.asarray(a), np.asarray(b)
    da = np.linalg.norm(xa[:, None, :] - xa[None, :, :], axis=-1)
    db = np.linalg.norm(xb[:, None, :] - xb[None, :, :], axis=-1)
    return bool(np.abs(da - db).max() < tol)


def _seconds(line: str) -> float | None:
    m = _TIME_RE.search(line)
    if not m:
        return None
    d, h, mi, s = m.groups()
    return int(d) * 86400 + int(h) * 3600 + int(mi) * 60 + float(s)


class GaussianParser(Parser):
    program = Program.GAUSSIAN

    def __init__(self, lines):
        super().__init__(lines)
        self.state.update(
            {
                "input_geoms": [],
                "std_geoms": [],
                "orb_alpha": [],
                "orb_beta": [],
                "orb_occ_alpha": 0,
                "orb_occ_beta": 0,
                "in_eigen": False,
                "job_index": 0,
                "pending_energy": None,
            }
        )

    # ------------------------------------------------------------------ header
    @section(lambda s: s.startswith(" Gaussian ") and ":" in s and "Rev" in s)
    def _version(self, line: str) -> None:
        m = _VERSION_RE.search(line)
        if m:
            major, build, date = m.groups()
            rev = build.split("Rev")[-1] if "Rev" in build else build
            self.calc.program.version = f"{major} Rev{rev}"
            self.calc.program.revision = build
            self.calc.program.date = date
            self.calc.program.version_tuple = (int(major),)
        self.lines.next()

    @section(lambda s: s.startswith(" %") and "=" in s)
    def _link0(self, line: str) -> None:
        key, _, val = line.strip().lstrip("%").partition("=")
        key = key.lower()
        if key in ("nprocshared", "nproc", "nprocs"):
            self.calc.input.nprocs = int(val)
        elif key == "mem":
            v = val.strip().upper()
            num = float(re.match(r"[\d.]+", v).group(0))
            unit = re.sub(r"[\d.]", "", v)
            factor = {"MB": 1, "GB": 1024, "MW": 8, "GW": 8192, "KB": 1 / 1024, "": 1 / 1024 / 1024}.get(unit, 1)
            self.calc.input.memory_mb = num * factor
        self.lines.next()

    @section(lambda s: s.startswith(" #") and not s.startswith(" ##"))
    def _route(self, line: str) -> None:
        lines = self.lines
        route = line.strip()
        lines.next()
        # Long routes wrap onto following lines until the dashed line.
        while not lines.done and lines.current.strip() and not lines.current.strip().startswith("---"):
            if re.match(r"^\s*\d+/", lines.current):
                break
            route += lines.next().strip()
        if self.state["job_index"] == 0 or self.calc.input.raw is None:
            self.calc.input.raw = route
        else:
            self.calc.input.extra.setdefault("additional_routes", []).append(route)
        self.state["job_index"] += 1
        self._parse_route(route)

    def _parse_route(self, route: str) -> None:
        inp = self.calc.input
        body = route.lstrip("#").strip()
        body = re.sub(r"^[pPnNtT]\s+", "", body)  # #p / #n / #t
        # "Opt=(TS, CalcFC)" -> one token; keep parenthesised groups together, drop spaces inside them
        body = re.sub(r"\(([^)]*)\)", lambda m: "(" + m.group(1).replace(" ", "") + ")", body)
        tokens = re.findall(r"[\w\-\+\(\)\*\.,=/]+", body)
        loose_candidates: list[str] = []
        for tok in tokens:
            low = tok.lower()
            key = low.split("=")[0].split("(")[0]
            if key in _ROUTE_JOBS and _ROUTE_JOBS[key] not in inp.job_types:
                inp.job_types.append(_ROUTE_JOBS[key])
            if "/" in tok and not low.startswith(("geom", "guess", "scrf", "iop", "pop", "density", "output")):
                method, _, basis = tok.partition("/")
                if inp.method is None:
                    self._set_method(method, basis.split("/")[0])
                continue
            m = _DISPERSION_RE.search(tok)
            if m:
                inp.dispersion = m.group(1).upper()
                continue
            if low.startswith("scrf"):
                m = _SCRF_RE.search(tok)
                inp.solvation = m.group(1) if m and m.group(1) else "PCM"
                continue
            if key in _METHOD_WITH_OPTIONS:
                loose_candidates.append(tok)
                continue
            if key in _ROUTE_JOBS or key in _ROUTE_OPTIONS or "=" in tok:
                continue
            loose_candidates.append(tok)
        if inp.method is None and loose_candidates:
            # "Opt EOMCCSD=(nstates=2) STO-3G" or "opt freq 6-31++g(d,p) pbe1pbe": method and basis as separate words
            basis = next((t for t in loose_candidates if _GAUSSIAN_BASIS_RE.match(t.lower())), None)
            rest = [t for t in loose_candidates if t is not basis]
            if rest:
                self._set_method(rest[0].split("=")[0].split("(")[0], basis)
            elif basis:
                inp.basis = basis
        if not inp.job_types:
            inp.job_types.append(JobType.ENERGY)
        if "ts" in [t.lower() for t in re.findall(r"opt\s*=\s*\(?([^)]*)\)?", body, re.I) for t in t.split(",")]:
            inp.job_types = [JobType.TS_OPTIMIZATION if j is JobType.OPTIMIZATION else j for j in inp.job_types]

    def _set_method(self, method: str, basis: str | None) -> None:
        inp = self.calc.input
        inp.method = method.upper()
        if basis:
            inp.basis = basis
        if (
            re.match(
                r"^[ru]?o?(hf|mp[2-5]|ccsd|ccd|qcisd|cis|eomccsd|cid|cisd|casscf|am1|pm[367]|mndo|zindo|huckel|uff|amber|dreiding|oniom)",
                method.lower(),
            )
            is None
        ):
            inp.functional = method.upper()
        if inp.method.startswith("RO"):
            inp.unrestricted = False
        elif inp.method.startswith("U") and inp.method not in ("UFF",):
            inp.unrestricted = True
        elif inp.method.startswith("R"):
            inp.unrestricted = False

    @section(" Charge = ")
    def _charge_mult(self, line: str) -> None:
        m = re.search(r"Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)", line)
        if m and self.calc.input.charge is None:
            self.calc.input.charge = int(m.group(1))
            self.calc.input.multiplicity = int(m.group(2))
            self.calc.molecule.charge = int(m.group(1))
            self.calc.molecule.multiplicity = int(m.group(2))
        self.lines.next()

    @section(" Standard basis:")
    def _basis(self, line: str) -> None:
        if self.calc.input.basis is None:
            self.calc.input.basis = line.split(":", 1)[1].split("(")[0].strip()
        self.lines.next()

    @section(lambda s: s.strip().startswith("NBasis=") or "basis functions," in s)
    def _nbasis(self, line: str) -> None:
        m = re.search(r"NBasis=\s*(\d+)", line) or re.match(r"\s*(\d+)\s+basis functions", line)
        if m:
            self.calc.molecule.n_basis_functions = int(m.group(1))
        self.lines.next()

    @section(lambda s: "alpha electrons" in s and "beta electrons" in s)
    def _nelec(self, line: str) -> None:
        vals = [int(v) for v in re.findall(r"(\d+)\s+(?:alpha|beta)", line)]
        if len(vals) == 2:
            self.calc.molecule.n_electrons = sum(vals)
        self.lines.next()

    @section(" Stoichiometry")
    def _stoich(self, line: str) -> None:
        self.calc.input.extra.setdefault("stoichiometry", line.split()[-1])
        self.lines.next()

    @section(lambda s: "Framework group" in s or "Full point group" in s)
    def _pg(self, line: str) -> None:
        m = re.search(r"(?:Framework group|Full point group)\s+(\S+)", line)
        if m:
            self.calc.molecule.point_group = m.group(1).split("[")[0]
        self.lines.next()

    # -------------------------------------------------------------- geometry
    @section(lambda s: s.strip() in ("Input orientation:", "Standard orientation:", "Z-Matrix orientation:"))
    def _orientation(self, line: str) -> None:
        kind = line.strip().split()[0].lower()
        lines = self.lines
        lines.next()
        # header: dashes, 2 title lines, dashes
        if not lines.find(lambda s: s.strip().startswith("-----"), limit=3):
            return
        lines.next()
        lines.find(lambda s: s.strip().startswith("-----"), limit=4)
        lines.next()
        zs: list[int] = []
        coords: list[list[float]] = []
        while not lines.done and not lines.current.strip().startswith("-----"):
            toks = lines.next().split()
            if len(toks) < 6:
                continue
            zs.append(int(toks[1]))
            coords.append([to_float(t) for t in toks[-3:]])
        lines.next()
        geo = Geometry(coordinates=coords, label=kind)
        if zs and (not self.calc.molecule.atomic_numbers or kind == "standard"):
            # Dummy atoms have Z = -1 / 0; drop them from the molecule but keep index alignment.
            self.calc.molecule.atomic_numbers = zs
            self.calc.molecule.symbols = [z_to_symbol(z) if z > 0 else ("Bq" if z == 0 else "X") for z in zs]
        if kind == "standard":
            self.state["std_geoms"].append(geo)
        elif kind == "input":
            self.state["input_geoms"].append(geo)
        else:
            self.state["input_geoms"].append(geo)

    # -------------------------------------------------------------------- SCF
    @section(" SCF Done:")
    def _scf_done(self, line: str) -> None:
        m = _SCF_DONE_RE.search(line)
        lines = self.lines
        if m:
            label, e, n = m.group(1), to_float(m.group(2)), int(m.group(3))
            cyc = SCFCycle(energy=e, converged=True, n_iterations=n)
            nxt = lines.peek()
            # "NFock= 7  Conv=0.28D-09  -V/T= 2.0060"
            self.calc.scf_cycles.append(cyc)
            self.calc.energies.scf = e
            self.calc.energies.total = e
            self.calc.input.extra["scf_label"] = label
            self.state["pending_energy"] = e
            if self.calc.input.unrestricted is None:
                self.calc.input.unrestricted = label.startswith("U")
            if "-V/T" in nxt:
                lines.next()
        else:
            lines.next()
            return
        lines.next()

    @section("Convergence failure -- run terminated")
    def _scf_fail(self, line: str) -> None:
        self.calc.scf_cycles.append(SCFCycle(converged=False))
        self.calc.termination = Termination.ERROR
        self.calc.error_message = "SCF convergence failure"
        self.lines.next()

    @section(lambda s: "<S**2>=" in s and "Annihilation" not in s and "Excited State" not in s)
    def _s2(self, line: str) -> None:
        m = re.search(r"<S\*\*2>=\s*(-?\d+\.\d+)", line)
        if m and self.calc.scf_cycles:
            self.calc.scf_cycles[-1].spin_squared = float(m.group(1))
            s = (self.calc.molecule.multiplicity or 1) - 1
            self.calc.scf_cycles[-1].spin_squared_ideal = (s / 2) * (s / 2 + 1)
        self.lines.next()

    @section(lambda s: " E2 =" in s and "EUMP2" in s)
    def _mp2(self, line: str) -> None:
        m = re.search(r"E2\s*=\s*(\S+)\s+EUMP2\s*=\s*(\S+)", line)
        if m:
            self.calc.energies.mp2_correlation = to_float(m.group(1))
            self.calc.energies.mp2 = to_float(m.group(2))
            self.calc.energies.total = self.calc.energies.mp2
        self.lines.next()

    @section(lambda s: s.strip().startswith("CCSD(T)=") or " CCSD(T)= " in s)
    def _ccsd_t(self, line: str) -> None:
        m = re.search(r"CCSD\(T\)=\s*(-?\d+\.\d+(?:[DdEe][-+]?\d+)?)", line)
        if m:
            self.calc.energies.ccsd_t = to_float(m.group(1))
            self.calc.energies.total = self.calc.energies.ccsd_t
        self.lines.next()

    @section(lambda s: "E(CORR)=" in s and "E(Corr)" not in s or "Wavefunction amplitudes converged. E(Corr)=" in s)
    def _ccsd(self, line: str) -> None:
        m = re.search(r"E\(C[Oo][Rr][Rr]\)=\s*(\S+)", line)
        if m:
            self.calc.energies.ccsd = to_float(m.group(1))
            if self.calc.energies.ccsd_t is None:
                self.calc.energies.total = self.calc.energies.ccsd
        self.lines.next()

    @section(lambda s: "Dispersion energy=" in s or "R6Disp:  Grimme-D" in s or "Empirical dispersion" in s)
    def _disp(self, line: str) -> None:
        m = re.search(r"Dispersion energy=\s*(\S+)\s+Hartrees", line)
        if m:
            self.calc.energies.dispersion = to_float(m.group(1))
        self.lines.next()

    # --------------------------------------------------------------- orbitals
    @section(lambda s: " eigenvalues --" in s and ("Alpha" in s or "Beta" in s))
    def _eigen(self, line: str) -> None:
        lines = self.lines
        # A fresh block resets the stored orbitals (frequency jobs reprint them).
        if not self.state["in_eigen"]:
            self.state.update(orb_alpha=[], orb_beta=[], orb_occ_alpha=0, orb_occ_beta=0, in_eigen=True)
        while not lines.done and " eigenvalues --" in lines.current:
            cur = lines.next()
            vals = [float(v) for v in _EIGEN_RE.findall(cur.split("--", 1)[1])]
            spin = "orb_alpha" if "Alpha" in cur else "orb_beta"
            self.state[spin].extend(vals)
            if " occ." in cur:
                self.state["orb_occ_alpha" if spin == "orb_alpha" else "orb_occ_beta"] += len(vals)
        self.state["in_eigen"] = False
        self._store_orbitals()

    def _store_orbitals(self) -> None:
        a, b = self.state["orb_alpha"], self.state["orb_beta"]
        if not a:
            return
        na, nb = self.state["orb_occ_alpha"], self.state["orb_occ_beta"]
        energies = [a] + ([b] if b else [])
        occs = [[2.0 if not b else 1.0] * na + [0.0] * (len(a) - na)]
        homo = [na - 1]
        if b:
            occs.append([1.0] * nb + [0.0] * (len(b) - nb))
            homo.append(nb - 1)
        orb = Orbitals(energies=energies, occupations=occs, homo=homo)
        orb.homo_energy = [e[h] for e, h in zip(energies, homo) if 0 <= h < len(e)]
        orb.lumo_energy = [e[h + 1] for e, h in zip(energies, homo) if 0 <= h + 1 < len(e)]
        # Guess=Fragment prints each fragment's eigenvalues after the molecule's: keep the largest block
        # (ties -> latest, which is what a frequency job after an optimisation wants).
        if self.calc.orbitals is None or len(a) >= len(self.calc.orbitals.energies[0]):
            self.calc.orbitals = orb

    # ---------------------------------------------------------------- charges
    def _read_charge_table(self, n_values: int) -> tuple[list[float], list[float] | None]:
        lines = self.lines
        charges: list[float] = []
        spins: list[float] = []
        while not lines.done:
            toks = lines.current.split()
            if len(toks) >= 2 + 1 and toks[0].isdigit() and not toks[1].replace(".", "").replace("-", "").isdigit():
                vals = floats_in(" ".join(toks[2:]))
                if not vals:
                    break
                charges.append(vals[0])
                if n_values > 1 and len(vals) > 1:
                    spins.append(vals[1])
                lines.next()
            else:
                break
        return charges, (spins if spins else None)

    @section(
        lambda s: s.strip() in ("Mulliken charges:", "Mulliken charges and spin densities:", "Mulliken atomic charges:")
    )
    def _mulliken(self, line: str) -> None:
        spin = "spin" in line
        self.lines.next()
        self.lines.next()  # column header line
        charges, spins = self._read_charge_table(2 if spin else 1)
        self._set_pop("mulliken", charges, spins)

    @section(lambda s: s.strip() in ("APT charges:", "APT atomic charges:", "APT charges and spin densities:"))
    def _apt(self, line: str) -> None:
        self.lines.next()
        self.lines.next()
        charges, spins = self._read_charge_table(2 if "spin" in line else 1)
        self._set_pop("apt", charges, spins)

    @section(lambda s: s.strip().startswith("ESP charges") or s.strip().startswith("Charges from ESP fit"))
    def _esp(self, line: str) -> None:
        lines = self.lines
        lines.next()
        # "Charges from ESP fit, RMS=..." is followed by "Charge=  0.0 Dipole=..." then header "1"
        lines.find(lambda s: s.strip() == "1" or re.match(r"^\s+\d+\s*$", s) is not None, limit=3)
        lines.next()
        charges, _ = self._read_charge_table(1)
        method = (
            "chelpg"
            if "chelpg" in (self.calc.input.raw or "").lower()
            else "mk"
            if "mk" in (self.calc.input.raw or "").lower()
            else "esp"
        )
        self._set_pop(method, charges, None)

    @section(lambda s: s.strip().startswith("Hirshfeld charges") and "CM5" in s)
    def _hirshfeld(self, line: str) -> None:
        lines = self.lines
        lines.next()
        lines.next()  # column names
        hirsh: list[float] = []
        spins: list[float] = []
        cm5: list[float] = []
        while not lines.done:
            toks = lines.current.split()
            if len(toks) >= 7 and toks[0].isdigit():
                vals = floats_in(" ".join(toks[2:]))
                hirsh.append(vals[0])
                spins.append(vals[1])
                cm5.append(vals[-1])
                lines.next()
            else:
                break
        has_spin = any(abs(s) > 0 for s in spins) or (self.calc.molecule.multiplicity or 1) > 1
        self._set_pop("hirshfeld", hirsh, spins if has_spin else None)
        self._set_pop("cm5", cm5, None)

    @section("Summary of Natural Population Analysis:")
    def _npa(self, line: str) -> None:
        lines = self.lines
        lines.next()
        if not lines.find(lambda s: s.strip().startswith("Atom") and "Charge" in s, limit=6):
            return
        lines.next()
        lines.next()  # dashes
        charges: list[float] = []
        while not lines.done and not lines.current.strip().startswith("==="):
            toks = lines.next().split()
            if len(toks) >= 3 and toks[1].isdigit():
                charges.append(to_float(toks[2]))
        self._set_pop("npa", charges, None)

    def _set_pop(self, method: str, charges: list[float], spins: list[float] | None) -> None:
        if not charges:
            return
        # Guess=Fragment and ONIOM print further charge tables for sub-systems, with zeros for atoms outside
        # the fragment; keep the table with the most non-zero entries (ties -> the latest one).
        existing = next((p for p in self.calc.populations if p.method == method), None)
        if existing is not None:
            nz_new = sum(1 for c in charges if abs(c) > 1e-6)
            nz_old = sum(1 for c in existing.charges if abs(c) > 1e-6)
            if nz_new < nz_old:
                return
        self.calc.populations = [p for p in self.calc.populations if p.method != method]
        self.calc.populations.append(PopulationAnalysis(method=method, charges=charges, spins=spins))

    @section("Dipole moment (field-independent basis, Debye)")
    def _dipole(self, line: str) -> None:
        self.lines.next()
        vals = floats_in(self.lines.next())
        if len(vals) >= 4:
            mp = self.calc.multipoles or Multipoles()
            mp.dipole = vals[:3]
            mp.dipole_magnitude = vals[3]
            self.calc.multipoles = mp
        if "Quadrupole moment" in self.lines.current:
            self.lines.next()
            q1 = floats_in(self.lines.next())
            q2 = floats_in(self.lines.next())
            if len(q1) == 3 and len(q2) == 3:
                self.calc.multipoles.quadrupole = q1 + q2

    @section("Isotropic polarizability for W=")
    def _polar(self, line: str) -> None:
        vals = floats_in(line.split("W=", 1)[1])
        if len(vals) >= 2:
            self.calc.input.extra["isotropic_polarizability_bohr3"] = vals[1]
        self.lines.next()

    # -------------------------------------------------------------- gradient
    @section(lambda s: "Forces (Hartrees/Bohr)" in s)
    def _forces(self, line: str) -> None:
        lines = self.lines
        lines.next()
        lines.next()  # Number Number X Y Z
        lines.next()  # dashes
        grad: list[list[float]] = []
        while not lines.done and not lines.current.strip().startswith("-----"):
            toks = lines.next().split()
            if len(toks) >= 5:
                grad.append([-to_float(t) for t in toks[2:5]])
        if grad:
            self.calc.gradient = grad
            if JobType.GRADIENT not in self.calc.input.job_types and not (
                {JobType.OPTIMIZATION, JobType.FREQUENCY, JobType.TS_OPTIMIZATION} & set(self.calc.input.job_types)
            ):
                self.calc.input.job_types.append(JobType.GRADIENT)

    # ---------------------------------------------------------- optimisation
    @section(lambda s: "Item               Value     Threshold  Converged?" in s)
    def _opt_criteria(self, line: str) -> None:
        lines = self.lines
        lines.next()
        opt = self.calc.optimization or Optimization()
        step = OptimizationStep(step=len(opt.steps) + 1, energy=self.state.get("pending_energy"))
        for _ in range(4):
            cur = lines.current
            m = re.match(r"\s*(Maximum|RMS)\s+(Force|Displacement)\s+(\S+)\s+(\S+)\s+(YES|NO)", cur)
            if not m:
                break
            lines.next()
            which, what, val, thr, ok = m.groups()
            key = f"{which.lower()}_{what.lower()}"
            v = to_float(val) if "*" not in val else float("nan")
            if what == "Force":
                if which == "Maximum":
                    step.max_gradient = v
                else:
                    step.rms_gradient = v
            else:
                if which == "Maximum":
                    step.max_step = v
                else:
                    step.rms_step = v
            step.converged_flags[key] = ok == "YES"
            opt.criteria[key] = to_float(thr)
        if opt.steps:
            prev = opt.steps[-1].energy
            if prev is not None and step.energy is not None:
                step.delta_energy = step.energy - prev
        opt.steps.append(step)
        opt.n_steps = len(opt.steps)
        self.calc.optimization = opt

    @section(lambda s: s.strip().startswith("Optimization completed"))
    def _opt_done(self, line: str) -> None:
        opt = self.calc.optimization or Optimization()
        opt.converged = True
        self.calc.optimization = opt
        self.lines.next()

    @section(lambda s: s.strip().startswith("Optimization stopped"))
    def _opt_stopped(self, line: str) -> None:
        opt = self.calc.optimization or Optimization()
        opt.converged = False
        self.calc.optimization = opt
        self.lines.next()

    # ------------------------------------------------------------ frequencies
    @section(lambda s: s.strip().startswith("Harmonic frequencies (cm**-1)"))
    def _frequencies(self, line: str) -> None:
        lines = self.lines
        lines.next()
        v = Vibrations()
        syms: list[str] = []
        ir: list[float] = []
        raman: list[float] = []
        rmass: list[float] = []
        frc: list[float] = []
        modes: list[list[list[float]]] = []
        n_atoms = len(self.calc.molecule.atomic_numbers)

        def is_normal_freq(s: str) -> bool:
            return " Frequencies --" in s and "Frequencies ---" not in s

        # Skip the prose header and, with freq=hpmodes, the whole high-precision block that precedes
        # the standard 3-modes-per-line block.
        if not lines.find(lambda s: is_normal_freq(s) or "Thermochemistry" in s or "Normal termination" in s):
            return
        if not is_normal_freq(lines.current):
            return
        while not lines.done and is_normal_freq(lines.current):
            cur = lines.current
            sym_line = lines.peek(-1)  # symmetry labels; mode numbers are one line above that
            freqs = floats_in(cur.split("--", 1)[1])
            n_here = len(freqs)
            labels = sym_line.split()
            syms.extend(labels if len(labels) == n_here else ["?"] * n_here)
            v.frequencies.extend(freqs)
            lines.next()
            # property rows: "Red. masses --", "IR Inten    --", "IRC Coupling -" (one dash), ...
            while (
                not lines.done
                and re.match(r"^\s*[A-Za-z][A-Za-z. ()]{2,20}\s-{1,3}\s", lines.current)
                and "Atom" not in lines.current
            ):
                row = lines.next()
                tag, _, rest = (
                    re.split(r"(\s-{1,3}\s)", row, maxsplit=1) if re.search(r"\s-{1,3}\s", row) else (row, "", "")
                )
                vals = floats_in(rest)
                tag = tag.strip()
                if tag.startswith("Red. masses"):
                    rmass.extend(vals)
                elif tag.startswith("Frc consts"):
                    frc.extend(vals)
                elif tag.startswith("IR Inten"):
                    ir.extend(vals)
                elif tag.startswith("Raman Activ"):
                    raman.extend(vals)
            if lines.current.split()[:2] == ["Atom", "AN"]:
                lines.next()
                block: list[list[list[float]]] = [[] for _ in range(n_here)]
                while not lines.done:
                    toks = lines.current.split()
                    if len(toks) < 2 + 3 * n_here or not toks[0].isdigit():
                        break
                    vals = [to_float(t) for t in toks[2:]]
                    for k in range(n_here):
                        block[k].append(vals[3 * k : 3 * k + 3])
                    lines.next()
                modes.extend(block)
            # Next chunk = mode-number line, symmetry line, Frequencies line.
            if is_normal_freq(lines.peek(2)):
                lines.skip(2)
            elif is_normal_freq(lines.peek(1)):
                lines.skip(1)
            else:
                break
        if v.frequencies:
            v.symmetries = syms if len(syms) == len(v.frequencies) else None
            v.ir_intensities = ir if len(ir) == len(v.frequencies) else None
            v.raman_activities = raman if len(raman) == len(v.frequencies) else None
            v.reduced_masses = rmass if len(rmass) == len(v.frequencies) else None
            v.force_constants = frc if len(frc) == len(v.frequencies) else None
            if modes and len(modes) == len(v.frequencies) and (n_atoms == 0 or all(len(m) == n_atoms for m in modes)):
                v.normal_modes = modes
            v.n_imaginary = sum(1 for f in v.frequencies if f < 0)
            self.calc.vibrations = v
            if JobType.FREQUENCY not in self.calc.input.job_types:
                self.calc.input.job_types.append(JobType.FREQUENCY)

    @section(lambda s: s.strip().startswith("Temperature") and "Kelvin" in s and "Pressure" in s)
    def _thermo_tp(self, line: str) -> None:
        vals = floats_in(line)
        t = self._thermo()
        t.temperature, t.pressure = vals[0], vals[1]
        self.lines.next()

    def _thermo(self) -> Thermochemistry:
        v = self.calc.vibrations or Vibrations()
        if v.thermochemistry is None:
            v.thermochemistry = Thermochemistry()
        self.calc.vibrations = v
        return v.thermochemistry

    @section(" Rotational symmetry number")
    def _rot_sym(self, line: str) -> None:
        m = re.search(r"number\s+(\d+)", line)
        if m:
            self._thermo().rotational_symmetry_number = int(m.group(1))
        self.lines.next()

    @section(" Zero-point correction=")
    def _thermo_block(self, line: str) -> None:
        lines = self.lines
        t = self._thermo()
        e = self.calc.energies
        mapping = {
            "Zero-point correction": ("zero_point_energy", "zero_point"),
            "Thermal correction to Energy": ("thermal_energy_correction", "thermal_energy_correction"),
            "Thermal correction to Enthalpy": ("enthalpy_correction", "enthalpy_correction"),
            "Thermal correction to Gibbs Free Energy": ("gibbs_correction", "gibbs_correction"),
            "Sum of electronic and thermal Enthalpies": ("enthalpy", "enthalpy"),
            "Sum of electronic and thermal Free Energies": ("gibbs_free_energy", "gibbs"),
        }
        while (
            not lines.done
            and "=" in lines.current
            and ("correction" in lines.current or "Sum of electronic" in lines.current)
        ):
            cur = lines.next()
            key, _, rest = cur.partition("=")
            vals = floats_in(rest)
            if not vals:
                continue
            attrs = mapping.get(key.strip())
            if attrs:
                setattr(t, attrs[0], vals[0])
                setattr(e, attrs[1], vals[0])
            elif "zero-point Energies" in key:
                e.extra["e_plus_zpe"] = vals[0]
            elif "thermal Energies" in key:
                e.extra["e_plus_thermal"] = vals[0]
        # Entropy table: "Total   17.079   5.968   45.281" (kcal/mol, cal/mol/K, cal/mol/K)
        if lines.find(lambda s: s.strip().startswith("Total") and len(floats_in(s)) == 3, limit=6):
            vals = floats_in(lines.next())
            t.entropy = vals[2] / 1000.0 / HARTREE_TO_KCAL

    # -------------------------------------------------------- excited states
    @section(" Excited State ")
    def _excited(self, line: str) -> None:
        lines = self.lines
        m = _EXC_RE.search(line)
        if not m:
            lines.next()
            return
        idx, label, ev, nm, f, s2 = m.groups()
        mult = None
        low = label.lower()
        if low.startswith("singlet"):
            mult = "singlet"
        elif low.startswith("triplet"):
            mult = "triplet"
        elif low.startswith("doublet"):
            mult = "doublet"
        sym = label.split("-", 1)[1] if "-" in label else None
        st = ExcitedState(
            index=int(idx),
            energy_ev=float(ev),
            wavelength_nm=float(nm),
            oscillator_strength=float(f),
            multiplicity=mult,
            symmetry=sym if sym not in ("?Sym",) else None,
            spin_squared=float(s2) if s2 else None,
        )
        lines.next()
        while not lines.done:
            tm = _TRANS_RE.match(lines.current)
            if not tm:
                break
            a, sa, arrow, b, sb, coef = tm.groups()
            if arrow == "->":
                st.transitions.append(
                    Transition(
                        from_orbital=int(a) - 1,
                        to_orbital=int(b) - 1,
                        coefficient=float(coef),
                        from_spin=sa.lower() or None,
                        to_spin=sb.lower() or None,
                    )
                )
            lines.next()
        # A TD re-print (e.g. after geometry optimisation) replaces earlier states with the same index.
        self.calc.excited_states = [s for s in self.calc.excited_states if s.index != st.index] + [st]
        self.calc.excited_states.sort(key=lambda s: s.index)
        if JobType.EXCITED_STATES not in self.calc.input.job_types:
            self.calc.input.job_types.append(JobType.EXCITED_STATES)

    # ------------------------------------------------------------------ NMR
    @section(lambda s: "Isotropic =" in s and "Anisotropy =" in s)
    def _nmr(self, line: str) -> None:
        from ..schema import NMRShielding

        m = re.match(r"\s*(\d+)\s+(\w+)\s+Isotropic =\s*(\S+)\s+Anisotropy =\s*(\S+)", line)
        if m:
            idx = int(m.group(1)) - 1
            self.calc.nmr = [n for n in self.calc.nmr if n.atom_index != idx]
            self.calc.nmr.append(
                NMRShielding(atom_index=idx, isotropic=to_float(m.group(3)), anisotropy=to_float(m.group(4)))
            )
            if JobType.NMR not in self.calc.input.job_types:
                self.calc.input.job_types.append(JobType.NMR)
        self.lines.next()

    # ----------------------------------------------------------- termination
    @section(" Normal termination of Gaussian")
    def _normal(self, line: str) -> None:
        self.calc.termination = Termination.NORMAL
        self.lines.next()

    @section(lambda s: s.strip().startswith("Error termination"))
    def _error(self, line: str) -> None:
        self.calc.termination = Termination.ERROR
        prev = [self.lines.peek(-k).strip() for k in range(1, 6)]
        msg = next((p for p in prev if p and not p.startswith("Error termination") and "l9999" not in p), line.strip())
        self.calc.error_message = msg
        self.lines.next()

    @section(" Job cpu time:")
    def _cpu(self, line: str) -> None:
        t = self.calc.timing or Timing()
        t.cpu_seconds = (t.cpu_seconds or 0.0) + (_seconds(line) or 0.0)
        self.calc.timing = t
        self.lines.next()

    @section(" Elapsed time:")
    def _wall(self, line: str) -> None:
        t = self.calc.timing or Timing()
        t.wall_seconds = (t.wall_seconds or 0.0) + (_seconds(line) or 0.0)
        self.calc.timing = t
        self.lines.next()

    # -------------------------------------------------------------- finalize
    def finalize(self) -> None:
        super().finalize()
        c = self.calc
        std, inp = self.state["std_geoms"], self.state["input_geoms"]
        geoms = std if std else inp
        # Attach SCF energies to geometries in order; Gaussian prints geometry then SCF Done.
        energies = [cyc.energy for cyc in c.scf_cycles]
        if geoms and energies:
            # A frequency job re-prints the last geometry without a new optimisation point; de-duplicate.
            if len(energies) < len(geoms):
                # map energies to the first len(energies) geometries; leftovers inherit last energy
                for g, e in zip(geoms, energies):
                    g.energy = e
                for g in geoms[len(energies) :]:
                    g.energy = energies[-1]
            else:
                for g, e in zip(geoms, energies):
                    g.energy = e
        # A frequency (or any follow-up) step re-prints the converged geometry, usually re-oriented: drop a
        # trailing geometry whose interatomic distances match the previous one.
        while len(geoms) >= 2 and _same_structure(geoms[-1].coordinates, geoms[-2].coordinates):
            dup = geoms.pop()
            if geoms[-1].energy is None:
                geoms[-1].energy = dup.energy
        c.geometries = geoms
        # Dummy atoms (X, Z = -1) are geometry scaffolding and are removed; ghost atoms (Bq, Z = 0) carry basis
        # functions and appear in population tables, so they stay.
        if c.molecule.atomic_numbers and any(z < 0 for z in c.molecule.atomic_numbers):
            keep = [i for i, z in enumerate(c.molecule.atomic_numbers) if z >= 0]
            c.molecule.atomic_numbers = [c.molecule.atomic_numbers[i] for i in keep]
            c.molecule.symbols = [c.molecule.symbols[i] for i in keep]
            for g in c.geometries:
                g.coordinates = [g.coordinates[i] for i in keep]
            c.provenance.warn("molecule", "dummy atoms removed from geometry")
        if c.timing and c.timing.n_procs is None:
            c.timing.n_procs = c.input.nprocs
        if c.termination is Termination.INCOMPLETE and c.scf_cycles:
            c.provenance.warn("termination", "no 'Normal termination' line: job still running or killed")
