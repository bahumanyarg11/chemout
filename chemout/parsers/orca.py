"""ORCA 4.x–6.x output parser.

ORCA prints a block per geometry-optimisation cycle (coordinates, SCF, energy, gradient), so most
handlers simply append/overwrite and `finalize` ties geometries to energies by order.
"""

from __future__ import annotations

import re

from ..schema import (
    AU_TO_DEBYE,
    ExcitedState,
    Geometry,
    JobType,
    Multipoles,
    NMRShielding,
    Optimization,
    OptimizationStep,
    Orbitals,
    PopulationAnalysis,
    Program,
    SCFCycle,
    SCFIteration,
    Termination,
    Thermochemistry,
    Timing,
    Transition,
    Vibrations,
)
from .base import Parser, floats_in, section, symbol_to_z, to_float, z_to_symbol

_VERSION_RE = re.compile(r"Program Version\s+(\S+)")
_SCF_HEADER_RE = re.compile(r"^\s*(ITER|Iteration)\s+Energy")
_INPUT_LINE_RE = re.compile(r"^\|\s*\d+>\s?(.*)$")
_RUNTIME_RE = re.compile(r"(\d+)\s+days?\s+(\d+)\s+hours?\s+(\d+)\s+minutes?\s+(\d+)\s+seconds?\s+(\d+)\s+msec")
_STATE_RE = re.compile(
    r"STATE\s+(\d+):\s+E=\s*(-?\d+\.\d+)\s*au\s+(-?\d+\.\d+)\s*eV\s+(-?\d+\.\d+)\s*cm\*\*-1(?:\s+<S\*\*2>\s*=\s*(-?\d+\.\d+))?(?:\s+Mult\s+(\d+))?"
)
_TRANS_RE = re.compile(r"^\s*(\d+)([ab]?)\s*->\s*(\d+)([ab]?)\s*:\s*(-?\d+\.\d+)")
_SPEC_ROW6_RE = re.compile(r"^\s*(\d+)-(\d+)([A-Za-z0-9'\"]*)\s*->\s*(\d+)-(\d+)([A-Za-z0-9'\"]*)\s+(.*)$")

_JOB_KEYWORDS = {
    "opt": JobType.OPTIMIZATION,
    "tightopt": JobType.OPTIMIZATION,
    "loosopt": JobType.OPTIMIZATION,
    "looseopt": JobType.OPTIMIZATION,
    "verytightopt": JobType.OPTIMIZATION,
    "normalopt": JobType.OPTIMIZATION,
    "copt": JobType.OPTIMIZATION,
    "zopt": JobType.OPTIMIZATION,
    "l-opt": JobType.OPTIMIZATION,
    "optts": JobType.TS_OPTIMIZATION,
    "freq": JobType.FREQUENCY,
    "numfreq": JobType.FREQUENCY,
    "anfreq": JobType.FREQUENCY,
    "nmr": JobType.NMR,
    "engrad": JobType.GRADIENT,
    "numgrad": JobType.GRADIENT,
    "md": JobType.MD,
    "irc": JobType.IRC,
    "goat": JobType.OTHER,
    "sp": JobType.ENERGY,
    "energy": JobType.ENERGY,
}
_DISPERSION = {"d3", "d3bj", "d3zero", "d4", "d2", "d3tz", "nl", "d4(eeq)"}
_BASIS_RE = re.compile(
    r"^(ma-|aug-|x2c-|zora-|dkh-|sarc-|old-|ano-|cc-|def2-|def-|pc-|pcseg-|sto-|[0-9]-[0-9]+\+*g|dz|tz|qz|svp|tzvp|tzvpp|qzvp|qzvpp|"
    r"6-31|6-311|3-21|4-31|cc-pv|cc-pcv|cc-pwcv|aug-cc|jun-cc|may-cc|apr-cc|def2|ma-def2|cc-pvdz|cc-pvtz|cc-pvqz|cc-pv5z|ano-rcc|"
    r"ahgbs|hgbs|pcj|pcs|pcx|lanl2dz|crenbl|sdd|sapporo|ugbs|wachters|partridge|epr-|iglo|dgauss|minix|minis|midi|"
    r"ano-pv|sarc2-|def2-mtzvp|def2-mtzvpp|def2-msvp|ecp)",
    re.I,
)
_AUX_RE = re.compile(r"^(def2/j|def2/jk|def2-\w+/c|autoaux|cc-pv\wz/c|aug-cc-pv\wz/c|sarc/j|x2c/j|/c$)", re.I)
_MISC_KEYWORDS = {
    "rijcosx",
    "rijk",
    "ri",
    "nori",
    "rijonx",
    "ri-jk",
    "rij",
    "cosx",
    "noiter",
    "miniprint",
    "largeprint",
    "normalprint",
    "smallprint",
    "printbasis",
    "printmos",
    "noprintmos",
    "tightscf",
    "verytightscf",
    "normalscf",
    "loosescf",
    "sloppyscf",
    "strongscf",
    "extremescf",
    "slowconv",
    "veryslowconv",
    "noautostart",
    "autostart",
    "moread",
    "uno",
    "unrestricted",
    "uks",
    "uhf",
    "rks",
    "rhf",
    "roks",
    "rohf",
    "cpcm",
    "cpcmc",
    "kdiis",
    "soscf",
    "nososcf",
    "diis",
    "nodiis",
    "damp",
    "nodamp",
    "pal",
    "pal2",
    "pal4",
    "pal8",
    "pal16",
    "grid4",
    "grid5",
    "grid6",
    "grid7",
    "gridx4",
    "gridx5",
    "gridx6",
    "defgrid1",
    "defgrid2",
    "defgrid3",
    "nofinalgrid",
    "finalgrid5",
    "finalgrid6",
    "bohrs",
    "angs",
    "frozencore",
    "nofrozencore",
    "nopop",
    "mulliken",
    "loewdin",
    "mayer",
    "npa",
    "nbo",
    "hirshfeld",
    "chelpg",
    "allpop",
    "reducedpop",
    "xyzfile",
    "pdbfile",
    "gcp(dft/tz)",
    "gcp(dft/svp)",
    "gcp(hf/minis)",
    "gcp(dft/sv(p))",
    "gcp(hf/sv)",
    "libxc",
    "conv",
    "direct",
    "semidirect",
    "nouseshark",
    "useshark",
    "noshark",
    "quick",
    "sym",
    "nosym",
    "usesym",
    "nousesym",
    "numgrad",
    "aim",
    "allowrhf",
    "noallowrhf",
    "nmrshielding",
    "nmrspinspin",
    "nmrss",
    "tda",
    "notda",
    "decontract",
    "decontractbas",
    "decontractaux",
    "readints",
    "keepints",
    "keepdens",
    "nokeepints",
    "nokeepdens",
    "normalconv",
    "zora",
    "x2c",
    "dkh",
    "dkh2",
    "nozora",
    "relativistic",
}
_WFN_METHODS = {
    "hf",
    "rhf",
    "uhf",
    "rohf",
    "hf-3c",
    "mp2",
    "ri-mp2",
    "dlpno-mp2",
    "scs-mp2",
    "ri-scs-mp2",
    "oo-ri-mp2",
    "mp2-f12",
    "ccsd",
    "ccsd(t)",
    "dlpno-ccsd",
    "dlpno-ccsd(t)",
    "dlpno-ccsd(t1)",
    "ccsd(t)-f12",
    "ri-ccsd(t)",
    "cepa/1",
    "qcisd(t)",
    "casscf",
    "nevpt2",
    "sc-nevpt2",
    "fic-nevpt2",
    "caspt2",
    "mrci",
    "cis",
    "cisd",
    "zindo/s",
    "zindo/1",
    "zindo/2",
    "indo/s",
    "am1",
    "pm3",
    "mndo",
    "xtb",
    "xtb1",
    "xtb2",
    "gfn-xtb",
    "gfn2-xtb",
    "gfn-ff",
    "xtbff",
    "hfld",
    "domp2",
    "eom-ccsd",
    "bt-pno-eom-ccsd",
    "stoeom-ccsd",
    "adc(2)",
    "ri-adc2",
    "cepa",
    "ncpf/1",
    "qcisd",
    "mp3",
    "mp4",
}


def _looks_like_basis(tok: str) -> bool:
    low = tok.lower()
    return bool(_BASIS_RE.match(low)) and "/" not in low[1:]


class OrcaParser(Parser):
    program = Program.ORCA

    def __init__(self, lines):
        super().__init__(lines)
        self.state.update(
            opt_cycle=0,
            scf_open=False,
            final_point_energies=[],
            spec_seen=False,
            geom_pending_energy=None,
        )

    # ------------------------------------------------------------------ header
    @section("Program Version")
    def _version(self, line: str) -> None:
        m = _VERSION_RE.search(line)
        if m and self.calc.program.version is None:
            v = m.group(1)
            self.calc.program.version = v
            nums = [int(x) for x in re.findall(r"\d+", v.split("-")[0])[:3]]
            self.calc.program.version_tuple = tuple(nums) or None
        self.lines.next()

    @section(lambda s: s.strip() == "INPUT FILE")
    def _input(self, line: str) -> None:
        lines = self.lines
        lines.next()
        raw: list[str] = []
        while not lines.done:
            cur = lines.current
            m = _INPUT_LINE_RE.match(cur)
            if m:
                txt = m.group(1)
                if "****END OF INPUT****" in txt:
                    lines.next()
                    break
                raw.append(txt.rstrip())
                lines.next()
            elif cur.startswith("NAME =") or cur.startswith("====") or not cur.strip():
                if cur.startswith("NAME ="):
                    self.calc.input.title = cur.split("=", 1)[1].strip()
                lines.next()
            else:
                break
        self.calc.input.raw = "\n".join(raw)
        self._parse_input(raw)

    def _parse_input(self, raw: list[str]) -> None:
        inp = self.calc.input
        keyword_lines = [l.strip()[1:] for l in raw if l.strip().startswith("!")]
        tokens: list[str] = []
        for kl in keyword_lines:
            # keep parenthesised groups together: CPCM(WATER), gCP(DFT/TZ)
            tokens.extend(re.findall(r"[^\s()]+(?:\([^)]*\))?", kl))
        method: str | None = None
        for tok in tokens:
            low = tok.lower()
            if low in _JOB_KEYWORDS:
                jt = _JOB_KEYWORDS[low]
                if jt not in inp.job_types:
                    inp.job_types.append(jt)
                continue
            if low in _DISPERSION:
                inp.dispersion = tok.upper()
                continue
            if low.startswith(("cpcm", "smd", "cosmo", "alpb", "ddcosmo")):
                inp.solvation = tok.upper() if "(" in tok else tok.upper()
                continue
            if _AUX_RE.match(low) or low.endswith("/c") or low.endswith("/j") or low.endswith("/jk"):
                inp.aux_basis.append(tok)
                continue
            if low in _MISC_KEYWORDS or low.startswith(
                ("pal", "grid", "defgrid", "finalgrid", "moread", "gcp(", "pmodel", "hcore", "hueckel")
            ):
                if low in ("uks", "uhf", "unrestricted"):
                    inp.unrestricted = True
                elif low in ("rks", "rhf"):
                    inp.unrestricted = False
                continue
            if _looks_like_basis(low) and method is not None:
                inp.basis = inp.basis or tok
                continue
            if _looks_like_basis(low) and low not in _WFN_METHODS:
                inp.basis = inp.basis or tok
                continue
            if method is None:
                method = tok
            else:
                inp.extra.setdefault("other_keywords", []).append(tok)
        if method is None:
            # "! ROHF cc-pVQZ": the HF-type flag is the method when nothing else names one
            hf_flag = next((t for t in tokens if t.lower() in ("rhf", "uhf", "rohf", "hf")), None)
            if hf_flag:
                method = hf_flag.upper()
        if method:
            inp.method = method
            if method.lower() not in _WFN_METHODS and not re.match(r"^(r|u|ro)?hf$", method.lower()):
                inp.functional = method
            if method.lower().endswith("-3c"):
                inp.extra["composite_method"] = True
        # charge / multiplicity from the coordinate block
        for l in raw:
            m = re.match(r"^\s*\*\s*(xyz|int|gzmt|xyzfile|pdbfile|intfile|gzmtfile)\s+(-?\d+)\s+(\d+)", l, re.I)
            if m:
                inp.charge = int(m.group(2))
                inp.multiplicity = int(m.group(3))
                break
        joined = "\n".join(raw).lower()
        m = re.search(r"%pal\s+nprocs\s+(\d+)", joined)
        if m:
            inp.nprocs = int(m.group(1))
        m = re.search(r"%maxcore\s+(\d+)", joined)
        if m:
            inp.memory_mb = float(m.group(1)) * (inp.nprocs or 1)
        if re.search(r"%tddft|%cis\b", joined) and JobType.EXCITED_STATES not in inp.job_types:
            inp.job_types.append(JobType.EXCITED_STATES)
        if re.search(r"%geom[\s\S]*?scan", joined):
            inp.job_types.append(JobType.SCAN)
        m = re.search(r"smdsolvent\s+\"?([^\"\s]+)", joined)
        if m:
            inp.solvation = f"SMD({m.group(1)})"
        if not inp.job_types:
            inp.job_types.append(JobType.ENERGY)

    @section("Your calculation utilizes the basis:")
    def _basis_line(self, line: str) -> None:
        if self.calc.input.basis is None:
            self.calc.input.basis = line.split(":", 1)[1].strip()
        self.lines.next()

    @section(lambda s: s.startswith(" Hartree-Fock type") and "HFTyp" in s)
    def _hftyp(self, line: str) -> None:
        typ = line.split("....")[-1].strip()
        self.calc.input.unrestricted = typ.startswith("U")
        self.calc.input.extra["hf_type"] = typ
        self.lines.next()

    @section(lambda s: s.startswith(" Total Charge") and "Charge" in s)
    def _charge(self, line: str) -> None:
        self.calc.molecule.charge = int(floats_in(line.split("....")[-1])[0])
        self.lines.next()

    @section(lambda s: s.startswith(" Multiplicity") and "Mult" in s)
    def _mult(self, line: str) -> None:
        self.calc.molecule.multiplicity = int(floats_in(line.split("....")[-1])[0])
        self.lines.next()

    @section(lambda s: s.startswith(" Number of Electrons") and "NEL" in s)
    def _nel(self, line: str) -> None:
        self.calc.molecule.n_electrons = int(floats_in(line.split("....")[-1])[0])
        self.lines.next()

    @section(lambda s: s.lstrip().startswith("Number of basis functions") and "..." in s)
    def _nbasis(self, line: str) -> None:
        self.calc.molecule.n_basis_functions = int(floats_in(line.split("...")[-1])[0])
        self.lines.next()

    @section(lambda s: s.startswith(" Exchange Functional") or s.startswith(" Correlation Functional"))
    def _functional_parts(self, line: str) -> None:
        key = "exchange_functional" if "Exchange" in line else "correlation_functional"
        self.calc.input.extra[key] = line.split("....")[-1].strip()
        self.lines.next()

    @section(lambda s: s.startswith(" Density Functional     Method") or s.startswith(" Ab initio Hamiltonian  Method"))
    def _hamiltonian(self, line: str) -> None:
        self.calc.input.extra["hamiltonian"] = line.split("....")[-1].strip()
        self.lines.next()

    @section("Program running with")
    def _mpi(self, line: str) -> None:
        m = re.search(r"with\s+(\d+)\s+parallel", line)
        if m:
            self.calc.input.nprocs = self.calc.input.nprocs or int(m.group(1))
        self.lines.next()

    # ------------------------------------------------------------- geometry
    @section("GEOMETRY OPTIMIZATION CYCLE")
    def _opt_cycle(self, line: str) -> None:
        m = re.search(r"CYCLE\s+(\d+)", line)
        if m:
            self.state["opt_cycle"] = int(m.group(1))
        if (
            JobType.OPTIMIZATION not in self.calc.input.job_types
            and JobType.TS_OPTIMIZATION not in self.calc.input.job_types
        ):
            self.calc.input.job_types.append(JobType.OPTIMIZATION)
        self.lines.next()

    @section(lambda s: s.strip() == "CARTESIAN COORDINATES (ANGSTROEM)")
    def _coords(self, line: str) -> None:
        lines = self.lines
        lines.next()
        if lines.current.strip().startswith("---"):
            lines.next()
        zs: list[int] = []
        coords: list[list[float]] = []
        while not lines.done and lines.current.strip():
            toks = lines.next().split()
            if len(toks) < 4:
                break
            sym = toks[0]
            z = symbol_to_z(sym.lstrip("-").split(":")[0])  # "C:" point charges / "-C" ghost atoms
            if z is None:
                continue
            zs.append(z)
            coords.append([to_float(t) for t in toks[1:4]])
        if not coords:
            return
        if not self.calc.molecule.atomic_numbers:
            self.calc.molecule.atomic_numbers = zs
            self.calc.molecule.symbols = [z_to_symbol(z) for z in zs]
        cyc = self.state["opt_cycle"]
        label = f"opt cycle {cyc}" if cyc else ("input" if not self.calc.geometries else "geometry")
        self.calc.geometries.append(Geometry(coordinates=coords, label=label, step=cyc or None))

    # ------------------------------------------------------------------- SCF
    @section(lambda s: _SCF_HEADER_RE.match(s) is not None)
    def _scf_table(self, line: str) -> None:
        lines = self.lines
        lines.next()
        # A new SCF starts unless the current cycle is still running (DIIS -> SOSCF prints a second header
        # inside one SCF). Geometry cycles without a TOTAL SCF ENERGY block leave the cycle "open" but
        # already converged, so the convergence flag decides, not the open flag.
        cur_cycle = self.calc.scf_cycles[-1] if self.calc.scf_cycles else None
        if cur_cycle is None or not self.state["scf_open"] or cur_cycle.converged is not None:
            self.calc.scf_cycles.append(SCFCycle())
            self.state["scf_open"] = True
        cyc = self.calc.scf_cycles[-1]
        while not lines.done:
            cur = lines.current
            toks = cur.split()
            if toks and toks[0].isdigit() and len(toks) >= 3:
                try:
                    it = SCFIteration(energy=to_float(toks[1]), delta_energy=to_float(toks[2]))
                    if len(toks) >= 5:
                        it.rms_density = to_float(toks[3])
                        it.max_density = to_float(toks[4])
                except ValueError:
                    break
                cyc.iterations.append(it)
                lines.next()
            elif _SCF_HEADER_RE.match(cur):
                lines.next()  # continuation header (SOSCF switch)
            elif not cur.strip() or cur.strip().startswith(("***", "---", "*")) or "S-O-S-C-F" in cur:
                if "CONVERGED" in cur or "NOT CONVERGED" in cur:
                    break
                lines.next()
            else:
                break

    @section(lambda s: re.match(r"^\s*!\s+ITERATION\s+\d+\s+!", s) is not None)
    def _verbose_iteration(self, line: str) -> None:
        # Large-print SCF: one "! ITERATION n !" block per iteration instead of a table.
        lines = self.lines
        lines.next()
        cur_cycle = self.calc.scf_cycles[-1] if self.calc.scf_cycles else None
        if cur_cycle is None or not self.state["scf_open"] or cur_cycle.converged is not None:
            if not (cur_cycle is not None and self.state["scf_open"] and cur_cycle.converged is None):
                self.calc.scf_cycles.append(SCFCycle())
                self.state["scf_open"] = True
        cyc = self.calc.scf_cycles[-1]
        it = SCFIteration(energy=float("nan"))
        for _ in range(14):
            if lines.done or lines.current.strip().startswith(("---", "***", "!")):
                break
            s = lines.current.strip()
            vals = floats_in(s.split(":", 1)[1]) if ":" in s else []
            if s.startswith("Total Energy") and vals:
                it.energy = vals[0]
            elif s.startswith("Energy Change") and vals:
                it.delta_energy = vals[0]
            elif s.startswith("MAX-DP") and vals:
                it.max_density = vals[0]
            elif s.startswith("RMS-DP") and vals:
                it.rms_density = vals[0]
            lines.next()
        if it.energy == it.energy:  # not NaN
            cyc.iterations.append(it)

    def _cycle_for(self, *, flag: bool = False, energy: bool = False) -> SCFCycle:
        """The SCF cycle a convergence flag / final energy belongs to: the last one if it still lacks it."""
        cycles = self.calc.scf_cycles
        if cycles:
            last = cycles[-1]
            if (flag and last.converged is None) or (energy and last.energy is None):
                return last
        cycles.append(SCFCycle())
        return cycles[-1]

    @section("SCF CONVERGED AFTER")
    def _scf_converged(self, line: str) -> None:
        cyc = self._cycle_for(flag=True)
        cyc.converged = True
        m = re.search(r"AFTER\s+(\d+)\s+CYCLES", line)
        if m:
            cyc.n_iterations = int(m.group(1))
        if cyc.energy is None and cyc.iterations:
            # Optimisation cycles may skip the TOTAL SCF ENERGY block; FINAL SINGLE POINT ENERGY or the last
            # iteration is then the best available SCF energy (overwritten if the block does appear).
            cyc.energy = cyc.iterations[-1].energy
        self.lines.next()

    @section("SCF NOT CONVERGED AFTER")
    def _scf_not_converged(self, line: str) -> None:
        cyc = self._cycle_for(flag=True)
        cyc.converged = False
        m = re.search(r"AFTER\s+(\d+)\s+CYCLES", line)
        if m:
            cyc.n_iterations = int(m.group(1))
        self.lines.next()

    @section(lambda s: s.strip() == "TOTAL SCF ENERGY")
    def _total_scf(self, line: str) -> None:
        lines = self.lines
        lines.next()
        # The block belongs to the SCF still open (its energy so far is only the last-iteration fallback).
        if self.calc.scf_cycles and self.state["scf_open"]:
            cyc = self.calc.scf_cycles[-1]
        else:
            cyc = self._cycle_for(energy=True)
        self.state["scf_open"] = False
        if lines.current.strip().startswith("----"):
            lines.next()  # the dashes under the header
        for _ in range(60):
            if lines.done:
                break
            cur = lines.current
            s = cur.strip()
            if s.startswith("----") or s.startswith("****") or (s.isupper() and ":" not in s and len(s) > 8):
                break  # next section header
            if ":" in cur and "Eh" in cur:
                key, _, rest = cur.partition(":")
                vals = floats_in(rest)
                key = key.strip()
                if vals:
                    if key == "Total Energy":
                        cyc.energy = vals[0]
                        self.calc.energies.scf = vals[0]
                    elif key == "Nuclear Repulsion":
                        cyc.nuclear_repulsion = vals[0]
                    elif key == "One Electron Energy":
                        cyc.one_electron = vals[0]
                    elif key == "Two Electron Energy":
                        cyc.two_electron = vals[0]
                    elif key == "E(XC)":
                        cyc.exchange_correlation = vals[0]
                    elif key == "Dispersion correction":
                        cyc.dispersion = vals[0]
                        self.calc.energies.dispersion = vals[0]
            elif cur.startswith("Dispersion correction"):
                vals = floats_in(cur)
                if vals:
                    cyc.dispersion = vals[0]
                    self.calc.energies.dispersion = vals[0]
            lines.next()

    @section(lambda s: s.startswith("Dispersion correction") and len(s.split()) == 3)
    def _dispersion_line(self, line: str) -> None:
        # Standalone line closing the DFT DISPERSION CORRECTION section; not part of "Total Energy".
        vals = floats_in(line)
        if vals:
            self.calc.energies.dispersion = vals[0]
            if self.calc.scf_cycles:
                self.calc.scf_cycles[-1].dispersion = vals[0]
        self.lines.next()

    @section(lambda s: s.startswith("gCP correction") and len(s.split()) == 3)
    def _gcp_line(self, line: str) -> None:
        vals = floats_in(line)
        if vals:
            self.calc.energies.gcp = vals[0]
        self.lines.next()

    @section("Expectation value of <S**2>")
    def _s2(self, line: str) -> None:
        if self.calc.scf_cycles:
            self.calc.scf_cycles[-1].spin_squared = floats_in(line.split(":")[-1])[0]
            nxt = self.lines.peek()
            if "Ideal value" in nxt:
                self.calc.scf_cycles[-1].spin_squared_ideal = floats_in(nxt.split(":")[-1])[0]
        self.lines.next()

    @section("FINAL SINGLE POINT ENERGY")
    def _final_sp(self, line: str) -> None:
        vals = floats_in(line.split("ENERGY", 1)[1])
        if vals:
            e = vals[0]
            self.calc.energies.total = e
            self.state["final_point_energies"].append(e)
            if self.calc.geometries:
                self.calc.geometries[-1].energy = e
        self.lines.next()

    @section(
        lambda s: (
            s.startswith("E(MP2)")
            or s.startswith("E(CCSD)")
            or s.startswith("E(CCSD(T))")
            or s.startswith("E(SCF)")
            or s.startswith("E(CORR)")
        )
    )
    def _corr_energies(self, line: str) -> None:
        key, _, rest = line.partition("=")
        vals = floats_in(rest)
        if vals:
            k = key.strip().lower()
            e = self.calc.energies
            if k == "e(mp2)":
                e.mp2 = vals[0]
            elif k == "e(ccsd)":
                e.ccsd = vals[0]
            elif k == "e(ccsd(t))":
                e.ccsd_t = vals[0]
            elif k == "e(corr)":
                e.extra["correlation"] = vals[0]
            elif k == "e(scf)" and e.scf is None:
                e.scf = vals[0]
        self.lines.next()

    @section(lambda s: s.startswith("Total MP2 energy") or s.startswith("MP2 TOTAL ENERGY"))
    def _mp2_total(self, line: str) -> None:
        vals = floats_in(line.split(":")[-1] if ":" in line else line.split("=")[-1])
        if vals:
            self.calc.energies.mp2 = vals[0]
        self.lines.next()

    @section(lambda s: s.startswith("Final correlation energy") or s.startswith("E(CORR)(corrected)"))
    def _corr(self, line: str) -> None:
        vals = floats_in(line.split("...")[-1])
        if vals:
            self.calc.energies.extra["correlation"] = vals[0]
        self.lines.next()

    # -------------------------------------------------------------- orbitals
    @section(lambda s: s.strip() == "ORBITAL ENERGIES")
    def _orbitals(self, line: str) -> None:
        lines = self.lines
        lines.next()
        if lines.current.strip().startswith("---"):
            lines.next()
        spins: list[tuple[list[float], list[float]]] = []
        truncated = False
        while not lines.done:
            cur = lines.current.strip()
            if cur in ("SPIN UP ORBITALS", "SPIN DOWN ORBITALS") or cur.startswith("NO   OCC"):
                if cur.startswith("NO   OCC"):
                    lines.next()
                else:
                    lines.next()
                    if lines.current.strip().startswith("NO   OCC"):
                        lines.next()
                energies: list[float] = []
                occs: list[float] = []
                while not lines.done:
                    toks = lines.current.split()
                    if len(toks) >= 3 and toks[0].isdigit():
                        occs.append(to_float(toks[1]))
                        energies.append(to_float(toks[2]))
                        lines.next()
                    else:
                        break
                spins.append((energies, occs))
                if lines.current.startswith("*Only the first"):
                    truncated = True
                    lines.next()
            elif not cur:
                lines.next()
                if lines.current.strip() not in (
                    "SPIN UP ORBITALS",
                    "SPIN DOWN ORBITALS",
                ) and not lines.current.strip().startswith("NO   OCC"):
                    break
            else:
                break
        if not spins:
            return
        orb = Orbitals(energies=[s[0] for s in spins], occupations=[s[1] for s in spins], complete=not truncated)
        homo: list[int] = []
        for _energies, occs in spins:
            occupied = [i for i, o in enumerate(occs) if o > 0.5]
            homo.append(occupied[-1] if occupied else -1)
        orb.homo = homo
        orb.homo_energy = [e[h] for e, h in zip(orb.energies, homo) if 0 <= h < len(e)]
        orb.lumo_energy = [e[h + 1] for e, h in zip(orb.energies, homo) if 0 <= h + 1 < len(e)]
        orb.n_basis = self.calc.molecule.n_basis_functions
        self.calc.orbitals = orb

    # --------------------------------------------------------------- charges
    @section(lambda s: s.startswith(("MULLIKEN ATOMIC CHARGES", "LOEWDIN ATOMIC CHARGES")) and "REDUCED" not in s)
    def _charges(self, line: str) -> None:
        method = "mulliken" if line.startswith("MULLIKEN") else "loewdin"
        has_spin = "SPIN" in line
        lines = self.lines
        lines.next()
        if lines.current.strip().startswith("---"):
            lines.next()
        charges: list[float] = []
        spins: list[float] = []
        while not lines.done:
            cur = lines.current
            if ":" not in cur or cur.startswith("Sum of"):
                break
            left, _, rest = cur.partition(":")
            if not left.split() or not left.split()[0].isdigit():
                break
            vals = floats_in(rest)
            if not vals:
                break
            charges.append(vals[0])
            if has_spin and len(vals) > 1:
                spins.append(vals[1])
            lines.next()
        if charges:
            self.calc.populations = [p for p in self.calc.populations if p.method != method]
            self.calc.populations.append(PopulationAnalysis(method=method, charges=charges, spins=spins or None))

    @section(lambda s: s.strip() == "HIRSHFELD ANALYSIS")
    def _hirshfeld(self, line: str) -> None:
        lines = self.lines
        lines.next()
        if not lines.find(lambda s: s.split()[:2] == ["ATOM", "CHARGE"], limit=12):
            return
        lines.next()
        charges: list[float] = []
        spins: list[float] = []
        while not lines.done:
            toks = lines.current.split()
            if len(toks) >= 3 and toks[0].isdigit():
                charges.append(to_float(toks[2]))
                if len(toks) >= 4:
                    spins.append(to_float(toks[3]))
                lines.next()
            else:
                break
        if charges:
            self.calc.populations = [p for p in self.calc.populations if p.method != "hirshfeld"]
            self.calc.populations.append(
                PopulationAnalysis(
                    method="hirshfeld", charges=charges, spins=spins if any(abs(s) > 0 for s in spins) else None
                )
            )

    @section(lambda s: s.strip() == "MAYER POPULATION ANALYSIS")
    def _mayer(self, line: str) -> None:
        lines = self.lines
        lines.next()
        if not lines.find(lambda s: s.split()[:3] == ["ATOM", "NA", "ZA"], limit=20):
            return
        lines.next()
        charges: list[float] = []
        while not lines.done:
            toks = lines.current.split()
            if len(toks) >= 5 and toks[0].isdigit():
                charges.append(to_float(toks[4]))
                lines.next()
            else:
                break
        if charges:
            self.calc.populations = [p for p in self.calc.populations if p.method != "mayer"]
            self.calc.populations.append(PopulationAnalysis(method="mayer", charges=charges))

    @section(lambda s: s.strip() == "CHELPG Charges")
    def _chelpg(self, line: str) -> None:
        lines = self.lines
        lines.next()
        if lines.current.strip().startswith("---"):
            lines.next()
        charges: list[float] = []
        while not lines.done:
            cur = lines.current
            toks = cur.split()
            if len(toks) >= 4 and toks[0].isdigit() and ":" in cur:
                charges.append(floats_in(cur.split(":")[-1])[0])
                lines.next()
            else:
                break
        if charges:
            self.calc.populations = [p for p in self.calc.populations if p.method != "chelpg"]
            self.calc.populations.append(PopulationAnalysis(method="chelpg", charges=charges))

    @section(lambda s: s.strip() == "DIPOLE MOMENT")
    def _dipole(self, line: str) -> None:
        lines = self.lines
        lines.next()
        mp = Multipoles()
        for _ in range(25):
            if lines.done:
                break
            cur = lines.current
            if cur.startswith("Total Dipole Moment"):
                vals = floats_in(cur.split(":")[-1])
                if len(vals) == 3:
                    mp.dipole = [v * AU_TO_DEBYE for v in vals]
            elif cur.startswith("Magnitude (Debye)"):
                mp.dipole_magnitude = floats_in(cur.split(":")[-1])[0]
                lines.next()
                break
            elif cur.strip().startswith("----") and "Rotational spectrum" in lines.peek(1):
                break
            lines.next()
        if mp.dipole is not None:
            if mp.dipole_magnitude is None:
                mp.dipole_magnitude = sum(v * v for v in mp.dipole) ** 0.5
            self.calc.multipoles = mp

    # -------------------------------------------------------------- gradient
    @section(lambda s: s.strip() == "CARTESIAN GRADIENT")
    def _gradient(self, line: str) -> None:
        lines = self.lines
        lines.next()
        if lines.current.strip().startswith("---"):
            lines.next()
        lines.skip_blank()
        grad: list[list[float]] = []
        while not lines.done and lines.current.strip():
            cur = lines.current
            if ":" not in cur:
                break
            vals = floats_in(cur.split(":")[-1])
            if len(vals) < 3:
                break
            grad.append(vals[:3])
            lines.next()
        if grad:
            self.calc.gradient = grad

    @section("Norm of the Cartesian gradient")
    def _gnorm(self, line: str) -> None:
        self.calc.energies.extra["gradient_norm"] = floats_in(line.split("...")[-1])[0]
        self.lines.next()

    @section("|Geometry convergence|")
    def _geom_conv(self, line: str) -> None:
        lines = self.lines
        lines.next()
        lines.next()  # Item value Tolerance Converged
        lines.next()  # dashes
        opt = self.calc.optimization or Optimization()
        step = OptimizationStep(step=self.state["opt_cycle"] or len(opt.steps) + 1)
        if self.state["final_point_energies"]:
            step.energy = self.state["final_point_energies"][-1]
            if len(self.state["final_point_energies"]) > 1:
                step.delta_energy = step.energy - self.state["final_point_energies"][-2]
        while not lines.done:
            cur = lines.current
            m = re.match(
                r"\s*(Energy change|RMS gradient|MAX gradient|RMS step|MAX step)\s+(-?\d+\.\d+)\s+(\d+\.\d+)"
                r"\s+(YES|NO)",
                cur,
            )
            if not m:
                break
            name, val, tol, ok = m.groups()
            key = name.lower().replace(" ", "_")
            v = float(val)
            if key == "rms_gradient":
                step.rms_gradient = v
            elif key == "max_gradient":
                step.max_gradient = v
            elif key == "rms_step":
                step.rms_step = v
            elif key == "max_step":
                step.max_step = v
            elif key == "energy_change":
                step.delta_energy = v
            step.converged_flags[key] = ok == "YES"
            opt.criteria[key] = float(tol)
            lines.next()
        opt.steps.append(step)
        opt.n_steps = len(opt.steps)
        self.calc.optimization = opt

    @section("THE OPTIMIZATION HAS CONVERGED")
    def _opt_converged(self, line: str) -> None:
        opt = self.calc.optimization or Optimization()
        opt.converged = True
        self.calc.optimization = opt
        self.lines.next()

    @section(lambda s: "The optimization did not converge" in s or "OPTIMIZATION RUN DONE WITHOUT CONVERGENCE" in s)
    def _opt_not_converged(self, line: str) -> None:
        opt = self.calc.optimization or Optimization()
        opt.converged = False
        self.calc.optimization = opt
        self.lines.next()

    @section("FINAL ENERGY EVALUATION AT THE STATIONARY POINT")
    def _stationary(self, line: str) -> None:
        self.state["opt_cycle"] = 0
        self.state["final_eval"] = True
        self.lines.next()

    # ------------------------------------------------------------ vibrations
    @section(lambda s: s.strip() == "VIBRATIONAL FREQUENCIES")
    def _frequencies(self, line: str) -> None:
        lines = self.lines
        lines.next()
        freqs: list[float] = []
        imag_flags: list[bool] = []
        for _ in range(8):  # skip dashes, blank, scaling-factor note
            if lines.done or re.match(r"^\s*\d+:\s+-?\d+\.\d+\s+cm\*\*-1", lines.current):
                break
            lines.next()
        while not lines.done:
            m = re.match(r"^\s*(\d+):\s+(-?\d+\.\d+)\s+cm\*\*-1(.*)$", lines.current)
            if not m:
                break
            f = float(m.group(2))
            freqs.append(f)
            imag_flags.append("imaginary" in m.group(3))
            lines.next()
        v = self.calc.vibrations or Vibrations()
        keep = [i for i, f in enumerate(freqs) if abs(f) > 1e-8]
        self.state["vib_keep"] = keep
        v.frequencies = [freqs[i] for i in keep]
        v.n_imaginary = sum(1 for i in keep if freqs[i] < 0 or imag_flags[i])
        self.calc.vibrations = v
        if JobType.FREQUENCY not in self.calc.input.job_types:
            self.calc.input.job_types.append(JobType.FREQUENCY)

    @section(lambda s: s.strip() == "NORMAL MODES")
    def _normal_modes(self, line: str) -> None:
        lines = self.lines
        lines.next()
        columns: dict[int, list[float]] = {}
        current_cols: list[int] = []
        started = False
        while not lines.done:
            cur = lines.current
            toks = cur.split()
            if not toks:
                if (
                    started
                    and lines.peek(1).split()
                    and not all(t.lstrip("-").isdigit() for t in lines.peek(1).split())
                ):
                    # blank followed by something that is not a column header or a row -> end
                    nxt = lines.peek(1).split()
                    if not (nxt and nxt[0].isdigit()):
                        lines.next()
                        break
                lines.next()
                continue
            if all(t.isdigit() for t in toks):
                current_cols = [int(t) for t in toks]
                for c in current_cols:
                    columns.setdefault(c, [])
                started = True
                lines.next()
                continue
            if toks[0].isdigit() and started and len(toks) == len(current_cols) + 1:
                vals = [to_float(t) for t in toks[1:]]
                for c, vv in zip(current_cols, vals):
                    columns[c].append(vv)
                lines.next()
                continue
            if cur.startswith("Point group") or (started and all(re.match(r"^\d+-\S+$", t) for t in toks)):
                lines.next()  # ORCA 6 prints the point group and an irrep label per mode column
                continue
            if started:
                break
            lines.next()
        n_atoms = len(self.calc.molecule.atomic_numbers)
        if not columns or n_atoms == 0:
            return
        keep = self.state.get("vib_keep") or sorted(columns)
        modes: list[list[list[float]]] = []
        for c in keep:
            col = columns.get(c)
            if col is None or len(col) != 3 * n_atoms:
                self.calc.provenance.warn(
                    "normal_modes", f"mode {c}: expected {3 * n_atoms} components, got {0 if col is None else len(col)}"
                )
                return
            modes.append([col[3 * a : 3 * a + 3] for a in range(n_atoms)])
        v = self.calc.vibrations or Vibrations()
        v.normal_modes = modes
        self.calc.vibrations = v

    @section(lambda s: s.strip() == "IR SPECTRUM")
    def _ir(self, line: str) -> None:
        lines = self.lines
        lines.next()
        if not lines.find(lambda s: s.strip().startswith("Mode") and "freq" in s, limit=6):
            return
        header = re.sub(r"\([^)]*\)", "", lines.next()).split()
        # ORCA>=5: Mode freq eps Int T**2 TX TY TZ ; ORCA 4: Mode freq (cm**-1) T**2 TX TY TZ
        if "Int" in header:
            int_col = header.index("Int")
        elif "T**2" in header:
            int_col = header.index("T**2")
        else:
            return
        # ORCA >= 5 adds a units line ("cm**-1  L/(mol*cm) km/mol  a.u.") before the dashes.
        for _ in range(3):
            if re.match(r"^\s*\d+:\s", lines.current):
                break
            lines.next()
        ints: dict[int, float] = {}
        while not lines.done:
            cur = lines.current
            m = re.match(r"^\s*(\d+):\s+(.*)$", cur)
            if not m:
                break
            vals = floats_in(m.group(2).replace("(", " ").replace(")", " "))
            # header has 'Mode' as col 0; values start at 'freq' => value index = int_col - 1
            idx = int_col - 1
            if len(vals) > idx:
                ints[int(m.group(1))] = vals[idx]
            lines.next()
        v = self.calc.vibrations or Vibrations()
        keep = self.state.get("vib_keep")
        if keep and all(k in ints for k in keep):
            v.ir_intensities = [ints[k] for k in keep]
        elif ints:
            v.ir_intensities = [ints[k] for k in sorted(ints)]
        self.calc.vibrations = v

    @section(lambda s: s.startswith("THERMOCHEMISTRY AT"))
    def _thermo(self, line: str) -> None:
        lines = self.lines
        lines.next()
        v = self.calc.vibrations or Vibrations()
        t = v.thermochemistry or Thermochemistry()
        e = self.calc.energies
        electronic = None

        def last(sv: str) -> float | None:
            vals = floats_in(sv.split("...")[-1])
            return vals[0] if vals else None

        for _ in range(120):
            if lines.done:
                break
            cur = lines.current
            s = cur.strip()
            if s.startswith("Temperature"):
                t.temperature = last(s)
            elif s.startswith("Pressure"):
                t.pressure = last(s)
            elif s.startswith("Electronic energy"):
                electronic = last(s)
            elif s.startswith("Zero point energy"):
                t.zero_point_energy = last(s)
                e.zero_point = t.zero_point_energy
            elif s.startswith("Total correction"):
                t.thermal_energy_correction = (floats_in(s) or [None])[0]
                e.thermal_energy_correction = t.thermal_energy_correction
            elif s.startswith("Total Enthalpy"):
                t.enthalpy = last(s)
                e.enthalpy = t.enthalpy
            elif s.startswith("Final entropy term") or s.startswith("Total entropy correction"):
                ts = last(s)
                if t.temperature and ts is not None:
                    t.entropy = ts / t.temperature
            elif s.startswith("Final Gibbs free energy") or s.startswith("Final Gibbs free enthalpy"):
                t.gibbs_free_energy = last(s)
                e.gibbs = t.gibbs_free_energy
            elif s.startswith("G-E(el)"):
                t.gibbs_correction = last(s)
                e.gibbs_correction = t.gibbs_correction
            elif s.startswith("Total free energy") and t.gibbs_free_energy is None and "Eh" in s:
                # ORCA <5 prints "Total free energy" meaning E + thermal corrections (not G)
                if last(s) is not None:
                    e.extra["total_free_energy_eh"] = last(s)
            elif s.startswith("Rotational symmetry number") or s.startswith("sn="):
                m = re.search(r"(\d+)", s.split("...")[-1])
                if m:
                    t.rotational_symmetry_number = int(m.group(1))
            elif (
                s.startswith("Timings for individual modules")
                or s.startswith("SUGGESTED CITATIONS")
                or s.startswith("****ORCA")
            ):
                break
            lines.next()
        if electronic is not None and t.enthalpy is not None:
            t.enthalpy_correction = t.enthalpy - electronic
            e.enthalpy_correction = t.enthalpy_correction
        v.thermochemistry = t
        self.calc.vibrations = v

    # -------------------------------------------------------- excited states
    @section(lambda s: "EXCITED STATES" in s and s.strip().endswith(")") and ("SINGLETS" in s or "TRIPLETS" in s))
    def _excited_block(self, line: str) -> None:
        # "TD-DFT/TDA EXCITED STATES (SINGLETS)" then "... (TRIPLETS)": both blocks restart at STATE 1, while
        # the spectrum numbers them consecutively. Offset the triplet indices so the two sets line up.
        if "TRIPLETS" in line:
            self.state["exc_offset"] = len(self.calc.excited_states)
            self.state["exc_mult"] = "triplet"
        else:
            self.calc.excited_states = []  # a new TD-DFT run (e.g. after an optimisation cycle)
            self.state["exc_offset"] = 0
            self.state["exc_mult"] = "singlet"
        self.lines.next()

    @section(lambda s: re.match(r"^STATE\s+\d+:\s+E=", s) is not None)
    def _state(self, line: str) -> None:
        lines = self.lines
        m = _STATE_RE.search(line)
        if not m:
            lines.next()
            return
        idx, e_au, e_ev, e_cm, s2, mult = m.groups()
        # the cm**-1 value carries one more significant digit than the eV one
        energy_ev = float(e_cm) / 8065.543937 if float(e_cm) else float(e_ev)
        multiplicity = (
            {"1": "singlet", "2": "doublet", "3": "triplet"}.get(mult) if mult else self.state.get("exc_mult")
        )
        st = ExcitedState(
            index=int(idx) + self.state.get("exc_offset", 0),
            energy_ev=energy_ev,
            wavelength_nm=1e7 / float(e_cm) if float(e_cm) else None,
            spin_squared=float(s2) if s2 else None,
            multiplicity=multiplicity,
        )
        lines.next()
        while not lines.done:
            tm = _TRANS_RE.match(lines.current)
            if not tm:
                break
            a, sa, b, sb, coef = tm.groups()
            st.transitions.append(
                Transition(
                    from_orbital=int(a),
                    to_orbital=int(b),
                    coefficient=float(coef),
                    from_spin=sa or None,
                    to_spin=sb or None,
                )
            )
            lines.next()
        # states are re-printed per optimisation cycle: latest wins
        self.calc.excited_states = [s for s in self.calc.excited_states if s.index != st.index] + [st]
        self.calc.excited_states.sort(key=lambda s: s.index)
        if JobType.EXCITED_STATES not in self.calc.input.job_types:
            self.calc.input.job_types.append(JobType.EXCITED_STATES)

    @section(lambda s: s.strip().startswith("ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"))
    def _absorption(self, line: str) -> None:
        lines = self.lines
        lines.next()
        if not lines.find(lambda s: s.strip().startswith(("Transition", "State")), limit=4):
            return
        header = lines.next().split()
        orca6 = header[0] == "Transition"
        # units line + dashes
        while not lines.done and (not lines.current.strip() or lines.current.strip().startswith(("(", "---"))):
            lines.next()
        by_index = {s.index: s for s in self.calc.excited_states}
        n = 0
        while not lines.done:
            cur = lines.current
            if not cur.strip():
                break
            forbidden = "spin forbidden" in cur  # "... spin forbidden (mult=3)" replaces fosc and the dipole
            mf = re.search(r"mult=(\d)", cur)
            if orca6:
                m = _SPEC_ROW6_RE.match(cur)
                if not m:
                    break
                idx = int(m.group(4))
                vals = floats_in(m.group(7).split("spin")[0])
                if len(vals) < 3:
                    break
                e_ev, e_cm, wl = vals[0], vals[1], vals[2]
                fosc = 0.0 if forbidden else (vals[3] if len(vals) > 3 else None)
                sym = m.group(6) or None
                mult_from_label = mf.group(1) if mf else m.group(5)
            else:
                toks = cur.split()
                if not toks[0].isdigit():
                    break
                vals = floats_in(" ".join(toks[1:]).split("spin")[0])
                if len(vals) < 2:
                    break
                idx = int(toks[0])
                e_cm, wl = vals[0], vals[1]
                fosc = 0.0 if forbidden else (vals[2] if len(vals) > 2 else None)
                e_ev = e_cm / 8065.543937
                sym = None
                mult_from_label = mf.group(1) if mf else None
            n += 1
            st = by_index.get(idx)
            if st is None:
                st = ExcitedState(index=idx, energy_ev=e_ev)
                self.calc.excited_states.append(st)
                by_index[idx] = st
            st.oscillator_strength = fosc
            st.wavelength_nm = wl
            if sym and st.symmetry is None:
                st.symmetry = sym
            if mult_from_label and st.multiplicity is None:
                st.multiplicity = {"1": "singlet", "2": "doublet", "3": "triplet"}.get(mult_from_label)
            lines.next()
        self.calc.excited_states.sort(key=lambda s: s.index)
        if n and JobType.EXCITED_STATES not in self.calc.input.job_types:
            self.calc.input.job_types.append(JobType.EXCITED_STATES)

    @section(lambda s: "CD SPECTRUM" in s and "VIA TRANSITION ELECTRIC DIPOLE" in s)
    def _cd(self, line: str) -> None:
        lines = self.lines
        lines.next()
        if not lines.find(lambda s: s.strip().startswith(("Transition", "State")), limit=4):
            return
        orca6 = lines.next().split()[0] == "Transition"
        while not lines.done and (not lines.current.strip() or lines.current.strip().startswith(("(", "---"))):
            lines.next()
        by_index = {s.index: s for s in self.calc.excited_states}
        while not lines.done and lines.current.strip():
            cur = lines.current
            if orca6:
                m = _SPEC_ROW6_RE.match(cur)
                if not m:
                    break
                idx, vals = int(m.group(4)), floats_in(m.group(7))
                r = vals[3] if len(vals) > 3 else None
            else:
                toks = cur.split()
                if not toks[0].isdigit():
                    break
                idx, vals = int(toks[0]), floats_in(" ".join(toks[1:]))
                r = vals[2] if len(vals) > 2 else None
            if idx in by_index and r is not None:
                by_index[idx].rotatory_strength = r
            lines.next()

    # ------------------------------------------------------------------- NMR
    @section("CHEMICAL SHIELDING SUMMARY (ppm)")
    def _nmr(self, line: str) -> None:
        lines = self.lines
        lines.next()
        if not lines.find(lambda s: "Nucleus" in s and "Isotropic" in s, limit=6):
            return
        lines.next()
        if lines.current.strip().startswith("---"):
            lines.next()
        out: list[NMRShielding] = []
        while not lines.done:
            toks = lines.current.split()
            if len(toks) >= 4 and toks[0].isdigit():
                out.append(
                    NMRShielding(atom_index=int(toks[0]), isotropic=to_float(toks[2]), anisotropy=to_float(toks[3]))
                )
                lines.next()
            else:
                break
        if out:
            self.calc.nmr = out
            if JobType.NMR not in self.calc.input.job_types:
                self.calc.input.job_types.append(JobType.NMR)

    # ----------------------------------------------------------- termination
    @section("****ORCA TERMINATED NORMALLY****")
    def _normal(self, line: str) -> None:
        self.calc.termination = Termination.NORMAL
        self.lines.next()

    @section(
        lambda s: (
            "ORCA finished by error termination" in s
            or "ORCA finished with error" in s
            or s.startswith("[file orca")
            and "Error" in s
            or "aborting the run" in s.lower()
        )
    )
    def _error(self, line: str) -> None:
        self.calc.termination = Termination.ERROR
        prev = [self.lines.peek(-k).strip() for k in range(1, 8)]
        detail = next((p for p in prev if p and not p.startswith(("---", "***", "[file"))), "")
        msg = line.strip()
        if detail and detail not in msg:
            msg = f"{msg} | {detail}"
        self.calc.error_message = msg if not self.calc.error_message else self.calc.error_message
        self.lines.next()

    @section("TOTAL RUN TIME:")
    def _runtime(self, line: str) -> None:
        m = _RUNTIME_RE.search(line)
        if m:
            d, h, mi, s, ms = (int(x) for x in m.groups())
            t = self.calc.timing or Timing()
            t.wall_seconds = d * 86400 + h * 3600 + mi * 60 + s + ms / 1000.0
            t.n_procs = self.calc.input.nprocs
            self.calc.timing = t
        self.lines.next()

    # -------------------------------------------------------------- finalize
    def finalize(self) -> None:
        super().finalize()
        c = self.calc
        if c.molecule.charge is None and c.input.charge is not None:
            c.molecule.charge = c.input.charge
        if c.molecule.multiplicity is None and c.input.multiplicity is not None:
            c.molecule.multiplicity = c.input.multiplicity
        if c.input.charge is None and c.molecule.charge is not None:
            c.input.charge = c.molecule.charge
        if c.input.multiplicity is None and c.molecule.multiplicity is not None:
            c.input.multiplicity = c.molecule.multiplicity
        if c.energies.total is None and c.energies.scf is not None:
            c.energies.total = c.energies.scf
        # Geometries that never received a FINAL SINGLE POINT ENERGY get the matching SCF energy.
        scf_energies = [cyc.energy for cyc in c.scf_cycles if cyc.energy is not None]
        if c.geometries and all(g.energy is None for g in c.geometries) and len(scf_energies) >= len(c.geometries):
            for g, e in zip(c.geometries, scf_energies):
                g.energy = e
        # Optimisations re-print the converged geometry: drop an exact duplicate at the end.
        if len(c.geometries) >= 2 and c.geometries[-1].coordinates == c.geometries[-2].coordinates:
            last = c.geometries.pop()
            if c.geometries[-1].energy is None:
                c.geometries[-1].energy = last.energy
            c.geometries[-1].label = "final"
        if c.termination is Termination.INCOMPLETE and c.scf_cycles:
            c.provenance.warn(
                "termination", "no 'ORCA TERMINATED NORMALLY' line: job still running, killed, or truncated"
            )
        if c.energies.total is None and c.energies.scf is None:
            c.provenance.warn("energies", "no total or SCF energy found in this file")
        if c.input.method is None:
            # No method on the "!" line (e.g. "! TightSCF STO-3G" with a %method block): use what ORCA ran.
            ham = c.input.extra.get("hamiltonian", "")
            if ham.startswith("Hartree-Fock"):
                c.input.method = "HF"
            elif ham.startswith("DFT"):
                c.input.method = "DFT"
                xc = (c.input.extra.get("exchange_functional"), c.input.extra.get("correlation_functional"))
                c.input.functional = "/".join(x for x in xc if x) or None
