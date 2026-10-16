"""xTB (Grimme group) output parser; tested on 6.5.x and 6.7.x, GFN1/GFN2/GFN-FF.

xtb prints almost no echo of its input: the starting geometry only appears for optimisations
(as ``final structure``) or with ``--verbose``. Everything else is parsed from the SETUP block,
the SCC table, the ``SUMMARY`` box and the property printout.
"""

from __future__ import annotations

import re

from ..schema import (
    Energies,
    Geometry,
    JobType,
    Multipoles,
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
    Vibrations,
)
from .base import Parser, floats_in, section, symbol_to_z, to_float, z_to_symbol

_VERSION_RE = re.compile(r"xtb version\s+(\S+)\s*(?:\(([0-9a-f]+)\))?")
_BOX_KV_RE = re.compile(r"::\s*(.*?)\s+(-?\d+\.\d+(?:[eE][-+]?\d+)?)\s+(\S+)\s*::")
_SETUP_RE = re.compile(r":\s+(.*?)\s{2,}(\S+)\s+(?:\S+\s+)?:$")
_TIME_RE = re.compile(r"(\d+)\s*d,\s*(\d+)\s*h,\s*(\d+)\s*min,\s*([\d.]+)\s*sec")


def _hms(line: str) -> float | None:
    m = _TIME_RE.search(line)
    if not m:
        return None
    d, h, mi, s = m.groups()
    return int(d) * 86400 + int(h) * 3600 + int(mi) * 60 + float(s)


class XtbParser(Parser):
    program = Program.XTB

    # ------------------------------------------------------------------ header
    @section("xtb version")
    def _version(self, line: str) -> None:
        m = _VERSION_RE.search(line)
        if m:
            self.calc.program.version = m.group(1)
            self.calc.program.revision = m.group(2)
            nums = re.findall(r"\d+", m.group(1))
            self.calc.program.version_tuple = tuple(int(n) for n in nums[:3]) or None
        self.lines.next()

    @section("program call")
    def _call(self, line: str) -> None:
        self.calc.input.raw = line.split(":", 1)[1].strip().strip('"')
        call = self.calc.input.raw.lower()
        jt = self.calc.input.job_types
        if "--opt" in call or "--ohess" in call:
            jt.append(JobType.OPTIMIZATION)
        if "--hess" in call or "--ohess" in call:
            jt.append(JobType.FREQUENCY)
        if "--md" in call or "--omd" in call:
            jt.append(JobType.MD)
        if "--grad" in call:
            jt.append(JobType.GRADIENT)
        if not jt:
            jt.append(JobType.ENERGY)
        m = re.search(r"--(alpb|gbsa)\s+(\S+)", call)
        if m:
            self.calc.input.solvation = f"{m.group(1).upper()}({m.group(2)})"
        m = re.search(r"--chrg\s+(-?\d+)", call)
        if m:
            self.calc.input.charge = int(m.group(1))
        m = re.search(r"--uhf\s+(\d+)", call)
        if m:
            self.calc.input.multiplicity = int(m.group(1)) + 1
        self.lines.next()

    @section(lambda s: "G F N" in s and "x T B" in s or "G F N - F F" in s)
    def _hamiltonian(self, line: str) -> None:
        name = re.sub(r"\s+", "", line.strip(" |")).replace("GFN-FF", "GFN-FF")
        if "FF" in name:
            self.calc.input.method = "GFN-FF"
        else:
            m = re.search(r"GFN(\d)", name)
            self.calc.input.method = f"GFN{m.group(1)}-xTB" if m else "GFN-xTB"
        self.lines.next()

    @section(lambda s: s.strip() == ":                      SETUP                      :")
    def _setup(self, line: str) -> None:
        lines = self.lines
        lines.next()
        while not lines.done and lines.current.strip().startswith(":"):
            cur = lines.next().strip()
            parts = [p for p in cur.strip(":").split("  ") if p.strip()]
            if len(parts) < 2:
                continue
            key = parts[0].strip().lower()
            val = parts[1].strip().split()[0] if parts[1].strip() else ""
            if key == "# basis functions":
                self.calc.molecule.n_basis_functions = int(val)
            elif key == "# electrons":
                self.calc.molecule.n_electrons = int(val)
            elif key == "net charge":
                self.calc.molecule.charge = int(val)
                if self.calc.input.charge is None:
                    self.calc.input.charge = int(val)
            elif key == "unpaired electrons":
                mult = int(val) + 1
                self.calc.molecule.multiplicity = mult
                if self.calc.input.multiplicity is None:
                    self.calc.input.multiplicity = mult
            elif key == "hamiltonian" and not self.calc.input.method:
                self.calc.input.method = val
            elif (
                key in ("gbsa solvation", "alpb solvation") and val.lower() == "true" and not self.calc.input.solvation
            ):
                self.calc.input.solvation = key.split()[0].upper()
            elif key == "electronic temp.":
                self.calc.input.extra["electronic_temperature_K"] = float(val)

    # --------------------------------------------------------------------- SCC
    @section(lambda s: s.split()[:3] == ["iter", "E", "dE"])
    def _scc(self, line: str) -> None:
        lines = self.lines
        lines.next()
        cycle = SCFCycle()
        while not lines.done:
            cur = lines.current
            toks = cur.split()
            if len(toks) >= 4 and toks[0].isdigit():
                cycle.iterations.append(
                    SCFIteration(
                        energy=to_float(toks[1]), delta_energy=to_float(toks[2]), rms_density=to_float(toks[3])
                    )
                )
                lines.next()
            elif not cur.strip():
                lines.next()
            elif "convergence criteria satisfied" in cur:
                cycle.converged = True
                m = re.search(r"after\s+(\d+)\s+iterations", cur)
                cycle.n_iterations = int(m.group(1)) if m else len(cycle.iterations)
                lines.next()
                break
            elif "convergence criteria cannot be satisfied" in cur or "SCC not converged" in cur or "***" in cur:
                cycle.converged = False
                cycle.n_iterations = len(cycle.iterations)
                lines.next()
                break
            else:
                break
        if cycle.iterations:
            cycle.energy = cycle.iterations[-1].energy
        self.calc.scf_cycles.append(cycle)

    # ------------------------------------------------------------------- orbitals
    @section(lambda s: s.split()[:2] == ["#", "Occupation"] and "Energy/Eh" in s)
    def _orbitals(self, line: str) -> None:
        lines = self.lines
        lines.next()  # header
        lines.next()  # dashes
        energies: list[float] = []
        occs: list[float] = []
        idx: list[int] = []
        homo = None
        while not lines.done:
            cur = lines.current
            if cur.strip().startswith("---"):
                lines.next()
                break
            toks = cur.replace("(HOMO)", "").replace("(LUMO)", "").split()
            if toks and toks[0].isdigit():
                n = int(toks[0])
                if len(toks) == 4:
                    occ, e_h = to_float(toks[1]), to_float(toks[2])
                else:
                    occ, e_h = 0.0, to_float(toks[1])
                if "(HOMO)" in cur:
                    homo = len(energies)
                idx.append(n)
                energies.append(e_h)
                occs.append(occ)
            lines.next()
        orb = Orbitals(energies=[energies], occupations=[occs], indices=[idx])
        if homo is not None:
            orb.homo = [homo]
            orb.homo_energy = [energies[homo]]
            if homo + 1 < len(energies):
                orb.lumo_energy = [energies[homo + 1]]
        # Fermi level / gap follow the table.
        if lines.find("HL-Gap", limit=3):
            self.calc.energies.extra["homo_lumo_gap_eh"] = floats_in(lines.current)[0]
            lines.next()
        self.calc.orbitals = orb  # last one printed (final structure) wins

    # ------------------------------------------------------------------- charges
    @section(lambda s: s.split()[:2] == ["#", "Z"] and "covCN" in s and " q " in s)
    def _charges(self, line: str) -> None:
        lines = self.lines
        lines.next()
        charges: list[float] = []
        zs: list[int] = []
        while not lines.done and lines.current.strip():
            toks = lines.next().split()
            if not toks[0].isdigit():
                break
            zs.append(int(toks[1]))
            # columns: # Z sym covCN q C6AA alpha ; GFN1 has no covCN column
            charges.append(to_float(toks[4]) if len(toks) >= 7 else to_float(toks[3]))
        if charges:
            self.calc.populations = [p for p in self.calc.populations if p.method != "mulliken"]
            self.calc.populations.append(PopulationAnalysis(method="mulliken", charges=charges))
            if not self.calc.molecule.atomic_numbers:
                self.calc.molecule.atomic_numbers = zs
                self.calc.molecule.symbols = [z_to_symbol(z) for z in zs]

    @section("molecular dipole:")
    def _dipole(self, line: str) -> None:
        lines = self.lines
        lines.next()
        if lines.find("full:", limit=4):
            vals = floats_in(lines.next().split("full:")[1])
            if len(vals) >= 4:
                self.calc.multipoles = self.calc.multipoles or Multipoles()
                self.calc.multipoles.dipole = vals[:3]
                self.calc.multipoles.dipole_magnitude = vals[3]

    # ------------------------------------------------------------------- summary
    @section(lambda s: s.strip() == "::                     SUMMARY                     ::")
    def _summary(self, line: str) -> None:
        lines = self.lines
        lines.next()
        e = self.calc.energies
        while not lines.done and lines.current.strip().startswith("::"):
            cur = lines.next()
            m = _BOX_KV_RE.search(cur)
            if not m:
                continue
            key, val = m.group(1).strip().lstrip("-> ").lower(), float(m.group(2))
            if key == "total energy":
                e.total = val
            elif key == "scc energy":
                e.scf = val
                e.extra["scc"] = val
            elif key == "dispersion":
                e.dispersion = val
            elif key == "repulsion energy":
                e.extra["repulsion"] = val
            elif key == "gradient norm":
                e.extra["gradient_norm"] = val
            elif key == "homo-lumo gap":
                e.extra["homo_lumo_gap_ev"] = val
            elif key in ("isotropic es", "anisotropic es", "anisotropic xc", "add. restraining", "total charge"):
                e.extra[key.replace(" ", "_").replace(".", "")] = val
            elif key in ("zero point energy", "total free energy", "g(rrho) contrib.", "g(rrho) w/o zpve"):
                self._thermo_kv(key, val)
            elif "solvation" in key or "gsolv" in key or key.startswith("g("):
                e.extra[key.replace(" ", "_")] = val
            else:
                e.extra[key.replace(" ", "_")] = val

    @section(lambda s: "::" in s and "THERMODYNAMIC" in s)
    def _thermo_box(self, line: str) -> None:
        lines = self.lines
        lines.next()
        while not lines.done and lines.current.strip().startswith("::"):
            cur = lines.next()
            m = _BOX_KV_RE.search(cur)
            if not m:
                continue
            key, val = m.group(1).strip().lower(), float(m.group(2))
            if key in ("zero point energy", "total free energy", "g(rrho) contrib.", "g(rrho) w/o zpve"):
                self._thermo_kv(key, val)
            elif key == "total energy":
                self.calc.energies.total = val
            else:
                self.calc.energies.extra[key.replace(" ", "_").replace(".", "")] = val

    @section(lambda s: s.split()[:1] == ["T/K"] and "G(T)/Eh" in s)
    def _thermo_table(self, line: str) -> None:
        # "T/K  H(0)-H(T)+PV  H(T)/Eh  T*S/Eh  G(T)/Eh" followed by one row per temperature
        lines = self.lines
        lines.next()
        rows = []
        while not lines.done:
            vals = floats_in(lines.current)
            if len(vals) < 5 or lines.current.strip().startswith(("-", ":")):
                break
            rows.append(vals)
            lines.next()
        if not rows:
            return
        temp, _h0, h_t, ts, g_t = rows[-1][:5]  # last row = requested temperature (298.15 by default)
        v = self.calc.vibrations or Vibrations()
        t = v.thermochemistry or Thermochemistry()
        t.temperature = temp
        t.enthalpy = h_t
        t.gibbs_free_energy = g_t
        t.entropy = ts / temp if temp else None
        v.thermochemistry = t
        self.calc.vibrations = v
        self.calc.energies.enthalpy = h_t
        self.calc.energies.gibbs = g_t

    @section(lambda s: "# imaginary freq." in s)
    def _n_imag(self, line: str) -> None:
        vals = floats_in(line.split("freq.")[-1])
        if vals and self.calc.vibrations:
            self.calc.vibrations.n_imaginary = int(vals[0])
        self.lines.next()

    def _thermo_kv(self, key: str, val: float) -> None:
        v = self.calc.vibrations or Vibrations()
        t = v.thermochemistry or Thermochemistry()
        if key == "zero point energy":
            t.zero_point_energy = val
            self.calc.energies.zero_point = val
        elif key == "total free energy":
            t.gibbs_free_energy = val
            self.calc.energies.gibbs = val
        elif key == "g(rrho) contrib.":
            t.gibbs_correction = val
            self.calc.energies.gibbs_correction = val
        v.thermochemistry = t
        self.calc.vibrations = v

    # -------------------------------------------------------------- optimisation
    @section(lambda s: s.strip().startswith("* total energy  :") and "change" in s)
    def _opt_step(self, line: str) -> None:
        opt = self.calc.optimization or Optimization()
        vals = floats_in(line.split(":", 1)[1])
        step = OptimizationStep(
            step=len(opt.steps) + 1, energy=vals[0], delta_energy=vals[1] if len(vals) > 1 else None
        )
        nxt = self.lines.peek()
        if "gradient norm" in nxt:
            step.rms_gradient = floats_in(nxt.split(":", 1)[1])[0]
        opt.steps.append(step)
        opt.n_steps = len(opt.steps)
        self.calc.optimization = opt
        self.lines.next()

    @section("GEOMETRY OPTIMIZATION CONVERGED AFTER")
    def _opt_converged(self, line: str) -> None:
        opt = self.calc.optimization or Optimization()
        opt.converged = True
        m = re.search(r"AFTER\s+(\d+)", line)
        if m:
            opt.n_steps = int(m.group(1))
        self.calc.optimization = opt
        self.lines.next()

    @section("FAILED TO CONVERGE GEOMETRY OPTIMIZATION")
    def _opt_failed(self, line: str) -> None:
        opt = self.calc.optimization or Optimization()
        opt.converged = False
        self.calc.optimization = opt
        self.lines.next()

    @section(lambda s: s.strip() == "final structure:")
    def _final_structure(self, line: str) -> None:
        lines = self.lines
        lines.next()
        if lines.current.startswith("==="):
            lines.next()
        zs: list[int] = []
        coords: list[list[float]] = []
        first = lines.current.split()
        if first and first[0].isdigit():  # xyz: natoms, comment, rows
            n = int(lines.next().split()[0])
            lines.next()
            for _ in range(n):
                toks = lines.next().split()
                z = symbol_to_z(toks[0])
                if z is None:
                    z = int(toks[0])
                zs.append(z)
                coords.append([to_float(t) for t in toks[1:4]])
        else:  # SDF/MOL V2000: 3 header lines, counts line, then "x y z Sym ..." rows
            lines.skip(3)
            counts = lines.next()
            n = int(counts[:3])
            for _ in range(n):
                toks = lines.next().split()
                z = symbol_to_z(toks[3])
                if z is None:
                    continue
                zs.append(z)
                coords.append([to_float(t) for t in toks[:3]])
        self.calc.molecule.atomic_numbers = zs
        self.calc.molecule.symbols = [z_to_symbol(z) for z in zs]
        self.calc.geometries.append(Geometry(coordinates=coords, label="final"))

    # --------------------------------------------------------------- frequencies
    @section(lambda s: "projected vibrational frequencies" in s or "vibrational frequencies (cm" in s)
    def _frequencies(self, line: str) -> None:
        lines = self.lines
        lines.next()
        freqs: list[float] = []
        while not lines.done and lines.current.strip().startswith("eigval"):
            freqs.extend(floats_in(lines.next().split(":", 1)[1]))
        v = self.calc.vibrations or Vibrations()
        # xtb lists the 6 (or 5) zero translational/rotational modes first; drop exact zeros.
        v.frequencies = [f for f in freqs if abs(f) > 1e-6]
        v.n_imaginary = sum(1 for f in v.frequencies if f < 0)
        self.calc.vibrations = v
        if JobType.FREQUENCY not in self.calc.input.job_types:
            self.calc.input.job_types.append(JobType.FREQUENCY)

    @section(lambda s: s.strip().startswith("reduced masses (amu)"))
    def _reduced_masses(self, line: str) -> None:
        self.lines.next()
        vals: list[float] = []
        while not self.lines.done and ":" in self.lines.current and self.lines.current.strip()[:1].isdigit():
            # "   1:  1.23    2:  4.56 ..." pairs
            toks = self.lines.next().replace(":", " ").split()
            vals.extend(to_float(t) for i, t in enumerate(toks) if i % 2 == 1)
        if self.calc.vibrations and vals:
            self.calc.vibrations.reduced_masses = vals[-len(self.calc.vibrations.frequencies) :] or None

    @section(lambda s: s.strip().startswith("IR intensities (km"))
    def _ir(self, line: str) -> None:
        self.lines.next()
        vals: list[float] = []
        while not self.lines.done and ":" in self.lines.current and self.lines.current.strip()[:1].isdigit():
            toks = self.lines.next().replace(":", " ").split()
            vals.extend(to_float(t) for i, t in enumerate(toks) if i % 2 == 1)
        if self.calc.vibrations and vals:
            self.calc.vibrations.ir_intensities = vals[-len(self.calc.vibrations.frequencies) :] or None

    @section(lambda s: s.strip().startswith(":") and "temperature" in s.lower() and "::" in s)
    def _thermo_box_temp(self, line: str) -> None:
        # ":: T/K  298.15 ..." style lines inside the thermo box
        vals = floats_in(line)
        if vals and self.calc.vibrations:
            t = self.calc.vibrations.thermochemistry or Thermochemistry()
            t.temperature = vals[0]
            self.calc.vibrations.thermochemistry = t
        self.lines.next()

    @section(lambda s: "| TOTAL ENERGY" in s)
    def _final_box(self, line: str) -> None:
        e = self.calc.energies
        e.total = floats_in(line)[0]
        lines = self.lines
        lines.next()
        while not lines.done and lines.current.strip().startswith("|"):
            cur = lines.next()
            vals = floats_in(cur)
            if not vals:
                continue
            if "TOTAL ENTHALPY" in cur:
                e.enthalpy = vals[0]
            elif "TOTAL FREE ENERGY" in cur:
                e.gibbs = vals[0]
            elif "GRADIENT NORM" in cur:
                e.extra["gradient_norm"] = vals[0]
            elif "HOMO-LUMO GAP" in cur:
                e.extra["homo_lumo_gap_ev"] = vals[0]

    # ---------------------------------------------------------------- termination
    @section("* finished run on")
    def _finished(self, line: str) -> None:
        self.calc.termination = Termination.NORMAL
        self.lines.next()

    @section(
        lambda s: (
            "abnormal termination of xtb" in s or s.strip().startswith("[ERROR]") or "######## ERROR ########" in s
        )
    )
    def _error(self, line: str) -> None:
        self.calc.termination = Termination.ERROR
        msg = line.strip()
        nxt = self.lines.peek().strip()
        if nxt and nxt != msg:
            msg = f"{msg} {nxt}"
        self.calc.error_message = (self.calc.error_message + " | " if self.calc.error_message else "") + msg
        self.lines.next()

    @section(lambda s: s.strip() == "total:")
    def _timing(self, line: str) -> None:
        self.lines.next()
        t = self.calc.timing or Timing()
        for _ in range(3):
            cur = self.lines.current
            if "wall-time" in cur:
                t.wall_seconds = _hms(cur)
            elif "cpu-time" in cur:
                t.cpu_seconds = _hms(cur)
            else:
                break
            self.lines.next()
        self.calc.timing = t

    def finalize(self) -> None:
        super().finalize()
        c = self.calc
        if c.energies.total is not None and c.geometries and c.geometries[-1].energy is None:
            c.geometries[-1].energy = c.energies.total
        if c.input.method and c.input.method.startswith("GFN") and c.input.functional is None:
            c.input.functional = None
        if c.scf_cycles and c.energies.scf is None:
            c.energies.scf = c.scf_cycles[-1].energy
        if not c.molecule.atomic_numbers:
            self.warn(
                "molecule", "xtb did not print atomic symbols (single point without --verbose); molecule is empty"
            )
        if c.optimization and JobType.OPTIMIZATION not in c.input.job_types:
            c.input.job_types.append(JobType.OPTIMIZATION)
        c.input.unrestricted = (c.molecule.multiplicity or 1) > 1
        e = Energies.model_validate(c.energies.model_dump())
        c.energies = e
