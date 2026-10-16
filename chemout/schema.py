"""Typed result schema shared by every parser.

Units are fixed and never depend on the program that produced the file:

* lengths in Angstrom, gradients in Hartree/Bohr
* energies in Hartree, except ``ExcitedState.energy_ev`` (eV)
* orbital energies in Hartree
* vibrational frequencies in cm^-1, IR intensities in km/mol
* dipole moments in Debye, temperatures in K, pressures in atm
* times in seconds

Every array field is a plain nested list in JSON and a numpy array via the
``*_array`` helpers, so the JSON round-trips without numpy.
"""

from __future__ import annotations

import datetime as _dt
from enum import Enum
from typing import Any

import numpy as np
from pydantic import BaseModel, ConfigDict, Field

UNITS: dict[str, str] = {
    "coordinates": "angstrom",
    "energy": "hartree",
    "orbital_energy": "hartree",
    "gradient": "hartree/bohr",
    "frequency": "cm^-1",
    "ir_intensity": "km/mol",
    "dipole": "debye",
    "excitation_energy": "eV",
    "wavelength": "nm",
    "temperature": "K",
    "pressure": "atm",
    "entropy": "hartree/K",
    "time": "s",
    "nmr_shielding": "ppm",
    "spin_squared": "hbar^2",
}

HARTREE_TO_EV = 27.211386245988
BOHR_TO_ANGSTROM = 0.529177210903
HARTREE_TO_KCAL = 627.5094740631
HARTREE_TO_KJ = 2625.4996394799
AU_TO_DEBYE = 2.541746473


class Program(str, Enum):
    ORCA = "orca"
    GAUSSIAN = "gaussian"
    XTB = "xtb"
    PSI4 = "psi4"
    UNKNOWN = "unknown"


class Termination(str, Enum):
    NORMAL = "normal"
    ERROR = "error"
    INCOMPLETE = "incomplete"  # file ends mid-run (job killed, still running, truncated)
    UNKNOWN = "unknown"


class _Model(BaseModel):
    model_config = ConfigDict(extra="forbid", validate_assignment=True)


class ParseWarning(_Model):
    """A sub-extractor failed or saw something it did not understand.

    Parsing never raises for these; it records them so the rest of the file still yields
    data and the user can see exactly which quantity is missing and why.
    """

    section: str
    message: str
    line: int | None = None


class ProgramInfo(_Model):
    name: Program = Program.UNKNOWN
    version: str | None = Field(None, description="Version string as printed, e.g. '6.0.1'")
    version_tuple: tuple[int, ...] | None = None
    revision: str | None = Field(None, description="Build / git revision if printed")
    date: str | None = Field(None, description="Build or run date string as printed")


class Provenance(_Model):
    path: str | None = None
    sha256: str | None = None
    size_bytes: int | None = None
    parsed_at: str = Field(default_factory=lambda: _dt.datetime.now(_dt.UTC).isoformat())
    chemout_version: str = "0.1.0"
    warnings: list[ParseWarning] = Field(default_factory=list)

    def warn(self, section: str, message: str, line: int | None = None) -> None:
        self.warnings.append(ParseWarning(section=section, message=message, line=line))


class JobType(str, Enum):
    ENERGY = "energy"
    GRADIENT = "gradient"
    OPTIMIZATION = "optimization"
    TS_OPTIMIZATION = "ts_optimization"
    FREQUENCY = "frequency"
    EXCITED_STATES = "excited_states"
    NMR = "nmr"
    SCAN = "scan"
    MD = "md"
    IRC = "irc"
    OTHER = "other"


class InputInfo(_Model):
    raw: str | None = Field(None, description="Route line / keyword block exactly as echoed")
    job_types: list[JobType] = Field(default_factory=list)
    method: str | None = Field(None, description="Normalised method name, e.g. 'B3LYP', 'DLPNO-CCSD(T)', 'GFN2-xTB'")
    functional: str | None = None
    dispersion: str | None = Field(None, description="e.g. 'D3BJ', 'D4'")
    basis: str | None = None
    aux_basis: list[str] = Field(default_factory=list)
    charge: int | None = None
    multiplicity: int | None = None
    unrestricted: bool | None = None
    solvation: str | None = Field(None, description="e.g. 'CPCM(water)', 'SMD(toluene)', 'ALPB(water)'")
    nprocs: int | None = None
    memory_mb: float | None = None
    title: str | None = None
    extra: dict[str, Any] = Field(default_factory=dict, description="Program-specific input details")


class Molecule(_Model):
    atomic_numbers: list[int] = Field(default_factory=list)
    symbols: list[str] = Field(default_factory=list)
    masses: list[float] | None = Field(None, description="amu, if printed")
    charge: int | None = None
    multiplicity: int | None = None
    n_electrons: int | None = None
    n_basis_functions: int | None = None
    point_group: str | None = None

    @property
    def n_atoms(self) -> int:
        return len(self.atomic_numbers)


class Geometry(_Model):
    coordinates: list[list[float]] = Field(description="n_atoms x 3, Angstrom")
    energy: float | None = Field(None, description="Total energy at this geometry, Hartree")
    step: int | None = None
    label: str | None = Field(None, description="e.g. 'input', 'opt step 3', 'final'")

    def coordinates_array(self) -> np.ndarray:
        return np.asarray(self.coordinates, dtype=float)


class SCFIteration(_Model):
    energy: float
    delta_energy: float | None = None
    rms_density: float | None = None
    max_density: float | None = None
    damping: float | None = None
    diis_error: float | None = None


class SCFCycle(_Model):
    """One SCF solve (one per geometry step, typically)."""

    energy: float | None = Field(None, description="Converged SCF energy, Hartree")
    converged: bool | None = None
    n_iterations: int | None = None
    iterations: list[SCFIteration] = Field(default_factory=list)
    nuclear_repulsion: float | None = None
    one_electron: float | None = None
    two_electron: float | None = None
    exchange_correlation: float | None = None
    dispersion: float | None = None
    spin_squared: float | None = Field(None, description="<S^2> as printed (unrestricted only)")
    spin_squared_ideal: float | None = None


class Energies(_Model):
    """All named total energies found in the file, Hartree. Each is the *final* value."""

    scf: float | None = None
    total: float | None = Field(None, description="Best available total energy (post-HF if present)")
    mp2: float | None = None
    mp2_correlation: float | None = None
    ccsd: float | None = None
    ccsd_t: float | None = None
    dispersion: float | None = None
    solvation: float | None = Field(None, description="Solvation free energy correction")
    gcp: float | None = Field(None, description="Geometrical counterpoise correction")
    zero_point: float | None = None
    thermal_energy_correction: float | None = None
    enthalpy_correction: float | None = None
    gibbs_correction: float | None = None
    enthalpy: float | None = Field(None, description="E + H correction (absolute)")
    gibbs: float | None = Field(None, description="E + G correction (absolute)")
    extra: dict[str, float] = Field(default_factory=dict, description="Program-specific energies")


class Orbitals(_Model):
    energies: list[list[float]] = Field(default_factory=list, description="[spin][mo] Hartree")
    occupations: list[list[float]] = Field(default_factory=list, description="[spin][mo]")
    symmetries: list[list[str]] = Field(default_factory=list)
    homo: list[int] = Field(default_factory=list, description="Index of HOMO per spin (0-based, into `energies`)")
    indices: list[list[int]] | None = Field(
        None,
        description="1-based orbital numbers as printed, only when the program printed a truncated window "
        "(xTB) so `energies` is not the complete set",
    )
    n_basis: int | None = None
    complete: bool = Field(True, description="False when the program printed only part of the virtual space")
    homo_energy: list[float] | None = Field(None, description="HOMO energy per spin, Hartree")
    lumo_energy: list[float] | None = Field(None, description="LUMO energy per spin, Hartree")

    def homo_lumo_gap(self, spin: int = 0) -> float | None:
        if not self.energies or spin >= len(self.energies) or spin >= len(self.homo):
            return None
        e, h = self.energies[spin], self.homo[spin]
        if h + 1 >= len(e):
            return None
        return e[h + 1] - e[h]


class PopulationAnalysis(_Model):
    method: str = Field(description="mulliken | loewdin | hirshfeld | cm5 | npa | chelpg | mayer | esp | ...")
    charges: list[float] = Field(default_factory=list)
    spins: list[float] | None = Field(None, description="Spin populations, if unrestricted")


class Multipoles(_Model):
    dipole: list[float] | None = Field(None, description="x,y,z in Debye")
    dipole_magnitude: float | None = Field(None, description="Debye")
    quadrupole: list[float] | None = Field(None, description="xx,yy,zz,xy,xz,yz in Debye*Angstrom")
    polarizability: list[list[float]] | None = Field(None, description="3x3 in bohr^3")


class Thermochemistry(_Model):
    temperature: float | None = None
    pressure: float | None = None
    zero_point_energy: float | None = None
    thermal_energy_correction: float | None = Field(None, description="Total thermal correction to E")
    enthalpy_correction: float | None = None
    gibbs_correction: float | None = None
    entropy: float | None = Field(None, description="Total entropy, Hartree/K")
    enthalpy: float | None = None
    gibbs_free_energy: float | None = None
    rotational_symmetry_number: int | None = None
    rotational_constants_ghz: list[float] | None = None


class Vibrations(_Model):
    frequencies: list[float] = Field(default_factory=list, description="cm^-1; negative = imaginary")
    ir_intensities: list[float] | None = Field(None, description="km/mol")
    raman_activities: list[float] | None = None
    reduced_masses: list[float] | None = None
    force_constants: list[float] | None = None
    symmetries: list[str] | None = None
    normal_modes: list[list[list[float]]] | None = Field(None, description="[mode][atom][xyz]")
    n_imaginary: int | None = None
    thermochemistry: Thermochemistry | None = None

    def normal_modes_array(self) -> np.ndarray | None:
        return None if self.normal_modes is None else np.asarray(self.normal_modes, dtype=float)


class Transition(_Model):
    from_orbital: int
    to_orbital: int
    coefficient: float
    from_spin: str | None = Field(None, description="'a' or 'b'")
    to_spin: str | None = None


class ExcitedState(_Model):
    index: int
    energy_ev: float
    wavelength_nm: float | None = None
    oscillator_strength: float | None = None
    rotatory_strength: float | None = None
    multiplicity: str | None = Field(None, description="'singlet', 'triplet', or None if unknown")
    symmetry: str | None = None
    spin_squared: float | None = None
    transitions: list[Transition] = Field(default_factory=list)
    transition_dipole: list[float] | None = None


class OptimizationStep(_Model):
    step: int
    energy: float | None = None
    delta_energy: float | None = None
    rms_gradient: float | None = None
    max_gradient: float | None = None
    rms_step: float | None = None
    max_step: float | None = None
    converged_flags: dict[str, bool] = Field(default_factory=dict)


class Optimization(_Model):
    converged: bool | None = None
    n_steps: int | None = None
    steps: list[OptimizationStep] = Field(default_factory=list)
    criteria: dict[str, float] = Field(default_factory=dict, description="Thresholds as printed")


class NMRShielding(_Model):
    atom_index: int
    isotropic: float
    anisotropy: float | None = None
    tensor: list[list[float]] | None = None


class Timing(_Model):
    wall_seconds: float | None = None
    cpu_seconds: float | None = None
    n_procs: int | None = None


class Calculation(_Model):
    """Everything parsed from one output file."""

    schema_version: str = "1"
    units: dict[str, str] = Field(default_factory=lambda: dict(UNITS))
    program: ProgramInfo = Field(default_factory=ProgramInfo)
    provenance: Provenance = Field(default_factory=Provenance)
    termination: Termination = Termination.UNKNOWN
    error_message: str | None = None
    input: InputInfo = Field(default_factory=InputInfo)
    molecule: Molecule = Field(default_factory=Molecule)
    geometries: list[Geometry] = Field(default_factory=list, description="In order of appearance")
    scf_cycles: list[SCFCycle] = Field(default_factory=list)
    energies: Energies = Field(default_factory=Energies)
    orbitals: Orbitals | None = None
    populations: list[PopulationAnalysis] = Field(default_factory=list)
    multipoles: Multipoles | None = None
    gradient: list[list[float]] | None = Field(None, description="Final gradient, n_atoms x 3, Hartree/Bohr")
    optimization: Optimization | None = None
    vibrations: Vibrations | None = None
    excited_states: list[ExcitedState] = Field(default_factory=list)
    nmr: list[NMRShielding] = Field(default_factory=list)
    timing: Timing | None = None

    # -- convenience -----------------------------------------------------------------
    @property
    def final_geometry(self) -> Geometry | None:
        return self.geometries[-1] if self.geometries else None

    @property
    def final_energy(self) -> float | None:
        return self.energies.total if self.energies.total is not None else self.energies.scf

    @property
    def success(self) -> bool:
        return self.termination is Termination.NORMAL

    def population(self, method: str) -> PopulationAnalysis | None:
        for p in self.populations:
            if p.method == method.lower():
                return p
        return None

    def to_dict(self) -> dict[str, Any]:
        return self.model_dump(mode="json")

    def to_json(self, indent: int | None = 2) -> str:
        return self.model_dump_json(indent=indent)

    @classmethod
    def from_json(cls, text: str) -> Calculation:
        return cls.model_validate_json(text)

    def to_hdf5(self, path: str, group: str | None = None) -> None:
        from .io.hdf5 import write_hdf5

        write_hdf5(self, path, group)
