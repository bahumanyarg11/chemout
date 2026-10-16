"""Hand-checked values from the in-repo corpus (numbers read directly from the output files)."""

from pathlib import Path

import pytest

import chemout
from chemout.schema import JobType, Termination

DATA = Path(__file__).parent / "data"


def _p(rel: str):
    return chemout.parse(DATA / rel)


class TestOrca6:
    def test_he_single_point(self):
        c = _p("cclib-data/ORCA6.0/commercial_version.out")
        assert c.program.version == "6.0.1-f.3"
        assert c.program.version_tuple == (6, 0, 1)
        assert c.termination is Termination.NORMAL
        assert c.input.method == "HF"
        assert c.input.basis == "STO-3G"
        assert c.input.charge == 0 and c.input.multiplicity == 1
        assert c.molecule.symbols == ["He"]
        assert c.energies.scf == pytest.approx(-2.80778395753546)
        assert c.energies.total == pytest.approx(-2.807783957535)
        assert c.scf_cycles[-1].converged is True and c.scf_cycles[-1].n_iterations == 2
        assert c.orbitals.energies[0][0] == pytest.approx(-0.876036)
        assert c.orbitals.occupations[0][0] == 2.0
        assert c.population("mulliken").charges == [0.0]
        assert c.population("loewdin").charges == [0.0]
        assert c.multipoles.dipole_magnitude == 0.0
        assert c.timing.wall_seconds == pytest.approx(0.223)

    def test_uhf_optimisation(self):
        c = _p("cclib-data/ORCA6.0/1618.out")
        assert c.program.version == "6.0.0"
        assert c.input.unrestricted is True
        assert c.input.multiplicity == 2
        assert c.input.nprocs == 36
        assert JobType.OPTIMIZATION in c.input.job_types
        assert c.molecule.symbols == ["H", "F", "F", "C"]
        assert len(c.geometries) == 5
        assert c.geometries[0].energy == pytest.approx(-237.198245745248)
        assert c.geometries[-1].label == "final"
        assert c.optimization.converged is True
        assert len(c.optimization.steps) == 4
        assert c.optimization.steps[0].rms_gradient == pytest.approx(0.0132738592)
        assert c.optimization.criteria["rms_gradient"] == pytest.approx(3e-5)
        assert c.scf_cycles[0].n_iterations == 14
        assert c.scf_cycles[0].iterations[0].energy == pytest.approx(-237.1543926738920618)
        assert c.scf_cycles[0].spin_squared == pytest.approx(0.755345)
        assert c.scf_cycles[0].spin_squared_ideal == pytest.approx(0.75)
        mull = c.population("mulliken")
        assert mull.charges[3] == pytest.approx(0.603628)
        assert mull.spins[3] == pytest.approx(0.964686)
        assert c.gradient[0] == pytest.approx([1.1827e-05, -5.801e-05, 4.7508e-05])
        assert c.multipoles.dipole_magnitude == pytest.approx(1.785720024)
        assert len(c.orbitals.energies) == 2  # spin up / spin down
        assert c.orbitals.energies[0][0] == pytest.approx(-26.316384)
        assert c.timing.wall_seconds == pytest.approx(9 * 60 + 49.114)

    def test_zindo_tddft_6_1(self):
        c = _p("cclib-data/ORCA6.1/1509.out")
        assert c.program.version == "6.1.1"
        assert c.input.method.lower() == "zindo/s"
        assert JobType.EXCITED_STATES in c.input.job_types
        assert len(c.excited_states) == 50
        s1 = c.excited_states[0]
        assert s1.energy_ev == pytest.approx(3.412305, abs=1e-3)
        assert s1.wavelength_nm == pytest.approx(363.3, abs=0.1)
        assert s1.oscillator_strength == 0.0
        assert s1.multiplicity == "singlet"
        assert s1.transitions[0].from_orbital == 1783 and s1.transitions[0].to_orbital == 1787
        assert s1.transitions[0].coefficient == pytest.approx(0.996871)
        assert c.energies.scf == pytest.approx(-13413.675475327)
        assert c.energies.total == pytest.approx(-13413.550075428942)
        assert c.molecule.n_basis_functions == 2761


class TestOrcaOlder:
    def test_mp2_freq_thermo_4_2(self):
        c = _p("cclib-data/ORCA4.2/MP2-gradient.out")
        v = c.vibrations
        assert v.frequencies == pytest.approx([1626.47, 3592.60, 3899.99])
        assert v.ir_intensities == pytest.approx([82.890396, 22.694685, 48.622069])
        assert v.n_imaginary == 0
        assert len(v.normal_modes) == 3 and len(v.normal_modes[0]) == 3
        t = v.thermochemistry
        assert t.temperature == 298.15 and t.pressure == 1.0
        assert t.zero_point_energy == pytest.approx(0.02077475)
        assert t.enthalpy == pytest.approx(-76.28802416)
        assert c.energies.zero_point == pytest.approx(0.02077475)

    def test_ccsd_t(self):
        c = _p("cclib-data/ORCA4.2/water_ccsd_t.out")
        assert c.energies.ccsd_t is not None or c.energies.total is not None
        assert c.energies.total == pytest.approx(-75.013556, abs=1e-5)


class TestGaussian16:
    def test_opt_freq_water(self):
        c = _p("cclib-data/Gaussian16/water_neutral_nbo_opt_freq.out")
        assert c.program.version == "16 RevA.03"
        assert c.termination is Termination.NORMAL
        assert c.input.raw.startswith("# hf/sto-3g opt freq pop=nbo")
        assert c.input.method == "HF" and c.input.basis == "sto-3g"
        assert JobType.OPTIMIZATION in c.input.job_types and JobType.FREQUENCY in c.input.job_types
        assert c.input.charge == 0 and c.input.multiplicity == 1
        assert c.molecule.atomic_numbers == [8, 1, 1]
        assert c.geometries[0].coordinates[0] == pytest.approx([0.0, 0.0, 0.119159])
        assert c.scf_cycles[0].energy == pytest.approx(-74.9643287914)
        assert c.scf_cycles[0].n_iterations == 7
        assert c.optimization.converged is True
        assert c.optimization.steps[0].max_gradient == pytest.approx(0.029587)
        assert c.optimization.criteria["maximum_force"] == pytest.approx(0.00045)
        assert c.orbitals.energies[0][:2] == pytest.approx([-20.25155, -1.25760])
        assert c.orbitals.homo == [4]
        assert c.population("mulliken").charges[0] == pytest.approx(-0.339215, abs=0.01)
        assert c.population("npa") is not None
        assert c.multipoles.dipole_magnitude == pytest.approx(1.6735, abs=0.1)
        v = c.vibrations
        assert v.frequencies == pytest.approx([2169.7613, 4141.3837, 4392.5759])
        assert v.ir_intensities == pytest.approx([7.2483, 44.2724, 29.9428])
        assert v.reduced_masses == pytest.approx([1.0785, 1.0491, 1.0774])
        assert v.symmetries == ["A1", "A1", "B2"]
        assert len(v.normal_modes) == 3 and len(v.normal_modes[0]) == 3
        t = v.thermochemistry
        assert t.temperature == 298.15
        assert t.zero_point_energy == pytest.approx(0.024385)
        assert t.gibbs_free_energy == pytest.approx(-74.959254)
        assert c.energies.gibbs == pytest.approx(-74.959254)
        assert t.entropy == pytest.approx(45.281 / 1000 / 627.5094740631)
        assert c.timing.cpu_seconds == pytest.approx(3.7 + 1.7)

    def test_eom_ccsd_excited_states(self):
        c = _p("cclib-data/Gaussian16/Pyridine_opt-eomccsd.log")
        assert len(c.excited_states) == 3
        assert all(s.energy_ev > 0 for s in c.excited_states)
        assert c.excited_states[0].energy_ev == pytest.approx(4.3331)
        assert c.excited_states[0].oscillator_strength == pytest.approx(0.0051)
        assert c.excited_states[0].multiplicity == "singlet"
        assert c.input.method == "EOMCCSD" and c.input.basis == "STO-3G"

    def test_hirshfeld_cm5(self):
        c = _p("cclib-data/Gaussian16/dvb_hirshfeld.out")
        assert c.population("hirshfeld") is not None
        assert c.population("cm5") is not None
        assert len(c.population("cm5").charges) == c.molecule.n_atoms

    def test_ccsd_t_archive_line_does_not_break(self):
        c = _p("cclib-data/Gaussian16/1615.log")
        assert c.termination is Termination.NORMAL
        assert c.energies.ccsd_t == pytest.approx(-94.7445478, abs=1e-6)


class TestXtb:
    def test_opt_6_7(self):
        c = _p("cclib-data/xTB/1768_xyz.out")
        assert c.program.version == "6.7.1pre"
        assert c.input.method == "GFN2-xTB"
        assert c.termination is Termination.NORMAL
        assert c.molecule.n_atoms == 21 and c.molecule.symbols[0] == "C"
        assert c.energies.total == pytest.approx(-39.629079572158)
        assert c.energies.extra["homo_lumo_gap_ev"] == pytest.approx(3.604684056817)
        assert c.optimization.converged is True and c.optimization.n_steps == 37
        assert c.orbitals.homo_energy[0] == pytest.approx(-0.4160635)
        assert c.orbitals.lumo_energy[0] == pytest.approx(-0.2835938)
        assert c.population("mulliken").charges[0] == pytest.approx(-0.149)
        assert c.multipoles.dipole_magnitude == pytest.approx(2.769)
        assert c.multipoles.dipole == pytest.approx([0.817, 0.186, -0.696])
        assert c.timing.wall_seconds == pytest.approx(0.181)
        assert c.timing.cpu_seconds == pytest.approx(5.344)

    def test_sdf_final_structure(self):
        c = _p("cclib-data/xTB/1768_sdf.out")
        assert c.molecule.n_atoms == 21
        assert c.geometries[-1].coordinates[0] == pytest.approx([-2.7248, -2.2295, -0.0449])

    def test_generated_ohess(self):
        c = _p("xtb_generated/water_ohess/water_ohess.out")
        assert c.termination is Termination.NORMAL
        assert JobType.OPTIMIZATION in c.input.job_types and JobType.FREQUENCY in c.input.job_types
        assert c.vibrations is not None and len(c.vibrations.frequencies) == 3
        assert all(f > 0 for f in c.vibrations.frequencies)
        assert c.energies.gibbs is not None
        assert c.energies.zero_point is not None

    def test_generated_charged_open_shell(self):
        c = _p("xtb_generated/methanol_cation_uhf/methanol_cation_uhf.out")
        assert c.input.charge == 1 and c.input.multiplicity == 2
        assert c.molecule.charge == 1
        c = _p("xtb_generated/acetate_alpb/acetate_alpb.out")
        assert c.input.charge == -1
        assert c.input.solvation and c.input.solvation.startswith("ALPB")
