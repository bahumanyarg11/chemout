# chemout

**One typed record for ORCA 6, Gaussian 16 and xTB output files.**

```bash
pip install chemout            # numpy + pydantic only
chemout info run.out           # what happened in this file?
chemout parse run.out -o run.json
chemout batch ./calcs --csv summary.csv -o all.h5   # every file under a directory, in parallel
```

```python
import chemout
calc = chemout.parse("water_opt_freq.out")        # program detected from the file
calc.program.version                               # '6.1.1'
calc.energies.total                                # Hartree, best available (post-HF if present)
calc.energies.gibbs                                # Hartree, from the thermochemistry block
calc.geometries[-1].coordinates_array()            # (n_atoms, 3) Angstrom, final geometry
calc.vibrations.frequencies                        # cm^-1, imaginary modes negative
calc.population("mulliken").charges                # e
calc.excited_states[0].energy_ev, calc.excited_states[0].oscillator_strength
calc.termination                                   # Termination.NORMAL / ERROR / INCOMPLETE
calc.provenance.warnings                           # what the parser could not read, and where
```

The same attribute names and units come back whichever program wrote the file. Everything is a
pydantic model, so `calc.to_json()`, `Calculation.from_json()`, `calc.to_hdf5("x.h5")` and IDE
autocompletion work out of the box.

## Why another parser

[cclib](https://cclib.github.io/) is the standard here and chemout borrows its scope. Three things
pushed us to write a new one:

1. **Version drift.** The latest cclib release (1.8.1, 2024) predates ORCA 6; the fixes live on an
   unreleased 2.0 branch. Groups running ORCA 6.0/6.1 today either patch cclib or write regexes.
   chemout is tested against **109 public ORCA 6.0–6.1 output files** collected from GitHub, the
   cclib test corpus (ORCA 4.1–6.1, Gaussian 03–16, xTB 6.5–6.7), and xTB runs we generate ourselves.
2. **One bad table should not lose the file.** Every section handler runs inside a guard. When a
   program changes a table, chemout records a `ParseWarning` (section, message, line number) and
   keeps going, so you still get the energy, geometry and termination status. Parsing never raises
   on file content.
3. **A typed schema, not a bag of attributes.** Fixed units, documented fields, `None` instead of
   missing attributes, JSON/HDF5 round trips, and a `batch` command that turns a directory of
   calculations into one table.

A `cclib` compatibility shim (`chemout.compat.cclib.to_ccdata_dict`) returns cclib's attribute names
and units (eV, Angstrom, cm^-1) if you have downstream code written for cclib.

## Conformance against cclib

`scripts/compare_cclib.py` runs cclib (in its own interpreter, with a timeout) and chemout on the same
files and compares final SCF energy, number of geometries, final coordinates, charge, multiplicity,
HOMO energy, Mulliken charges, vibrational frequencies, excitation energies and the success flag.
Results, with every disagreement listed by file, are in
[docs/CONFORMANCE_cclib_1.8.1.md](docs/CONFORMANCE_cclib_1.8.1.md) and
[docs/CONFORMANCE_cclib_master.md](docs/CONFORMANCE_cclib_master.md). The definitions behind
the remaining differences are in [docs/CONVENTIONS.md](docs/CONVENTIONS.md).

<!-- CONFORMANCE_SUMMARY -->
| | cclib 1.8.1 (latest release) | cclib master (2.0 dev) |
|---|---|---|
| files tried | 246 | 246 |
| **cclib raised / returned nothing** | **119** | **12** |
| chemout returned a record | 246 | 246 |

Where both produced a record, per quantity (agree / differ by < 1e-4 Eh / disagree / only one side has it):

| quantity | cclib 1.8.1 (latest release) | cclib master (2.0 dev) |
|---|---|---|
| scf_energy | 109 / 0 / 4 / 2 | 212 / 0 / 4 / 2 |
| n_geometries | 105 / 0 / 13 / 6 | 207 / 0 / 17 / 6 |
| final_coordinates | 113 / 0 / 5 / 6 | 219 / 0 / 5 / 6 |
| charge | 124 / 0 / 0 / 3 | 231 / 0 / 0 / 3 |
| multiplicity | 113 / 0 / 0 / 12 | 216 / 0 / 0 / 16 |
| homo_energy | 99 / 0 / 9 / 10 | 198 / 0 / 12 / 14 |
| mulliken_charges | 116 / 0 / 1 / 4 | 216 / 0 / 1 / 4 |
| vibfreqs | 20 / 0 / 3 / 0 | 54 / 0 / 3 / 0 |
| etenergies | 10 / 0 / 0 / 2 | 24 / 0 / 2 / 3 |
| success_flag | 125 / 0 / 2 / 0 | 228 / 0 / 0 / 6 |

Corpus: cclib-data ORCA 4.2–6.1, Gaussian 09/16, xTB; 109 public ORCA 6.0–6.1 outputs from GitHub; 9 xTB runs generated here. Every disagreement is listed by file in the linked reports, with the convention behind it where one exists (e.g. cclib master reports the last SCF-iteration energy rather than ORCA's converged `Total Energy`; cclib keeps xTB's six zero-frequency modes; chemout drops Gaussian dummy atoms but keeps ghost atoms).
<!-- /CONFORMANCE_SUMMARY -->

## What is parsed

| Group | Fields |
|---|---|
| `program` | name, version, version tuple, revision/date |
| `input` | raw route/keyword block, method, functional, dispersion, basis, auxiliary bases, charge, multiplicity, unrestricted, solvation model, nprocs, memory, job types |
| `molecule` | atomic numbers, symbols, charge, multiplicity, electrons, basis functions, point group |
| `geometries` | every geometry in order (optimisation trajectory), each with its total energy |
| `scf_cycles` | per-SCF energy, convergence, iteration table (energy, ΔE, density changes), energy components, dispersion, ⟨S²⟩ |
| `energies` | scf, total, mp2, ccsd, ccsd_t, dispersion, zero-point, thermal/enthalpy/Gibbs corrections, enthalpy, Gibbs, program-specific extras |
| `orbitals` | energies and occupations per spin, HOMO index/energy, LUMO energy, symmetries, `complete` flag |
| `populations` | mulliken, loewdin, hirshfeld, cm5, npa, apt, esp/mk/chelpg, mayer — charges and spin populations |
| `multipoles` | dipole vector and magnitude (Debye), quadrupole |
| `gradient` | final Cartesian gradient (Hartree/Bohr) |
| `optimization` | converged flag, per-step energies/gradients/steps with YES/NO flags, thresholds |
| `vibrations` | frequencies, IR intensities, Raman activities, reduced masses, force constants, symmetries, normal modes, imaginary count, thermochemistry (T, p, ZPE, H, G, S) |
| `excited_states` | energy (eV), wavelength, oscillator and rotatory strengths, multiplicity, symmetry, ⟨S²⟩, orbital transitions with coefficients |
| `nmr` | isotropic shielding and anisotropy per atom |
| `timing` | wall and CPU seconds, process count |
| `termination` | normal / error (with message) / incomplete |

Units never depend on the program: Angstrom, Hartree, cm^-1, Debye, Hartree/Bohr; excitation
energies in eV. `calc.units` carries the table.

## Programs and versions

| Program | Versions in the regression corpus | Notes |
|---|---|---|
| ORCA | 4.1, 4.2, 5.0, 6.0.0, 6.0.1, 6.1.0, 6.1.1 | SP, Opt, Freq, TD-DFT/CIS/EOM/ADC spectra, CASSCF, NMR, NEB; ORCA-DOCKER/GOAT/VPT2 module output not yet |
| Gaussian | 09 (Rev D/E), 16 (Rev A/B/C) | multi-step `Link1` files are one record (last termination wins) |
| xTB | 6.5.1, 6.7.1 | GFN1, GFN2, GFN-FF, `--opt/--ohess/--hess`, ALPB/GBSA, open shell; single points print no coordinates (xTB limitation, use `--verbose`) |

## Install from source

```bash
git clone https://github.com/bahumanyarg11/chemout && cd chemout
python -m venv .venv && . .venv/bin/activate
pip install -e ".[dev]"
pytest -q                               # 64 tests, 42 golden files
python scripts/make_goldens.py --check  # byte-for-byte regression against tests/golden
```

Goldens are regenerated deliberately with `--update`; the diff is reviewed because it is a change in
what users receive.

## Limitations

- Only ORCA, Gaussian and xTB. The handler model makes adding a program a few hundred lines; Psi4 is next.
- Gaussian `freq=hpmodes` high-precision displacement blocks are skipped in favour of the standard block.
- Gaussian route lines are parsed heuristically; `method/basis` is reliable, exotic `method basis` orderings fall back to the `Standard basis:` line.
- ORCA `.property.txt` / `.property.json` are not read; the text output is.

## License

BSD 3-Clause. Test files under `tests/data/cclib-data` come from the cclib-data repository (BSD 3-Clause).
