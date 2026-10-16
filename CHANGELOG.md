# Changelog

All notable changes to `chemout` are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- Typed, program-independent result schema (`chemout.schema.Calculation`, pydantic v2) with
  fixed units: Angstrom, Hartree, cm^-1, Debye, eV for excitation energies.
- Parsers for ORCA 4.1–6.1, Gaussian 09/16 and xTB 6.5–6.7 (GFN1, GFN2, GFN-FF): geometry
  trajectories, SCF cycles and iterations, post-HF energies, orbital energies and occupations,
  Mulliken/Löwdin/Hirshfeld/CM5/NPA/APT/ESP/Mayer/CHELPG charges and spin populations, dipoles,
  gradients, optimisation convergence tables, frequencies with IR intensities and normal modes,
  thermochemistry, TD-DFT/CIS/EOM excited states with transitions, NMR shieldings, timings,
  termination status and error messages.
- Every section handler runs in a guard: a format change in one table is recorded as a
  `ParseWarning` and the rest of the file is still parsed. Parsing never raises on content.
- `chemout` CLI: `parse` (JSON / JSONL / HDF5), `info`, `batch` (parallel, CSV summary), `sniff`.
- HDF5 export with the full JSON record embedded.
- `chemout.compat.cclib.to_ccdata_dict` / `to_ccdata`: cclib attribute names and units.
- Golden-file regression suite over 42 in-repo outputs (`tests/golden`) plus a conformance
  script against cclib (`scripts/compare_cclib.py`) and the resulting reports in `docs/`.
- Corpus-driven fixes found while building the conformance report: ORCA 6 `NORMAL MODES` point-group
  and irrep lines, ORCA 5/6 units line under `IR SPECTRUM`, ORCA post-HF keywords (`MP2`, `CCSD(T)`,
  `ROHF`) as methods, singlet/triplet TD-DFT blocks that both restart at `STATE 1`, "spin forbidden"
  spectrum rows, SCF cycles that omit the `TOTAL SCF ENERGY` block or print one block per iteration,
  Gaussian `freq=hpmodes`, IRC `IRC Coupling -` rows, `Guess=Fragment`/ONIOM extra charge tables,
  ghost atoms (Bq) kept and dummy atoms (X) dropped, re-oriented duplicate geometries from a
  following frequency job, xTB SDF-format final structures and the `THERMODYNAMIC` box.
