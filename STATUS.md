# chemout build status

Resume: read this, `git log --oneline | head`, `.venv/bin/python -m pytest -q`, then continue with the first
unchecked item. Commit after each item.

Everything lives on the external disk: repo `/Volumes/Seagate Expansion Drive/TOPRESEARCH/chemout`
(symlinked from `TOPRESEARCH/09_PROJECTS/chemout`), corpora in `../dl/` (cclib-data clone, `orca6_github/` = 109
public ORCA 6 outputs with `manifest.tsv` giving repo + path + URL; **not redistributable**, external only),
cclib venvs in `../envs/cclib_release` (1.8.1) and `../envs/cclib_master` (2.0 dev, installed with a setuptools
override), xtb in `../envs/xtb`. Package caches: `UV_CACHE_DIR=../cache/uv`. Do not use the A6000 for this project.

## Done
- [x] Schema (pydantic v2), parser base with guarded section handlers, program sniffing
- [x] ORCA 4.1–6.1, Gaussian 09/16, xTB 6.5–6.7 parsers; 0 crashes on 63 + 103 + 108 corpus files
- [x] CLI (parse/info/batch/sniff), HDF5 export, cclib compat shim
- [x] 42 in-repo golden files + value tests (64 tests), ruff clean, CI workflow, LICENSE, CHANGELOG
- [x] Fresh xTB goldens generated locally (GFN1/GFN2/GFN-FF, opt/ohess/hess, ALPB, open shell)

## Queue
1. [x] Conformance report vs cclib 1.8.1 and cclib master (`scripts/compare_cclib.py`) → `docs/CONFORMANCE_*.md`
2. [x] Fix chemout-side disagreements (see CHANGELOG); remaining ones are conventions, listed in the reports
3. [x] README numbers via `scripts/readme_conformance.py`, commit
4. [x] Pushed to https://github.com/bahumanyarg11/chemout (public), tag v0.1.0. PyPI upload needs the user's token (`uv publish dist/*`)
5. [ ] Later: Psi4 parser (corpus exists in cclib-data), ORCA `.property.txt`/JSON reader, ORCA 6 DOCKER/GOAT/VPT2 modules, Gaussian `freq=hpmodes` high-precision block, ONIOM extrapolated energies, cclib issues to file: 1.8.1 IndexError on ORCA 6 (109/115 files), master uses last-iteration SCF energy, master xTB homos=0

## Log
- 2026-10-16 10:50 chose this project over MLIP harness (saturated: MLIPAudit, mlipx, MLIP Arena, LAMBench) and IBP reducer (months). cclib master already has ORCA 6 branches; the honest gap is: released cclib (1.8.1) vs master, crash tolerance, typed schema, HDF5/JSON, CLI.
- 11:30 first full smoke: 0 crashes across all corpora. Bugs fixed from corpus: Gaussian freq header prose, archive-line CCSD(T), SDF final structures in xtb, ORCA IR header `(cm**-1)` token, ORCA 6 units line under IR SPECTRUM, ORCA 6 point-group/irrep lines in NORMAL MODES, post-HF keywords wrongly in ORCA option list, xtb THERMODYNAMIC box.
- 11:55 final conformance: cclib 1.8.1 fails on 119/246 files (109 of 115 ORCA 6), master on 12; chemout 246/246. Remaining disagreements are conventions (docs/CONVENTIONS.md). Repo public, v0.1.0 tagged.
