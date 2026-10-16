# Conventions, and where they differ from cclib

The conformance reports (`CONFORMANCE_cclib_1.8.1.md`, `CONFORMANCE_cclib_master.md`) count
disagreements per quantity. Most of the remaining ones are not parsing errors on either side but a
different definition of the quantity. This page states chemout's definition and the cclib behaviour it
was compared against (cclib 1.8.1 release and the 2.0 development branch as of 2026-10-16).

## Geometries (`geometries`, cclib `atomcoords`)

**chemout:** one entry per distinct structure at which the program evaluated something, in file order,
each carrying its total energy. A structure re-printed without a new evaluation is dropped when its
interatomic distances match the previous entry to 1e-4 Å (Gaussian re-prints the converged geometry
after `Optimization completed` and again at the start of a following `Freq` job, re-oriented; ORCA
re-prints it under `FINAL ENERGY EVALUATION AT THE STATIONARY POINT`).
**cclib:** keeps Gaussian's re-prints, so an opt+freq file has one more geometry than SCF energies.
Gaussian dummy atoms (`X`, atomic number −1) are removed by chemout and kept by cclib; ghost atoms
(`Bq`, atomic number 0) are kept by both. ORCA dummy centres (`-`) are removed by chemout, kept by
cclib master (`1038.out`: 6 vs 9 "atoms").

## SCF energy (`energies.scf`, `scf_cycles[i].energy`; cclib `scfenergies`)

**chemout:** `scf_cycles[i].energy` is the converged SCF total as the program prints it (ORCA
`TOTAL SCF ENERGY → Total Energy`; Gaussian `SCF Done`; xTB `SCC energy`). `energies.total` is the
best total energy in the file (ORCA `FINAL SINGLE POINT ENERGY`, which includes dispersion and any
post-HF correlation; Gaussian post-HF totals when present). `energies.dispersion` is kept separately.
**cclib:** `scfenergies` is SCF-level: ORCA `Total Energy` + the printed `Dispersion correction`
(added unconditionally, although ORCA ≥ 5 already includes it in `Total Energy`); xTB per-step
`* total energy`. cclib master takes the *last SCF-iteration* energy when the `TOTAL SCF ENERGY`
block exists but is parsed from the table, which differs from the converged total by up to 1e-5 Eh —
those cases are the "differ < 1e-4 Eh" column. The `chemout.compat.cclib` shim follows cclib's
definitions (dispersion added only for ORCA < 5).

## HOMO (`orbitals.homo`, `orbitals.homo_energy`; cclib `homos`, `moenergies`)

**chemout:** `orbitals.energies` holds the orbitals actually printed. xTB prints only a window around
the gap; `orbitals.indices` then gives the true orbital numbers and `orbitals.complete` is False.
`homo` indexes `energies`. **cclib master:** reports `homos = [0]` for xTB files with a 60-element
`moenergies`, which is not the HOMO. For Gaussian `Guess=Fragment`, cclib and chemout both report the
molecule's orbitals (chemout keeps the largest eigenvalue block, the fragments' blocks are smaller).
ORCA fractional-occupation atoms: chemout takes the last orbital with occupation > 0.5, cclib master
reports orbital 0.

## Vibrational frequencies

**chemout:** translational/rotational zero modes are removed for every program (ORCA lists them as
`0.00 cm**-1`, xTB as `±0.00`). **cclib:** removes them for ORCA and Gaussian but keeps xTB's six.

## Excited states

**chemout:** `excited_states` is ordered by state index; singlet and triplet TD-DFT blocks (which
both restart at `STATE 1`) are numbered consecutively as in ORCA's spectrum table; "spin forbidden"
spectrum rows have `oscillator_strength = 0.0`. Energies come from the cm⁻¹ column (one more
significant digit than the eV column). **cclib:** `etenergies` in cm⁻¹ sorted by energy. The
conformance script compares sorted lists.

## Termination

**chemout:** `termination` is `normal` only if the program's own normal-termination line is present;
a multi-step Gaussian file is one record and the last job's status wins. **cclib:** `metadata.success`;
differs on `BH3CO_fragment_guess.log` (Gaussian 09, `guess=only`), where cclib reports failure for a
file that ends in `Normal termination`.

## Units

Fixed per quantity and independent of the program: Å, Hartree, cm⁻¹, Debye, Hartree/Bohr, eV for
excitation energies; see `Calculation.units`. cclib uses eV for all energies and Å for lengths; the shim
converts.
