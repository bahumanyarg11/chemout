# chemout vs cclib master (2.0 development branch, cloned 2026-10-16)

Files: 247 (1 skipped: not ORCA/Gaussian/xTB). cclib failed on **12** files (exception, timeout > 300s, or `ccread` returned None); chemout returned a record for every file.

| quantity | agree | differ < 1e-4 Eh | disagree | chemout only | cclib only | neither |
|---|---|---|---|---|---|---|
| scf_energy | 181 | 1 | 34 | 2 | 0 | 16 |
| n_geometries | 207 | 0 | 17 | 2 | 4 | 4 |
| final_coordinates | 219 | 0 | 5 | 2 | 4 | 4 |
| charge | 231 | 0 | 0 | 1 | 2 | 0 |
| multiplicity | 216 | 0 | 0 | 16 | 0 | 2 |
| homo_energy | 198 | 0 | 12 | 14 | 0 | 10 |
| mulliken_charges | 216 | 0 | 1 | 3 | 1 | 13 |
| vibfreqs | 54 | 0 | 3 | 0 | 0 | 177 |
| etenergies | 24 | 0 | 2 | 1 | 2 | 205 |
| success_flag | 228 | 0 | 0 | 6 | 0 | 0 |

Counts are over the files where cclib produced a record. 'differ < 1e-4 Eh' marks energy differences
below 1e-4 Hartree, i.e. one implementation reporting the last SCF-iteration energy and the other the
converged total: not a parsing error on either side.

## cclib failures

- `JEFF7712__cctop__orca.out` (orca 6.0.0): ?
- `ReactionMechanismGenerator__ARC__neb_res.out` (orca 6.0.0): ?
- `avcopan__project-cyclopentene-oh-addition__sp_orca.out` (orca 6.1.1): ?
- `avcopan__project-dimethyl-oxirane-cl-abstraction__opt.log` (orca 6.1.0): ?
- `chpc-uofu__chpc-myjobs-templates__water.out` (orca 6.1.1): ?
- `grimme-lab__q-vSZPs__orca.out` (orca 6.1.1): ?
- `ushcodex__Vanadium-zeolite-y-dealumination__V-PRC-T1-V03.out` (orca 6.1.0): ?
- `ushcodex__Vanadium-zeolite-y-dealumination__V-TS2-T1-V03.out` (orca 6.1.0): ?
- `ushcodex__Vanadium-zeolite-y-dealumination__irc_backward_opt.out` (orca 6.1.0): ?
- `ushcodex__Vanadium-zeolite-y-dealumination__irc_forward_opt.out` (orca 6.1.0): ?
- `g98.out` (gaussian None): ?
- `g98.out` (gaussian None): ?

## Disagreements

- `ho-h_cas_aDZ_scanSA.out` (orca 6.0.0): etenergies
- `1038.out` (orca 5.0.4): final_coordinates
- `947.out` (orca 4.2.1): scf_energy
- `cope.out` (orca 4.2.1): scf_energy
- `dvb_dispersion_bp86_d3zero.out` (orca 4.2.0): scf_energy
- `Benzene_anion.log` (gaussian 16 RevC.01): n_geometries
- `Pyridine_opt-eomccsd.log` (gaussian 16 RevC.01): etenergies
- `co_pbe1pbe_631ppGss.log` (gaussian 16 RevApple): n_geometries
- `dol-1-pen-5-pen-trip.out` (gaussian 16 RevC.01): n_geometries, final_coordinates
- `issue1153.log` (gaussian 16 RevC.01): n_geometries
- `issue962.log` (gaussian 16 RevA.03): scf_energy, n_geometries, final_coordinates, mulliken_charges
- `water.log` (gaussian 16 RevC.01): n_geometries
- `water_neutral_nbo_opt_freq.out` (gaussian 16 RevA.03): n_geometries
- `1080.out` (gaussian 09 RevC.01): n_geometries
- `1504.log` (gaussian 09 RevA.02): n_geometries
- `25DMF_HRANH.log` (gaussian 09 RevB.01): n_geometries
- `BSL_opt_freq_DFT.out` (gaussian 09 RevD.01): n_geometries
- `OPT_oniom.log` (gaussian 09 RevD.01): scf_energy
- `ansaBPFP2_acetylene_ircf_LQA.out` (gaussian 09 RevE.01): n_geometries, final_coordinates
- `encaps1.log` (gaussian 09 RevC.01): scf_energy
- `issue_460.log` (gaussian 09 RevD.01): n_geometries, final_coordinates
- `oniom_IR_intensity.log` (gaussian 09 RevC.01): scf_energy, vibfreqs
- `test0200.log` (gaussian 09 RevD.01): n_geometries
- `1448.out` (xtb 6.5.1): homo_energy
- `1768_sdf.out` (xtb 6.7.1pre): homo_energy
- `1768_xyz.out` (xtb 6.7.1pre): homo_energy
- `Svatunek-Lab__teaching__H2O.out` (orca 6.0.0): scf_energy
- `Xingyao-Wu__Chebyshev-Tree-Tensor-Network-for-Vibrational-Spectra__opt.out` (orca 6.1.0): n_geometries
- `aligfellow__graphRC__dihedral.out` (orca 6.0.1): scf_energy
- `avcopan__project-cyclopentene-oh-addition__opt_orca.out` (orca 6.1.1): scf_energy
- `chpc-uofu__chpc-myjobs-templates__cadi.out` (orca 6.1.1): n_geometries
- `chpc-uofu__cluster-validation__cadi.out` (orca 6.0.1): n_geometries
- `dhwan0112__dhk.github.io__water.out` (orca 6.1.1): scf_energy
- `geoffreyweal__ORCA_Mechanism_Procedure__output.out` (orca 6.1.0-f.0): scf_energy
- `geoffreyweal__viewORCA__output.out` (orca 6.1.0-f.0): scf_energy
- `jhoanlond__QCXMS2__orca.out` (orca 6.0.1): scf_energy
- `jhoanlond__QCXMS2__orca_atom6.out` (orca 6.0.1): homo_energy
- `ll23002__Orca__Etanol.out` (orca 6.1.1): scf_energy
- `ll23002__Orca__agua.out` (orca 6.1.1): scf_energy
- `ll23002__Orca__water.out` (orca 6.1.0-f.0): scf_energy
- `lucasandradechem-lab__chg2itp__dft.out` (orca 6.1.1): scf_energy
- `lucasgarayy__Machine-Learned-Interatomic-Potentials-for-Photophysics__orca.out` (orca 6.0.0): scf_energy
- `madsr2d2__masterThesis__benzaldehyde_tddft.out` (orca 6.0.1): scf_energy
- `madsr2d2__masterThesis__benzoic_acid_tddft.out` (orca 6.0.1): scf_energy
- `mersalas__qm_cmd__h2_frag.out` (orca 6.0.1): scf_energy
- `mersalas__qm_cmd__h3_frag.out` (orca 6.0.1): scf_energy
- `mersalas__qm_cmd__hh-pg_cplx.out` (orca 6.0.1): scf_energy
- `mersalas__qm_cmd__k2_frag.out` (orca 6.0.1): scf_energy
- `mersalas__qm_cmd__kh-pg_cplx.out` (orca 6.0.1): scf_energy
- `mersalas__qm_cmd__pg1_frag.out` (orca 6.0.1): scf_energy
- `mersalas__qm_cmd__r2_frag.out` (orca 6.0.1): scf_energy
- `mersalas__qm_cmd__rh-pg_cplx.out` (orca 6.0.1): scf_energy
- `patonlab__pyDFTD3__ibuprofen_b3lyp_orca.out` (orca 6.1.0): scf_energy
- `patonlab__pyDFTD3__ibuprofen_pbe_orca.out` (orca 6.1.0): scf_energy
- `rgora__photochemistry__uracil_s0.out` (orca 6.0.1): scf_energy
- `rgora__photochemistry__uracil_s0_vee.out` (orca 6.0.1): scf_energy
- `rgora__photochemistry__uracil_s0_vg.out` (orca 6.0.1): scf_energy
- `acetate_alpb.out` (xtb 6.7.1): homo_energy
- `ethanol_gfn1_opt.out` (xtb 6.7.1): homo_energy
- `ethanol_grad.out` (xtb 6.7.1): homo_energy
- `methanol_cation_uhf.out` (xtb 6.7.1): homo_energy
- `water_hess_verbose.out` (xtb 6.7.1): homo_energy, vibfreqs
- `water_ohess.out` (xtb 6.7.1): homo_energy, vibfreqs
- `water_opt.out` (xtb 6.7.1): homo_energy
- `water_sp.out` (xtb 6.7.1): homo_energy

## Only chemout

- `1618.out`: scf_energy, n_geometries, final_coordinates, multiplicity, homo_energy, mulliken_charges
- `casscf_beryllium_atom_nosym.out`: homo_energy
- `casscf_beryllium_atom_sym.out`: homo_energy
- `ho-h_cas_aDZ_scanSA.out`: homo_energy
- `1509.out`: multiplicity
- `1177.out`: multiplicity
- `casscf_beryllium_atom_nosym.out`: homo_energy
- `casscf_beryllium_atom_sym.out`: homo_energy
- `casscf.out`: homo_energy
- `casscf_beryllium_atom_nosym.out`: homo_energy
- `casscf_beryllium_atom_sym.out`: homo_energy
- `BH3CO_fragment_guess.log`: success_flag
- `BH3CO_fragment_guess.pop_minimal.log`: success_flag
- `OPT_oniom.log`: homo_energy
- `ansaBPFP2_acetylene_ircf_LQA.out`: success_flag
- `encaps1.log`: homo_energy, mulliken_charges
- `issue_460.log`: success_flag
- `oniom_IR_intensity.log`: homo_energy
- `stopiter_gaussian.out`: success_flag
- `1768_sdf.out`: multiplicity
- `1768_xyz.out`: multiplicity
- `ChemistryCourses__NPTEL_Computational_Quantum_Chemistry__docker.out`: charge, multiplicity
- `Xingyao-Wu__Chebyshev-Tree-Tensor-Network-for-Vibrational-Spectra__vib.out`: multiplicity
- `cclib__cclib-data__1618.out`: scf_energy, n_geometries, final_coordinates, multiplicity, homo_energy, mulliken_charges
- `jhonapontee__methyl-dimer-spin-coupling__caspt2.out`: homo_energy
- `lucasgarayy__Machine-Learned-Interatomic-Potentials-for-Photophysics__orca.out`: etenergies
- `madsr2d2__comp_chem__BnOH_H2O2.out`: multiplicity
- `yamashun95__orca-work__tsopt.log`: success_flag
- `acetate_alpb.out`: multiplicity
- `ethanol_gfn1_opt.out`: multiplicity
- `ethanol_grad.out`: multiplicity
- `water_hess_verbose.out`: multiplicity
- `water_ohess.out`: multiplicity
- `water_opt.out`: multiplicity
- `water_sp.out`: multiplicity

## Only cclib

- `steom_orbs.log`: etenergies
- `casscf.out`: etenergies
- `1448.out`: charge
- `acetate_alpb.out`: n_geometries, final_coordinates
- `ethanol_gfn1_opt.out`: mulliken_charges
- `ethanol_gfnff_opt.out`: charge
- `ethanol_grad.out`: n_geometries, final_coordinates
- `water_hess_verbose.out`: n_geometries, final_coordinates
- `water_sp.out`: n_geometries, final_coordinates
