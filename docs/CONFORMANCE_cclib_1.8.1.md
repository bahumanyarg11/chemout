# chemout vs cclib 1.8.1 (latest PyPI release, 2024)

Files: 247 (1 skipped: not ORCA/Gaussian/xTB). cclib failed on **119** files (exception, timeout > 300s, or `ccread` returned None); chemout returned a record for every file.

| quantity | agree | disagree | chemout only | cclib only | neither |
|---|---|---|---|---|---|
| scf_energy | 96 | 16 | 2 | 1 | 12 |
| n_geometries | 92 | 26 | 2 | 4 | 3 |
| final_coordinates | 113 | 5 | 2 | 4 | 3 |
| charge | 124 | 0 | 1 | 2 | 0 |
| multiplicity | 113 | 0 | 12 | 0 | 2 |
| homo_energy | 98 | 10 | 10 | 0 | 9 |
| mulliken_charges | 114 | 3 | 3 | 1 | 6 |
| vibfreqs | 20 | 3 | 0 | 0 | 104 |
| etenergies | 4 | 6 | 0 | 2 | 115 |
| success_flag | 125 | 2 | 0 | 0 | 0 |

Counts are over the files where cclib produced a record.

## cclib failures

- `casscf_beryllium_atom_nosym.out` (orca 6.0.0): ?
- `casscf_beryllium_atom_sym.out` (orca 6.0.0): ?
- `commercial_version.out` (orca 6.0.1-f.3): ?
- `ho-h_cas_aDZ_scanSA.out` (orca 6.0.0): ?
- `1509.out` (orca 6.1.1): ?
- `1038.out` (orca 5.0.4): ?
- `1177.out` (orca 5.0.3): ?
- `a_DLPNO.txt` (orca 4.0.1.2): ?
- `Pyridine_opt-eomccsd.log` (gaussian 16 RevC.01): ?
- `co_pbe1pbe_631ppGss.log` (gaussian 16 RevApple): ?
- `mp2_avdz_freq_tight.log` (gaussian 16 RevC.01): ?
- `1768_sdf.out` (xtb 6.7.1pre): ?
- `1768_xyz.out` (xtb 6.7.1pre): ?
- `FoleyLab__nbr_data__ccsdpt_meta.out` (orca 6.1.1): ?
- `FoleyLab__nbr_data__ccsdpt_ortho.out` (orca 6.1.1): ?
- `FoleyLab__nbr_data__ccsdpt_para.out` (orca 6.1.1): ?
- `GitNindjapatrick__Quantum_chemistry_ORCA-Simulations__water.out` (orca 6.0.0): ?
- `HabershonLab__cde__pes_orca.out` (orca 6.1.0): ?
- `JEFF7712__cctop__orca.out` (orca 6.0.0): ?
- `Kev-Park__CBE423__n2orca.out` (orca 6.1.1): ?
- `Kev-Park__CBE423__n3orca.out` (orca 6.1.1): ?
- `Konjkov__pycasino__orca.out` (orca 6.1.0): ?
- `OpenChemistry__avogadrodata__camphor.out` (orca 6.1.0): ?
- `OpenChemistry__avogadrodata__formaldehyde.out` (orca 6.0.1): ?
- `OpenChemistry__avogadrodata__nmr.out` (orca 6.1.0): ?
- `OpenChemistry__avogadrodata__raman.out` (orca 6.1.0): ?
- `PoloTier__PSiNad__orca-QM.log` (orca 6.0.0): ?
- `ReactionMechanismGenerator__ARC__neb_res.out` (orca 6.0.0): ?
- `Spin-Chemistry-Labs__radicalpy__NH2_A.out` (orca 6.0.0): ?
- `Svatunek-Lab__teaching__H2O.out` (orca 6.0.0): ?
- `Tovila01__urea-derivate__adenine_opt_freq.out` (orca 6.1.1): ?
- `Tovila01__urea-derivate__biuret_opt_freq.out` (orca 6.1.1): ?
- `Tovila01__urea-derivate__cytosine_opt_freq.out` (orca 6.1.1): ?
- `Tovila01__urea-derivate__dimethylurea_opt_freq.out` (orca 6.1.1): ?
- `Tovila01__urea-derivate__formamide_opt_freq.out` (orca 6.1.1): ?
- `Tovila01__urea-derivate__guanine_opt_freq.out` (orca 6.1.1): ?
- `Tovila01__urea-derivate__methylurea_opt_freq.out` (orca 6.1.1): ?
- `Tovila01__urea-derivate__thiourea_opt_freq.out` (orca 6.1.1): ?
- `Xingyao-Wu__Chebyshev-Tree-Tensor-Network-for-Vibrational-Spectra__opt.out` (orca 6.1.0): ?
- `adole2008__compchem__watergeoopt.out` (orca 6.1.1): ?
- `ahmetaltunfatih__LEDAW__mono2.mpi16.out` (orca 6.1.x): ?
- `aligfellow__graphRC__dihedral.out` (orca 6.0.1): ?
- `avcopan__project-cyclopentene-oh-addition__opt.log` (orca 6.1.0): ?
- `avcopan__project-cyclopentene-oh-addition__opt_orca.out` (orca 6.1.1): ?
- `avcopan__project-cyclopentene-oh-addition__sp_orca.out` (orca 6.1.1): ?
- `avcopan__project-dimethyl-oxirane-cl-abstraction__opt.log` (orca 6.1.0): ?
- `cclib__cclib__dvb_adc2.log` (orca 6.0.1): ?
- `cclib__cclib__dvb_eom_ccsd.log` (orca 6.0.1): ?
- `cclib__cclib__dvb_ir.out` (orca 6.0.1): ?
- `cclib__cclib__dvb_pno_eom_ccsd.log` (orca 6.0.1): ?
- `cclib__cclib__dvb_raman.out` (orca 6.0.1): ?
- `cclib__cclib__dvb_rocis.out` (orca 6.0.1): ?
- `cclib__cclib__dvb_steom_ccsd.log` (orca 6.0.1): ?
- `cclib__cclib__dvb_steom_dlpno_ccsd.log` (orca 6.0.1): ?
- `chpc-uofu__chpc-myjobs-templates__cadi.out` (orca 6.1.1): ?
- `chpc-uofu__chpc-myjobs-templates__methane.out` (orca 6.1.1): ?
- `chpc-uofu__chpc-myjobs-templates__water.out` (orca 6.1.1): ?
- `chpc-uofu__cluster-validation__cadi.out` (orca 6.0.1): ?
- `dhwan0112__dhk.github.io__water.out` (orca 6.1.1): ?
- `dydtkddl__DFT_HONOKIOL__DMC.out` (orca 6.0.1): ?
- `dydtkddl__DFT_HONOKIOL__EC.out` (orca 6.0.1): ?
- `dydtkddl__DFT_HONOKIOL__EMC.out` (orca 6.0.1): ?
- `dydtkddl__DFT_HONOKIOL__FEC.out` (orca 6.0.1): ?
- `gentle1999__MolOP__water_mp2.out` (orca 6.0.1): ?
- `geoffreyweal__ORCA_Mechanism_Procedure__output.out` (orca 6.1.0-f.0): ?
- `geoffreyweal__viewORCA__output.out` (orca 6.1.0-f.0): ?
- `grimme-lab__q-vSZPs__orca.out` (orca 6.1.1): ?
- `grimme-lab__q-vSZPs__orca_atom8.out` (orca 6.1.1): ?
- `hbar-team__SPyCCI__water_0_1_orca_BP86_def2-TZVP_vacuum_opt.out` (orca 6.1.0): ?
- `jhoanlond__QCXMS2__orca.out` (orca 6.0.1): ?
- `jhoanlond__QCXMS2__orca_atom6.out` (orca 6.0.1): ?
- `jhoanlond__QCXMS2__orca_atom8.out` (orca 6.0.1): ?
- `jhonapontee__methyl-dimer-spin-coupling__caspt2.out` (orca 6.1.0): ?
- `jhonapontee__methyl-dimer-spin-coupling__ch3.out` (orca 6.1.0): ?
- `jhonapontee__methyl-dimer-spin-coupling__d3h.out` (orca 6.1.0): ?
- `leonardopagliochini__Energy_Additivity_Evaluation__AR2LD_2.8.log` (orca 6.0.0): ?
- `leonardopagliochini__Energy_Additivity_Evaluation__AR2LD_3.0.log` (orca 6.0.0): ?
- `leonardopagliochini__Energy_Additivity_Evaluation__AR2LD_3.2.log` (orca 6.0.0): ?
- `leonardopagliochini__Energy_Additivity_Evaluation__AR2LD_3.4.log` (orca 6.0.0): ?
- `leonardopagliochini__Energy_Additivity_Evaluation__AR2LD_3.6.log` (orca 6.0.0): ?
- `leonardopagliochini__Energy_Additivity_Evaluation__AR2LD_3.8.log` (orca 6.0.0): ?
- `leonardopagliochini__Energy_Additivity_Evaluation__AR2LD_4.0.log` (orca 6.0.0): ?
- `leonardopagliochini__Energy_Additivity_Evaluation__AR2LD_4.2.log` (orca 6.0.0): ?
- `ll23002__Orca__Etanol.out` (orca 6.1.1): ?
- `ll23002__Orca__agua.out` (orca 6.1.1): ?
- `ll23002__Orca__water.out` (orca 6.1.0-f.0): ?
- `lucasandradechem-lab__chg2itp__dft.out` (orca 6.1.1): ?
- `lucasgarayy__Machine-Learned-Interatomic-Potentials-for-Photophysics__orca.out` (orca 6.0.0): ?
- `madsr2d2__masterThesis__benzaldehyde_tddft.out` (orca 6.0.1): ?
- `madsr2d2__masterThesis__benzoic_acid_tddft.out` (orca 6.0.1): ?
- `madsr2d2__masterThesis__job.out` (orca 6.1.1): ?
- `mersalas__qm_cmd__h2_frag.out` (orca 6.0.1): ?
- `mersalas__qm_cmd__h3_frag.out` (orca 6.0.1): ?
- `mersalas__qm_cmd__hh-pg_cplx.out` (orca 6.0.1): ?
- `mersalas__qm_cmd__k2_frag.out` (orca 6.0.1): ?
- `mersalas__qm_cmd__kh-pg_cplx.out` (orca 6.0.1): ?
- `mersalas__qm_cmd__pg1_frag.out` (orca 6.0.1): ?
- `mersalas__qm_cmd__r2_frag.out` (orca 6.0.1): ?
- `mersalas__qm_cmd__rh-pg_cplx.out` (orca 6.0.1): ?
- `patonlab__pyDFTD3__ibuprofen_b3lyp_orca.out` (orca 6.1.0): ?
- `patonlab__pyDFTD3__ibuprofen_pbe_orca.out` (orca 6.1.0): ?
- `patonlab__pyQRC__acetaldehyde.out` (orca 6.0.0): ?
- `patonlab__yadda__ch3_rad.out` (orca 6.1.0): ?
- `patonlab__yadda__propene.out` (orca 6.1.0): ?
- `physicien__spectroplot__FRQ_tungsten_hexacarbonyl_f.out` (orca 6.0.1): ?
- `re-gr__rgmol__CH3Cl.out` (orca 6.0.1): ?
- `rgora__photochemistry__uracil_s0.out` (orca 6.0.1): ?
- `rgora__photochemistry__uracil_s0_vee.out` (orca 6.0.1): ?
- `rgora__photochemistry__uracil_s0_vg.out` (orca 6.0.1): ?
- `ushcodex__Vanadium-zeolite-y-dealumination__V-PRC-T1-V03.out` (orca 6.1.0): ?
- `ushcodex__Vanadium-zeolite-y-dealumination__V-TS2-T1-V03.out` (orca 6.1.0): ?
- `ushcodex__Vanadium-zeolite-y-dealumination__irc_backward_opt.out` (orca 6.1.0): ?
- `ushcodex__Vanadium-zeolite-y-dealumination__irc_forward_opt.out` (orca 6.1.0): ?
- `yamashun95__orca-work__tsopt.log` (orca 6.1.1): ?
- `yamashun95__orca-work__wateropt.log` (orca 6.1.1): ?
- `yamashun95__orca-work__waterpcm.log` (orca 6.1.1): ?
- `zwang20__honours-data-7__censo.out` (orca 6.0.0): ?
- `g98.out` (gaussian None): ?
- `g98.out` (gaussian None): ?

## Disagreements

- `ADBNA_Me_Mes_MesCz.log` (orca 5.0.3): etenergies
- `Benzene_opt_etsyms.log` (orca 5.0.3): etenergies
- `dvb_soc.log` (orca 5.0.3): etenergies
- `947.out` (orca 4.2.1): scf_energy
- `cope.out` (orca 4.2.1): scf_energy
- `dvb_dispersion_bp86_d3zero.out` (orca 4.2.0): scf_energy
- `dvb_rocis.out` (orca 4.2.0): etenergies
- `dvb_td.out` (orca 4.2.0): etenergies
- `ligando-30-SRM1-S_ZINDO.out` (orca 4.2.1): etenergies
- `Benzene_anion.log` (gaussian 16 RevC.01): n_geometries
- `dol-1-pen-5-pen-trip.out` (gaussian 16 RevC.01): n_geometries
- `issue1153.log` (gaussian 16 RevC.01): n_geometries
- `issue962.log` (gaussian 16 RevA.03): scf_energy, n_geometries, final_coordinates
- `water.log` (gaussian 16 RevC.01): n_geometries
- `water_cation_nbo_opt.out` (gaussian 16 RevA.03): n_geometries
- `water_neutral_nbo_opt.out` (gaussian 16 RevA.03): n_geometries
- `water_neutral_nbo_opt_freq.out` (gaussian 16 RevA.03): n_geometries
- `1080.out` (gaussian 09 RevC.01): n_geometries
- `1504.log` (gaussian 09 RevA.02): n_geometries
- `25DMF_HRANH.log` (gaussian 09 RevB.01): n_geometries
- `2D-PES-all-converged.log` (gaussian 09 RevD.01): n_geometries
- `2D-PES-one-unconverged.log` (gaussian 09 RevD.01): n_geometries
- `534.out` (gaussian 09 RevA.02): n_geometries
- `BH3CO_fragment_guess.log` (gaussian 09 RevD.01): success_flag
- `BH3CO_fragment_guess.pop_minimal.log` (gaussian 09 RevD.01): homo_energy, mulliken_charges, success_flag
- `BSL_opt_freq_DFT.out` (gaussian 09 RevD.01): n_geometries
- `OPT_oniom.log` (gaussian 09 RevD.01): scf_energy, n_geometries, mulliken_charges
- `ansaBPFP2_acetylene_ircf_LQA.out` (gaussian 09 RevE.01): n_geometries, final_coordinates
- `benzene_excited_states_optimization_issue889.log` (gaussian 09 RevD.01): n_geometries
- `dvb_gopt_revA.02.out` (gaussian 09 RevA.02): n_geometries
- `dvb_gopt_unconverged.log` (gaussian 09 RevD.01): n_geometries
- `dvb_scan_revA.02.log` (gaussian 09 RevA.02): n_geometries
- `encaps1.log` (gaussian 09 RevC.01): scf_energy, n_geometries
- `irc_point.log` (gaussian 09 RevD.01): vibfreqs
- `issue1150.log` (gaussian 09 RevE.01): final_coordinates
- `issue_460.log` (gaussian 09 RevD.01): n_geometries, final_coordinates
- `jos.teunissen.out` (gaussian 09 RevD.01): n_geometries
- `oniom_IR_intensity.log` (gaussian 09 RevC.01): scf_energy, mulliken_charges
- `relaxed_PES_testCO2.log` (gaussian 09 RevD.01): n_geometries
- `relaxed_PES_testH2.log` (gaussian 09 RevD.01): n_geometries
- `1448.out` (xtb 6.5.1): scf_energy, final_coordinates, homo_energy
- `acetate_alpb.out` (xtb 6.7.1): scf_energy, homo_energy
- `ethanol_gfn1_opt.out` (xtb 6.7.1): scf_energy, homo_energy
- `ethanol_grad.out` (xtb 6.7.1): scf_energy, homo_energy
- `methanol_cation_uhf.out` (xtb 6.7.1): scf_energy, homo_energy
- `water_hess_verbose.out` (xtb 6.7.1): scf_energy, homo_energy, vibfreqs
- `water_ohess.out` (xtb 6.7.1): scf_energy, homo_energy, vibfreqs
- `water_opt.out` (xtb 6.7.1): scf_energy, homo_energy
- `water_sp.out` (xtb 6.7.1): scf_energy, homo_energy

## Only chemout

- `1618.out`: scf_energy, n_geometries, final_coordinates, multiplicity, homo_energy, mulliken_charges
- `casscf_beryllium_atom_nosym.out`: homo_energy
- `casscf_beryllium_atom_sym.out`: homo_energy
- `casscf.out`: homo_energy
- `casscf_beryllium_atom_nosym.out`: homo_energy
- `casscf_beryllium_atom_sym.out`: homo_energy
- `OPT_oniom.log`: homo_energy
- `encaps1.log`: homo_energy, mulliken_charges
- `oniom_IR_intensity.log`: homo_energy
- `ChemistryCourses__NPTEL_Computational_Quantum_Chemistry__docker.out`: charge, multiplicity
- `Xingyao-Wu__Chebyshev-Tree-Tensor-Network-for-Vibrational-Spectra__vib.out`: multiplicity
- `cclib__cclib-data__1618.out`: scf_energy, n_geometries, final_coordinates, multiplicity, homo_energy, mulliken_charges
- `madsr2d2__comp_chem__BnOH_H2O2.out`: multiplicity
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
- `ethanol_gfnff_opt.out`: scf_energy, charge
- `ethanol_grad.out`: n_geometries, final_coordinates
- `water_hess_verbose.out`: n_geometries, final_coordinates
- `water_sp.out`: n_geometries, final_coordinates
