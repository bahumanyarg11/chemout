"""Map a Calculation onto cclib's attribute names and units.

cclib conventions (https://cclib.github.io/data.html): coordinates in Angstrom, energies in eV,
orbital energies in eV, frequencies in cm^-1, excitation energies in cm^-1, charges in e.
``to_ccdata_dict`` returns a plain dict; ``to_ccdata`` builds a real ``cclib.parser.data.ccData``
when cclib is importable.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from ..schema import HARTREE_TO_EV, Calculation

_EV_TO_CM = 8065.543937


def to_ccdata_dict(calc: Calculation) -> dict[str, Any]:
    d: dict[str, Any] = {}
    mol = calc.molecule
    if mol.atomic_numbers:
        d["natom"] = len(mol.atomic_numbers)
        d["atomnos"] = np.asarray(mol.atomic_numbers, dtype=int)
    if mol.masses:
        d["atommasses"] = np.asarray(mol.masses)
    if calc.input.charge is not None:
        d["charge"] = calc.input.charge
    if calc.input.multiplicity is not None:
        d["mult"] = calc.input.multiplicity
    if mol.n_basis_functions is not None:
        d["nbasis"] = mol.n_basis_functions
    if calc.geometries:
        d["atomcoords"] = np.asarray([g.coordinates for g in calc.geometries], dtype=float)
    # cclib's scfenergies: one SCF-level energy per SCF (never a post-HF total).
    #  - Gaussian: "SCF Done" energies.
    #  - ORCA: the SCF block's "Total Energy" plus the dispersion correction, which every ORCA version prints
    #    separately and leaves out of "Total Energy" (FINAL SINGLE POINT ENERGY = Total + dispersion + gCP).
    #  - xTB: the "* total energy" of every optimisation step (the SCC energy alone is not a total energy).
    prog = calc.program.name.value
    if prog == "xtb":
        steps = [s.energy for s in (calc.optimization.steps if calc.optimization else []) if s.energy is not None]
        if calc.energies.total is not None:
            steps = steps[:-1] + [calc.energies.total] if steps else [calc.energies.total]
        if steps:
            d["scfenergies"] = np.asarray(steps) * HARTREE_TO_EV
    else:
        scf = [
            c.energy + (c.dispersion if prog == "orca" and c.dispersion is not None else 0.0)
            for c in calc.scf_cycles
            if c.energy is not None
        ]
        if scf:
            d["scfenergies"] = np.asarray(scf) * HARTREE_TO_EV
    if calc.energies.mp2 is not None:
        d["mpenergies"] = np.asarray([[calc.energies.mp2 * HARTREE_TO_EV]])
    cc = calc.energies.ccsd_t if calc.energies.ccsd_t is not None else calc.energies.ccsd
    if cc is not None:
        d["ccenergies"] = np.asarray([cc * HARTREE_TO_EV])
    if calc.energies.dispersion is not None:
        d["dispersionenergies"] = np.asarray([calc.energies.dispersion * HARTREE_TO_EV])
    if calc.orbitals and calc.orbitals.energies:
        orb = calc.orbitals
        if orb.indices:
            # Only a window of orbitals was printed (xTB): expand to the full orbital count with NaN gaps so
            # that moenergies[homos] is still the HOMO.
            n_mo = max(mol.n_basis_functions or 0, max(max(ix) for ix in orb.indices))
            full = []
            homos = []
            for spin, (e, ix) in enumerate(zip(orb.energies, orb.indices)):
                arr = np.full(n_mo, np.nan)
                for val, i in zip(e, ix):
                    arr[i - 1] = val * HARTREE_TO_EV
                full.append(arr)
                h = orb.homo[spin] if spin < len(orb.homo) else -1
                homos.append(ix[h] - 1 if 0 <= h < len(ix) else -1)
            d["moenergies"] = full
            d["homos"] = np.asarray(homos, dtype=int)
        else:
            d["moenergies"] = [np.asarray(e) * HARTREE_TO_EV for e in orb.energies]
            d["homos"] = np.asarray(orb.homo, dtype=int)
        if calc.orbitals.symmetries:
            d["mosyms"] = calc.orbitals.symmetries
    if calc.populations:
        d["atomcharges"] = {p.method: np.asarray(p.charges) for p in calc.populations}
        spins = {p.method: np.asarray(p.spins) for p in calc.populations if p.spins is not None}
        if spins:
            d["atomspins"] = spins
    if calc.multipoles and calc.multipoles.dipole is not None:
        moments = [np.zeros(3), np.asarray(calc.multipoles.dipole)]
        if calc.multipoles.quadrupole is not None:
            moments.append(np.asarray(calc.multipoles.quadrupole))
        d["moments"] = moments
    if calc.gradient is not None:
        d["grads"] = np.asarray([calc.gradient], dtype=float)
    if calc.optimization:
        d["optdone"] = bool(calc.optimization.converged)
        if calc.optimization.steps:
            d["optstatus"] = np.asarray(
                [1] + [0] * (len(calc.optimization.steps) - 2) + [4 if calc.optimization.converged else 0]
                if len(calc.optimization.steps) > 1
                else [5 if calc.optimization.converged else 1],
                dtype=int,
            )
    v = calc.vibrations
    if v and v.frequencies:
        d["vibfreqs"] = np.asarray(v.frequencies)
        if v.ir_intensities is not None:
            d["vibirs"] = np.asarray(v.ir_intensities)
        if v.raman_activities is not None:
            d["vibramans"] = np.asarray(v.raman_activities)
        if v.reduced_masses is not None:
            d["vibrmasses"] = np.asarray(v.reduced_masses)
        if v.force_constants is not None:
            d["vibfconsts"] = np.asarray(v.force_constants)
        if v.symmetries is not None:
            d["vibsyms"] = list(v.symmetries)
        if v.normal_modes is not None:
            d["vibdisps"] = np.asarray(v.normal_modes)
        t = v.thermochemistry
        if t:
            if t.temperature is not None:
                d["temperature"] = t.temperature
            if t.pressure is not None:
                d["pressure"] = t.pressure
            if t.zero_point_energy is not None:
                d["zpve"] = t.zero_point_energy
            if t.enthalpy is not None:
                d["enthalpy"] = t.enthalpy
            if t.gibbs_free_energy is not None:
                d["freeenergy"] = t.gibbs_free_energy
            if t.entropy is not None:
                d["entropy"] = t.entropy
    if calc.excited_states:
        d["etenergies"] = np.asarray([s.energy_ev * _EV_TO_CM for s in calc.excited_states])
        d["etoscs"] = np.asarray(
            [s.oscillator_strength if s.oscillator_strength is not None else np.nan for s in calc.excited_states]
        )
        syms = []
        for s in calc.excited_states:
            m = (s.multiplicity or "").capitalize()
            syms.append(f"{m}-{s.symmetry}" if m and s.symmetry else (m or s.symmetry or "?"))
        d["etsyms"] = syms
        d["etsecs"] = [
            [
                (
                    (t.from_orbital, 0 if (t.from_spin or "a") == "a" else 1),
                    (t.to_orbital, 0 if (t.to_spin or "a") == "a" else 1),
                    t.coefficient,
                )
                for t in s.transitions
            ]
            for s in calc.excited_states
        ]
        if any(s.rotatory_strength is not None for s in calc.excited_states):
            d["etrotats"] = np.asarray(
                [s.rotatory_strength if s.rotatory_strength is not None else np.nan for s in calc.excited_states]
            )
    if calc.nmr:
        d["nmrtensors"] = {n.atom_index: {"isotropic": n.isotropic, "anisotropy": n.anisotropy} for n in calc.nmr}
    d["metadata"] = {
        "package": {"orca": "ORCA", "gaussian": "Gaussian", "xtb": "xTB", "psi4": "Psi4"}.get(
            calc.program.name.value, calc.program.name.value
        ),
        "package_version": calc.program.version,
        "success": calc.success,
        "methods": [calc.input.method] if calc.input.method else [],
        "basis_set": calc.input.basis,
        "functional": calc.input.functional,
        "num_cpu": calc.input.nprocs,
        "wall_time": calc.timing.wall_seconds if calc.timing else None,
        "cpu_time": calc.timing.cpu_seconds if calc.timing else None,
    }
    return d


def to_ccdata(calc: Calculation):
    """Return a ``cclib.parser.data.ccData`` populated from ``calc`` (requires cclib)."""
    from cclib.parser.data import ccData

    attrs = to_ccdata_dict(calc)
    data = ccData()
    for k, v in attrs.items():
        if k in ccData._attrlist or k == "metadata":
            setattr(data, k, v)
    return data
