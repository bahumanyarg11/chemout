"""HDF5 export: one group per calculation, arrays as datasets, scalars/strings as attributes.

Layout (``/`` or the requested group)::

    attrs: program, version, termination, path, sha256, chemout_version, json (full record)
    molecule/atomic_numbers          int   (n_atoms,)
    geometries/coordinates           float (n_geom, n_atoms, 3)   Angstrom
                                     (or coordinates_0000, ... when the atom count changes within the file)
    geometries/energies              float (n_geom,)              Hartree (nan when unknown)
    scf/energies                     float (n_scf,)
    orbitals/energies                float (n_spin, n_mo)
    orbitals/occupations             float (n_spin, n_mo)
    populations/<method>/charges     float (n_atoms,)
    gradient                         float (n_atoms, 3)
    vibrations/frequencies, ir_intensities, normal_modes
    excited_states/energies_ev, oscillator_strengths

The complete JSON record is stored in the ``json`` attribute so nothing is lost; the datasets are
the fast path for analysis across thousands of files.
"""

from __future__ import annotations

import math

import numpy as np

from ..schema import Calculation


def _pad(rows: list[list[float]]) -> np.ndarray:
    """Stack per-spin rows into a 2-D array, padding a shorter spin block with NaN."""
    width = max(len(r) for r in rows)
    out = np.full((len(rows), width), np.nan)
    for i, r in enumerate(rows):
        out[i, : len(r)] = r
    return out


def write_hdf5(calc: Calculation, path: str, group: str | None = None, mode: str = "a") -> None:
    try:
        import h5py
    except ImportError as exc:  # pragma: no cover
        raise ImportError("HDF5 export needs h5py: pip install 'chemout[hdf5]'") from exc

    with h5py.File(path, mode) as f:
        g = f.require_group(group) if group else f
        for k in list(g.keys()):
            del g[k]
        g.attrs["program"] = calc.program.name.value
        g.attrs["version"] = calc.program.version or ""
        g.attrs["termination"] = calc.termination.value
        g.attrs["path"] = calc.provenance.path or ""
        g.attrs["sha256"] = calc.provenance.sha256 or ""
        g.attrs["chemout_version"] = calc.provenance.chemout_version
        g.attrs["json"] = calc.to_json(indent=None)

        if calc.molecule.atomic_numbers:
            g.create_dataset("molecule/atomic_numbers", data=np.asarray(calc.molecule.atomic_numbers, dtype=np.int32))
        if calc.geometries:
            energies = np.asarray(
                [geo.energy if geo.energy is not None else math.nan for geo in calc.geometries], dtype=float
            )
            g.create_dataset("geometries/energies", data=energies)
            sizes = {len(geo.coordinates) for geo in calc.geometries}
            if len(sizes) == 1:
                g.create_dataset(
                    "geometries/coordinates", data=np.asarray([geo.coordinates for geo in calc.geometries], dtype=float)
                )
            else:
                # Ragged (atom count changes within the file, e.g. fragments or scans): one dataset per geometry.
                for i, geo in enumerate(calc.geometries):
                    g.create_dataset(f"geometries/coordinates_{i:04d}", data=np.asarray(geo.coordinates, dtype=float))
        if calc.scf_cycles:
            g.create_dataset(
                "scf/energies",
                data=np.asarray([c.energy if c.energy is not None else math.nan for c in calc.scf_cycles]),
            )
        if calc.orbitals and calc.orbitals.energies:
            g.create_dataset("orbitals/energies", data=_pad(calc.orbitals.energies))
            if calc.orbitals.occupations:
                g.create_dataset("orbitals/occupations", data=_pad(calc.orbitals.occupations))
        for pop in calc.populations:
            g.create_dataset(f"populations/{pop.method}/charges", data=np.asarray(pop.charges, dtype=float))
            if pop.spins is not None:
                g.create_dataset(f"populations/{pop.method}/spins", data=np.asarray(pop.spins, dtype=float))
        if calc.gradient is not None:
            g.create_dataset("gradient", data=np.asarray(calc.gradient, dtype=float))
        if calc.vibrations:
            v = calc.vibrations
            g.create_dataset("vibrations/frequencies", data=np.asarray(v.frequencies, dtype=float))
            if v.ir_intensities is not None:
                g.create_dataset("vibrations/ir_intensities", data=np.asarray(v.ir_intensities, dtype=float))
            if v.normal_modes is not None:
                g.create_dataset("vibrations/normal_modes", data=np.asarray(v.normal_modes, dtype=float))
        if calc.excited_states:
            g.create_dataset("excited_states/energies_ev", data=np.asarray([s.energy_ev for s in calc.excited_states]))
            g.create_dataset(
                "excited_states/oscillator_strengths",
                data=np.asarray(
                    [
                        s.oscillator_strength if s.oscillator_strength is not None else math.nan
                        for s in calc.excited_states
                    ]
                ),
            )
        for key, val in calc.energies.model_dump().items():
            if isinstance(val, (int, float)) and val is not None:
                g.attrs[f"energy_{key}"] = float(val)


def read_hdf5(path: str, group: str | None = None) -> Calculation:
    import h5py

    with h5py.File(path, "r") as f:
        g = f[group] if group else f
        return Calculation.from_json(g.attrs["json"])
