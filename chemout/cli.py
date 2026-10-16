"""Command line interface.

chemout parse FILE [FILE...] [-o out.json | -o out.h5 | -o out.jsonl]
chemout info FILE                 human-readable summary
chemout batch DIR [-o out.jsonl|out.h5] [--jobs N] [--csv summary.csv]
chemout sniff FILE                print the detected program
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from . import __version__, parse, sniff_program
from .parse import read_text
from .schema import HARTREE_TO_EV, Calculation

EXTS = {".out", ".log", ".txt", ".output", ".orcaout", ".gout"}


def _write(calc: Calculation, out: str | None, append: bool = False) -> None:
    if out is None or out == "-":
        sys.stdout.write(calc.to_json() + "\n")
        return
    suffix = Path(out).suffix.lower()
    if suffix in (".h5", ".hdf5"):
        group = Path(calc.provenance.path or "calc").stem if append else None
        calc.to_hdf5(out, group=group)
    elif suffix == ".jsonl":
        with open(out, "a" if append else "w") as fh:
            fh.write(calc.to_json(indent=None) + "\n")
    else:
        Path(out).write_text(calc.to_json() + "\n")


def cmd_parse(args: argparse.Namespace) -> int:
    rc = 0
    multi = len(args.files) > 1
    for i, f in enumerate(args.files):
        calc = parse(f, program=args.program)
        if calc.provenance.warnings and not args.quiet:
            for w in calc.provenance.warnings:
                print(f"warning: {f}: [{w.section}] {w.message}", file=sys.stderr)
        if args.strict and calc.provenance.warnings:
            rc = 2
        _write(
            calc,
            args.output,
            append=multi
            and i > 0
            or (multi and args.output is not None and i == 0 and Path(args.output).suffix in (".h5", ".hdf5")),
        )
    return rc


def summary_row(calc: Calculation) -> dict[str, object]:
    e = calc.energies
    v = calc.vibrations
    return {
        "path": calc.provenance.path,
        "program": calc.program.name.value,
        "version": calc.program.version,
        "termination": calc.termination.value,
        "method": calc.input.method,
        "basis": calc.input.basis,
        "charge": calc.input.charge,
        "multiplicity": calc.input.multiplicity,
        "n_atoms": calc.molecule.n_atoms,
        "n_geometries": len(calc.geometries),
        "scf_energy_hartree": e.scf,
        "total_energy_hartree": e.total,
        "gibbs_hartree": e.gibbs,
        "homo_lumo_gap_ev": (calc.orbitals.homo_lumo_gap() or 0) * HARTREE_TO_EV
        if calc.orbitals and calc.orbitals.homo_lumo_gap() is not None
        else None,
        "dipole_debye": calc.multipoles.dipole_magnitude if calc.multipoles else None,
        "n_imaginary": v.n_imaginary if v else None,
        "lowest_frequency_cm": min(v.frequencies) if v and v.frequencies else None,
        "n_excited_states": len(calc.excited_states),
        "wall_seconds": calc.timing.wall_seconds if calc.timing else None,
        "n_warnings": len(calc.provenance.warnings),
    }


def cmd_info(args: argparse.Namespace) -> int:
    calc = parse(args.file, program=args.program)
    row = summary_row(calc)
    width = max(len(k) for k in row)
    for k, val in row.items():
        if val is None:
            continue
        if isinstance(val, float):
            val = f"{val:.6f}" if abs(val) < 1e6 else f"{val:.6g}"
        print(f"{k:<{width}}  {val}")
    if calc.input.job_types:
        print(f"{'job_types':<{width}}  {', '.join(j.value for j in calc.input.job_types)}")
    if calc.populations:
        print(f"{'populations':<{width}}  {', '.join(p.method for p in calc.populations)}")
    if calc.optimization:
        print(f"{'opt_converged':<{width}}  {calc.optimization.converged} ({calc.optimization.n_steps} steps)")
    if calc.error_message:
        print(f"{'error':<{width}}  {calc.error_message}")
    for w in calc.provenance.warnings:
        print(f"{'warning':<{width}}  [{w.section}] {w.message}" + (f" (line {w.line})" if w.line else ""))
    return 0


def _parse_one(path: str) -> tuple[str, str | None, str | None]:
    try:
        calc = parse(path)
        return path, calc.to_json(indent=None), None
    except Exception as exc:  # noqa: BLE001 - a crash in one file must not kill the batch
        return path, None, f"{type(exc).__name__}: {exc}"


def cmd_batch(args: argparse.Namespace) -> int:
    root = Path(args.directory)
    files = sorted(str(p) for p in root.rglob("*") if p.suffix.lower() in EXTS and p.is_file())
    if not files:
        print(f"no output files under {root}", file=sys.stderr)
        return 1
    out = args.output
    if out and Path(out).exists() and not args.append:
        os.remove(out)
    rows: list[dict[str, object]] = []
    done = 0
    errors = 0
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(_parse_one, f): f for f in files}
        for fut in as_completed(futures):
            path, js, err = fut.result()
            done += 1
            if err:
                errors += 1
                print(f"ERROR {path}: {err}", file=sys.stderr)
                continue
            calc = Calculation.from_json(js)
            rows.append(summary_row(calc))
            if out:
                try:
                    _write(calc, out, append=True)
                except Exception as exc:  # noqa: BLE001 - one unwritable record must not stop the batch
                    errors += 1
                    print(f"ERROR writing {path}: {type(exc).__name__}: {exc}", file=sys.stderr)
            if not args.quiet and done % 50 == 0:
                print(f"{done}/{len(files)}", file=sys.stderr)
    if args.csv:
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else ["path"])
            w.writeheader()
            w.writerows(rows)
    normal = sum(1 for r in rows if r["termination"] == "normal")
    print(f"{len(files)} files: {normal} normal, {len(rows) - normal} other, {errors} parser errors", file=sys.stderr)
    return 1 if errors else 0


def cmd_sniff(args: argparse.Namespace) -> int:
    print(sniff_program(read_text(args.file)).value)
    return 0


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="chemout", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--version", action="version", version=f"chemout {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("parse", help="parse one or more files to JSON / JSONL / HDF5")
    p.add_argument("files", nargs="+")
    p.add_argument("-o", "--output", help="output path (.json, .jsonl or .h5); default stdout")
    p.add_argument("--program", choices=["orca", "gaussian", "xtb"], help="override program detection")
    p.add_argument("--strict", action="store_true", help="exit 2 if any parse warning was recorded")
    p.add_argument("-q", "--quiet", action="store_true")
    p.set_defaults(func=cmd_parse)

    p = sub.add_parser("info", help="print a human-readable summary")
    p.add_argument("file")
    p.add_argument("--program", choices=["orca", "gaussian", "xtb"])
    p.set_defaults(func=cmd_info)

    p = sub.add_parser("batch", help="parse every output file under a directory in parallel")
    p.add_argument("directory")
    p.add_argument("-o", "--output", help=".jsonl or .h5 collecting every record")
    p.add_argument("--csv", help="write a one-row-per-file summary table")
    p.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    p.add_argument("--append", action="store_true", help="append to an existing output instead of replacing it")
    p.add_argument("-q", "--quiet", action="store_true")
    p.set_defaults(func=cmd_batch)

    p = sub.add_parser("sniff", help="print which program wrote a file")
    p.add_argument("file")
    p.set_defaults(func=cmd_sniff)
    return ap


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
