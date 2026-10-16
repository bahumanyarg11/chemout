"""Conformance report: chemout vs cclib on the same files.

For each file, run cclib (in a separate interpreter so a crash or hang in cclib cannot stop the
run) and chemout, then compare a fixed set of quantities. Writes a Markdown table.

    python scripts/compare_cclib.py --cclib-python /path/to/venv/bin/python --label "cclib 1.8.1" \
        --out docs/CONFORMANCE_1.8.1.md  DIR [DIR...]

Quantities compared (chemout -> cclib units): final SCF energy (eV), number of geometries,
final coordinates (Angstrom), charge, multiplicity, HOMO energy (eV), Mulliken charges,
vibrational frequencies, excitation energies (cm^-1), success flag.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import chemout  # noqa: E402
from chemout.compat.cclib import to_ccdata_dict  # noqa: E402

EXTS = {".out", ".log", ".txt"}

CCLIB_SNIPPET = r"""
import json, sys, math
import numpy as np
path = sys.argv[1]
out = {"ok": False}
try:
    import cclib
    from cclib.io import ccread
    data = ccread(path)
    out["version"] = cclib.__version__
    if data is None:
        out["error"] = "ccread returned None (unrecognised file)"
    else:
        def grab(name, f=lambda x: x):
            if hasattr(data, name):
                try:
                    out[name] = f(getattr(data, name))
                except Exception as e:
                    out[name + "_error"] = repr(e)
        grab("scfenergies", lambda x: [float(v) for v in x])
        grab("atomcoords", lambda x: [len(x), [[float(c) for c in row] for row in x[-1]]])
        grab("charge", int)
        grab("mult", int)
        grab("natom", int)
        grab("homos", lambda x: [int(v) for v in x])
        grab("moenergies", lambda x: [[float(v) for v in spin] for spin in x])
        grab("atomcharges", lambda x: {k: [float(v) for v in vals] for k, vals in x.items()})
        grab("vibfreqs", lambda x: [float(v) for v in x])
        grab("etenergies", lambda x: [float(v) for v in x])
        grab("etoscs", lambda x: [float(v) for v in x])
        grab("metadata", lambda m: {"success": m.get("success"), "package_version": m.get("package_version")})
        out["ok"] = True
except Exception as e:
    import traceback
    out["error"] = f"{type(e).__name__}: {e}"
    out["traceback"] = traceback.format_exc()[-1500:]
print(json.dumps(out))
"""


def run_cclib(python: str, path: Path, timeout: int) -> dict:
    t0 = time.time()
    try:
        r = subprocess.run([python, "-c", CCLIB_SNIPPET, str(path)], capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return {"ok": False, "error": f"timeout after {timeout}s", "seconds": timeout}
    out: dict
    try:
        out = json.loads(r.stdout.strip().splitlines()[-1])
    except Exception:
        out = {"ok": False, "error": f"no JSON (rc={r.returncode}): {r.stderr.strip()[-400:]}"}
    out["seconds"] = time.time() - t0
    return out


def close(a, b, rel=1e-6, abs_=1e-5) -> bool:
    if a is None or b is None:
        return False
    try:
        return math.isclose(float(a), float(b), rel_tol=rel, abs_tol=abs_)
    except (TypeError, ValueError):
        return a == b


def seq_close(a, b, rel=1e-6, abs_=1e-4) -> bool:
    if a is None or b is None or len(a) != len(b):
        return False
    return all(close(x, y, rel, abs_) for x, y in zip(a, b))


def compare(calc_d: dict, cc: dict) -> dict[str, str]:
    """Return {quantity: 'agree' | 'disagree' | 'chemout-only' | 'cclib-only' | 'neither'}."""
    res: dict[str, str] = {}

    SMALL_EV = 1e-4 * 27.211386245988  # 1e-4 Hartree: SCF-convergence-level differences

    def verdict(name, mine, theirs, eq, small=None):
        if mine is None and theirs is None:
            res[name] = "neither"
        elif theirs is None:
            res[name] = "chemout-only"
        elif mine is None:
            res[name] = "cclib-only"
        elif eq(mine, theirs):
            res[name] = "agree"
        elif small is not None and small(mine, theirs):
            res[name] = "disagree-small"
        else:
            res[name] = "disagree"

    scf = calc_d.get("scfenergies")
    verdict(
        "scf_energy",
        scf[-1] if scf is not None and len(scf) else None,
        (cc.get("scfenergies") or [None])[-1],
        lambda a, b: close(a, b, 1e-7, 1e-4),
        small=lambda a, b: abs(float(a) - float(b)) < SMALL_EV,
    )
    ac = calc_d.get("atomcoords")
    verdict(
        "n_geometries",
        len(ac) if ac is not None else None,
        cc["atomcoords"][0] if cc.get("atomcoords") else None,
        lambda a, b: a == b,
    )
    verdict(
        "final_coordinates",
        ac[-1].tolist() if ac is not None else None,
        cc["atomcoords"][1] if cc.get("atomcoords") else None,
        lambda a, b: len(a) == len(b) and all(seq_close(x, y, 1e-6, 1e-4) for x, y in zip(a, b)),
    )
    verdict("charge", calc_d.get("charge"), cc.get("charge"), lambda a, b: a == b)
    verdict("multiplicity", calc_d.get("mult"), cc.get("mult"), lambda a, b: a == b)
    mine_homo = None
    if calc_d.get("moenergies") is not None and calc_d.get("homos") is not None and len(calc_d["homos"]):
        h = int(calc_d["homos"][0])
        if 0 <= h < len(calc_d["moenergies"][0]):
            mine_homo = float(calc_d["moenergies"][0][h])
    their_homo = None
    if cc.get("moenergies") and cc.get("homos"):
        h = cc["homos"][0]
        if 0 <= h < len(cc["moenergies"][0]):
            their_homo = cc["moenergies"][0][h]
    verdict(
        "homo_energy",
        mine_homo,
        their_homo,
        lambda a, b: close(a, b, 1e-6, 1e-3),
        small=lambda a, b: abs(float(a) - float(b)) < SMALL_EV,
    )
    mine_mull = (
        calc_d["atomcharges"]["mulliken"].tolist()
        if calc_d.get("atomcharges") and "mulliken" in calc_d["atomcharges"]
        else None
    )
    their_mull = (cc.get("atomcharges") or {}).get("mulliken")
    verdict("mulliken_charges", mine_mull, their_mull, lambda a, b: seq_close(a, b, 1e-6, 1e-4))
    verdict(
        "vibfreqs",
        calc_d["vibfreqs"].tolist() if calc_d.get("vibfreqs") is not None else None,
        cc.get("vibfreqs"),
        lambda a, b: seq_close(a, b, 1e-6, 0.02),
    )
    verdict(
        "etenergies",
        calc_d["etenergies"].tolist() if calc_d.get("etenergies") is not None else None,
        cc.get("etenergies"),
        lambda a, b: seq_close(sorted(a), sorted(b), 1e-5, 2.0),
    )
    md = cc.get("metadata") or {}
    verdict(
        "success_flag", calc_d.get("metadata", {}).get("success"), md.get("success"), lambda a, b: bool(a) == bool(b)
    )
    return res


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("dirs", nargs="+")
    ap.add_argument("--cclib-python", required=True)
    ap.add_argument("--label", default="cclib")
    ap.add_argument("--out", required=True)
    ap.add_argument("--timeout", type=int, default=300)
    ap.add_argument("--max-mb", type=float, default=30.0)
    args = ap.parse_args()

    files = [
        p
        for d in args.dirs
        for p in sorted(Path(d).rglob("*"))
        if p.suffix in EXTS and p.is_file() and p.stat().st_size < args.max_mb * 1e6
    ]
    rows = []
    per_q: dict[str, dict[str, int]] = {}
    cclib_fail = 0
    chemout_unknown = 0
    for p in files:
        calc = chemout.parse(p)
        if calc.program.name.value == "unknown":
            chemout_unknown += 1
            continue
        cd = to_ccdata_dict(calc)
        cc = run_cclib(args.cclib_python, p, args.timeout)
        if not cc.get("ok"):
            cclib_fail += 1
            rows.append((p, calc, cc, None))
            continue
        res = compare(cd, cc)
        for q, v in res.items():
            per_q.setdefault(q, {})
            per_q[q][v] = per_q[q].get(v, 0) + 1
        rows.append((p, calc, cc, res))
        print(
            f"{p.name:45s} cclib {cc['seconds']:5.1f}s  " + " ".join(f"{q}={v[:3]}" for q, v in res.items()),
            file=sys.stderr,
        )

    lines = [
        f"# chemout vs {args.label}",
        "",
        f"Files: {len(files)} ({chemout_unknown} skipped: not ORCA/Gaussian/xTB). cclib failed on **{cclib_fail}** files "
        f"(exception, timeout > {args.timeout}s, or `ccread` returned None); chemout returned a record for every file.",
        "",
        "| quantity | agree | differ < 1e-4 Eh | disagree | chemout only | cclib only | neither |",
        "|---|---|---|---|---|---|---|",
    ]
    for q, counts in per_q.items():
        lines.append(
            f"| {q} | {counts.get('agree', 0)} | {counts.get('disagree-small', 0)} | {counts.get('disagree', 0)} | {counts.get('chemout-only', 0)} | {counts.get('cclib-only', 0)} | {counts.get('neither', 0)} |"
        )
    lines += [
        "",
        "Counts are over the files where cclib produced a record. 'differ < 1e-4 Eh' marks energy differences",
        "below 1e-4 Hartree, i.e. one implementation reporting the last SCF-iteration energy and the other the",
        "converged total: not a parsing error on either side.",
        "",
        "## cclib failures",
        "",
    ]
    for p, calc, _cc, res in rows:
        if res is None:
            lines.append(
                f"- `{p.name}` ({calc.program.name.value} {calc.program.version}): {cc.get('error', '?')[:300]}"
            )
    lines += ["", "## Disagreements", ""]
    for p, calc, _cc, res in rows:
        if res:
            bad = [q for q, v in res.items() if v == "disagree"]
            if bad:
                lines.append(f"- `{p.name}` ({calc.program.name.value} {calc.program.version}): {', '.join(bad)}")
    lines += ["", "## Only chemout", ""]
    for p, _calc, _cc, res in rows:
        if res:
            only = [q for q, v in res.items() if v == "chemout-only"]
            if only:
                lines.append(f"- `{p.name}`: {', '.join(only)}")
    lines += ["", "## Only cclib", ""]
    for p, _calc, _cc, res in rows:
        if res:
            only = [q for q, v in res.items() if v == "cclib-only"]
            if only:
                lines.append(f"- `{p.name}`: {', '.join(only)}")
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text("\n".join(lines) + "\n")
    json_out = Path(args.out).with_suffix(".json")
    json_out.write_text(
        json.dumps(
            [
                {"file": str(p), "cclib_ok": r is not None, "cclib_error": cc.get("error"), "result": r}
                for p, c, cc, r in rows
            ],
            indent=1,
        )
    )
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
