"""Generate or check golden JSON files for the regression corpus.

Usage:
    python scripts/make_goldens.py --update            # (re)write tests/golden/**.json
    python scripts/make_goldens.py --check             # compare and exit 1 on any difference

Corpus roots are read from tests/corpus.txt (one "label<TAB>path" per line); paths may be
absolute (external corpora such as cclib-data) or relative to the repo. Each golden file stores
the complete Calculation JSON minus volatile provenance fields.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import chemout  # noqa: E402

GOLDEN = ROOT / "tests" / "golden"
CORPUS = ROOT / "tests" / "corpus.txt"
EXTS = {".out", ".log", ".txt"}
ATOL = 1e-6


def corpus_files() -> list[tuple[str, Path, Path]]:
    out = []
    for line in CORPUS.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        label, _, path = line.partition("\t")
        root = Path(path) if Path(path).is_absolute() else ROOT / path
        if not root.exists():
            print(f"skip missing corpus root {root}", file=sys.stderr)
            continue
        files = [root] if root.is_file() else sorted(p for p in root.rglob("*") if p.suffix in EXTS and p.is_file())
        for f in files:
            rel = f.relative_to(root.parent if root.is_file() else root)
            out.append((label, f, GOLDEN / label / rel.with_suffix(rel.suffix + ".json")))
    return out


def strip_volatile(d: dict) -> dict:
    prov = d.get("provenance", {})
    for k in ("path", "parsed_at", "sha256", "size_bytes"):
        prov.pop(k, None)
    return d


def diff(a, b, path="") -> list[str]:
    out: list[str] = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if k not in a or k not in b:
                out.append(f"{path}/{k}: only in {'golden' if k in a else 'new'}")
            else:
                out.extend(diff(a[k], b[k], f"{path}/{k}"))
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            out.append(f"{path}: length {len(a)} -> {len(b)}")
        else:
            for i, (x, y) in enumerate(zip(a, b)):
                out.extend(diff(x, y, f"{path}[{i}]"))
                if len(out) > 20:
                    break
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool):
        if not (math.isclose(a, b, rel_tol=1e-9, abs_tol=ATOL) or (math.isnan(a) and math.isnan(b))):
            out.append(f"{path}: {a} -> {b}")
    elif a != b:
        out.append(f"{path}: {a!r} -> {b!r}")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--update", action="store_true")
    g.add_argument("--check", action="store_true")
    ap.add_argument("--only", help="substring filter on file path")
    args = ap.parse_args()

    failures = 0
    n = 0
    for _label, src, gold in corpus_files():
        if args.only and args.only not in str(src):
            continue
        n += 1
        calc = chemout.parse(src)
        new = strip_volatile(json.loads(calc.to_json(indent=None)))
        if args.update:
            gold.parent.mkdir(parents=True, exist_ok=True)
            gold.write_text(json.dumps(new, indent=1, sort_keys=True) + "\n")
            continue
        if not gold.exists():
            print(f"MISSING golden for {src}")
            failures += 1
            continue
        old = json.loads(gold.read_text())
        d = diff(old, new)
        if d:
            failures += 1
            print(f"DIFF {src.name} ({len(d)} changes)")
            for line in d[:12]:
                print("   ", line)
    print(f"{n} files, {failures} failures" if args.check else f"{n} goldens written")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
