"""Fill the <!-- CONFORMANCE_SUMMARY --> block in README.md from docs/CONFORMANCE_*.json."""

from __future__ import annotations

import json
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
README = ROOT / "README.md"
REPORTS = [
    ("cclib 1.8.1 (latest release)", "docs/CONFORMANCE_cclib_1.8.1.json", "docs/CONFORMANCE_cclib_1.8.1.md"),
    ("cclib master (2.0 dev)", "docs/CONFORMANCE_cclib_master.json", "docs/CONFORMANCE_cclib_master.md"),
]
ORDER = [
    "scf_energy",
    "n_geometries",
    "final_coordinates",
    "charge",
    "multiplicity",
    "homo_energy",
    "mulliken_charges",
    "vibfreqs",
    "etenergies",
    "success_flag",
]


def summarise(path: Path) -> tuple[int, int, dict[str, Counter]]:
    rows = json.loads(path.read_text())
    n = len(rows)
    failed = sum(1 for r in rows if not r["cclib_ok"])
    per: dict[str, Counter] = defaultdict(Counter)
    for r in rows:
        if r["result"]:
            for q, v in r["result"].items():
                per[q][v] += 1
    return n, failed, per


def main() -> None:
    lines = []
    head = ["| | " + " | ".join(label for label, _, _ in REPORTS) + " |", "|---|" + "---|" * len(REPORTS)]
    stats = [summarise(ROOT / js) for _, js, _ in REPORTS]
    head.append("| files tried | " + " | ".join(str(n) for n, _, _ in stats) + " |")
    head.append("| **cclib raised / returned nothing** | " + " | ".join(f"**{f}**" for _, f, _ in stats) + " |")
    head.append("| chemout returned a record | " + " | ".join(str(n) for n, _, _ in stats) + " |")
    lines += head
    lines.append("")
    lines.append(
        "Where both produced a record, per quantity (agree / differ by < 1e-4 Eh / disagree / only one side has it):"
    )
    lines.append("")
    lines.append("| quantity | " + " | ".join(label for label, _, _ in REPORTS) + " |")
    lines.append("|---|" + "---|" * len(REPORTS))
    for q in ORDER:
        cells = []
        for _, _, per in stats:
            c = per.get(q, Counter())
            only = c.get("chemout-only", 0) + c.get("cclib-only", 0)
            cells.append(f"{c.get('agree', 0)} / {c.get('disagree-small', 0)} / {c.get('disagree', 0)} / {only}")
        lines.append(f"| {q} | " + " | ".join(cells) + " |")
    lines.append("")
    lines.append(
        "Corpus: cclib-data ORCA 4.2–6.1, Gaussian 09/16, xTB; 109 public ORCA 6.0–6.1 outputs from GitHub; 9 xTB runs generated "
        "here. Every disagreement is listed by file in the linked reports, with the convention behind it where one exists "
        "(e.g. cclib master reports the last SCF-iteration energy rather than ORCA's converged `Total Energy`; cclib keeps "
        "xTB's six zero-frequency modes; chemout drops Gaussian dummy atoms but keeps ghost atoms)."
    )
    block = "<!-- CONFORMANCE_SUMMARY -->\n" + "\n".join(lines) + "\n<!-- /CONFORMANCE_SUMMARY -->"
    text = README.read_text()
    if "<!-- /CONFORMANCE_SUMMARY -->" in text:
        text = re.sub(r"<!-- CONFORMANCE_SUMMARY -->[\s\S]*?<!-- /CONFORMANCE_SUMMARY -->", block, text)
    else:
        text = text.replace("<!-- CONFORMANCE_SUMMARY -->", block)
    README.write_text(text)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
