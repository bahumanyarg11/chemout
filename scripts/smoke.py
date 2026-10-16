"""Parse every file under the given dirs; print a one-line summary per file and crash/warning stats."""

import sys
import time
import traceback
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import chemout

exts = {".out", ".log", ".txt"}
files = [p for d in sys.argv[1:] for p in sorted(Path(d).rglob("*")) if p.suffix in exts and p.is_file()]
crashes, warn_sections = 0, Counter()
t0 = time.time()
for p in files:
    try:
        c = chemout.parse(p)
    except Exception:
        crashes += 1
        print(f"CRASH {p}\n{traceback.format_exc()}")
        continue
    w = c.provenance.warnings
    for x in w:
        warn_sections[x.section] += 1
    e = c.final_energy
    print(
        f"{c.program.name.value:8s} {c.program.version or '?':14s} {c.termination.value:10s} geo={len(c.geometries):3d} scf={len(c.scf_cycles):3d} "
        f"E={f'{e:.6f}' if e is not None else '   None   ':>14s} pop={len(c.populations)} orb={'y' if c.orbitals else '-'} "
        f"vib={len(c.vibrations.frequencies) if c.vibrations else 0:3d} exc={len(c.excited_states):3d} warn={len(w):2d} {p.name}"
        + (f"  | {w[0].section}: {w[0].message[:90]}" if w else "")
    )
print(
    f"\n{len(files)} files, {crashes} crashes, {time.time() - t0:.1f}s; warnings by section: {dict(warn_sections.most_common())}"
)
