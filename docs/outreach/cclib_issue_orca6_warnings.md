Title: ORCA 6: a WARNINGS block without a closing "====" line makes the parser consume the whole file

**File:** `data/regression/ORCA/ORCA6.0/1618.out` (already in cclib-data; a UHF geometry optimisation).

**What happens** (master):

```python
>>> d = ccread("data/regression/ORCA/ORCA6.0/1618.out")
>>> sorted(k for k in vars(d) if not k.startswith("_"))
['charge', 'metadata']
```

No exception, no warning: `atomcoords`, `scfenergies`, `mult`, `natom` … are all missing, and `metadata["success"]` is `True`.

**Cause:** ORCA ≤ 5 closes the WARNINGS section with a line of `=`:

```
================================================================================
                                        WARNINGS
                       Please study these warnings very carefully!
================================================================================

WARNING: ...
================================================================================
```

ORCA 6 prints no closing line; after the warnings come blank lines and the next banner (`* Geometry Optimization Run *` here, lines 241–262 of `1618.out`). The loop in `orcaparser.py`

```python
line = next(inputfile)
while line[0] != "=":
    ...
```

therefore runs to the end of the file, `StopIteration` ends `parse()` quietly, and the record is empty.

**Fix:** also stop at the next banner (`*`) or the `INPUT FILE` header. PR follows; with it `1618.out` yields 6 geometries, 6 SCF energies, `mult == 2`, and `metadata["warnings"]` holds the two warnings. Every ORCA 6 output that prints a WARNINGS block (open-shell switch, AutoStart, basis-set notices, …) is affected.
