Title: xTB parser: homos is always orbital 0 and moenergies are 27.2× too large

**Files:** `data/regression/XTB/basicXTB6.7.1/1768_xyz.out`, `1768_sdf.out`, `XTB/basicXTB6.5.1/1448.out` (all already in cclib-data).

**What happens** (master, `cclib.__version__` 2.0.0.dev):

```python
>>> from cclib.io import ccread
>>> d = ccread("data/regression/XTB/basicXTB6.7.1/1768_xyz.out")
>>> d.homos
array([ 0, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33])
>>> d.moenergies[0][33]
-308.079138120585
```

The file says orbital 34 (1-based) is the HOMO at −0.4160635 Eh = −11.32 eV:

```
        34        2.0000           -0.4124250             -11.2227 (HOMO)
        35                         -0.2928383              -7.9685 (LUMO)
```

**Cause** (`cclib/parser/xtbparser.py`, `_extract_orbitals`):

1. `bool("(HOMO) in line")` is `bool` of a non-empty string literal, so it is always `True`; every occupied orbital is appended to `homos`, and `homos[0]` is orbital 0.
2. The function returns `float(line_split[3])` (and `line_split[2]` for virtuals), which is the **eV** column. The parser then converts Hartree → eV once more, so `moenergies` are 27.21 times too large.

**Fix:** check `"(HOMO)" in line` and return the Hartree column. PR with a regression test on `1768_xyz.out` follows. Expected after the fix: `homos == [33]`, `moenergies[0][33] == -11.3217 eV`, `moenergies[0][34] == -7.7170 eV`.
