Comment on cclib/cclib#1832 (Release v1.9):

A data point in favour of releasing soon: the latest PyPI release (1.8.1) cannot read ORCA 6 at all, while master can.

I ran cclib 1.8.1 and master over 246 output files — the cclib-data ORCA 4.2–6.1, Gaussian 09/16 and xTB sets plus 109 public ORCA 6.0.0–6.1.1 outputs collected from GitHub repositories. Results:

| | cclib 1.8.1 | cclib master (76df867) |
|---|---|---|
| files on which `ccread` raised or returned `None` | 119 / 246 | 12 / 246 |
| of which ORCA 6 | 109 / 115 (96 `IndexError: list index out of range`, 11 `AssertionError`) | 3 / 115 |

So anyone who `pip install cclib` today and feeds it an ORCA 6 file gets a traceback in ~95% of cases. The 12 remaining master failures are a separate matter (two I have filed with fixes: #xtb-homos, #orca6-warnings; one is #1728).

Happy to share the per-file table if useful.
