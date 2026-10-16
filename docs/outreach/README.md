# Upstream contributions from the conformance work

| What | Where | Status (2026-10-16) |
|---|---|---|
| xTB parser: `homos` always orbital 0, `moenergies` 27.2× too large | [cclib/cclib#1934](https://github.com/cclib/cclib/issues/1934) → PR [#1936](https://github.com/cclib/cclib/pull/1936) | open |
| ORCA 6 WARNINGS block without closing `====` empties the record | [cclib/cclib#1935](https://github.com/cclib/cclib/issues/1935) → PR [#1937](https://github.com/cclib/cclib/pull/1937) | open |
| 1.8.1 release cannot read ORCA 6 (109/115 files): data point on the release issue | [cclib/cclib#1832 comment](https://github.com/cclib/cclib/issues/1832#issuecomment-6209845536) | posted |

Retracted: the earlier "cclib master reports the last SCF-iteration energy" finding was a unit-constant difference (cclib uses the CODATA 2010 Hartree→eV value); see CONVENTIONS.md.
