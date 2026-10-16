"""Regression test: every in-repo corpus file must parse exactly as its golden JSON records.

Regenerate goldens deliberately with ``python scripts/make_goldens.py --update`` and review the
diff; a change here is a change in what users get.
"""

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from make_goldens import corpus_files, diff, strip_volatile  # noqa: E402

import chemout  # noqa: E402

CASES = [(src, gold) for _label, src, gold in corpus_files() if gold.exists()]


@pytest.mark.parametrize("src,gold", CASES, ids=[c[0].name for c in CASES])
def test_golden(src, gold):
    calc = chemout.parse(src)
    new = strip_volatile(json.loads(calc.to_json(indent=None)))
    old = json.loads(gold.read_text())
    d = diff(old, new)
    assert not d, "\n".join(d[:20])


def test_every_corpus_file_has_a_golden():
    missing = [src for _l, src, gold in corpus_files() if not gold.exists()]
    assert not missing, f"run scripts/make_goldens.py --update; missing: {missing}"
