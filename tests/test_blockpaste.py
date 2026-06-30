"""Regression tests for the Ma-Zheng block-pasting Stage-4 candidate generator.

Fast tests guard the library adapter and encoding; the `slow` test guards the
end-to-end count on a validated anchor (type 6 = census P2 = 49 polytopes).
Full-census block-paste coverage is run via run_blockpaste_parallel.py, not here.
"""
import json

import numpy as np
import pytest

from pipeline.utils import mazheng_lib as ml
from pipeline import stage4_blockpaste as bp


# --- fast: library adapter + encoding --------------------------------------

def test_library_sizes():
    assert len(ml.S(3)) == 31
    assert len(ml.S(4)) == 242
    assert len(ml.E(3)) == 10
    assert len(ml.E(4)) == 27
    assert len(ml.L4()) == 392
    # every library tuple is in the {2..7} alphabet (7 = wildcard for >= 7)
    for t in ml.S(4):
        assert all(2 <= m <= 7 for m in t)


def test_facet_pair_columns_order():
    assert ml.facet_pair_columns([0, 1, 3, 6]) == [(0, 1), (0, 3), (0, 6), (1, 3), (1, 6), (3, 6)]


def test_collapse7():
    assert ml.collapse7((2, 8, 12, 3, 7, 5)) == (2, 7, 7, 3, 7, 5)


def test_encode_roundtrip_membership():
    lib = frozenset({(2, 3, 4), (5, 6, 7)})
    codes = bp._lib_codes(lib, 3)
    rows = np.array([[2, 3, 4], [5, 6, 7], [2, 2, 2]], dtype=np.int8)
    mask = np.isin(bp._encode(rows), codes)
    assert mask.tolist() == [True, True, False]


def test_expand_wildcards():
    ordinary = [(0, 1), (0, 2)]
    cands = np.array([[3, 7]], dtype=np.int8)  # one wildcard -> 5 expansions
    out = list(bp.expand_label_assignments(cands, ordinary))
    assert len(out) == len(bp.HIGH_LABELS)
    assert {la[(0, 2)] for la in out} == set(bp.HIGH_LABELS)
    assert all(la[(0, 1)] == 3 for la in out)


# --- slow: end-to-end anchor count -----------------------------------------

@pytest.mark.slow
def test_type6_blockpaste_count():
    """Block-paste + screen + solve reproduces census P2 = 49 for d=4 type 6,
    with exactly the brute label-assignment count (166)."""
    from run_blockpaste_parallel import process_type
    from pathlib import Path
    import tempfile

    survivors = {t["type_id"]: t
                 for t in json.load(open("runs/d4_n8/stage2/types.json"))}
    with tempfile.TemporaryDirectory() as d:
        res = process_type(6, survivors, nproc=4, outdir=Path(d), p=1)
    assert res["label_assigns"] == 166
    assert res["distinct"] == 49
