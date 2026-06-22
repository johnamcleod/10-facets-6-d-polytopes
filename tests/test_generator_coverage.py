"""Full-census coverage gates for the Stage 1/2 generator (marked `slow`).

The generator must reproduce the published combinatorial census exactly.  These
re-run Stage 2 over the AAK order-type database and compare to the Ma–Zheng
ground truth.  Slow (d=4 ~40s, d=5 ~30-50min); run with:  pytest -m slow

The oracle files (data/ground_truth/{4d8m,5d9m}.txt) are the full simple
(d,d+4)-polytope candidate sets; the generator gate is Stage-2 coverage, and for
d=5 the comparison is against the k>=2 (compact-relevant) subset (our p>=2
pre-filter correctly excludes the k<2 types, which are non-hyperbolic by
Felikson–Tumarkin).
"""
from pathlib import Path
from collections import defaultdict

import pytest

from pipeline.stage2_gale import run_stage2_from_chi
from pipeline.validate_coverage import load_ground_truth, type_key, k_of

pytestmark = pytest.mark.slow

CHI = {4: "data/aak/otypes08.chi", 5: "data/aak/otypes09.chi"}
TRUTH = {4: "data/ground_truth/4d8m.txt", 5: "data/ground_truth/5d9m.txt"}


def _coverage(d, tmp_path):
    n = d + 4
    types = run_stage2_from_chi(d, Path(CHI[d]), output_dir=str(tmp_path / "stage2"),
                                exact_dedup=True, verbose=False)
    ours = {type_key([frozenset(m) for m in t["missing_faces"]], n) for t in types}
    truth = load_ground_truth(TRUTH[d], n)
    tk = defaultdict(list)
    for t in truth:
        tk[type_key(t["mnf"], n)].append(t)
    truth_k2 = {key for key in tk if k_of(tk[key][0]["mnf"]) >= 2}
    covered = truth_k2 & ours
    spurious = ours - set(tk)
    return covered, truth_k2, spurious


@pytest.mark.skipif(not Path(CHI[4]).exists(), reason="AAK n=8 db missing")
def test_d4_generator_30_of_30(tmp_path):
    covered, truth_k2, spurious = _coverage(4, tmp_path)
    assert len(truth_k2) == 30
    assert len(covered) == 30, f"covered {len(covered)}/30"
    assert len(spurious) == 0, f"{len(spurious)} spurious types"


@pytest.mark.skipif(not Path(CHI[5]).exists(), reason="AAK n=9 db missing")
def test_d5_generator_109_of_109(tmp_path):
    covered, truth_k2, spurious = _coverage(5, tmp_path)
    assert len(truth_k2) == 109
    assert len(covered) == 109, f"covered {len(covered)}/109"
    assert len(spurious) == 0, f"{len(spurious)} spurious types"
