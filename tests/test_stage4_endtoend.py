"""End-to-end Stage-4 gate (marked `slow`, ~8min): the FULL process_type_stage4
(label enumeration + field-agnostic solver) must recover P^B6 from its
combinatorial type alone.  Run with:  pytest -m slow tests/test_stage4_endtoend.py
"""
import math
from itertools import combinations

import pytest

np = pytest.importorskip("numpy")
from pipeline.stage4_gram import process_type_stage4

pytestmark = pytest.mark.slow


def _pb6_type():
    """Reconstruct P^B6's combinatorial type (missing faces + vertices) from its
    known Gram matrix, so Stage 4 must rediscover the labels and weights."""
    n, d = 10, 6
    ordin = {(0, 1): 5, (1, 2): 3, (2, 3): 3, (3, 4): 3, (4, 5): 3,
             (8, 9): 5, (1, 6): 4, (1, 7): 4, (4, 8): 3, (5, 9): 5}
    x56 = 2 * math.sqrt(2) + math.sqrt(10)
    x67 = 17 + 8 * math.sqrt(5)
    G = np.eye(n)
    for (i, j), m in ordin.items():
        G[i, j] = G[j, i] = -math.cos(math.pi / m)
    for (i, j), x in {(5, 6): x56, (6, 7): x67, (7, 8): x56}.items():
        G[i, j] = G[j, i] = -x

    def is_pd(S):
        return bool(np.all(np.linalg.eigvalsh(G[np.ix_(S, S)]) > 1e-9))

    nf, nf_set = [], set()
    for size in range(2, 6):
        for S in combinations(range(n), size):
            if any(frozenset(sub) in nf_set
                   for ss in range(2, size) for sub in combinations(S, ss)):
                continue
            if not is_pd(list(S)):
                nf.append(frozenset(S)); nf_set.add(frozenset(S))
    return {
        "type_id": 379,
        "missing_faces": [sorted(m) for m in nf],
        "vertex_sets": [sorted(S) for S in combinations(range(n), d) if is_pd(list(S))],
        "p_count": sum(1 for m in nf if len(m) == 2),
    }


def test_pb6_recovered_end_to_end():
    t = _pb6_type()
    assert t["p_count"] == 3
    res = process_type_stage4(t, 6, max_assignments=2_000_000,
                              enum_timeout=600.0, solve_timeout=300.0,
                              label_indices=(0, 1, 2, 3), verbose=False)
    assert len(res) >= 1, "Stage 4 failed to recover P^B6"
    # Every recovered config must carry P^B6's exact dotted weights.
    for r in res:
        polys = {kk: tuple(v["minpoly"]) for kk, v in r["dot_values"].items()}
        vals = sorted(tuple(p) for p in polys.values())
        assert (1, -2, -1, 1) not in vals          # not the wrong field
        assert (4, 0, -36, 0, 1) in vals           # 2√2+√10 appears
        assert (-31, -34, 1) in vals               # 17+8√5 appears
