"""Fast regression tests for the exact affine-Gale criterion (gale_exact).

Guards the 2026-06-22 generator rewrite that fixed the broken float/C criteria
(which reproduced 0/30 of the d=4 census).  Full-census coverage (30/30, 109/109)
lives in test_generator_coverage.py (marked `slow`); these are quick structural
checks on the criterion itself.
"""
from pipeline.utils.gale_exact import AffineGale, _strictly_inside

D, N = 4, 8
LBT = (D - 1) * N - (D + 1) * (D - 2)   # = 14: min vertices of a simple 4-polytope, 8 facets


def test_strictly_inside_basic():
    # A point clearly inside a square, and one outside.
    sq = [(0, 0), (10, 0), (10, 10), (0, 10)]
    assert _strictly_inside((5, 5), sq)
    assert not _strictly_inside((20, 5), sq)
    assert not _strictly_inside((0, 0), sq)        # boundary -> not strict


def test_affine_gale_known_polytope():
    """A concrete order-type + positive pair that yields a valid d=4 type
    (extracted from the validated d=4 Stage-2 output)."""
    pts = [(72, 71), (5, 62), (6, 74), (250, 188), (239, 170),
           (233, 150), (132, 93), (137, 84)]
    pos = (4, 6)
    ag = AffineGale(pts, pos, D)
    assert ag.is_polytope()
    assert ag.f0() == 14                       # vertices
    mf = ag.missing_faces(max_size=5)
    sizes = sorted(len(m) for m in mf)
    assert sizes == [2, 2, 2, 2, 2, 2, 4, 4, 4]
    assert sum(1 for m in mf if len(m) == 2) == 6     # k = 6 disjoint pairs


def test_polytopality_gate_enforces_lbt():
    """is_polytope() must reject a config whose vertex count is below the Lower
    Bound Theorem minimum (the gate that removed the old spurious high-k garbage)."""
    pts = [(72, 71), (5, 62), (6, 74), (250, 188), (239, 170),
           (233, 150), (132, 93), (137, 84)]
    ag = AffineGale(pts, (4, 6), D)
    assert ag.is_polytope()                    # genuine polytope (f0=14) passes
    # Force the vertex count below LBT -> the gate must reject.
    ag.f0 = lambda: LBT - 1
    assert not ag.is_polytope()


def test_odd_vertex_counts_reachable():
    """The old criterion could only produce EVEN f0; a correct criterion must
    reach odd vertex counts (the d=4 census has types with 15 and 17 vertices)."""
    from pipeline.stage1_order_types import stream_chi_file
    from pathlib import Path
    from itertools import combinations
    seen_odd = False
    count = 0
    for chi, pts_chunk, s, e, T in stream_chi_file(Path("data/aak/otypes08.chi"), N, 50000):
        for rec in range(min(len(pts_chunk), 400)):   # a small slice is enough
            pts = [tuple(map(int, p)) for p in pts_chunk[rec]]
            for u, v in combinations(range(N), 2):
                neg = [pts[i] for i in range(N) if i not in (u, v)]
                if not (_strictly_inside(pts[u], neg) and _strictly_inside(pts[v], neg)):
                    continue
                ag = AffineGale(pts, (u, v), D)
                if ag.is_polytope() and ag.f0() % 2 == 1:
                    seen_odd = True
                    break
            if seen_odd:
                break
        break   # only the first chunk's slice
    assert seen_odd, "no odd-f0 polytope found — criterion may have the even-only bug"
