"""Tests for the tolerance-free forward-checking gates.

These two predicates decide all but one of the searched types on their own, so a
wrong verdict from either would silently lose polytopes.  The tests check them
three ways: against the floating-point tests they replace, against a direct
determinant computation, and against Lanner's published classification.
"""
import itertools
import math

import numpy as np
import pytest

from pipeline.utils import exact_gates as G

ALPHABET = [2, 3, 4, 5, 6, 7]     # 7 is the wildcard proxy


def gram(k, labs):
    M = np.eye(k)
    for (a, b), m in zip(G._pairs(k), labs):
        M[a, b] = M[b, a] = -math.cos(math.pi / m)
    return M


@pytest.mark.parametrize("k", [3, 4])
def test_elliptic_agrees_with_float_exhaustively(k):
    """Ellipticity is combinatorial; it must still agree with the eigenvalue test
    on every label tuple where that test is reliable."""
    n = 0
    for labs in itertools.product(ALPHABET, repeat=len(G._pairs(k))):
        ev = np.linalg.eigvalsh(gram(k, labs))
        assert G.elliptic_core(k, labs) == bool(np.all(ev > 1e-8)), labs
        n += 1
    assert n == len(ALPHABET) ** len(G._pairs(k))


@pytest.mark.parametrize("k", [3, 4])
def test_exact_lanner_implies_one_negative_eigenvalue(k):
    """Lanner membership is strictly stronger than the necessary condition a float
    gate would test, so it must imply it everywhere -- that is what makes swapping
    the tests sound rather than merely different."""
    for labs in itertools.product(ALPHABET, repeat=len(G._pairs(k))):
        if G.lanner_core(k, labs):
            ev = np.linalg.eigvalsh(gram(k, labs))
            assert int((ev < -1e-8).sum()) == 1, labs
            assert int((abs(ev) <= 1e-8).sum()) == 0, labs


def test_triangle_invariant_is_the_determinant_sign():
    """On three nodes the gate uses 1/p+1/q+1/r-1 in integer arithmetic in place of
    a determinant.  Checked against the determinant over a wide label range,
    including labels far beyond the wildcard."""
    for labs in itertools.product([2, 3, 4, 5, 6, 7, 8, 10, 20, 100, 1000],
                                  repeat=3):
        det = float(np.linalg.det(gram(3, labs)))
        s = sum(1 / m for m in labs) - 1
        assert (det < -1e-11) == (s < -1e-11), (labs, det, s)


def test_wildcard_verdict_is_uniform_on_three_nodes():
    """Evaluating a wild edge at m = 7 must give the verdict it gives for every
    m >= 7, or the substitution would lose configurations."""
    for p in range(2, 8):
        for q in range(2, 8):
            v7 = G.lanner_core(3, (p, q, 7))
            for r in (8, 9, 12, 50, 1000, 10 ** 6):
                assert G.lanner_core(3, (p, q, r)) == v7, (p, q, r)


def test_no_high_label_in_lanner_of_order_at_least_four():
    """Lemma 4.12(c): a Lanner diagram of order >= 4 has no edge of label >= 6.
    The gate relies on this to reject wild edges outright, which is also what keeps
    cos(pi/7) -- cubic over Q, hence outside the field the gate computes in -- from
    ever arising.  Verified from the definition at 4 nodes."""
    for labs in itertools.product([2, 3, 4, 5, 6, 7], repeat=6):
        if max(labs) >= 6:
            assert not G.lanner_core(4, labs), labs


def test_recovers_lanners_classification():
    """Counting the diagrams the exact test calls Lanner must reproduce Lanner's
    lists: 9 of order 4 and 5 of order 5, up to relabelling."""
    for k, expected in ((4, 9), (5, 5)):
        seen = set()
        for labs in itertools.product([2, 3, 4, 5], repeat=len(G._pairs(k))):
            if not G.lanner_core(k, labs):
                continue
            lab = {p: m for p, m in zip(G._pairs(k), labs)}
            best = None
            for perm in itertools.permutations(range(k)):
                img = tuple(lab[tuple(sorted((perm[a], perm[b])))]
                            for a, b in G._pairs(k))
                best = img if best is None or img < best else best
            seen.add(best)
        assert len(seen) == expected, (k, len(seen))


def test_gate_shapes_match_the_core():
    """The wrappers the enumerator calls must agree with the memoized cores."""
    nodes = (2, 5, 7, 9)
    assignment = {(2, 5): 3, (2, 7): 2, (2, 9): 4, (5, 7): 5, (5, 9): 2, (7, 9): 3}
    labs = G._labs_from_assignment(nodes, assignment)
    assert G.vertex_pd_exact(nodes, assignment) == G.elliptic_core(4, labs)
    assert G.lanner_exact(nodes, assignment, 4) == G.lanner_core(4, labs)
    # a pair absent from the assignment is orthogonal, i.e. label 2
    assert G._labs_from_assignment((0, 1), {}) == (2,)
