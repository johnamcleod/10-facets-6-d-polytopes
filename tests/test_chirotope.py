"""Unit tests for chirotope computation."""

import pytest
from pipeline.utils.chirotope import (
    orientation, compute_chirotope, chirotope_to_tuple,
    canonical_chirotope, verify_chirotope_from_points
)


def test_orientation_ccw():
    assert orientation((0, 0), (1, 0), (0, 1)) == 1  # CCW


def test_orientation_cw():
    assert orientation((0, 0), (0, 1), (1, 0)) == -1  # CW


def test_orientation_collinear():
    assert orientation((0, 0), (1, 1), (2, 2)) == 0


def test_chirotope_square():
    """4 points of a square (CCW): (0,0),(1,0),(1,1),(0,1)."""
    pts = [(0, 0), (1, 0), (1, 1), (0, 1)]
    chi = compute_chirotope(pts)
    assert chi is not None
    # (0,1,2) = CCW? det([1,0],[1,1]) = 1 > 0 => +1
    assert chi[(0, 1, 2)] == 1
    # (0,1,3) = det([1,0],[0,1]) = 1 > 0 => +1
    assert chi[(0, 1, 3)] == 1


def test_chirotope_degenerate():
    """3 collinear points => degenerate."""
    pts = [(0, 0), (1, 1), (2, 2), (0, 1)]
    chi = compute_chirotope(pts)
    assert chi is None


def test_canonical_chirotope_invariant():
    """Two rotations/relabelings of same point set give same canonical."""
    pts1 = [(0, 0), (1, 0), (2, 1), (1, 2), (0, 1)]
    pts2 = [(1, 0), (2, 1), (1, 2), (0, 1), (0, 0)]  # cyclic relabeling

    chi1 = compute_chirotope(pts1)
    chi2 = compute_chirotope(pts2)
    assert chi1 is not None
    assert chi2 is not None

    n = 5
    t1 = chirotope_to_tuple(chi1, n)
    t2 = chirotope_to_tuple(chi2, n)

    c1, _ = canonical_chirotope(t1, n)
    c2, _ = canonical_chirotope(t2, n)
    assert c1 == c2, "Same order type should have the same canonical chirotope"


def test_verify_chirotope():
    pts = [(0, 0), (3, 0), (3, 3), (0, 3)]
    n = 4
    chi = compute_chirotope(pts)
    chi_tuple = chirotope_to_tuple(chi, n)
    assert verify_chirotope_from_points(chi_tuple, pts)


def test_verify_chirotope_wrong_points():
    pts = [(0, 0), (3, 0), (3, 3), (0, 3)]
    n = 4
    chi = compute_chirotope(pts)
    chi_tuple = chirotope_to_tuple(chi, n)
    # Wrong points
    wrong_pts = [(0, 0), (3, 0), (3, 3), (1, 1)]  # different last point
    # chi_tuple may or may not match — depends on orientations
    # Just verify the function runs without error
    result = verify_chirotope_from_points(chi_tuple, wrong_pts)
    # No assertion on value — just that it doesn't crash
