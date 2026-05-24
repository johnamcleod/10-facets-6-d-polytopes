"""
Chirotope (order type) computation and canonicalization for planar point sets.

A chirotope χ: C(n,3) → {-1, 0, +1} records the signed area (orientation)
of every triple of points. For points in general position (no 3 collinear),
χ takes values in {-1, +1} only — these are the "uniform" chirotopes
corresponding to realizable order types.

Convention: χ(i,j,k) = sign(det([p_j-p_i, p_k-p_i]))
  = +1 if (p_i, p_j, p_k) are in counter-clockwise orientation
  = -1 if clockwise
"""

from fractions import Fraction
from itertools import combinations, permutations
from functools import lru_cache
import math


def orientation(p0, p1, p2):
    """Exact integer orientation test.

    Returns +1 (CCW), -1 (CW), or 0 (collinear).
    Inputs must be integer or Fraction coordinates.
    """
    # det([p1-p0, p2-p0]) = (p1.x-p0.x)*(p2.y-p0.y) - (p1.y-p0.y)*(p2.x-p0.x)
    ax, ay = p1[0] - p0[0], p1[1] - p0[1]
    bx, by = p2[0] - p0[0], p2[1] - p0[1]
    det = ax * by - ay * bx
    if det > 0:
        return 1
    elif det < 0:
        return -1
    else:
        return 0


def compute_chirotope(points):
    """Compute the chirotope of a point configuration.

    Args:
        points: list of (x, y) pairs (integer or Fraction coordinates).

    Returns:
        dict mapping sorted triple (i,j,k) with i<j<k to orientation ±1.
        Returns None if any triple is collinear (degenerate configuration).
    """
    n = len(points)
    chi = {}
    for i, j, k in combinations(range(n), 3):
        o = orientation(points[i], points[j], points[k])
        if o == 0:
            return None  # degenerate
        chi[(i, j, k)] = o
    return chi


def chirotope_to_tuple(chi, n):
    """Flatten chirotope to a tuple in lexicographic triple order."""
    return tuple(chi[(i, j, k)] for i, j, k in combinations(range(n), 3))


def relabel_chirotope(chi_tuple, n, perm):
    """Apply a relabeling permutation to a chirotope tuple.

    perm[i] = new label of old label i.
    Returns the chirotope tuple under the relabeled point indices
    (reindexed in sorted canonical form).
    """
    # Build map: new sorted triple -> sign (adjusted for perm parity)
    new_chi = {}
    for idx, (i, j, k) in enumerate(combinations(range(n), 3)):
        pi, pj, pk = perm[i], perm[j], perm[k]
        sign = chi_tuple[idx]
        # Account for the sign flip when sorting the permuted triple
        triple = [pi, pj, pk]
        inv = _inversion_parity(triple)
        triple_sorted = tuple(sorted(triple))
        new_chi[triple_sorted] = sign * inv
    return tuple(new_chi[(i, j, k)] for i, j, k in combinations(range(n), 3))


def _inversion_parity(lst):
    """Return +1 or -1 for the parity of the permutation sorting lst."""
    lst = list(lst)
    n = len(lst)
    parity = 1
    for i in range(n):
        while lst[i] != sorted(lst)[i]:
            j = lst.index(sorted(lst)[i], i)
            lst[i], lst[j] = lst[j], lst[i]
            parity *= -1
    return parity


def canonical_chirotope(chi_tuple, n):
    """Compute the lexicographically minimal relabeling of a chirotope.

    This is the canonical form for the order type equivalence class.
    Returns (canonical_tuple, representative_permutation).
    """
    best = None
    best_perm = None
    # Iterate over all n! permutations — feasible for n ≤ 10
    from itertools import permutations as perms
    for perm in perms(range(n)):
        relabeled = relabel_chirotope(chi_tuple, n, perm)
        if best is None or relabeled < best:
            best = relabeled
            best_perm = perm
    return best, best_perm


def points_to_canonical_chirotope(points):
    """Full pipeline: points -> canonical chirotope tuple.

    Returns (canonical_chi_tuple, n) or None if degenerate.
    """
    chi = compute_chirotope(points)
    if chi is None:
        return None
    n = len(points)
    chi_tuple = chirotope_to_tuple(chi, n)
    canonical, _ = canonical_chirotope(chi_tuple, n)
    return canonical, n


def verify_chirotope_from_points(chi_tuple, points):
    """Verify that chi_tuple matches the chirotope of points."""
    n = len(points)
    chi = compute_chirotope(points)
    if chi is None:
        return False
    expected = chirotope_to_tuple(chi, n)
    return expected == chi_tuple
