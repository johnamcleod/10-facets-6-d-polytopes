"""Affine Gale diagram: face lattice and missing-face computation.

Setup (d=6, n=10, k=2):
  n signed points (a_i, sigma_i) in R^2, sigma_i in {+1,-1}.
  Exactly 2 positive points lie in interior of conv(negative points).
  Encodes a simple d-polytope with n facets.

Face criterion (Ziegler, Lectures on Polytopes Ch. 6 / Burcroff Sec 3):
  A subset F of facets is a FACE of P iff in T = [n]\\F every positive point
  of T lies in conv(negative points of T), i.e. conv(pos(T)) subset conv(neg(T)).

A MISSING FACE is a minimal subset M subset [n] that is NOT a face.
"""

from itertools import combinations
from fractions import Fraction

try:
    import numpy as np
    from scipy.spatial import ConvexHull, Delaunay
    _SCIPY = True
except ImportError:
    _SCIPY = False


# ---------------------------------------------------------------------------
# Exact convex-hull / point-in-convex-hull tests in R^2
# ---------------------------------------------------------------------------

def _sign(x):
    if x > 0:
        return 1
    elif x < 0:
        return -1
    return 0


def _cross2(o, a, b):
    """2D cross product of vectors (a-o) and (b-o). Exact."""
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def point_in_convex_hull_2d(pt, hull_pts, fast=True):
    """Test whether pt is in conv(hull_pts) in R^2.

    When fast=True (default), uses scipy/numpy for speed. This is reliable
    for generic integer coordinates but may mis-classify boundary cases.
    Set fast=False for exact arithmetic (slower).
    """
    m = len(hull_pts)
    if m == 0:
        return False
    if m == 1:
        return pt == hull_pts[0]
    if m == 2:
        return _on_segment(pt, hull_pts[0], hull_pts[1])

    if fast and _SCIPY:
        return _point_in_hull_numpy(pt, hull_pts)

    # Exact path
    hull = _convex_hull_2d(hull_pts)
    if hull is None:
        return _point_on_segment_set(pt, hull_pts)
    h = len(hull)
    if h == 1:
        return pt == hull[0]
    if h == 2:
        return _on_segment(pt, hull[0], hull[1])
    for i in range(h):
        a = hull[i]
        b = hull[(i + 1) % h]
        if _cross2(a, b, pt) < 0:
            return False
    return True


def _point_in_hull_numpy(pt, hull_pts):
    """Fast point-in-convex-hull for 2D using numpy cross products.

    Uses Andrew's monotone chain (pure numpy) then cross-product test.
    Reliable for generic integer/float coordinates with no near-degeneracy.
    """
    pts = np.array(hull_pts, dtype=float)
    p = np.array(pt, dtype=float)
    m = len(pts)
    if m == 1:
        return np.allclose(pts[0], p)
    if m == 2:
        # On segment?
        ab = pts[1] - pts[0]
        ap = p - pts[0]
        cross = ab[0] * ap[1] - ab[1] * ap[0]
        if abs(cross) > 1e-9:
            return False
        t = np.dot(ap, ab) / np.dot(ab, ab)
        return 0 <= t <= 1

    # Monotone-chain convex hull
    idx = np.lexsort((pts[:, 1], pts[:, 0]))
    pts_sorted = pts[idx]

    def cross_np(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower = []
    for q in pts_sorted:
        while len(lower) >= 2 and cross_np(lower[-2], lower[-1], q) <= 0:
            lower.pop()
        lower.append(q)
    upper = []
    for q in reversed(pts_sorted):
        while len(upper) >= 2 and cross_np(upper[-2], upper[-1], q) <= 0:
            upper.pop()
        upper.append(q)
    hull = np.array(lower[:-1] + upper[:-1])

    if len(hull) <= 2:
        # Collinear hull — use segment check
        return bool(np.all(
            np.all((pts_sorted[0] <= p) & (p <= pts_sorted[-1])) or
            np.all((pts_sorted[-1] <= p) & (p <= pts_sorted[0]))
        ))

    h = len(hull)
    for i in range(h):
        a, b = hull[i], hull[(i + 1) % h]
        if cross_np(a, b, p) < -1e-9:
            return False
    return True


def _point_strictly_in_hull_numpy(pt, hull_pts):
    """Strict interior test (same as _point_in_hull_numpy but strict inequalities)."""
    pts = np.array(hull_pts, dtype=float)
    p = np.array(pt, dtype=float)
    m = len(pts)
    if m <= 2:
        return False

    def cross_np(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    idx = np.lexsort((pts[:, 1], pts[:, 0]))
    pts_sorted = pts[idx]
    lower, upper = [], []
    for q in pts_sorted:
        while len(lower) >= 2 and cross_np(lower[-2], lower[-1], q) <= 0:
            lower.pop()
        lower.append(q)
    for q in reversed(pts_sorted):
        while len(upper) >= 2 and cross_np(upper[-2], upper[-1], q) <= 0:
            upper.pop()
        upper.append(q)
    hull = np.array(lower[:-1] + upper[:-1])
    if len(hull) <= 2:
        return False
    h = len(hull)
    for i in range(h):
        a, b = hull[i], hull[(i + 1) % h]
        if cross_np(a, b, p) <= 1e-9:
            return False
    return True


def _on_segment(pt, a, b):
    """Is pt on segment [a,b]?"""
    # cross product must be 0 and pt within bounding box
    cross = _cross2(a, b, pt)
    if cross != 0:
        return False
    # Check bounding box
    if not (min(a[0], b[0]) <= pt[0] <= max(a[0], b[0])):
        return False
    if not (min(a[1], b[1]) <= pt[1] <= max(a[1], b[1])):
        return False
    return True


def _point_on_segment_set(pt, pts):
    """Is pt on any segment between consecutive pts (all collinear)?"""
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            if _on_segment(pt, pts[i], pts[j]):
                return True
    return False


def _convex_hull_2d(pts):
    """Compute convex hull of pts in CCW order using Andrew's monotone chain.

    Returns list of hull points in CCW order, or None if all collinear.
    Exact arithmetic (integer/Fraction).
    """
    pts = sorted(set(map(tuple, pts)))
    n = len(pts)
    if n < 2:
        return pts

    lower = []
    for p in pts:
        while len(lower) >= 2 and _cross2(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)

    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and _cross2(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)

    hull = lower[:-1] + upper[:-1]

    # Check if all points are collinear
    if len(hull) <= 2:
        return None  # collinear
    return hull


def _segments_cross_2d(a, b, c, d):
    """True iff segments [a,b] and [c,d] intersect (including endpoints)."""
    d1 = _cross2(c, d, a)
    d2 = _cross2(c, d, b)
    d3 = _cross2(a, b, c)
    d4 = _cross2(a, b, d)
    if _sign(d1) * _sign(d2) < 0 and _sign(d3) * _sign(d4) < 0:
        return True
    if d1 == 0 and _on_segment(a, c, d):
        return True
    if d2 == 0 and _on_segment(b, c, d):
        return True
    if d3 == 0 and _on_segment(c, a, b):
        return True
    if d4 == 0 and _on_segment(d, a, b):
        return True
    return False


def _segment_meets_conv_hull_2d(a, b, hull_pts):
    """True iff segment [a,b] intersects conv(hull_pts) in R^2 (exact arithmetic)."""
    if point_in_convex_hull_2d(a, hull_pts, fast=False):
        return True
    if point_in_convex_hull_2d(b, hull_pts, fast=False):
        return True
    hull = _convex_hull_2d(hull_pts)
    if hull is None or len(hull) <= 2:
        return False
    h = len(hull)
    return any(
        _segments_cross_2d(a, b, hull[i], hull[(i + 1) % h])
        for i in range(h)
    )


def point_strictly_in_convex_hull_2d(pt, hull_pts, fast=True):
    """Test whether pt is in the INTERIOR of conv(hull_pts) in R^2."""
    m = len(hull_pts)
    if m == 0:
        return False
    if m <= 2:
        return False

    if fast and _SCIPY:
        return _point_strictly_in_hull_numpy(pt, hull_pts)

    # Exact fallback
    hull = _convex_hull_2d(hull_pts)
    if hull is None or len(hull) <= 2:
        return False
    h = len(hull)
    for i in range(h):
        a = hull[i]
        b = hull[(i + 1) % h]
        if _cross2(a, b, pt) <= 0:
            return False
    return True


# ---------------------------------------------------------------------------
# Affine Gale diagram: face / coface / missing-face computation
# ---------------------------------------------------------------------------

class GaleDiagram:
    """An affine Gale diagram for a simple d-polytope with n facets.

    Args:
        points: list of n (x,y) coordinates (integer/Fraction).
        positive_indices: set/list of indices with sigma=+1.
            The rest are sigma=-1.
        d: ambient polytope dimension.
    """

    def __init__(self, points, positive_indices, d):
        self.n = len(points)
        self.points = [tuple(p) for p in points]
        self.positive = frozenset(positive_indices)
        self.negative = frozenset(range(self.n)) - self.positive
        self.d = d
        self._face_cache = {}

    def is_face(self, S):
        """Test whether S (subset of facet indices) is a face of P.

        A subset F of facets is a face iff in T = [n]\\F every positive
        point lies in conv(negative points of T).
        """
        S = frozenset(S)
        if S in self._face_cache:
            return self._face_cache[S]

        T = frozenset(range(self.n)) - S
        pos_T = [self.points[i] for i in sorted(T) if i in self.positive]
        neg_T = [self.points[i] for i in sorted(T) if i in self.negative]

        if not T:
            # Empty complement (S = [n]): the empty face — vacuously a face
            result = True
        elif not pos_T:
            # T non-empty but no positive points: in the 3D Gale space all
            # vectors in T have the same sign (z = -1), so 0 ∉ relint(conv(T)).
            # Equivalently the two positive facets cannot both lie on a common face.
            result = False
        elif not neg_T:
            # Positive points but no negative => cannot be captured
            result = False
        else:
            # Affine Gale criterion: conv(pos_T) ∩ conv(neg_T) ≠ ∅.
            # For 1 positive point: p₀ ∈ conv(neg_T).
            # For 2 positive points: segment [p₀,p₁] meets conv(neg_T)
            # (weaker than "both inside"; correct per Ziegler §6 since
            # the conic-hull condition translates to segment intersection).
            if len(pos_T) == 1:
                result = point_in_convex_hull_2d(pos_T[0], neg_T, fast=False)
            else:
                result = _segment_meets_conv_hull_2d(pos_T[0], pos_T[1], neg_T)

        self._face_cache[S] = result
        return result

    def valid_gale_diagram(self):
        """Check that the diagram is valid for a bounded polytope.

        The empty face (S=∅) must be a face strictly: conv(pos) meets
        conv(neg) in the interior.  For the 2-positive-point case, the
        segment [p₀,p₁] must have a point strictly inside conv(neg).
        """
        pos_pts = [self.points[i] for i in self.positive]
        neg_pts = [self.points[i] for i in self.negative]
        # Each positive pt strictly inside conv(neg) is sufficient
        # (implies segment meets interior of hull).
        return all(
            point_strictly_in_convex_hull_2d(p, neg_pts) for p in pos_pts
        )

    def compute_missing_faces(self):
        """Compute all missing faces (minimal non-faces).

        Returns sorted list of frozensets, each a missing face.
        """
        # A missing face has size >= 2 (size 1 = singleton, always a face
        # since removing one facet still leaves room for many vertices)
        # Strategy: enumerate all subsets and find minimal non-faces.

        non_faces = []
        for size in range(2, self.n + 1):
            for S in combinations(range(self.n), size):
                S_frozen = frozenset(S)
                if not self.is_face(S_frozen):
                    # Check if it's minimal (no proper subset is a non-face)
                    is_minimal = True
                    for sub_size in range(2, size):
                        for sub in combinations(S, sub_size):
                            if frozenset(sub) in non_faces_set:
                                is_minimal = False
                                break
                        if not is_minimal:
                            break
                    if is_minimal:
                        non_faces.append(S_frozen)
            # Build set for quick lookup after each size
            if size == 2:
                non_faces_set = set(non_faces)
            else:
                non_faces_set = set(non_faces)

        # Actually redo with proper tracking
        return self._compute_missing_faces_correct()

    def _compute_missing_faces_correct(self, max_size=None):
        """Correct missing-face computation via upward closure.

        max_size: stop searching beyond this size (default: n).
        For Coxeter polytopes, Lannér diagrams require size <= 5 (F2),
        so pass max_size=5 to skip larger subsets.
        """
        if max_size is None:
            max_size = self.n
        non_faces_set = set()
        missing_faces = []

        for size in range(2, max_size + 1):
            for S in combinations(range(self.n), size):
                S_frozen = frozenset(S)
                # Skip if already covered by a sub-non-face
                if any(frozenset(sub) in non_faces_set
                       for sub_size in range(2, size)
                       for sub in combinations(S, sub_size)):
                    continue
                if not self.is_face(S_frozen):
                    missing_faces.append(S_frozen)
                    non_faces_set.add(S_frozen)

        return missing_faces

    def disjoint_pairs(self):
        """Return the set of pairs of non-intersecting facets (size-2 missing faces)."""
        mf = self._compute_missing_faces_correct()
        return [m for m in mf if len(m) == 2]

    def p_count(self):
        """Number of disjoint facet pairs."""
        return len(self.disjoint_pairs())

    def vertex_sets(self):
        """Return all maximal faces of size d (vertices of the simple polytope).

        For a simple d-polytope each vertex is on exactly d facets.
        Returns list of frozensets of size d.
        """
        vertices = []
        for S in combinations(range(self.n), self.d):
            S_frozen = frozenset(S)
            if self.is_face(S_frozen):
                # Check it's a vertex (maximal face of size d)
                # For simple polytopes this is automatic for size-d faces
                vertices.append(S_frozen)
        return vertices

    def is_simple(self):
        """Check simplicity: every vertex is on exactly d facets."""
        # For a simple d-polytope, all vertices have degree d.
        # We verify by checking the f-vector conditions.
        # Quick check: every (d-1)-face is contained in exactly 2 facets.
        # Here we just count: every maximal face of size d should exist.
        # Actually for Gale diagrams of simple polytopes this is automatic
        # when the configuration is in general position.
        # We check: no set of d+1 facets is a face (that would give a
        # non-simple vertex).
        for S in combinations(range(self.n), self.d + 1):
            if self.is_face(frozenset(S)):
                return False
        return True

    def combinatorial_type_key(self):
        """A hashable key representing the combinatorial type.

        Uses the sorted tuple of missing faces (as sorted tuples).
        This is a combinatorial invariant but not a canonical form;
        use canonical_form() for deduplication.
        """
        mf = self._compute_missing_faces_correct()
        return tuple(sorted(tuple(sorted(m)) for m in mf))
