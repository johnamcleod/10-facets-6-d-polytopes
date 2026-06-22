"""Exact affine-Gale-diagram face criterion (integer arithmetic).

Reference implementation, written 2026-06-22 to replace the buggy float
(`gale.py`) and chirotope (`c/stage2_filter.c`) criteria that reproduced 0/30
of the d=4 ground-truth combinatorial types.  See memory `stage2-generator-bugs`.

Setup.  A simple d-polytope P with n = d+4 facets has an affine Gale diagram of
n labelled points in R^2 (= R^{n-d-2}), each with a sign sigma_i in {+1,-1}.
The Gale vectors live in R^3 (= R^{n-d-1}): positive points lift to (a_i, +1),
negative points to (a_i, -1).

Face criterion (Ziegler, Lectures on Polytopes, Thm 6.19 + the affine reduction
of §6.4).  F subset [n] indexes a face of P (the facets in F share a common
vertex/face)  iff  0 in relint conv{ Gale vectors of the complement T=[n]\\F }.
Translating the lift to the plane (derivation: write the nonneg combination
giving 0, split by sign, divide by the common mass):

    0 in relint conv{ v_j : j in T }
      <=>  relint conv(POS_T) intersect relint conv(NEG_T) != empty

where POS_T / NEG_T are the positive / negative *plane* points indexed by T.
BOTH relative interiors -- this is the fix; the old code used closed hulls.

We assume the AAK input is in general position (no 3 collinear), so the only
lower-dimensional relints come from |POS_T| or |NEG_T| <= 2.
"""

from itertools import combinations


# ---------------------------------------------------------------------------
# Exact 2D primitives (integer / any exact numeric tuples)
# ---------------------------------------------------------------------------

def _orient(o, a, b):
    """Twice the signed area of triangle (o,a,b); sign = orientation. Exact."""
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def _sgn(x):
    return (x > 0) - (x < 0)


def _hull_ccw(pts):
    """Convex hull (CCW vertex list) of distinct points via monotone chain.

    Returns the list of hull vertices in CCW order.  For <=2 distinct points
    returns them as-is (caller treats that as a degenerate/lower-dim hull).
    Assumes general position (no 3 collinear) so no collinear hull vertices.
    """
    P = sorted(set(map(tuple, pts)))
    if len(P) <= 2:
        return P
    lower = []
    for p in P:
        while len(lower) >= 2 and _orient(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(P):
        while len(upper) >= 2 and _orient(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def _strictly_inside(p, pts):
    """Is p in the relative interior (= interior, full-dim) of conv(pts)?

    Requires >=3 points spanning 2D.  Strict: returns False on the boundary.
    """
    hull = _hull_ccw(pts)
    if len(hull) < 3:
        return False           # lower-dimensional hull: empty interior in R^2
    h = len(hull)
    for i in range(h):
        if _orient(hull[i], hull[(i + 1) % h], p) <= 0:
            return False       # on/outside an edge (CCW hull => interior is > 0)
    return True


def _proper_seg_cross(a, b, c, d):
    """Do OPEN segments (a,b) and (c,d) cross at an interior point of both?

    Strict (proper) crossing: endpoints touching / collinear overlap -> False.
    """
    o1 = _orient(a, b, c)
    o2 = _orient(a, b, d)
    o3 = _orient(c, d, a)
    o4 = _orient(c, d, b)
    return (_sgn(o1) * _sgn(o2) < 0) and (_sgn(o3) * _sgn(o4) < 0)


def _seg_meets_relint_hull(a, b, pts):
    """Does the OPEN segment (a,b) meet the relative interior of conv(pts)?"""
    m = len(set(map(tuple, pts)))
    if m <= 1:
        return False
    if m == 2:
        u, v = sorted(set(map(tuple, pts)))
        return _proper_seg_cross(a, b, u, v)
    # m >= 3: full-dim polygon interior (general position).
    if _strictly_inside(a, pts) or _strictly_inside(b, pts):
        return True
    # Both endpoints outside/boundary: open segment meets the open interior iff
    # it properly crosses the boundary at >=2 places (enters and exits).
    hull = _hull_ccw(pts)
    h = len(hull)
    crossings = 0
    for i in range(h):
        if _proper_seg_cross(a, b, hull[i], hull[(i + 1) % h]):
            crossings += 1
            if crossings >= 2:
                return True
    return False


# ---------------------------------------------------------------------------
# Affine Gale diagram
# ---------------------------------------------------------------------------

class AffineGale:
    """Exact affine Gale diagram of a candidate simple d-polytope, n=d+4 facets.

    points: list of n integer (x,y) tuples.
    positive: iterable of the indices with sign +1 (the rest are -1).
    """

    def __init__(self, points, positive, d):
        self.n = len(points)
        self.pts = [tuple(p) for p in points]
        self.pos = frozenset(positive)
        self.positive = self.pos          # alias for callers expecting `.positive`
        self.neg = frozenset(range(self.n)) - self.pos
        self.d = d
        self._fc = {}

    def is_face(self, F):
        F = frozenset(F)
        c = self._fc.get(F)
        if c is not None:
            return c
        T = frozenset(range(self.n)) - F
        POS = [self.pts[i] for i in T if i in self.pos]
        NEG = [self.pts[i] for i in T if i in self.neg]
        if not POS or not NEG:
            r = False
        elif len(POS) == 1:
            r = _strictly_inside(POS[0], NEG)
        elif len(POS) == 2:
            r = _seg_meets_relint_hull(POS[0], POS[1], NEG)
        else:
            # >2 positive points: relint conv(POS) is a 2D region; it meets
            # relint conv(NEG) iff the two open polygons overlap.  Not needed
            # for the 2-positive setup, but handle defensively.
            r = self._relint_polys_overlap(POS, NEG)
        self._fc[F] = r
        return r

    def _relint_polys_overlap(self, A, B):
        # General fallback: open convex polygons overlap iff some vertex of one
        # is strictly inside the other, or their boundaries properly cross.
        if any(_strictly_inside(p, B) for p in A):
            return True
        if any(_strictly_inside(p, A) for p in B):
            return True
        ha, hb = _hull_ccw(A), _hull_ccw(B)
        for i in range(len(ha)):
            for j in range(len(hb)):
                if _proper_seg_cross(ha[i], ha[(i + 1) % len(ha)],
                                     hb[j], hb[(j + 1) % len(hb)]):
                    return True
        return False

    def missing_faces(self, max_size=None):
        """Minimal non-faces (sizes 2..max_size), via upward closure."""
        if max_size is None:
            max_size = self.n
        nf = set()
        out = []
        for size in range(2, max_size + 1):
            for S in combinations(range(self.n), size):
                Sf = frozenset(S)
                if size > 2 and any(
                    frozenset(sub) in nf
                    for ss in range(2, size)
                    for sub in combinations(S, ss)
                ):
                    continue
                if not self.is_face(Sf):
                    nf.add(Sf)
                    out.append(Sf)
        return out

    def vertex_sets(self):
        """All vertices of the simple d-polytope = the size-d faces.

        (Compact/simple => every vertex lies on exactly d facets, so vertices
        are exactly the d-subsets that are faces.)
        """
        return [frozenset(S) for S in combinations(range(self.n), self.d)
                if self.is_face(frozenset(S))]

    def f0(self):
        """Number of vertices = number of size-d faces."""
        return len(self.vertex_sets())

    def is_polytope(self):
        """Necessary structural checks that the diagram is a simple d-polytope.

        - every singleton {i} is a face (each facet is real);
        - no (d+1)-subset is a face (simplicity: no vertex on > d facets);
        - Lower Bound Theorem on f0 for the dual simplicial d-polytope on n verts.
        """
        if not all(self.is_face(frozenset([i])) for i in range(self.n)):
            return False
        for S in combinations(range(self.n), self.d + 1):
            if self.is_face(frozenset(S)):
                return False
        d, n = self.d, self.n
        lbt = (d - 1) * n - (d + 1) * (d - 2)   # min facets of simplicial d-poly, n verts
        return self.f0() >= lbt
