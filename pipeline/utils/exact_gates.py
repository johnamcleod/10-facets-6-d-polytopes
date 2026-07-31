"""Tolerance-free forward-checking gates for the label enumeration.

The enumerator prunes a partial labelling with two predicates on a labelled
diagram: a face's subdiagram must be ELLIPTIC (positive definite), and a missing
face's subdiagram must be LANNER.  Deciding either by a floating-point eigenvalue
test leaves an emptiness verdict resting on a tolerance, and for the 6-node vertex
condition -- the one that does most of the pruning -- the tolerance cannot be
audited by enumeration, there being 6^15 label tuples.

Neither predicate needs a tolerance.

ELLIPTIC is a combinatorial condition: a diagram is elliptic if and only if every
connected component is one of the finite Coxeter types.  No arithmetic at all.
Delegated to pipeline/utils/coxeter_exact.py.

LANNER is decided here as: connected, not elliptic, every proper subdiagram
elliptic, and Gram determinant < 0.  The last clause is what separates Lanner from
PARABOLIC, and it is exact:

  * on three nodes the sign of the determinant is the sign of
    1/p + 1/q + 1/r - 1, where a non-edge counts as r = 2 (the classical
    spherical/Euclidean/hyperbolic trichotomy for the triangle group).  That is
    INTEGER arithmetic, and it holds for every label including m >= 7, so the
    wildcard needs no special case: raising a label only decreases the sum, and
    the sum is already below 1 at m = 7 for every connected three-node diagram,
    so the verdict is uniform over m >= 7.

  * on four or five nodes every proper subdiagram is elliptic, hence the leading
    (k-1)-minor is positive, so by Jacobi's rule the inertia is fixed by the sign
    of the determinant alone.  With all labels at most 6 the entries -cos(pi/m)
    lie in K = Q(sqrt2, sqrt3, sqrt5) and the determinant is computed exactly
    there (pipeline/utils/exact_field.py).  A label >= 7 cannot occur in a Lanner
    diagram of order >= 4 at all -- for any third node w of such a diagram the
    subdiagram on {u, v, w} is elliptic, hence w is orthogonal to both endpoints
    of the heavy edge, so that edge would be a connected component -- so those are
    rejected outright and cos(pi/7), which is cubic over Q, never arises.

SOUNDNESS AND STRENGTH.  The float gate tested the weaker necessary condition
"exactly one negative eigenvalue" in place of Lanner membership.  Requiring
membership is sound, since the subdiagram of a missing face IS Lanner, and it
prunes at least as much, so it can only shrink the search.  Both gates are
memoized on the label tuple: the enumerator asks the same questions about the same
small diagrams many millions of times.
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache

from pipeline.utils import coxeter_exact as CX
from pipeline.utils import exact_field as EF

# -cos(pi/m) as exact elements of K, for the labels that can occur in a Lanner
# diagram of order >= 4.  Built by radicand, never by basis position.
_ENTRY = {
    2: EF.rat(0),
    3: EF.rat(Fraction(-1, 2)),
    4: EF.scale(EF.sqrt_basis(2), Fraction(-1, 2)),
    5: EF.add(EF.rat(Fraction(-1, 4)),
              EF.scale(EF.sqrt_basis(5), Fraction(-1, 4))),
    6: EF.scale(EF.sqrt_basis(3), Fraction(-1, 2)),
}


def _pairs(k):
    return [(a, b) for a in range(k) for b in range(a + 1, k)]


@lru_cache(maxsize=1 << 20)
def elliptic_core(k, labs):
    """Exact: is the diagram with these labels (lex pair order) elliptic?"""
    lab = {p: m for p, m in zip(_pairs(k), labs)}
    return CX.is_elliptic(k, lab)


def _connected(k, lab):
    edges = {(a, b) for (a, b), m in lab.items() if m >= 3}
    return len(CX._components(k, edges)) == 1


def _sub_labels(k, lab, drop):
    keep = [u for u in range(k) if u != drop]
    idx = {u: i for i, u in enumerate(keep)}
    sub = {}
    for (a, b), m in lab.items():
        if a in idx and b in idx:
            sub[tuple(sorted((idx[a], idx[b])))] = m
    return tuple(sub.get(p, 2) for p in _pairs(k - 1))


def _det_negative_in_K(k, lab):
    """Exact sign test det < 0 over K.  All labels must be at most 6."""
    rows = []
    for i in range(k):
        row = []
        for j in range(k):
            if i == j:
                row.append(EF.ONE)
            else:
                m = lab.get((i, j) if i < j else (j, i), 2)
                row.append(_ENTRY[m])
        rows.append(row)
    det = EF.det(rows)
    if EF.is_zero(det):
        return False                      # parabolic, not Lanner
    iv = EF.alg_iv(det)
    if iv.hi < 0:
        return True
    if iv.lo > 0:
        return False
    # The enclosure straddles zero although the element is nonzero: decide by the
    # exact sign of the norm times a rational multiple.  This has not been observed
    # to arise (the enclosures are tight to ~1e-30), so it is a guard rather than a
    # code path, and it refuses to guess.
    raise ArithmeticError(f"cannot decide sign of determinant exactly: {det}")


@lru_cache(maxsize=1 << 20)
def lanner_core(k, labs):
    """Exact: is the diagram with these labels (lex pair order) Lanner?"""
    lab = {p: m for p, m in zip(_pairs(k), labs)}
    if not _connected(k, lab):
        return False
    if k == 2:
        # a single edge: Lanner iff the two facets diverge, which is a dashed edge
        # and not an ordinary label; an ordinary edge of finite label is elliptic
        return False
    if k == 3:
        # sign(det) = sign(1/p + 1/q + 1/r - 1), integer arithmetic, valid for
        # every label including the wildcard proxy
        s = sum(Fraction(1, m) for m in labs)
        return s < 1
    if any(m >= 7 for m in labs):
        return False                      # no Lanner diagram of order >= 4 has one
    if elliptic_core(k, labs):
        return False
    for drop in range(k):
        if not elliptic_core(k - 1, _sub_labels(k, lab, drop)):
            return False
    return _det_negative_in_K(k, lab)


# ---------------------------------------------------------------------------
# The shapes the enumerator calls with
# ---------------------------------------------------------------------------

def _labs_from_assignment(nodes, assignment):
    nodes = list(nodes)
    pos = {u: i for i, u in enumerate(nodes)}
    lab = {}
    for a in range(len(nodes)):
        for b in range(a + 1, len(nodes)):
            p = (nodes[a], nodes[b]) if nodes[a] < nodes[b] else (nodes[b], nodes[a])
            lab[(a, b)] = assignment.get(p, 2)
    return tuple(lab[p] for p in _pairs(len(nodes)))


def vertex_pd_exact(v_sorted, assignment):
    """Exact replacement for the on-the-fly float PD test on a face."""
    return elliptic_core(len(v_sorted), _labs_from_assignment(v_sorted, assignment))


def lanner_exact(face_sorted, assignment, k):
    """Exact replacement for the on-the-fly float Lanner test on a missing face."""
    return lanner_core(k, _labs_from_assignment(face_sorted, assignment))


def cache_info():
    return {"elliptic": elliptic_core.cache_info(),
            "lanner": lanner_core.cache_info()}
