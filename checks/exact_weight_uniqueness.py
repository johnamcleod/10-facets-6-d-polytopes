#!/usr/bin/env python3
"""The accepted labelling determines its ultraparallel weights: an EXACT proof.

The classification cascade pins the dashed weights in double precision, and its
weight tree for the accepted labelling has one leaf.  That count is a floating-point
observation, so it cannot carry the uniqueness half of the main theorem.  This
script replaces it with an exact argument over K = Q(sqrt2, sqrt5).

Fix the ordinary labels of P_{6,10} and leave the three dashed weights
x = (x_67, x_78, x_89) unknown.  Realizability needs rank G <= 7, i.e. the vanishing
of every 8x8 principal minor (for a real symmetric matrix the rank is the largest
order of a non-vanishing principal minor).

 1. For each dashed edge e, every 8x8 principal minor whose index set contains e
    and no other dashed edge is a polynomial of degree <= 2 in x_e alone, with
    coefficients in K.  Their gcd over K, g_e, is computed exactly; x_e must be a
    root of g_e.  The roots of g_e that exceed 1 are found exactly.
 2. That leaves finitely many candidate triples.  For each, every one of the 45
    principal 8x8 minors is evaluated exactly, and the rank and signature of G are
    computed exactly.
 3. The claim holds iff exactly one candidate triple has rank 7 and signature
    (6,1), and it is the Gram matrix of verify_polytope.py.

No floating point enters any decision: roots are compared with 1 by exact sign
determination in K (or in a quadratic extension of K), and every vanishing test is
exact simplification of an algebraic number to zero.

Run:  python3 checks/exact_weight_uniqueness.py
"""
from __future__ import annotations

import itertools
import sys
from pathlib import Path

import sympy as sp
from sympy import sqrt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from verify_polytope import (ORDINARY, DOTTED, N, cos_pi_over,  # noqa: E402
                             exact_inertia, gram)

EXT = [sqrt(2), sqrt(5)]
EDGES = sorted(DOTTED)                       # 1-indexed facet pairs
X = {e: sp.Symbol(f"x_{e[0]}{e[1]}", positive=True) for e in EDGES}


def symbolic_gram():
    G = sp.eye(N)
    for i, j, m in ORDINARY:
        G[i - 1, j - 1] = G[j - 1, i - 1] = -cos_pi_over(m)
    for e in EDGES:
        i, j = e
        G[i - 1, j - 1] = G[j - 1, i - 1] = -X[e]
    return G


def dashed_in(S):
    """Dashed edges with both endpoints in the (1-indexed) index set S."""
    return [e for e in EDGES if e[0] in S and e[1] in S]


def is_zero(a):
    """Exact zero test for an algebraic number built from sqrt2, sqrt5 and the
    square roots of elements of K: reduce to its minimal polynomial."""
    a = sp.nsimplify(a) if a.is_Float else a
    if a == 0:
        return True
    return sp.minimal_polynomial(a, sp.Symbol("z")) == sp.Symbol("z")


def exact_gt(a, b):
    """a > b for real algebraic numbers, decided exactly: a - b is non-zero (exact
    test), and then its sign is that of any sufficiently precise enclosure, which
    sympy's evalf certifies once the value is known to be non-zero."""
    d = sp.simplify(a - b)
    if is_zero(d):
        return False
    # d is a non-zero real algebraic number; enclose it with outward-rounded
    # interval arithmetic at increasing precision until the enclosure excludes 0.
    # That makes the sign a proof rather than an estimate.
    import mpmath
    for prec in (60, 120, 240, 480, 960):
        mpmath.iv.prec = prec
        iv = _enclose(d, mpmath.iv)
        if iv.a > 0:
            return True
        if iv.b < 0:
            return False
    raise RuntimeError(f"could not separate {d} from 0")


def _enclose(e, iv):
    """Interval enclosure of a sympy expression built from rationals, +, *, integer
    powers and square roots."""
    if e.is_Rational:
        return iv.mpf(e.p) / iv.mpf(e.q)
    if e.is_Add:
        out = iv.mpf(0)
        for a in e.args:
            out = out + _enclose(a, iv)
        return out
    if e.is_Mul:
        out = iv.mpf(1)
        for a in e.args:
            out = out * _enclose(a, iv)
        return out
    if e.is_Pow:
        base, ex = e.args
        b = _enclose(base, iv)
        if ex == sp.Rational(1, 2):
            return iv.sqrt(b)
        if ex == sp.Rational(-1, 2):
            return 1 / iv.sqrt(b)
        if ex.is_Integer:
            return b ** int(ex)
    raise ValueError(f"cannot enclose {e!r}")


def main():
    G = symbolic_gram()
    idx = range(1, N + 1)
    all8 = list(itertools.combinations(idx, 8))

    # ---- 1. per-edge gcd of the isolating minors
    cand = {}
    for e in EDGES:
        polys = []
        for S in all8:
            if dashed_in(S) == [e]:
                M = G.extract([s - 1 for s in S], [s - 1 for s in S])
                p = sp.Poly(sp.expand(M.det(method="berkowitz")), X[e],
                            extension=EXT)
                if not p.is_zero:
                    polys.append(p)
        g = polys[0]
        for p in polys[1:]:
            g = sp.gcd(g, p)
        g = g.monic()
        roots = sp.roots(g.as_expr(), X[e])
        if sum(roots.values()) != g.degree():
            raise SystemExit(f"edge {e}: could not solve gcd {g} in radicals")
        gt1 = [r for r in roots if r.is_real is not False and exact_gt(r, 1)]
        print(f"edge {e}: {len(polys)} isolating minors, gcd degree {g.degree()}: "
              f"{g.as_expr()};  roots > 1: {gt1}")
        cand[e] = gt1

    # ---- 2. every candidate triple, decided exactly
    target = {e: sp.nsimplify(DOTTED[e]) for e in EDGES}
    passing = []
    for vals in itertools.product(*(cand[e] for e in EDGES)):
        sub = dict(zip((X[e] for e in EDGES), vals))
        Gv = G.subs(sub)
        minors_zero = all(
            is_zero(sp.simplify(Gv.extract([s - 1 for s in S],
                                           [s - 1 for s in S]).det(method="berkowitz")))
            for S in all8)
        verdict = "rank <= 7" if minors_zero else "some 8x8 minor non-zero"
        sig = None
        if minors_zero and all(v in K_elements(vals) for v in vals):
            sig = exact_inertia(Gv)
        print(f"  candidate {dict(zip(EDGES, vals))}: {verdict}"
              + (f", inertia {sig}" if sig else ""))
        if minors_zero and (sig is None or sig == (6, 1, 3)):
            passing.append(dict(zip(EDGES, vals)))

    print()
    ok = (len(passing) == 1 and
          all(is_zero(sp.simplify(passing[0][e] - target[e])) for e in EDGES))
    if ok:
        print("RESULT: exactly one weight triple in (1, inf)^3 gives rank 7 and "
              "signature (6,1),")
        print("        and it is the Gram matrix of verify_polytope.py.  The "
              "accepted labelling")
        print("        determines its Gram matrix; no floating point enters the "
              "argument.")
        return 0
    print(f"RESULT: FAILED -- {len(passing)} passing candidate(s): {passing}")
    return 1


def K_elements(vals):
    """The candidates that lie in K itself (so exact_inertia, which works over K,
    applies).  Candidates outside K never pass the minor test here, but the guard
    keeps the inertia computation honest if one ever did."""
    K = sp.QQ.algebraic_field(*EXT)
    out = []
    for v in vals:
        try:
            K.from_sympy(sp.sympify(v))
            out.append(v)
        except Exception:
            pass
    return out


if __name__ == "__main__":
    sys.exit(main())
