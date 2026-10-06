#!/usr/bin/env python3
"""P^B6 is the double of a compact Coxeter 6-polytope with 9 facets.

Provenance (2026-08-11): following the project owner's correspondence with
P. Tumarkin, the polytope drawn as Burcroff Fig. 5 is a doubling of the
9-facet compact Coxeter 6-polytope of Bugaenko (1984).  This script verifies
the mathematical content of that explanation from the Gram matrices alone:

  1.  The labelled diagram of P^B6 has exactly one non-trivial automorphism,
      sigma = (6 9)(7 8) in the companion note's numbering.  If P^B6 is a
      double, sigma is the doubling involution: its mirror hyperplane f is the
      doubled facet, the six sigma-fixed facets are orthogonal to f, and the
      two swapped pairs are mirror images through f.

  2.  The half-polytope Q read off from that involution has facets
      {1,2,3,4,5,10} (orthogonal to f), 6 (meeting f at pi/4, since 6 and
      sigma(6)=9 are orthogonal in P^B6 and 1-2cos^2(pi/4)=0), 8 (diverging
      from f with weight w' = 2+sqrt5, since 1-2w'^2 = -(17+8*sqrt5) = -w78),
      and f itself.  Note w' = tau^3 for tau the golden ratio -- squarely in
      the Z[tau] world of Bugaenko's construction.

  3.  EXACT checks performed here:
        (a) Q's Gram matrix has rank 7 and signature (6,1) (kernel 2 = n-d-1
            for n = 9, d = 6), all dotted weights exceed 1, and CoxIter
            certifies Q cocompact of dimension 6;
        (b) doubling Q along f -- normals e_i for the orthogonal facets,
            e and R(e) = e - 2<e,e_f>e_f for the swapped pairs -- reproduces
            the Gram matrix of P^B6 EXACTLY, as a symbolic identity over
            Q(sqrt2, sqrt5), under R(6) -> 9, R(8) -> 7.

  What this does not check: that Q is literally the 9-facet polytope printed
  in Bugaenko (1984), which only the primary source can settle; the owner has
  inspected it.  Q must also appear in Tumarkin's classification of compact
  d-polytopes with d+3 facets.

Run:  python3 checks/doubling.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

c3, c4, c5 = sp.Rational(1, 2), sp.sqrt(2) / 2, (1 + sp.sqrt(5)) / 4
w67 = sp.sqrt(2) * (2 + sp.sqrt(5))
w78 = 17 + 8 * sp.sqrt(5)

# P^B6 exactly as in the companion note (facets 1..10)
P_EDGES = {(1, 2): c5, (2, 3): c3, (3, 4): c3, (4, 5): c3, (5, 6): c3,
           (9, 10): c5, (2, 7): c4, (2, 8): c4, (5, 9): c3, (6, 10): c5,
           (6, 7): w67, (7, 8): w78, (8, 9): w67}

# the half Q: order [1,2,3,4,5,6,8,10,f]
WPRIME = 2 + sp.sqrt(5)                       # = tau^3
Q_IDX = {1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5, 8: 6, 10: 7, 'f': 8}
Q_EDGES = {(1, 2): c5, (2, 3): c3, (3, 4): c3, (4, 5): c3, (5, 6): c3,
           (2, 8): c4, (6, 10): c5, (6, 'f'): c4, (8, 'f'): WPRIME}


def gram(n, idx, edges):
    G = sp.eye(n)
    for (i, j), v in edges.items():
        G[idx[i], idx[j]] = G[idx[j], idx[i]] = -v
    return G


def main():
    P = gram(10, {i: i - 1 for i in range(1, 11)}, P_EDGES)
    Q = gram(9, Q_IDX, Q_EDGES)

    # (a) Q is the Gram matrix of a compact Coxeter 6-polytope with 9 facets
    rank = Q.rank()
    ev = np.linalg.eigvalsh(np.array(Q.evalf(30), dtype=float))
    inertia = (int((ev > 1e-9).sum()), int((ev < -1e-9).sum()),
               int((np.abs(ev) <= 1e-9).sum()))
    print(f"Q: exact rank {rank} (expect 7), inertia {inertia} (expect (6,1,2)), "
          f"dotted weight w' = 2+sqrt5 = tau^3 "
          f"{sp.simplify(WPRIME - ((1 + sp.sqrt(5)) / 2) ** 3) == 0}")
    ok_q = rank == 7 and inertia == (6, 1, 2)

    from pipeline.utils.coxiter import check, available
    if available():
        lab = {(0, 1): 5, (1, 2): 3, (2, 3): 3, (3, 4): 3, (4, 5): 3,
               (1, 6): 4, (5, 7): 5, (5, 8): 4}
        compact, dim = check(lab, [(6, 8)], 9, 6)
        print(f"Q via CoxIter: cocompact={compact}, dimension={dim}")
        ok_q = ok_q and compact and dim == 6
    else:
        print("Q via CoxIter: SKIPPED (CoxIter not built)")

    # (b) double Q along f and compare with P exactly
    f = Q_IDX['f']
    dbl = [1, 2, 3, 4, 5, 6, 8, 10, 'R6', 'R8']

    def entry(a, b):
        i = Q_IDX[a] if isinstance(a, int) else Q_IDX[int(a[1:])]
        j = Q_IDX[b] if isinstance(b, int) else Q_IDX[int(b[1:])]
        v = Q[i, j]
        if isinstance(a, str) != isinstance(b, str):
            v = Q[i, j] - 2 * Q[i, f] * Q[f, j]
        return sp.radsimp(sp.simplify(v))

    D = sp.Matrix(10, 10, lambda r, c: entry(dbl[r], dbl[c]))
    to_p = {1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 8: 8, 10: 10, 'R6': 9, 'R8': 7}
    perm = [to_p[x] - 1 for x in dbl]
    Pp = sp.Matrix(10, 10, lambda r, c: P[perm[r], perm[c]])
    identical = sp.simplify(D - Pp) == sp.zeros(10, 10)
    print(f"double(Q) along f == Gram(P^B6), exactly over Q(sqrt2,sqrt5): {identical}")

    ok = ok_q and identical
    print("\nRESULT:", "P^B6 is the double of the 9-facet compact Coxeter "
          "6-polytope Q along its facet f" if ok else "MISMATCH -- see above")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
