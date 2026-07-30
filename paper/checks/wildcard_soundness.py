#!/usr/bin/env python3
"""Numerical confirmation of Lemma "wildcard" in the paper (Section 4.1).

The enumerator carries an ordinary edge of label m >= 7 as a WILDCARD and
evaluates the two forward-checking predicates at cos(pi/7):

    PD      : the Gram submatrix is positive definite   (elliptic subdiagram)
    ONE-NEG : the Gram submatrix has exactly one negative eigenvalue
              (necessary condition for a Lanner subdiagram)

The paper proves the substitution is lossless.  This script confirms it
numerically: for every label tuple on 3 and 4 nodes over {2,3,4,5,6,*} with at
least one wildcard (and for a random sample on 5 nodes), it compares the verdict
at m = 7 with the verdict at every m in a wide list.  Expected output: zero
discrepancies.

Note: the list deliberately stops at m = 1000.  At m = 10^6 the PD verdict
appears to flip for the disconnected tuples (2,2,...,*), but only because
1 - cos(pi/10^6) ~ 5e-12 falls below the 1e-9 eigenvalue threshold -- a
float64 tolerance artefact of this script, not of the enumerator, which never
instantiates m beyond 100.

Run:  python3 paper/checks/wildcard_soundness.py
"""
import itertools
import math
import sys

import numpy as np

FIXED = [2, 3, 4, 5, 6]
COSV = {m: -math.cos(math.pi / m) for m in FIXED}
CMIN = math.cos(math.pi / 7)
MS = [7, 8, 9, 10, 12, 15, 20, 30, 50, 100, 300, 1000]
CS = [math.cos(math.pi / m) for m in MS]
TOL = 1e-9


def check(k, sample=None, seed=0):
    pairs = [(a, b) for a in range(k) for b in range(a + 1, k)]
    npair = len(pairs)
    alpha = FIXED + ["W"]
    if sample:
        rng = np.random.default_rng(seed)
        combos = [tuple(alpha[i] for i in rng.integers(0, 6, npair))
                  for _ in range(sample)]
    else:
        combos = list(itertools.product(alpha, repeat=npair))
    combos = [c for c in combos if "W" in c]

    bad_pd = bad_ln = 0
    evals = 0
    for combo in combos:
        widx = [j for j, x in enumerate(combo) if x == "W"]
        w = len(widx)
        if w <= 2:
            grids = list(itertools.product(CS, repeat=w))
        else:                       # diagonal sample, keeps the cost bounded
            grids = [tuple([c] * w) for c in CS]
        M = np.zeros((len(grids) + 1, k, k))
        M[:, range(k), range(k)] = 1.0
        for j, (a, b) in enumerate(pairs):
            if j in widx:
                col = widx.index(j)
                vals = np.array([-CMIN] + [-g[col] for g in grids])
            else:
                vals = np.full(len(grids) + 1, COSV[combo[j]])
            M[:, a, b] = vals
            M[:, b, a] = vals
        ev = np.linalg.eigvalsh(M)
        neg = (ev < -TOL).sum(1)
        zer = (np.abs(ev) <= TOL).sum(1)
        pd = (neg == 0) & (zer == 0)
        ln = neg == 1
        evals += len(grids)
        if (pd[1:] != pd[0]).any():
            bad_pd += 1
        if (ln[1:] != ln[0]).any():
            bad_ln += 1

    print(f"k={k}: {len(combos):,} wildcard-bearing tuples, "
          f"{evals:,} (tuple, m) evaluations")
    print(f"    PD verdict differing from m=7      : {bad_pd}")
    print(f"    ONE-NEG verdict differing from m=7 : {bad_ln}")
    return bad_pd + bad_ln


def main():
    total = 0
    total += check(3)
    total += check(4)
    total += check(5, sample=8000)
    print()
    print("RESULT:", "no discrepancies -- proxy is lossless on the tested range"
          if total == 0 else f"{total} DISCREPANCIES")
    return 0 if total == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
