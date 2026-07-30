#!/usr/bin/env python3
"""Exhaustive margin audit of the forward-checking predicates (Section 8, item 4).

Most labellings never reach the cascade: they are pruned during enumeration by
one of two float64 tests on a Gram submatrix,

    PD      : positive definite            (elliptic subdiagram, vertex/face groups)
    ONE-NEG : exactly one negative eigenvalue (necessary for a Lanner subdiagram)

with eigenvalue thresholds of 1e-8 to 1e-10.  Unlike the cascade, these are
decided on a FINITE set: the label alphabet in wildcard mode is
{2,3,4,5,6,*} and the groups have 3 to 6 nodes.  For 3 and 4 nodes the whole set
is small enough to audit exhaustively, which is what this script does.

WHICH ERRORS MATTER.  Both predicates are used to PRUNE, so only a false
rejection is unsound:

  * PD.  We prune when the submatrix is not PD.  A parabolic submatrix has an
    eigenvalue exactly 0, is genuinely not PD, and is correctly pruned however
    the float lands.  The only unsound case is a submatrix that is GENUINELY
    positive definite with its smallest eigenvalue below the threshold.
  * ONE-NEG.  We prune when the negative count is not 1.  Again the risk is a
    genuine eigenvalue too close to 0 to be counted correctly.

So the quantity to measure is the smallest NONZERO |eigenvalue| occurring
anywhere, compared with the threshold.  Exact zeros (the parabolic boundary --
e.g. affine G_2-tilde with labels 3,6) are not a hazard; they are reported
separately so the two populations are not confused.

Expected output: the smallest nonzero |eigenvalue| is many orders of magnitude
above 1e-8, so no forward-checking decision on 3 or 4 nodes is anywhere near its
threshold.

Group sizes 5 and 6 exceed the code's table cap and are tested on the fly; they
are sampled here rather than enumerated.

Run:  python3 paper/checks/forward_check_margins.py
"""
import itertools
import math
import sys

import numpy as np

ALPHABET = [2, 3, 4, 5, 6, 7]     # 7 == the wildcard, evaluated at cos(pi/7)
THRESH = 1e-8                     # the code's eigenvalue tolerance
ZERO = 1e-12                      # below this we call an eigenvalue exactly zero


def batch_eigs(k, grid):
    pairs = [(a, b) for a in range(k) for b in range(a + 1, k)]
    M = np.zeros((len(grid), k, k))
    M[:, range(k), range(k)] = 1.0
    v = -np.cos(np.pi / grid)
    for j, (a, b) in enumerate(pairs):
        M[:, a, b] = v[:, j]
        M[:, b, a] = v[:, j]
    return np.linalg.eigvalsh(M)


def audit(k, sample=None, seed=0):
    npair = k * (k - 1) // 2
    if sample:
        rng = np.random.default_rng(seed)
        grid = rng.choice(ALPHABET, size=(sample, npair))
        how = f"{sample:,} random tuples"
    else:
        grid = np.array(list(itertools.product(ALPHABET, repeat=npair)))
        how = f"all {len(grid):,} tuples (exhaustive)"

    out = []
    for start in range(0, len(grid), 200_000):
        out.append(batch_eigs(k, grid[start:start + 200_000]))
    ev = np.concatenate(out)

    a = np.abs(ev)
    zeros = int((a <= ZERO).sum())
    nz = a[a > ZERO]
    smallest = nz.min()
    print(f"nodes={k}: {how}, {ev.size:,} eigenvalues")
    print(f"    exactly-zero eigenvalues (parabolic boundary, harmless): "
          f"{zeros:,}")
    print(f"    smallest NONZERO |eigenvalue|: {smallest:.6e}   "
          f"threshold {THRESH:.0e}   margin factor {smallest/THRESH:.2e}")
    return smallest


def main():
    print("Exhaustive where feasible; the code tabulates 3- and 4-node groups and")
    print("falls back to on-the-fly tests for 5 and 6 nodes.\n")
    worst = min(audit(3), audit(4))
    print()
    print("Group sizes beyond the table cap (sampled, not exhaustive):")
    worst = min(worst, audit(5, sample=300_000), audit(6, sample=300_000))
    print()
    ok = worst > 100 * THRESH
    print(f"RESULT: smallest nonzero |eigenvalue| anywhere is {worst:.3e}, "
          f"a factor {worst/THRESH:.2e} above the {THRESH:.0e} threshold."
          if ok else
          f"RESULT: MARGIN TOO TIGHT -- smallest nonzero |eigenvalue| {worst:.3e}")
    if ok:
        print("        No forward-checking decision is near its tolerance, so no")
        print("        labelling can have been pruned by a sign error there.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
