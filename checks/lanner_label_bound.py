#!/usr/bin/env python3
"""Verify the classification fact that the wildcard lemma rests on:

    no Lanner diagram of order >= 4 contains an edge of label >= 6.

A Lanner diagram of order k is a connected Coxeter diagram on k nodes that is
not elliptic but all of whose proper subdiagrams are elliptic.  The paper cites
Lanner's classification for this; here we re-derive it, so that nothing depends
on reading a published figure correctly.

The derivation is in three steps.

STEP 1 (verified exhaustively below).  There is no CONNECTED elliptic diagram on
3 nodes with an edge of label >= 6.  Equivalently, every elliptic rank-3 diagram
containing an edge of label m >= 6 is A_1 x I_2(m): the third node is orthogonal
to both endpoints.  (This is the classification of rank-3 finite Coxeter groups
-- A_3, B_3, H_3 all have labels <= 5.)  We check it by enumerating all label
triples over {2,...,MMAX}.  Raising a label only decreases the 3x3 Gram
determinant, monotonically, so MMAX = 200 is far past sufficient.

STEP 2 (immediate consequence, no enumeration needed).  Let L be a Lanner
diagram of order k >= 4 and suppose the edge uv has label m >= 6.  For any third
node w of L, the subdiagram on {u, v, w} is a proper subdiagram of L, hence
elliptic; by Step 1 it must be A_1 x I_2(m), so w is orthogonal to both u and v.
As this holds for every w outside {u, v}, and k >= 4 guarantees such a w exists,
{u, v} is a connected component of L -- contradicting the connectedness of a
Lanner diagram.  Hence no edge of L has label >= 6.

STEP 3 (sanity check, verified below).  With labels therefore confined to
{2,...,5}, enumerate the order-4 and order-5 Lanner diagrams directly from the
definition and report how many there are and the largest label occurring.

Expected output: no violation in Step 1; Step 3 finds Lanner diagrams in both
orders, with largest edge label 5.

Run:  python3 checks/lanner_label_bound.py
"""
import itertools
import math
import sys

import numpy as np

MMAX = 200

# Positive-definiteness tolerance.  It is genuinely needed: the connected rank-3
# diagram with labels (3, 6) -- affine G_2-tilde -- has Gram determinant exactly
# 1 - cos^2(pi/3) - cos^2(pi/6) = 1 - 1/4 - 3/4 = 0, so it is PARABOLIC, not
# elliptic.  Without a tolerance float64 renders that zero eigenvalue as a tiny
# positive number and the diagram is misread as elliptic, producing a spurious
# "violation" of Step 1.  Parabolic subdiagrams are excluded from Lanner diagrams
# by definition (all proper subdiagrams must be positive DEFINITE), so screening
# them out here is correct, not a fudge.
PD_TOL = 1e-10


def cosv(m):
    return -math.cos(math.pi / m)


def is_pd(M):
    try:
        np.linalg.cholesky(M)
        return True
    except np.linalg.LinAlgError:
        return False


def gram(k, pairs, labels):
    M = np.eye(k)
    for (a, b), m in zip(pairs, labels):
        M[a, b] = M[b, a] = cosv(m)
    return M


def connected(k, pairs, labels):
    adj = {i: set() for i in range(k)}
    for (a, b), m in zip(pairs, labels):
        if m >= 3:                      # m == 2 means no edge
            adj[a].add(b)
            adj[b].add(a)
    seen, stack = {0}, [0]
    while stack:
        for y in adj[stack.pop()]:
            if y not in seen:
                seen.add(y)
                stack.append(y)
    return len(seen) == k


def step1(mmax=MMAX):
    """No CONNECTED elliptic rank-3 diagram has an edge of label >= 6."""
    pairs = [(0, 1), (0, 2), (1, 2)]
    ms = np.arange(2, mmax + 1)
    grid = np.array(np.meshgrid(ms, ms, ms, indexing="ij")).reshape(3, -1).T
    grid = grid[grid.max(1) >= 6]
    tested = len(grid)
    conn = np.array([connected(3, pairs, tuple(g)) for g in grid])
    grid = grid[conn]
    M = np.zeros((len(grid), 3, 3))
    M[:, range(3), range(3)] = 1.0
    v = -np.cos(np.pi / grid)
    for j, (a, b) in enumerate(pairs):
        M[:, a, b] = v[:, j]
        M[:, b, a] = v[:, j]
    ev = np.linalg.eigvalsh(M)
    violations = [tuple(g) for g in grid[(ev > PD_TOL).all(1)]]
    parabolic = int(((np.abs(ev) <= PD_TOL).any(1)
                     & (ev > -PD_TOL).all(1)).sum())
    print(f"STEP 1: scanned m in 2..{mmax}; {tested:,} label triples having some "
          f"label >= 6")
    print(f"        connected elliptic rank-3 diagrams with an edge of label "
          f">= 6: {len(violations)}"
          + (f"   e.g. {violations[:3]}" if violations else "   (none)  OK"))
    print(f"        (of these, {parabolic} are parabolic -- e.g. affine G_2-tilde, "
          f"labels (3,6);\n         parabolic diagrams are not elliptic and are "
          f"excluded from Lanner diagrams by definition)")
    return len(violations) == 0


def step3():
    """Enumerate order-4 and order-5 Lanner diagrams over labels {2,...,5}."""
    ok = True
    for k in (4, 5):
        pairs = [(a, b) for a in range(k) for b in range(a + 1, k)]
        grid = np.array([g for g in itertools.product(range(2, 6),
                                                       repeat=len(pairs))
                         if connected(k, pairs, g)])
        M = np.zeros((len(grid), k, k))
        M[:, range(k), range(k)] = 1.0
        v = -np.cos(np.pi / grid)
        for j, (a, b) in enumerate(pairs):
            M[:, a, b] = v[:, j]
            M[:, b, a] = v[:, j]
        # Lanner = connected, HYPERBOLIC (exactly one negative eigenvalue), all
        # proper subdiagrams elliptic.  Requiring merely "not elliptic" would also
        # admit the affine (parabolic) diagrams, whose proper subdiagrams are all
        # elliptic too -- those are not Lanner.  (Note that once every proper
        # principal submatrix is positive definite, Cauchy interlacing already
        # forces at most one non-positive eigenvalue, so the test below is
        # exactly "the last eigenvalue is negative rather than zero".)
        keep = (np.linalg.eigvalsh(M) < -PD_TOL).sum(1) == 1
        for r in range(1, k):                            # every proper sub elliptic
            for sub in itertools.combinations(range(k), r):
                idx = np.ix_(np.arange(len(grid)), sub, sub)
                keep &= (np.linalg.eigvalsh(M[idx]) > PD_TOL).all(1)
        found = [tuple(g) for g in grid[keep]]
        # collapse node relabellings to count ABSTRACT diagrams, so the total can
        # be compared with Lanner's published classification (9 of order 4,
        # 5 of order 5)
        pos = {p: i for i, p in enumerate(pairs)}
        classes = set()
        for g in found:
            best = None
            for perm in itertools.permutations(range(k)):
                img = tuple(g[pos[tuple(sorted((perm[a], perm[b])))]]
                            for (a, b) in pairs)
                if best is None or img < best:
                    best = img
            classes.add(best)
        if not found:
            print(f"STEP 3: order {k}: none found -- unexpected")
            ok = False
            continue
        maxlab = max(max(d) for d in found)
        expect = {4: 9, 5: 5}[k]
        print(f"STEP 3: order {k}: {len(found)} labelled tuples = "
              f"{len(classes)} diagrams up to relabelling "
              f"(Lanner's classification: {expect})"
              f"{'  OK' if len(classes) == expect else '  MISMATCH'}")
        print(f"        largest edge label occurring = {maxlab}")
        ok = ok and maxlab <= 5 and len(classes) == expect
    return ok


def main():
    a = step1()
    print()
    print("STEP 2: immediate from Step 1 (see the module docstring). In a Lanner")
    print("        diagram of order k >= 4, any triple containing a label->=6 edge")
    print("        is a proper subdiagram, hence elliptic, hence splits that edge")
    print("        off; so its endpoints form a component, contradicting")
    print("        connectedness.")
    print()
    b = step3()
    print()
    print("RESULT:",
          "no Lanner diagram of order 4 or 5 has an edge of label >= 6"
          if (a and b) else "CLAIM VIOLATED")
    return 0 if (a and b) else 1


if __name__ == "__main__":
    sys.exit(main())
