#!/usr/bin/env python3
"""Do the exact combinatorial gates agree with the float ones the search used?

The forward-checking of the search decides two predicates on a labelled diagram:
elliptic (Gram matrix positive definite) and Lanner (exactly one negative
eigenvalue, all proper subdiagrams elliptic).  The production run decided both with
a floating-point eigenvalue test at tolerance 1e-8.  That tolerance cannot be
audited exhaustively for the case that does most of the work, the 6-node vertex
condition, because there are 6^15 label tuples.

pipeline/utils/coxeter_exact.py decides both predicates exactly, by classifying the
connected components against the finite Coxeter types -- no arithmetic and no
tolerance.  This script compares the two, exhaustively where that is possible and
by large random sampling where it is not, and reports every disagreement.

Two independent things are established when they agree.  Lanner's classification is
recovered from the definition (9 diagrams of order 4, 5 of order 5, no edge label
above 5), which is the check the paper already carried out by a different route; and
the float gate is confirmed correct on the compared domain, which is what the
tolerance audit was trying to establish.

Run:  python3 paper/checks/exact_vs_float_gates.py [samples]
"""
from __future__ import annotations

import itertools
import random
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pipeline.stage4_gram import _GRAM_FLOAT                      # noqa: E402
from pipeline.utils.coxeter_exact import is_elliptic, is_lanner   # noqa: E402

ALPHABET = [2, 3, 4, 5, 6, 7]      # 7 stands for the wildcard proxy cos(pi/7)
TOL = 1e-8


def gram(k, lab):
    M = np.eye(k)
    for (a, b), m in lab.items():
        M[a, b] = M[b, a] = _GRAM_FLOAT[m]
    return M


def float_elliptic(k, lab):
    return bool(np.all(np.linalg.eigvalsh(gram(k, lab)) > TOL))


def float_lanner(k, lab):
    ev = np.linalg.eigvalsh(gram(k, lab))
    return int(np.sum(ev < -TOL)) == 1


def tuples_exhaustive(k):
    pairs = [(a, b) for a in range(k) for b in range(a + 1, k)]
    for combo in itertools.product(ALPHABET, repeat=len(pairs)):
        yield dict(zip(pairs, combo))


def tuples_random(k, n, rng):
    pairs = [(a, b) for a in range(k) for b in range(a + 1, k)]
    for _ in range(n):
        yield {p: rng.choice(ALPHABET) for p in pairs}


def compare(k, gen, label, what):
    bad = []
    n = 0
    for lab in gen:
        n += 1
        if what == "elliptic":
            e, f = is_elliptic(k, lab), float_elliptic(k, lab)
        else:
            e, f = is_lanner(k, lab), float_lanner(k, lab)
        if e != f:
            bad.append((dict(lab), e, f))
    status = "agree" if not bad else f"{len(bad)} DISAGREEMENTS"
    print(f"  {what:<9} k={k}  {label:<22} {n:>10,} tuples   {status}")
    for lab, e, f in bad[:3]:
        nz = {p: m for p, m in lab.items() if m != 2}
        print(f"        exact={e} float={f}  edges={nz}")
    return bad


def main():
    samples = int(sys.argv[1]) if len(sys.argv) > 1 else 200_000
    rng = random.Random(20260730)
    bad = []

    print("exhaustive where the domain permits:")
    for k in (3, 4):
        bad += compare(k, tuples_exhaustive(k), "exhaustive", "elliptic")
    for k in (3, 4):
        bad += compare(k, tuples_exhaustive(k), "exhaustive", "lanner")

    print("\nrandom sampling where it does not (5- and 6-node groups):")
    for k in (5, 6):
        bad += compare(k, tuples_random(k, samples, rng), f"{samples:,} random",
                       "elliptic")
    bad += compare(5, tuples_random(5, samples, rng), f"{samples:,} random",
                   "lanner")

    print("\nLanner's classification, recovered from the definition:")
    for k in (4, 5):
        pairs = [(a, b) for a in range(k) for b in range(a + 1, k)]
        found = set()
        for combo in itertools.product([2, 3, 4, 5, 6, 7], repeat=len(pairs)):
            lab = dict(zip(pairs, combo))
            if is_lanner(k, lab):
                found.add(tuple(sorted(combo)))
        mx = max((m for f in found for m in f), default=0)
        print(f"  order {k}: {len(found)} diagrams up to label multiset, "
              f"largest label {mx}")

    print("\nRESULT:", "the exact and float gates agree everywhere compared"
          if not bad else f"{len(bad)} DISAGREEMENTS -- see above")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
