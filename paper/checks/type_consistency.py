#!/usr/bin/env python3
"""Consistency of the generated combinatorial types, checked exhaustively.

The generator produces two representations of each type -- a missing-face
hypergraph and a vertex list -- and the whole pipeline relies on them agreeing, so
the agreement is worth checking rather than assuming.  This script checks every generated
type in every dimension available, and reports three distinct things, which must
not be conflated.

  (A) ERRORS.  A vertex that is not a 6-subset, or a vertex that CONTAINS a
      recorded missing face.  Either would mean the two representations
      contradict each other, since a subset of a face is a face and therefore no
      missing face can lie inside one.  Any hit here is a bug.

  (B) A THEOREM-KILLABLE TYPE, not an error.  A 6-subset that contains no recorded
      missing face but is not a vertex.  The generator records minimal non-faces up
      to size 5 only (`missing_faces(max_size=5)`), so such a subset is a minimal
      non-face of size 6 that simply was not written down.  A type possessing one
      cannot carry a compact Coxeter structure at all: by Vinberg's correspondence
      its Coxeter diagram would need a Lannér subdiagram of order 6, and Lannér
      diagrams exist only up to order 5.  These types are therefore eliminable by
      theorem, and the paper's claim that the size-<=5 bound "eliminates nothing"
      is false -- it eliminates them.

  (C) GENERAL POSITION of the stored order-type coordinates: no three of the ten
      points collinear, tested by exact integer determinants.  The affine Gale
      criterion is stated for configurations in general position, so this is a
      hypothesis of the construction rather than a conclusion.

Run:  python3 paper/checks/type_consistency.py
"""
from __future__ import annotations

import itertools
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pipeline.stage4_blockpaste import _setup   # noqa: E402

CASES = [(6, 10), (5, 9), (4, 8)]


def vertices(t):
    V, *_ = _setup(dict(t))
    return {frozenset(v) for v in V}


def collinear_triples(pts):
    """Number of collinear triples, by exact integer determinant."""
    bad = 0
    for a, b, c in itertools.combinations(pts, 3):
        det = ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))
        if det == 0:
            bad += 1
    return bad


def check(d, n):
    path = ROOT / f"runs/d{d}_n{n}/stage2/types.json"
    if not path.exists():
        return None
    types = {t["type_id"]: t for t in json.load(open(path))}
    errors, killable, degen = [], [], []
    for tid, t in types.items():
        mf = [frozenset(m) for m in t["missing_faces"]]
        V = vertices(t)
        # (A) errors
        for v in V:
            if len(v) != d:
                errors.append((tid, "vertex not of size d", sorted(v)))
                break
        for v in V:
            for m in mf:
                if m <= v:
                    errors.append((tid, "missing face inside a vertex",
                                   (sorted(m), sorted(v))))
                    break
            else:
                continue
            break
        # (B) 6-subsets avoiding every recorded missing face that are not vertices
        pred = {frozenset(c) for c in itertools.combinations(range(n), d)
                if not any(m <= frozenset(c) for m in mf)}
        if pred - V:
            killable.append((tid, len(pred - V)))
        # (C) general position
        pts = t.get("example_points")
        if pts:
            k = collinear_triples([tuple(p) for p in pts])
            if k:
                degen.append((tid, k))
    return types, errors, killable, degen


def main():
    rc = 0
    for d, n in CASES:
        r = check(d, n)
        if r is None:
            continue
        types, errors, killable, degen = r
        print(f"=== d={d}, n={n}: {len(types)} generated types ===")

        print(f"  (A) errors (representations contradict each other): {len(errors)}")
        for tid, kind, detail in errors:
            print(f"        tid {tid}: {kind}: {detail}")
            rc = 1

        print(f"  (B) types with an unrecorded minimal non-face of size {d}: "
              f"{len(killable)}")
        if killable:
            surv_p = ROOT / f"runs/d{d}_n{n}/facet_profile_survivors.json"
            surv = set(json.load(open(surv_p))["survivors"]) if surv_p.exists() else set()
            ks = {t for t, _ in killable}
            print(f"        these cannot be compact Coxeter polytopes: a Lannér")
            print(f"        subdiagram of order {d} would be required, and none exists.")
            if surv:
                print(f"        of the {len(surv)} facet-filter survivors, "
                      f"{len(ks & surv)} are of this kind:")
                print(f"          {sorted(ks & surv)}")

        print(f"  (C) types whose stored coordinates are not in general position: "
              f"{len(degen)}")
        for tid, k in degen[:5]:
            print(f"        tid {tid}: {k} collinear triple(s)")
            rc = 1
        print()

    print("RESULT:", "no contradictions between the two representations"
          if rc == 0 else "ERRORS FOUND (see (A)/(C) above)")
    return rc


if __name__ == "__main__":
    sys.exit(main())
