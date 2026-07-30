#!/usr/bin/env python3
"""Reconcile our d=4 count against the published census, type by type.

Our pipeline finds 338 compact hyperbolic Coxeter 4-polytopes with 8 facets; the
published number is 348.  This script localises the difference by matching our
combinatorial types to Burcroff's G-types (exact missing-face hypergraph
isomorphism) and comparing per-type counts.

The published per-type counts are transcribed below from the first table of
[Burcroff 2024, Appendix A] ("List of Combinatorial Types", pages 48-50 of
arXiv:2201.03437), whose 34 rows sum to exactly 348 and so agree with her
abstract and with [Ma-Zheng 2024, Thm 1.1].

Result: the whole 10-polytope difference sits in exactly two types, and 27 types
agree exactly, including the three largest (130, 115, 49).

Run:  python3 paper/checks/d4_census_reconciliation.py
"""
import glob
import itertools
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
N = 8

# [Burcroff 2024, Appendix A], first table.  Only the realizing types are needed
# here; the other 20 rows carry 0 and are omitted.  G4's twelve are attributed in
# her text to Jacquemet and Tschantz -- they are the compact hyperbolic Coxeter
# 4-cubes, and G4's missing-face list (four disjoint pairs, nothing else) is
# exactly the 3-free / cube type.
BURCROFF = {
    "G1":  ("01, 02, 03, 12, 14, 25, 3467, 3567, 4567", 130),
    "G2":  ("01, 02, 03, 12, 14, 2567, 34, 3567, 4567", 115),
    "G3":  ("01, 02, 03, 1245, 1456, 27, 36, 37, 4567", 49),
    "G4":  ("01, 23, 45, 67", 12),
    "G5":  ("01, 02, 13, 24, 34, 567", 3),
    "G6":  ("01, 02, 03, 124, 145, 26, 357, 367, 4567", 2),
    "G7":  ("01, 02, 03, 12, 14, 256, 347, 3567, 4567", 2),
    "G8":  ("01, 02, 03, 1245, 146, 257, 36, 37, 4567", 1),
    "G9":  ("01, 02, 03, 14, 25, 367, 4567", 15),
    "G10": ("01, 02, 034, 134, 156, 256, 27, 347, 567", 1),
    "G11": ("01, 02, 13, 245, 345, 67", 8),
    "G12": ("01, 02, 034, 134, 15, 256, 267, 347, 567", 4),
    "G13": ("01, 02, 034, 15, 26, 347, 567", 4),
    "G14": ("01, 023, 145, 236, 237, 46, 57", 2),
}
PUBLISHED_TOTAL = 348


def parse(mf):
    return [frozenset(int(c) for c in tok.strip()) for tok in mf.split(",")]


def canon(faces, n=N):
    """Canonical form of a missing-face hypergraph under relabelling of facets."""
    deg = [0] * n
    for f in faces:
        for x in f:
            deg[x] += 1
    sig = [(deg[i], tuple(sorted(len(f) for f in faces if i in f))) for i in range(n)]
    groups = {}
    for i, s in enumerate(sig):
        groups.setdefault(s, []).append(i)
    slots = [groups[k] for k in sorted(groups)]
    best = [None]

    def rec(k, mp):
        if k == len(slots):
            enc = tuple(sorted(tuple(sorted(mp[x] for x in f)) for f in faces))
            if best[0] is None or enc < best[0]:
                best[0] = enc
            return
        base = sum(len(slots[j]) for j in range(k))
        for perm in itertools.permutations(slots[k]):
            m2 = dict(mp)
            for idx, node in enumerate(perm):
                m2[node] = base + idx
            rec(k + 1, m2)

    rec(0, {})
    return best[0]


def main():
    # The d=4 census of record is runs/d4_n8/survivors_discfix, produced on the
    # code with all three corrections of paper section 7.3.  The earlier
    # blockpaste/ and survivors_final/ directories are read too, so that the
    # per-type effect of each correction is visible rather than asserted.
    import collections

    def from_jsonl(d):
        p = ROOT / f"runs/d4_n8/{d}/state.jsonl"
        if not p.exists():
            return {}
        k = collections.defaultdict(set)
        for line in p.read_text().splitlines():
            line = line.strip()
            if line:
                r = json.loads(line)
                k[int(r["key"].split("|")[0])].update(r.get("keys", []))
        return {t: len(v) for t, v in k.items()}

    ours_counts = from_jsonl("survivors_discfix")
    prior = {}
    for f in glob.glob(str(ROOT / "runs/d4_n8/blockpaste/type_*.json")):
        tid = int(os.path.basename(f)[5:-5])
        x = json.load(open(f))
        prior[tid] = (x.get("distinct") if isinstance(x, dict) and "distinct" in x
                      else (len(x.get("keys", [])) if isinstance(x, dict)
                            else len(x))) or 0
    mid = from_jsonl("survivors_final")
    if not ours_counts:
        print("runs/d4_n8/survivors_discfix not present -- nothing to reconcile")
        return 0

    types = {t["type_id"]: t for t in
             json.load(open(ROOT / "runs/d4_n8/stage2/types.json"))}
    ourc = {tid: canon([frozenset(m) for m in t["missing_faces"]])
            for tid, t in types.items()}
    burc = {g: canon(parse(mf)) for g, (mf, _) in BURCROFF.items()}
    rev = {v: k for k, v in burc.items()}

    print(f"published per-type counts sum to "
          f"{sum(c for _, c in BURCROFF.values())} (expected {PUBLISHED_TOTAL})")
    print(f"our per-type counts sum to {sum(ours_counts.values())}\n")
    print(f"{'our tid':>8} {'G-type':>7} {'ours':>6} {'published':>10} {'delta':>6}")

    tot_o = tot_b = 0
    diffs = []
    for tid in sorted(ourc):
        g = rev.get(ourc[tid])
        if g is None:
            continue
        o, b = ours_counts.get(tid, 0), BURCROFF[g][1]
        tot_o += o
        tot_b += b
        if o or b:
            mark = "   <-- DIFFERS" if o != b else ""
            print(f"{tid:>8} {g:>7} {o:>6} {b:>10} {b - o:>6}{mark}")
            if o != b:
                diffs.append((tid, g, o, b))

    print(f"\nmatched types: ours {tot_o}, published {tot_b}, gap {tot_b - tot_o}")
    print("\nthe gap is confined to:")
    for tid, g, o, b in diffs:
        note = ""
        sizes = {len(m) for m in types[tid]["missing_faces"]}
        if sizes == {2}:
            note = ("  <-- 3-free: this is the 4-CUBE type. Its published count is "
                    "attributed\n         to Jacquemet-Tschantz's classification of "
                    "hyperbolic Coxeter cubes,\n         not to generic label "
                    "enumeration. It has NO d=6 analogue: a 3-free\n         "
                    "compact Coxeter 6-polytope would need >= 12 facets (Esselmann).")
        print(f"   our tid {tid} = {g}: {o} vs {b}, short by {b - o}{note}")

    if prior or mid:
        print("\nper-type effect of the three corrections "
              "(block-pasting -> +pair/refine -> +tangency -> published):")
        for tid in sorted(ourc):
            g = rev.get(ourc[tid])
            if g is None:
                continue
            row = (prior.get(tid, 0), mid.get(tid, 0), ours_counts.get(tid, 0),
                   BURCROFF[g][1])
            if len(set(row)) > 1:
                print(f"   tid {tid:>3} = {g:<4} {row[0]:>4} -> {row[1]:>4} -> "
                      f"{row[2]:>4}   published {row[3]}")
        print(f"   totals: {sum(prior.values()) or '-'} -> {sum(mid.values()) or '-'}"
              f" -> {sum(ours_counts.values())}   published {PUBLISHED_TOTAL}")

    ok = (tot_b == PUBLISHED_TOTAL) and (tot_o == PUBLISHED_TOTAL) and not diffs
    msg = ("our census agrees with the published one type by type"
           if ok else f"MISMATCH -- {len(diffs)} type(s) differ")
    print(f"\nRESULT: {msg}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
