#!/usr/bin/env python3
"""Facet-realizability filter for the d=6 / 10-facet census (sound, combinatorial).

Ma-Zheng (arXiv:2203.16049, p.20 + Tables 13/14) prove that EXACTLY six simple
5-polytopes with 9 facets admit a compact hyperbolic Coxeter structure.  Their
missing-face profiles (all have 9 missing faces):

    (2,2,2,2,2,2,5,5,5)  (2,2,2,2,2,3,4,5,5)  (2,2,2,2,3,3,4,4,5)  (2,2,2,3,3,3,4,4,4)

Every facet of a compact Coxeter 6-polytope is itself a compact Coxeter 5-polytope
(CLAUDE.md §2.5).  A facet F_i with 9 neighbors is a 5-polytope with 9 facets, so its
missing-face profile must be one of the four above — else the d=6 type cannot realize.

Combinatorics (simple polytope): facet i's polytope has vertex sets
{v \\ {i} : v in vertex_sets(P), i in v}; a set S of its facets is a face iff S is
contained in some vertex set; missing faces = minimal non-faces.  Facet i has 9
neighbors iff i occurs in no size-2 missing face of P (no dotted partner).

This regenerates (as a committed, reproducible artifact) the 54-survivor list that
run_survivors_fulllabel.py / run_survivors_rigorous.py hardcode.
"""
import json, sys
from itertools import combinations
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline.stage4_blockpaste import _setup

REALIZING_N9_PROFILES = {
    (2, 2, 2, 2, 2, 2, 5, 5, 5),
    (2, 2, 2, 2, 2, 3, 4, 5, 5),
    (2, 2, 2, 2, 3, 3, 4, 4, 5),
    (2, 2, 2, 3, 3, 3, 4, 4, 4),
}


def missing_face_profile(vertex_sets, facets):
    """Sorted sizes of minimal non-faces of the polytope with the given
    vertex-facet incidences (vertex_sets: list of frozensets over `facets`)."""
    vsets = [frozenset(v) for v in vertex_sets]

    def is_face(S):
        return any(S <= v for v in vsets)

    missing = []
    for size in range(2, len(facets) + 1):
        for S in combinations(facets, size):
            fs = frozenset(S)
            if is_face(fs):
                continue
            # minimal: every proper subset is a face
            if all(is_face(fs - {x}) for x in fs):
                missing.append(fs)
    return tuple(sorted(len(m) for m in missing))


def facet_polytope_profile(P_vertex_sets, i):
    """Missing-face profile of facet i's polytope (i assumed to meet all others)."""
    sub = [frozenset(v) - {i} for v in P_vertex_sets if i in v]
    facets = sorted(set().union(*sub)) if sub else []
    return missing_face_profile(sub, facets)


def main():
    types = json.load(open('runs/d6_n10/stage2/types.json'))
    survivors, killed = [], []
    for t in types:
        tid = t['type_id']
        V, *_ = _setup(dict(t))
        vsets = [frozenset(v) for v in V]
        in_dotted = set()
        for m in t['missing_faces']:
            if len(m) == 2:
                in_dotted.update(m)
        n9_facets = [i for i in range(10) if i not in in_dotted]
        bad = None
        for i in n9_facets:
            prof = facet_polytope_profile(vsets, i)
            if prof not in REALIZING_N9_PROFILES:
                bad = (i, prof)
                break
        if bad is None:
            survivors.append(tid)
        else:
            killed.append((tid, bad))
        if tid % 50 == 0:
            print(f"  ...through tid {tid}", flush=True)

    print(f"\nsurvivors ({len(survivors)}): {survivors}")
    print(f"killed: {len(killed)}")
    Path('runs/d6_n10/facet_profile_survivors.json').write_text(
        json.dumps({"survivors": survivors,
                    "killed": {str(t): [b[0], list(b[1])] for t, b in killed}}, indent=1))

    # cross-check against the list the survivor drivers hardcode
    SURV = [8, 12, 17, 34, 36, 38, 40, 51, 55, 59, 60, 61, 69, 70, 92, 103, 120, 127,
            132, 140, 154, 159, 162, 168, 173, 206, 214, 218, 220, 229, 234, 239, 255,
            265, 273, 284, 286, 287, 295, 297, 308, 315, 317, 320, 329, 332, 344, 352,
            354, 356, 360, 378, 379, 382]
    match = sorted(survivors) == sorted(SURV)
    print(f"matches hardcoded 54-list: {match}")
    if not match:
        print("  only-here:", sorted(set(survivors) - set(SURV)))
        print("  only-hardcoded:", sorted(set(SURV) - set(survivors)))


if __name__ == '__main__':
    main()
