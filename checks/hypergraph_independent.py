#!/usr/bin/env python3
"""Recompute every type's missing-face hypergraph from its stored coordinates.

Independent of the pipeline: this script imports no
pipeline code.  For each of the 387 generated types it takes the stored integer
realization of its order type (`example_points`, ten points in the plane) and the
two positive points (`example_positive`), and:

  1. checks the conditions of Burcroff's Theorem 3.5: the ten points are in general
     position (no three collinear) and both positive points lie strictly inside the
     convex hull of the other eight;
  2. decides, for every subset S of the ten facets, whether S is a face, by the
     affine Gale criterion: S is a face iff the relative interiors of conv(T+) and
     conv(T-) meet, T the complement of S.  All tests are exact: orientation
     determinants of integer points, and rational parameters when an open segment is
     clipped against an open polygon;
  3. computes ALL minimal non-faces, of every size;
  4. compares with the pipeline: the recorded missing faces must be exactly the
     computed ones of size at most 5 (the generator records no larger ones), and the
     type must have a minimal non-face of size 6 iff runs/d6_n10/type_flags.json
     says so; no minimal non-face may have size 7 or more;
  5. checks that the 387 full hypergraphs are pairwise non-isomorphic, so the list
     has no duplicates.

What this does NOT establish is completeness of the type list, which rests on the
order-type database and on Burcroff's Theorem 3.5, as the paper states.

Run:  python3 checks/hypergraph_independent.py
Artifact: runs/d6_n10/hypergraph_independent.json
"""
from __future__ import annotations

import itertools
import json
import sys
from fractions import Fraction
from pathlib import Path

import networkx as nx
from networkx.algorithms import isomorphism as iso

ROOT = Path(__file__).resolve().parents[1]
N = 10


def orient(a, b, c):
    """Twice the signed area of triangle abc (exact, integers)."""
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def hull(pts):
    """Convex hull, counter-clockwise, of points in general position."""
    pts = sorted(pts)
    if len(pts) <= 2:
        return pts
    lower, upper = [], []
    for p in pts:
        while len(lower) >= 2 and orient(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    for p in reversed(pts):
        while len(upper) >= 2 and orient(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def strictly_inside(p, poly):
    """p in the interior of the convex polygon poly (counter-clockwise, >= 3)."""
    return all(orient(poly[i], poly[(i + 1) % len(poly)], p) > 0
               for i in range(len(poly)))


def open_segment_meets_open_polygon(u, v, poly):
    """Some point u + t(v - u), 0 < t < 1, strictly inside poly."""
    lo, hi = Fraction(0), Fraction(1)
    for i in range(len(poly)):
        a, b = poly[i], poly[(i + 1) % len(poly)]
        c0, c1 = orient(a, b, u), orient(a, b, v)
        slope = c1 - c0                       # f(t) = c0 + t * slope must be > 0
        if slope == 0:
            if c0 <= 0:
                return False
        elif slope > 0:
            lo = max(lo, Fraction(-c0, slope))
        else:
            hi = min(hi, Fraction(-c0, slope))
        if lo >= hi:
            return False
    return lo < hi


def relints_meet(A, B):
    """Do the relative interiors of conv(A) and conv(B) meet?  A has at most two
    points (the positive ones); all points are in general position and distinct."""
    if not A or not B:
        return False
    if len(A) == 1:
        p = A[0]
        if len(B) <= 2:
            return False                      # a point on a point or open segment
        return strictly_inside(p, hull(B))    # would need three collinear points
    u, v = A
    if len(B) == 1:
        return False
    if len(B) == 2:
        a, b = B
        return (orient(u, v, a) * orient(u, v, b) < 0 and
                orient(a, b, u) * orient(a, b, v) < 0)
    return open_segment_meets_open_polygon(u, v, hull(B))


def hypergraph(t):
    pts = [tuple(p) for p in t["example_points"]]
    pos = set(t["example_positive"])
    assert len(pts) == N and len(pos) == 2
    # 1. Theorem 3.5's conditions
    assert all(orient(a, b, c) != 0 for a, b, c in itertools.combinations(pts, 3)), \
        "three collinear points"
    negs = [pts[i] for i in range(N) if i not in pos]
    H = hull(negs)
    assert all(strictly_inside(pts[i], H) for i in pos), "a positive point not interior"
    # 2. faces
    face = {}
    for r in range(N + 1):
        for S in itertools.combinations(range(N), r):
            T = [i for i in range(N) if i not in S]
            A = [pts[i] for i in T if i in pos]
            B = [pts[i] for i in T if i not in pos]
            face[frozenset(S)] = relints_meet(A, B)
    # sanity: faces form a simplicial complex (closed under subsets)
    for S, f in face.items():
        if f:
            assert all(face[S - {i}] for i in S), "faces not closed under subsets"
    # 3. minimal non-faces
    mnf = [S for S, f in face.items()
           if not f and all(face[S - {i}] for i in S)]
    vertices = [S for S, f in face.items() if f and len(S) == 6]
    return mnf, vertices


def main():
    types = json.load(open(ROOT / "runs/d6_n10/stage2/types.json"))
    flags = json.load(open(ROOT / "runs/d6_n10/type_flags.json"))["d6_n10"]
    bad, full, big = [], {}, 0
    for t in types:
        tid = t["type_id"]
        mnf, V = hypergraph(t)
        sizes = [len(S) for S in mnf]
        assert max(sizes) <= 6, f"type {tid}: a minimal non-face of size {max(sizes)}"
        rec = {frozenset(m) for m in t["missing_faces"]}
        small = {S for S in mnf if len(S) <= 5}
        has6 = any(len(S) == 6 for S in mnf)
        big += has6
        ok = (small == rec and has6 == flags[str(tid)]["has_missing_face_of_size_d"])
        if not ok:
            bad.append(tid)
        full[tid] = mnf
    print(f"{len(types)} types recomputed from their coordinates: "
          f"{len(types) - len(bad)} agree with the pipeline, disagreements: {bad or 'none'}")
    print(f"types with a minimal non-face of size 6: {big}")

    # 5. pairwise non-isomorphic
    def inc(h):
        g = nx.Graph()
        g.add_nodes_from((("v", i) for i in range(N)), side=0)
        for k, m in enumerate(h):
            g.add_node(("e", k), side=1)
            g.add_edges_from((("v", i), ("e", k)) for i in m)
        return g
    graphs = {tid: inc(h) for tid, h in full.items()}
    hashes = {tid: nx.weisfeiler_lehman_graph_hash(g, node_attr=None, iterations=4)
              + str(sorted(len(m) for m in full[tid])) for tid, g in graphs.items()}
    groups = {}
    for tid, h in hashes.items():
        groups.setdefault(h, []).append(tid)
    dup = []
    nm = iso.categorical_node_match("side", None)
    for g in groups.values():
        for a, b in itertools.combinations(g, 2):
            if nx.is_isomorphic(graphs[a], graphs[b], node_match=nm):
                dup.append((a, b))
    print(f"pairwise isomorphic types: {dup or 'none'}")
    art = ROOT / "runs/d6_n10/hypergraph_independent.json"
    art.write_text(json.dumps({"types": len(types), "agree": len(types) - len(bad),
                               "disagree": bad, "size6_types": big,
                               "isomorphic_pairs": dup}, indent=1) + "\n")
    print(f"artifact -> {art.relative_to(ROOT)}")
    return 0 if not bad and not dup else 1


if __name__ == "__main__":
    sys.exit(main())
