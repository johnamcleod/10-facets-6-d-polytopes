#!/usr/bin/env python3
"""Exact, tolerance-free recognition of elliptic and Lanner Coxeter diagrams.

The forward-checking of the search asks two questions of a labelled diagram on k
nodes: is it elliptic (positive definite Gram matrix), and is it Lanner (exactly
one negative eigenvalue, every proper subdiagram elliptic).  Both were answered by
a floating-point eigenvalue test with a tolerance of 1e-8, which for the 6-node
vertex condition cannot be audited exhaustively -- there are 6^15 label tuples.

Neither question needs arithmetic.  A diagram is elliptic if and only if every
connected component is one of the finite Coxeter types, and Lanner if and only if
it appears in Lanner's classification.  Both are finite combinatorial conditions on
a labelled graph, so they can be decided exactly, with no tolerance anywhere.

Conventions.  A diagram on nodes 0..k-1 is given as a dict {(a,b): m} with a<b and
m >= 2; m = 2 means the nodes are NOT joined (orthogonal), m >= 3 means an edge of
label m.  Absent pairs are taken to be 2.

The finite (spherical) Coxeter types, following the standard classification:

    A_n  (n>=1)   path, all labels 3
    B_n  (n>=2)   path, labels 3 except one end edge 4
    D_n  (n>=4)   fork: a path with two extra nodes attached to one end, labels 3
    E_6,E_7,E_8   path of length n-1 with one node attached to the third
    F_4           path with labels 3,4,3
    H_3,H_4       path with labels 5,3(,3)
    I_2(m) (m>=5) single edge of label m  (m=3,4 are A_2,B_2)

Lanner diagrams (compact hyperbolic simplices) exist only for orders 2..5 and are
recognised here as the connected diagrams all of whose proper subdiagrams are
elliptic and which are not themselves elliptic.  That characterisation is the
definition, and with the exact elliptic test above it is itself exact.
"""
from __future__ import annotations

import itertools
from functools import lru_cache


def _components(k, edges):
    """Connected components of the graph on 0..k-1 with the given edge set."""
    adj = {i: set() for i in range(k)}
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    seen, out = set(), []
    for s in range(k):
        if s in seen:
            continue
        stack, comp = [s], []
        seen.add(s)
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in adj[u]:
                if v not in seen:
                    seen.add(v)
                    stack.append(v)
        out.append(sorted(comp))
    return out


def _shape(nodes, lab):
    """Degree sequence and label multiset of the induced subgraph, as a signature."""
    deg = {u: 0 for u in nodes}
    labels = []
    for (a, b), m in lab.items():
        if a in deg and b in deg:
            deg[a] += 1
            deg[b] += 1
            labels.append(m)
    return sorted(deg.values()), sorted(labels)


def _is_finite_component(nodes, lab):
    """Is this CONNECTED labelled component a finite Coxeter type?"""
    n = len(nodes)
    if n == 1:
        return True
    # only pairs with m >= 3 are EDGES; m == 2 means the nodes are not joined, and
    # including those was a bug: it made every tree fail the edge-count test.
    edges = {(a, b): m for (a, b), m in lab.items()
             if a in nodes and b in nodes and m >= 3}
    if n == 2:
        # I_2(m) is finite for every finite m (A_2 = I_2(3), B_2 = I_2(4))
        return len(edges) == 1
    degs, labels = _shape(nodes, edges)
    n_edges = len(edges)
    # a finite type of rank >= 3 is a tree
    if n_edges != n - 1:
        return False
    high = [m for m in labels if m >= 4]
    if len(high) > 1:
        return False
    maxdeg = degs[-1]
    if maxdeg > 3:
        return False
    n_deg3 = sum(1 for d in degs if d == 3)

    if not high:                      # all labels 3
        if maxdeg <= 2:
            return True               # A_n
        if n_deg3 == 1:
            # D_n: the branch node has two leaves attached; E_6/7/8: branch at
            # distance >= 2 from every leaf on one arm.  Both are finite for the
            # ranks that occur; distinguish by arm lengths.
            arms = _arm_lengths(nodes, edges)
            if arms is None:
                return False
            arms = sorted(arms)
            if arms[:2] == [1, 1]:
                return True                          # D_n
            if arms == [1, 2, 2] and n == 6:
                return True                          # E_6
            if arms == [1, 2, 3] and n == 7:
                return True                          # E_7
            if arms == [1, 2, 4] and n == 8:
                return True                          # E_8
            return False
        return False
    m = high[0]
    if maxdeg > 2:
        return False                  # a branch node forces all labels 3
    # a path with exactly one label m >= 4
    pos = _edge_position_on_path(nodes, edges, m)
    if pos is None:
        return False
    if m == 4:
        if pos == 0 or pos == n - 2:
            return True               # B_n
        if n == 4 and pos == 1:
            return True               # F_4
        return False
    if m == 5:
        if n in (3, 4) and (pos == 0 or pos == n - 2):
            return True               # H_3, H_4
        return False
    # m >= 6 on a path of >= 3 nodes is never finite
    return False


def _arm_lengths(nodes, edges):
    """For a tree with a single degree-3 node, the three arm lengths."""
    adj = {u: [] for u in nodes}
    for (a, b) in edges:
        adj[a].append(b)
        adj[b].append(a)
    centres = [u for u in nodes if len(adj[u]) == 3]
    if len(centres) != 1:
        return None
    c = centres[0]
    arms = []
    for nb in adj[c]:
        length, prev, cur = 1, c, nb
        while True:
            nxt = [w for w in adj[cur] if w != prev]
            if len(nxt) != 1:
                break
            prev, cur = cur, nxt[0]
            length += 1
        arms.append(length)
    return arms


def _edge_position_on_path(nodes, edges, m):
    """Index of the unique label-m edge along the path, or None if not a path."""
    adj = {u: [] for u in nodes}
    for (a, b) in edges:
        adj[a].append(b)
        adj[b].append(a)
    ends = [u for u in nodes if len(adj[u]) == 1]
    if len(ends) != 2:
        return None
    order, prev, cur = [ends[0]], None, ends[0]
    while True:
        nxt = [w for w in adj[cur] if w != prev]
        if not nxt:
            break
        prev, cur = cur, nxt[0]
        order.append(cur)
    if len(order) != len(nodes):
        return None
    for i in range(len(order) - 1):
        a, b = sorted((order[i], order[i + 1]))
        if edges.get((a, b)) == m:
            return i
    return None


def is_elliptic(k, lab):
    """Exact: is the diagram on 0..k-1 elliptic (finite)?  No arithmetic."""
    edges = {(a, b) for (a, b), m in lab.items() if m >= 3}
    for comp in _components(k, edges):
        if not _is_finite_component(set(comp), lab):
            return False
    return True


def _exact_inertia(k, lab):
    """Inertia of the Gram matrix by exact symmetric congruence elimination.

    Entries are -cos(pi/m).  For m in {2,3,4,5,6} these lie in Q(sqrt2, sqrt3,
    sqrt5) and sympy decides pivot signs exactly, so the inertia is exact with no
    tolerance.  For m >= 7 the entries are algebraic of higher degree and sympy
    cannot in general decide the sign symbolically; those cases raise
    NotImplementedError rather than fall back to floating point, so a caller can
    never mistake a numerical decision for an exact one.

    This is needed only for the Lanner test: separating hyperbolic from PARABOLIC is
    a statement about a determinant, not about graph shape.  The affine diagram
    G~_2 -- a path with labels 3 and 6 -- is connected, not elliptic, and has all
    proper subdiagrams elliptic, yet its Gram determinant is exactly zero, so it is
    parabolic and not Lanner.
    """
    import sympy as sp
    if any(m >= 7 for m in lab.values()):
        raise NotImplementedError("exact inertia not available for labels >= 7")
    M = sp.eye(k)
    for (a, b), m in lab.items():
        if m >= 3:
            v = -sp.cos(sp.pi / sp.Integer(m))
            M[a, b] = v
            M[b, a] = v
    M = sp.Matrix(sp.nsimplify(M, rational=False))
    pos = neg = zero = 0
    idx = list(range(k))
    while idx:
        # find a nonzero diagonal pivot; if none, the block is hyperbolic-degenerate
        piv = None
        for t in idx:
            if sp.simplify(M[t, t]) != 0:
                piv = t
                break
        if piv is None:
            # all diagonal entries zero: each nonzero off-diagonal contributes one
            # positive and one negative to the inertia
            rest = len(idx)
            offs = any(sp.simplify(M[a, b]) != 0
                       for a in idx for b in idx if a < b)
            if offs:
                pos += 1
                neg += 1
                zero += rest - 2
            else:
                zero += rest
            break
        d = sp.simplify(M[piv, piv])
        if d > 0:
            pos += 1
        elif d < 0:
            neg += 1
        else:
            zero += 1
        for a in idx:
            if a == piv:
                continue
            f = sp.simplify(M[a, piv] / d)
            for b in idx:
                if b == piv:
                    continue
                M[a, b] = sp.simplify(M[a, b] - f * M[piv, b])
        idx = [t for t in idx if t != piv]
    return pos, neg, zero


def is_lanner(k, lab):
    """Exact: is the diagram a Lanner diagram?

    Connected, every proper subdiagram elliptic, and signature (k-1, 1).  The last
    condition is what separates Lanner from parabolic and is decided exactly by
    _exact_inertia.
    """
    edges = {(a, b) for (a, b), m in lab.items() if m >= 3}
    if len(_components(k, edges)) != 1:
        return False
    if is_elliptic(k, lab):
        return False
    for drop in range(k):
        keep = [u for u in range(k) if u != drop]
        idx = {u: i for i, u in enumerate(keep)}
        sub = {tuple(sorted((idx[a], idx[b]))): m
               for (a, b), m in lab.items() if a in idx and b in idx}
        if not is_elliptic(k - 1, sub):
            return False
    pos, neg, zero = _exact_inertia(k, lab)
    return neg == 1 and zero == 0
