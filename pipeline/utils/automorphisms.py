"""Combinatorial automorphism group of a facet-vertex incidence structure, and
sound symmetry-breaking support for the Stage-4 label backtracker.

Aut(type) = permutations of the n facet indices that map the vertex-set list
to itself (equivalently: preserve the missing-face structure, since faces are
exactly the subsets of some vertex's facet-set).  Computed via VF2 graph
automorphism on the bipartite facet/vertex incidence graph (networkx).

Soundness of the pruning below (see run_survivors_rigorous session notes):
for any sigma in Aut(type) whose induced permutation of ordinary-pair
POSITIONS (under a fixed enumeration order ``pairs_list``) maps the prefix
domain {0,...,k-1} to itself (as a set), comparing the length-k assignment
vector against its sigma-image is a genuine same-orbit comparison — if the
image is lexicographically smaller, the current prefix (and everything
below it in the search) is dominated by an equivalent branch elsewhere in
the same search and can be pruned with no loss of solutions.  This holds
for ARBITRARY sigma (not just single transpositions) because it is a full
vector comparison, not an independent per-pair inequality — the latter can
silently drop whole orbits for compound permutations, which is why this
module only ever compares full prefix vectors.
"""
from __future__ import annotations

import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher


def compute_aut_group(vertex_sets, n):
    """Aut(type): permutations (tuples of length n) of facet indices 0..n-1
    that map ``vertex_sets`` (list of frozensets/tuples of facet indices) to
    itself.  Always includes the identity."""
    G = nx.Graph()
    for i in range(n):
        G.add_node(("f", i), bipartite=0)
    for vi, v in enumerate(vertex_sets):
        G.add_node(("v", vi), bipartite=1)
        for f in v:
            G.add_edge(("f", f), ("v", vi))
    gm = GraphMatcher(G, G, node_match=lambda a, b: a["bipartite"] == b["bipartite"])
    perms = set()
    for iso in gm.isomorphisms_iter():
        perms.add(tuple(iso[("f", i)][1] for i in range(n)))
    return sorted(perms)


def orbit_aligned_order(base_order, aut_group):
    """Reorder ``base_order`` (a list of ordinary facet-pairs) so that every
    orbit of Aut(type) acting on the pairs is contiguous — this guarantees
    that at every orbit-block boundary depth, the FULL group stabilizes the
    prefix domain (orbits are always globally invariant sets), maximising
    where the runtime canonicity check in ``symmetry_prefix_checks`` can
    actually fire.  Orbit blocks are ordered by the earliest position their
    members occupy in ``base_order`` (preserves as much of the caller's
    heuristic ordering as possible); pairs within a block keep their
    relative base_order position.
    """
    pos = {p: i for i, p in enumerate(base_order)}
    seen = set()
    orbits = []
    for p in base_order:
        if p in seen:
            continue
        orb = set()
        for sigma in aut_group:
            q = tuple(sorted((sigma[p[0]], sigma[p[1]])))
            orb.add(q)
        orbits.append(orb)
        seen |= orb
    orbits.sort(key=lambda orb: min(pos[p] for p in orb))
    new_order = []
    for orb in orbits:
        new_order.extend(sorted(orb, key=lambda p: pos[p]))
    return new_order


def symmetry_prefix_checks(pairs_list, aut_group):
    """Precompute, for every prefix length k = 1..len(pairs_list), the list of
    inverse-position arrays (each a tuple of length k) for every non-identity
    sigma in aut_group whose induced pair-position permutation stabilizes the
    prefix domain {0,...,k-1} as a set.  Depends only on structure (pairs_list
    order + the group), never on label values — computed once per type/order
    and reused across every node of the backtracking search.

    Returned as a list ``checks`` of length len(pairs_list)+1; checks[k] is a
    (possibly empty) list of tuples ``inv`` of length k such that, for an
    assignment vector v (label-index values) of length >= k, the sigma-image
    of v restricted to the first k positions is ``[v[j] for j in inv]``.
    """
    n_pairs = len(pairs_list)
    pos = {p: i for i, p in enumerate(pairs_list)}
    pis = []
    for sigma in aut_group:
        pi = [None] * n_pairs
        ok = True
        for m, p in enumerate(pairs_list):
            q = tuple(sorted((sigma[p[0]], sigma[p[1]])))
            if q not in pos:
                ok = False
                break
            pi[m] = pos[q]
        if ok and pi != list(range(n_pairs)):
            pis.append(pi)

    checks = [None] * (n_pairs + 1)
    checks[0] = []
    for k in range(1, n_pairs + 1):
        dom = set(range(k))
        lst = []
        for pi in pis:
            if set(pi[:k]) == dom:
                inv = [0] * k
                for m in range(k):
                    inv[pi[m]] = m
                lst.append(tuple(inv))
        checks[k] = lst
    return checks
