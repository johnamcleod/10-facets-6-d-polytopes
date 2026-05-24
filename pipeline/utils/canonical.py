"""Canonical forms for missing-face hypergraphs.

Fast approach: use a structural invariant (degree sequence + local structure)
to assign canonical node labels, then express the missing faces with those labels.

Two missing-face hypergraphs are isomorphic iff they have the same canonical form.
Our canonical form is exact for most practical cases (small n, few missing faces).
For cases where the fast form isn't discriminating enough, an exact n!-search
is available but only run on the small set of surviving types.
"""

from itertools import permutations


def _node_invariant(node, mf_frozen, n):
    """Compute a structural invariant for node in the hypergraph.

    Returns a hashable tuple that is the same for isomorphic nodes.
    """
    # Faces containing this node, sorted by size
    containing = sorted(len(m) for m in mf_frozen if node in m)
    # For each face containing this node, the sorted sizes of co-faces
    # (other faces that share at least one other node)
    coface_structure = []
    for m in mf_frozen:
        if node not in m:
            continue
        # Other faces that share a neighbor with m
        neighbors = m - {node}
        neighbor_faces = tuple(sorted(
            len(f) for f in mf_frozen if f != m and f & neighbors
        ))
        coface_structure.append((len(m), neighbor_faces))
    coface_structure.sort()
    return (tuple(containing), tuple(coface_structure))


def canonical_missing_face_hypergraph(mf_list, n, exact=False):
    """Compute canonical form of a missing-face hypergraph on n nodes.

    Algorithm:
    1. Compute a structural invariant for each node.
    2. Sort nodes by invariant (ties broken by index for stability).
    3. Assign new labels 0..n-1 based on this sorted order.
    4. Express all missing faces in new labels and sort.

    This is correct when all node invariants are distinct (the common case).
    When there are ties, it's a heuristic — good enough for our dedup needs.
    exact=True falls back to restricted permutation search within tied classes.

    Returns a canonical tuple (used as deduplication key).
    """
    mf_frozen = [frozenset(m) for m in mf_list]

    # Compute invariant for each node
    inv = {i: _node_invariant(i, mf_frozen, n) for i in range(n)}

    # Sort nodes by invariant (deterministic: break ties by index)
    sorted_nodes = sorted(range(n), key=lambda i: (inv[i], i))
    # perm[old_index] = new_label
    perm = {sorted_nodes[new_idx]: new_idx for new_idx in range(n)}

    # Apply permutation
    canonical = tuple(sorted(
        tuple(sorted(perm[i] for i in m))
        for m in mf_frozen
    ))

    if not exact:
        return canonical

    # Exact: search within equivalence classes of tied nodes
    from collections import defaultdict
    classes = defaultdict(list)
    for i in range(n):
        classes[inv[i]].append(i)

    # If all invariants unique, we're done
    if all(len(v) == 1 for v in classes.values()):
        return canonical

    # Try permutations within each tied class
    tied_classes = [sorted(v) for v in classes.values() if len(v) > 1]
    from itertools import product as iproduct
    best = canonical
    for class_perms in iproduct(*[permutations(cls) for cls in tied_classes]):
        # Build full perm: fixed nodes from singleton classes + permuted tied nodes
        full_perm = dict(perm)  # start with the default perm
        for cls, cls_perm in zip(tied_classes, class_perms):
            # old nodes in cls get remapped
            for old_node, new_node in zip(cls, cls_perm):
                full_perm[old_node] = perm[new_node]
        relabeled = tuple(sorted(
            tuple(sorted(full_perm[i] for i in m))
            for m in mf_frozen
        ))
        if relabeled < best:
            best = relabeled

    return best
