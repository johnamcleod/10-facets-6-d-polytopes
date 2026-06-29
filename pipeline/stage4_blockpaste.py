"""Block-pasting Stage-4 candidate generator (Ma-Zheng method), d=4 / n=8.

Replaces the brute-force integer-label enumeration of ``stage4_gram.process_type_stage4``
for the *candidate generation* step.  Instead of enumerating every label assignment over
the ordinary edges (~10^24 for low-k types) and screening each, we build the Gram matrix
VERTEX BY VERTEX and prune any partial candidate whose sub-angle tuples violate the
precomputed spherical / Euclidean / Lanner libraries (``pipeline/utils/mazheng_lib``).
This keeps the candidate set bounded and makes per-type enumeration
exhaustive-by-construction with no timeout wall.

The constraint logic is a faithful port of Ma-Zheng's ``scratchpad/HCPdm/pyFile/chcp48.py``
(``run48``):

  SEED   each vertex's 6 ordinary angles from the spherical-4 library (``S4``); a vertex is
         therefore spherical by construction.
  KILL   (drop rows whose sub-tuple IS in the library):
           * s3  -- a non-face triple must NOT be spherical (else it would force a 2-face).
           * s4  -- a non-vertex 4-subset must NOT be spherical (else a phantom vertex).
           * e3,e4,e5,e6,e7 -- no present sub-config may be Euclidean (compact => no
             parabolic subdiagram).
           * se5,se6,se7 -- no rank-5+ sub-config may be spherical (impossible in d=4).
           * i4 -- a 4-subset with a "matching" pair of dotted edges may not have all four
             cross angles = pi/2 (degenerate/reducible).
  SAVE   (keep only rows whose sub-tuple IS in the library):
           * l4 -- every size-4 missing face must be a Lanner tetrahedron.
           * l4_basis -- the ridges around a simplex (tetrahedral) facet must be pi/2
             (prism-end orthogonality).

Size-3 missing faces need no explicit Lanner saver: ``s3`` (not spherical) + ``e3`` (not
Euclidean) together force a rank-3 non-face to be Lanner.

Alphabet is ``{2,3,4,5,6,7}`` with ``7`` a wildcard for "label >= 7"; the join stays in
this alphabet.  Concrete labels >= 7 are introduced (7 -> {7,8,9,10,12}) only at the solver
hand-off and resolved by the exact signature solve.  We do NOT apply Ma-Zheng's symmetry
quotient -- dedup happens post-solve via ``run_d4.canonical_key``.
"""
from __future__ import annotations

import itertools

import numpy as np
import pandas as pd

from pipeline.utils import mazheng_lib as ml

N, D = 8, 4
_UNSET = 0  # sentinel for an unassigned ordinary column
USE_L4_BASIS = False  # Ma-Zheng prism-end orthogonality saver (over-prunes -- see below)


# ---------------------------------------------------------------------------
# uint64 encoding of small-alphabet tuples (values 0..7 -> 3 bits each).
# width <= 21 fits in uint64 (21*3 = 63 bits).
# ---------------------------------------------------------------------------
def _encode(arr: np.ndarray) -> np.ndarray:
    """Encode each row of ``arr`` (values 0..7) as a uint64 via 3-bit packing."""
    arr = np.ascontiguousarray(arr, dtype=np.uint64)
    out = np.zeros(arr.shape[0], dtype=np.uint64)
    for p in range(arr.shape[1]):
        out |= arr[:, p] << np.uint64(3 * p)
    return out


def _lib_codes(lib: frozenset, width: int) -> np.ndarray:
    rows = np.array(sorted(lib), dtype=np.uint64)
    if rows.shape[1] != width:
        raise ValueError(f"library width {rows.shape[1]} != {width}")
    return _encode(rows)


# ---------------------------------------------------------------------------
# Per-type setup: ordinary columns, vertices, constraints.
# ---------------------------------------------------------------------------
def _setup(t):
    V = [tuple(sorted(v)) for v in t["vertex_sets"]]
    missing = [tuple(sorted(m)) for m in t["missing_faces"]]
    dotted = {m for m in missing if len(m) == 2}
    mf4 = [m for m in missing if len(m) == 4]

    allpairs = [(i, j) for i in range(N) for j in range(i + 1, N)]
    ordinary = [p for p in allpairs if p not in dotted]
    op_idx = {p: c for c, p in enumerate(ordinary)}
    Nr = len(ordinary)

    def cols_of(subset):
        """Ordinary column indices for ``subset``'s pairs in combinations order, or None
        if any pair is dotted (subset straddles a missing face -> excluded)."""
        pcs = ml.facet_pair_columns(subset)
        if any(p in dotted for p in pcs):
            return None
        return [op_idx[p] for p in pcs]

    # every vertex must be all-ordinary (compact => simple => no missing face inside a vertex)
    for v in V:
        assert cols_of(v) is not None, f"vertex {v} contains a dotted pair (malformed type)"

    return V, missing, dotted, mf4, ordinary, op_idx, Nr, cols_of


def _faces(V):
    """The set of 2-faces (mutually-meeting facet triples) = triples appearing in a vertex."""
    edges = set()
    for v in V:
        for tri in itertools.combinations(v, 3):
            edges.add(tri)
    return edges


def _build_constraints(t):
    """Return (constraints, l4_basis_cols, ctx) where each constraint is
    (cols_tuple, lib_codes, keep_in: bool).  ``keep_in`` True = saver (keep rows IN lib),
    False = killer (keep rows NOT in lib)."""
    V, missing, dotted, mf4, ordinary, op_idx, Nr, cols_of = _setup(t)
    faces = _faces(V)
    vert_set = set(V)

    cons = []  # (tuple(cols), lib_codes uint64 array, keep_in)

    def add(subset, lib, keep_in):
        cols = cols_of(subset)
        if cols is None:
            return
        width = len(cols)
        cons.append((tuple(cols), _lib_codes(lib, width), keep_in))

    all4 = list(itertools.combinations(range(N), 4))
    all3 = list(itertools.combinations(range(N), 3))

    # killers: non-vertex 4-subsets not spherical (s4); all 4-subsets not Euclidean (e4)
    for q in all4:
        if q not in vert_set:
            add(q, ml.S(4), keep_in=False)      # s4
        add(q, ml.E(4), keep_in=False)          # e4
    # killers: non-face triples not spherical (s3); all triples not Euclidean (e3)
    for tri in all3:
        if tri not in faces:
            add(tri, ml.S(3), keep_in=False)    # s3
        add(tri, ml.E(3), keep_in=False)        # e3
    # killers: rank-5/6/7 sub-configs neither spherical nor Euclidean (se5/6/7)
    for k in (5, 6, 7):
        for sub in itertools.combinations(range(N), k):
            add(sub, ml.S(k), keep_in=False)
            add(sub, ml.E(k), keep_in=False)
    # savers: size-4 missing faces are Lanner (l4)
    for m in mf4:
        add(m, ml.L4(), keep_in=True)

    # i4 killer: 4-subset with a "matching" pair of dotted edges -> the 4 cross angles
    # may not be all == 2.  Port of chcp48.infty2 (positions in combinations order:
    # [0]=(a,b),[1]=(a,c),[2]=(a,d),[3]=(b,c),[4]=(b,d),[5]=(c,d); opposite pairs are
    # {0,5},{1,4},{2,3}).
    OPP = [((0, 5), (1, 2, 3, 4)), ((1, 4), (0, 2, 3, 5)), ((2, 3), (0, 1, 4, 5))]
    i4_codes = _lib_codes(frozenset({(2, 2, 2, 2)}), 4)
    for q in all4:
        pcs = ml.facet_pair_columns(q)  # 6 pairs in combinations order
        for opp, rest in OPP:
            if all(pcs[i] in dotted for i in opp) and all(pcs[i] not in dotted for i in rest):
                cross_cols = [op_idx[pcs[i]] for i in rest]
                cons.append((tuple(cross_cols), i4_codes, False))

    # l4_basis: Ma-Zheng force the ridges around a simplex (tetrahedral) facet to pi/2.
    # DISABLED: empirically this over-prunes -- verified compact type-6 (=P2) polytopes
    # from the brute run have tetrahedral-facet ridges (e.g. (1,3),(3,7),(5,6)) that are
    # NOT pi/2, yet pass every other library constraint.  So orthogonality is not a valid
    # universal pruning here (it is an efficiency saver in their cluster pipeline, not a
    # soundness condition).  Leaving it off keeps the candidate set complete; the exact
    # signature+isolation solver remains the final arbiter.  Set USE_L4_BASIS=True to study
    # the Ma-Zheng behaviour.
    l4_basis_cols = set()
    if USE_L4_BASIS:
        facet_nverts = [sum(1 for v in V if f in v) for f in range(N)]
        for f in range(N):
            if facet_nverts[f] != 4:
                continue
            nbrs = set()
            for v in V:
                if f in v:
                    nbrs |= set(v)
            nbrs.discard(f)
            for g in nbrs:
                p = (min(f, g), max(f, g))
                if p not in dotted:
                    l4_basis_cols.add(op_idx[p])

    ctx = dict(V=V, ordinary=ordinary, op_idx=op_idx, Nr=Nr, dotted=dotted,
               cols_of=cols_of, mf4=mf4)
    return cons, sorted(l4_basis_cols), ctx


# ---------------------------------------------------------------------------
# The block-paste join.
# ---------------------------------------------------------------------------
_S4_ARR = None


def _s4_array():
    global _S4_ARR
    if _S4_ARR is None:
        _S4_ARR = np.array(sorted(ml.S(4)), dtype=np.int8)  # (242, 6) combinations order
    return _S4_ARR


def _apply_constraints(data, cons, covered, l4_basis_cols, forced=None):
    """Apply every constraint whose columns are all covered.  Returns filtered ``data``
    and the list of still-pending constraints.  ``forced`` is an optional {col: value} map
    (partition pins + l4_basis); each is enforced as soon as its column is covered."""
    pending = []
    for cols, codes, keep_in in cons:
        if not all(c in covered for c in cols):
            pending.append((cols, codes, keep_in))
            continue
        sub = data[:, list(cols)]
        code = _encode(sub)
        inlib = np.isin(code, codes)
        data = data[inlib if keep_in else ~inlib]
        if data.shape[0] == 0:
            return data, pending
    # l4_basis: force covered basis columns to 2
    for c in l4_basis_cols:
        if c in covered:
            data = data[data[:, c] == 2]
            if data.shape[0] == 0:
                return data, pending
    # partition pins
    if forced:
        for c, v in forced.items():
            if c in covered:
                data = data[data[:, c] == v]
                if data.shape[0] == 0:
                    return data, pending
    return data, pending


def _order_vertices(V, cols_of):
    """Greedy paste order maximising overlap with the growing cover (keeps joins small)."""
    remaining = list(range(len(V)))
    order = [remaining.pop(0)]
    cover = set(cols_of(V[order[0]]))
    while remaining:
        best, best_ov = None, -1
        for idx in remaining:
            ov = len(cover & set(cols_of(V[idx])))
            if ov > best_ov:
                best, best_ov = idx, ov
        order.append(best)
        cover |= set(cols_of(V[best]))
        remaining.remove(best)
    return order


def paste_candidates(t, max_candidates=8_000_000, verbose=False, forced=None):
    """Generate all candidate ordinary-label vectors (alphabet {2..7}, 7=wildcard) for a
    combinatorial type ``t`` via block-pasting.  Returns (cands, ordinary) where ``cands``
    is an int8 array of shape (n_candidates, Nr) and ``ordinary`` the column->pair list.

    ``forced`` is an optional {column_index: value} map that pins those ordinary columns --
    used to PARTITION the search across cores (call once per value of a seed column; the
    partitions are disjoint and their union is the full set).  Raises ``OverflowError`` if
    the candidate table exceeds ``max_candidates`` (caller should partition further)."""
    cons, l4_basis_cols, ctx = _build_constraints(t)
    V, cols_of, Nr = ctx["V"], ctx["cols_of"], ctx["Nr"]
    s4 = _s4_array()
    forced = dict(forced or {})

    order = _order_vertices(V, cols_of)
    # seed from first vertex
    v0 = V[order[0]]
    c0 = cols_of(v0)
    data = np.zeros((s4.shape[0], Nr), dtype=np.int8)
    data[:, c0] = s4
    covered = set(c0)
    data, cons = _apply_constraints(data, cons, covered, l4_basis_cols, forced)
    if verbose:
        print(f"  seed v{order[0]}: {data.shape[0]} rows")

    for step, vi in enumerate(order[1:], 1):
        vcols = cols_of(V[vi])
        link = [c for c in vcols if c in covered]
        new = [c for c in vcols if c not in covered]
        # map link/new ordinary columns back to S4 library column positions
        link_k = [vcols.index(c) for c in link]
        new_k = [vcols.index(c) for c in new]

        if link:
            ld = _encode(data[:, link])
            ls = _encode(s4[:, link_k])
        else:  # no shared columns -> full cartesian product
            ld = np.zeros(data.shape[0], dtype=np.uint64)
            ls = np.zeros(s4.shape[0], dtype=np.uint64)
        dfd = pd.DataFrame({"key": ld, "ridx": np.arange(data.shape[0])})
        dfs = pd.DataFrame({"key": ls, "sidx": np.arange(s4.shape[0])})
        m = dfd.merge(dfs, on="key", sort=False)
        if m.shape[0] == 0:
            return np.zeros((0, Nr), dtype=np.int8), ctx["ordinary"]
        out = data[m["ridx"].to_numpy()]
        if new:
            out[:, new] = s4[:, new_k][m["sidx"].to_numpy()]
        data = out
        covered |= set(vcols)
        data, cons = _apply_constraints(data, cons, covered, l4_basis_cols, forced)
        # dedupe to bound memory
        data = np.unique(data, axis=0)
        if verbose:
            print(f"  paste v{vi} (step {step}): {data.shape[0]} rows")
        if data.shape[0] > max_candidates:
            raise OverflowError(f"candidate table {data.shape[0]} > {max_candidates}")

    # final safety: apply any still-pending constraints (all columns covered now)
    data, leftover = _apply_constraints(data, cons, covered, l4_basis_cols, forced)
    assert not leftover, f"{len(leftover)} constraints never covered"
    data = np.unique(data, axis=0)
    return data, ctx["ordinary"]


def paste_candidates_partitioned(t, part_col=None, max_candidates=8_000_000, verbose=False):
    """Run :func:`paste_candidates` partitioned over the 6 possible labels {2..7} of a
    single seed column ``part_col`` and concatenate the (disjoint) results.  Bounds peak
    memory to ~1/6 of the unpartitioned run.  If ``part_col`` is None, the first column of
    the greedily-ordered seed vertex is used."""
    if part_col is None:
        _, _, ctx = (None, None, None), None, None
        cons, l4b, ctx = _build_constraints(t)
        order = _order_vertices(ctx["V"], ctx["cols_of"])
        part_col = ctx["cols_of"](ctx["V"][order[0]])[0]
    parts = []
    ordinary = None
    for v in (2, 3, 4, 5, 6, 7):
        cands, ordinary = paste_candidates(t, max_candidates=max_candidates,
                                           verbose=verbose, forced={part_col: v})
        if cands.shape[0]:
            parts.append(cands)
        if verbose:
            print(f" partition col{part_col}={v}: {cands.shape[0]} rows")
    if not parts:
        return np.zeros((0, len(ordinary)), dtype=np.int8), ordinary
    return np.unique(np.concatenate(parts, axis=0), axis=0), ordinary


# ---------------------------------------------------------------------------
# Wildcard expansion at the solver hand-off.
# ---------------------------------------------------------------------------
HIGH_LABELS = (7, 8, 9, 10, 12)  # concrete values that the wildcard 7 stands for


def expand_label_assignments(cands, ordinary):
    """Yield concrete label-assignment dicts {(i,j): m} from candidate rows, expanding each
    wildcard ``7`` to every value in :data:`HIGH_LABELS`.  A row with ``w`` wildcards yields
    ``len(HIGH_LABELS)**w`` assignments."""
    ordinary = list(ordinary)
    for row in cands:
        wild = [c for c, m in enumerate(row) if m == 7]
        base = {ordinary[c]: int(m) for c, m in enumerate(row) if m != 7}
        if not wild:
            yield dict(base)
            continue
        for combo in itertools.product(HIGH_LABELS, repeat=len(wild)):
            la = dict(base)
            for c, val in zip(wild, combo):
                la[ordinary[c]] = val
            yield la
