"""Rank-propagation Stage-4 candidate generation (solve labels, don't enumerate them).

The block-paste generator enumerates integer dihedral labels over the ordinary edges; for
low-k types (~30 ordinary edges, no prism pruning) that is billions of candidates.  This
module instead PROPAGATES ordinary labels from the rank condition:

  A compact Coxeter d-polytope's Gram matrix has rank d+1, so EVERY (d+2)x(d+2) minor
  vanishes.  A symmetric minor determinant is degree <= 2 in any single off-diagonal entry
  (symmetry -> the entry sits in exactly two positions).  So once a pure-ordinary (d+2)-subset
  has all but one of its C(d+2,2) pairs pinned, the last ordinary entry g is a root of a
  quadratic; it is only a valid label if g == -cos(pi/m) for some m in VALID_LABELS.

The search therefore SEEDS vertices from the spherical-d library (each vertex is a finite
Coxeter group; its C(d,2) angles are one S(d) row) and, after each seed, PROPAGATES: any
pure-ordinary (d+2)-minor with a single unknown ordinary pair is solved and pinned (or the
branch pruned if the forced value is non-cyclotomic), and any fully-pinned pure-ordinary
minor whose determinant is not ~0 prunes the branch.  Every ordinary pair lies in some vertex
(adjacent facets share a ridge -> a common vertex), so once all vertices are consistently
labelled all ordinary entries are fixed; the DOTTED weights are then handed to the existing
exact solver, unchanged.

This is a research prototype: the payoff (branch nodes << block-paste candidates) is
measured empirically (see run_rankprop.py).
"""
from __future__ import annotations

import itertools

import numpy as np

from pipeline.stage4_blockpaste import _setup as _bp_setup
from pipeline.utils import mazheng_lib as ml
from pipeline.stage4_gram import _GRAM_FLOAT, VALID_LABELS

# label -> float Gram entry (-cos(pi/m)); and the reverse match table
_G = dict(_GRAM_FLOAT)
_ENTRY_TO_LABEL = sorted(((v, m) for m, v in _G.items()))   # (entry, m) ascending
HIGH = (7, 8, 9, 10, 12)                                    # concrete values of wildcard 7


def _match_label(g, tol=1e-6):
    """Return the integer label m with -cos(pi/m) == g (within tol), or None."""
    best_m, best_d = None, tol
    for m, v in _G.items():
        dd = abs(g - v)
        if dd < best_d:
            best_m, best_d = m, dd
    return best_m


def _expand_seed_rows(d):
    """Concrete S(d) vertex labelings: each library row with wildcard 7 expanded to HIGH."""
    rows = []
    for tup in sorted(ml.S(d)):
        wild = [i for i, m in enumerate(tup) if m == 7]
        if not wild:
            rows.append(tuple(int(m) for m in tup))
        else:
            for combo in itertools.product(HIGH, repeat=len(wild)):
                r = list(int(m) for m in tup)
                for i, val in zip(wild, combo):
                    r[i] = val
                rows.append(tuple(r))
    # dedupe
    return sorted(set(rows))


def _pure_ordinary_minors(n, d, dotted):
    """(d+2)-subsets whose C(d+2,2) pairs are all ordinary (contain no dotted pair), with
    their pair list in a fixed order."""
    out = []
    for S in itertools.combinations(range(n), d + 2):
        pairs = [(a, b) for a, b in itertools.combinations(S, 2)]
        if any(p in dotted for p in pairs):
            continue
        out.append((S, pairs))
    return out


def _det_entry_roots(sub, ii, jj):
    """Real roots g of det(M)=0 where M=sub but M[ii,jj]=M[jj,ii]=g.  det is quadratic in g
    (sampled at g in {0,-0.5,-1})."""
    def det_at(g):
        M = sub.copy(); M[ii, jj] = M[jj, ii] = g
        return np.linalg.det(M)
    c = det_at(0.0)
    f1 = det_at(-0.5)
    f2 = det_at(-1.0)
    # c + b*g + a*g^2 ; g=-0.5 -> c -0.5b +0.25a ; g=-1 -> c -b +a
    a = 2.0 * (f2 - 2.0 * f1 + c)
    b = -(f2 - c - a)
    if abs(a) < 1e-12:
        if abs(b) < 1e-12:
            return []
        return [-c / b]
    disc = b * b - 4 * a * c
    if disc < 0:
        return []
    r = np.sqrt(disc)
    return [(-b + r) / (2 * a), (-b - r) / (2 * a)]


class _Prop:
    def __init__(self, t):
        (self.V, self.missing, self.dotted, self.ordinary, self.op_idx,
         self.Nr, self.cols_of, self.n, self.d) = _bp_setup(t)
        self.Vsorted = [tuple(sorted(v)) for v in self.V]
        self.dotted = set(self.dotted)
        self.ordinary = set(self.ordinary)
        self.seed_rows = _expand_seed_rows(self.d)
        self.minors = _pure_ordinary_minors(self.n, self.d, self.dotted)
        # greedy vertex order (maximise overlap with covered pairs)
        self.order = self._order_vertices()
        self.nodes = 0
        self.max_nodes = 10 ** 12
        self.solutions = []   # list of {ordinary pair: m}

    def _order_vertices(self):
        rem = list(range(len(self.Vsorted)))
        order = [rem.pop(0)]
        cov = set(ml.facet_pair_columns(self.Vsorted[order[0]]))
        while rem:
            best, bo = None, -1
            for i in rem:
                ov = len(cov & set(ml.facet_pair_columns(self.Vsorted[i])))
                if ov > bo:
                    best, bo = i, ov
            order.append(best); cov |= set(ml.facet_pair_columns(self.Vsorted[best]))
            rem.remove(best)
        return order

    def _build_sub(self, S, lab, unknown_pos=None):
        """float submatrix for subset S from pinned ordinary labels; unknown pair left 0."""
        k = len(S); M = np.eye(k)
        idx = {f: a for a, f in enumerate(S)}
        for a, b in itertools.combinations(S, 2):
            if (a, b) == unknown_pos:
                continue
            m = lab.get((a, b))
            g = _G[m] if m is not None else 0.0
            M[idx[a], idx[b]] = M[idx[b], idx[a]] = g
        return M, idx

    def _propagate(self, lab, tol=1e-6):
        """Pin every forced ordinary pair; return (ok, lab) or (False, None) if pruned."""
        lab = dict(lab)
        changed = True
        while changed:
            changed = False
            for S, pairs in self.minors:
                unk = [p for p in pairs if p not in lab]
                if len(unk) == 0:
                    M, _ = self._build_sub(S, lab)
                    if abs(np.linalg.det(M)) > 1e-6:
                        return False, None            # rank condition violated
                elif len(unk) == 1:
                    e = unk[0]
                    M, idx = self._build_sub(S, lab, unknown_pos=e)
                    ms = []
                    for g in _det_entry_roots(M, idx[e[0]], idx[e[1]]):
                        m = _match_label(g, tol)
                        if m is not None and m not in ms:
                            ms.append(m)
                    if not ms:
                        return False, None            # forced value non-cyclotomic
                    if len(ms) == 1:
                        lab[e] = ms[0]; changed = True
                    # multiple candidate labels -> leave for vertex branching
        return True, lab

    def _dfs(self, vi, lab):
        if self.nodes > self.max_nodes:
            return
        ok, lab = self._propagate(lab)
        if not ok:
            return
        if vi == len(self.order):
            if all(p in lab for p in self.ordinary):
                self.solutions.append({p: lab[p] for p in self.ordinary})
            return
        v = self.Vsorted[self.order[vi]]
        vpairs = ml.facet_pair_columns(v)                 # C(d,2) pairs, combinations order
        if all(p in lab for p in vpairs):
            self._dfs(vi + 1, lab); return
        for row in self.seed_rows:
            assign = {vpairs[k]: row[k] for k in range(len(vpairs))}
            if all(lab.get(p, assign[p]) == assign[p] for p in vpairs):
                self.nodes += 1
                self._dfs(vi + 1, {**lab, **assign})

    def run(self):
        self._dfs(0, {})
        return self.solutions, self.nodes


def rank_propagation_candidates(t):
    """Return (list of ordinary label-assignment dicts {(i,j): m}, branch_node_count) for a
    combinatorial type, via rank-propagation.  Dotted weights are solved downstream."""
    p = _Prop(t)
    return p.run()
