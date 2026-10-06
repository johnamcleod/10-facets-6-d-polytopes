#!/usr/bin/env python3
"""An independently checkable certificate for the search tree: SAT + DRAT.

The classification enumerator decides 260 of the 304 searched types, and every
pruned branch of the other 44, by forward checking.  Its verdict "exhausted" is
evidence produced by the pipeline's own code.  This script replaces that evidence
by a proof that an off-the-shelf checker can verify.

It shares NO code with the pipeline.  From a type's missing-face hypergraph alone it
builds a CNF whose satisfying assignments include every Coxeter labelling a compact
polytope of that type could carry:

  variables   x[e, l]  for each ordinary pair e (a non-dashed pair of facets) and
              each label l in {2,3,4,5,6,*}, * standing for any m >= 7;
  (1)  each ordinary pair carries exactly one label;
  (2)  every FACE carries an elliptic labelling (Proposition 2.2).  A diagram is
       non-elliptic iff it contains a minimal non-elliptic subdiagram, and a subset
       of a face is a face, so this is encoded by forbidding, on every face S of
       size 3..6, every minimal non-elliptic labelling of S;
  (3)  every MISSING FACE of size >= 3 carries a Lannér labelling (Proposition 2.2).
       Its proper subsets are faces, hence elliptic by (2); given that, the
       labelling is Lannér iff it is neither elliptic nor parabolic, so the allowed
       labellings are enumerated and selected by auxiliary variables;
  (4)  the low-weight caps of Burcroff Lemma 5.5(b).

Ellipticity is decided by the classification of finite Coxeter groups, implemented
here from scratch as a test on labelled graphs.  Lannér vs parabolic is decided for
3 nodes by the triangle trichotomy (integer arithmetic, * read as 7, which Lemma
4.6(b) justifies) and for 4 and 5 nodes, where no label exceeds 5, by the sign of the
Gram determinant computed exactly in Q(sqrt2, sqrt5) with SymPy.  The resulting
Lannér lists are cross-checked against Lannér's classification: 9 diagrams of order
4 and 5 of order 5 up to isomorphism.

  (5)  the diagram is connected (Proposition 2.1), one clause per split of the
       facets into two parts that no dashed edge crosses.

Omitted, soundly: all symmetry breaking.  Omitting a necessary condition can only
add solutions, never remove one.  All 304 types are certified with (1)-(5); types
57 and 81 add the lex-leader predicate of run_type_sb.

Verdicts.
  * For a type in which no complete labelling reaches a screen, the CNF must be
    UNSATISFIABLE; the solver's DRAT proof is written out and checked by drat-trim.
  * For the other types, every labelling the search screened is refuted by an exact
    certificate, except P_{6,10}.  Those labellings,
    and their images under Aut(T), are added as blocking clauses (an image of an
    unrealizable labelling is unrealizable), as is every solution the solver finds
    that is disconnected (Proposition 2.1).  The final CNF must then be
    unsatisfiable, with a DRAT proof checked in the same way.  What this certifies is
    that every labelling satisfying the local conditions is a relabelling of one the
    exact certificates refute, or of P_{6,10}.

Run:  python3 checks/tree_sat.py [tids=1,2,...|all] [drat=path/to/drat-trim]
      [out=DIR]
      python3 checks/tree_sat.py tids=57 split=57   (36-way case split, for a
      type whose blocking clauses are too many for one formula)
      python3 checks/tree_sat.py tids=57,81 sb=57,81   (lex-leader symmetry
      breaking in place of blocking every image; see run_type_sb)
      python3 checks/tree_sat.py --link   (match every blocked labelling to its
      exact certificate, or to P_{6,10}; adds the counts to the artifact)
      python3 checks/tree_sat.py merge=log1,log2,...   (rebuild the artifact
      from the JSON lines of runs made in parallel)
Artifact: runs/d6_n10/tree_sat.json, with the CNF and DRAT files under out=DIR.
"""
from __future__ import annotations

import itertools
import json
import subprocess
import sys
import time
from fractions import Fraction
from pathlib import Path

import networkx as nx
import sympy as sp
from networkx.algorithms import isomorphism as iso
from pysat.solvers import Solver

ROOT = Path(__file__).resolve().parents[1]
N, D = 10, 6
LABELS = (2, 3, 4, 5, 6, 7)          # 7 is the wildcard *: "any m >= 7"
ARG = {a.split("=", 1)[0]: a.split("=", 1)[1] for a in sys.argv[1:] if "=" in a}


# ------------------------------------------------------------ finite Coxeter types
def elliptic(nodes, lab):
    """True iff the labelled diagram on `nodes` is elliptic, i.e. every connected
    component is a finite Coxeter group: A_n, B_n, D_n, E_6,7,8, F_4, H_3, H_4 or a
    rank-2 I_2(m) (any finite m, including the wildcard).  lab(u, v) in LABELS, with
    2 meaning 'no edge'."""
    nodes = list(nodes)
    adj = {u: [v for v in nodes if v != u and lab(u, v) != 2] for u in nodes}
    seen = set()
    for s in nodes:
        if s in seen:
            continue
        comp, stack = [], [s]
        seen.add(s)
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in adj[u]:
                if v not in seen:
                    seen.add(v)
                    stack.append(v)
        if not _finite_connected(comp, adj, lab):
            return False
    return True


def _finite_connected(comp, adj, lab):
    k = len(comp)
    if k <= 2:
        return True                                   # A_1, I_2(m)
    edges = [(u, v) for u, v in itertools.combinations(comp, 2) if lab(u, v) != 2]
    if len(edges) != k - 1:
        return False                                  # not a tree
    labs = [lab(u, v) for u, v in edges]
    if any(m > 5 for m in labs):
        return False
    deg = {u: len(adj[u]) for u in comp}
    if max(deg.values()) > 3 or sum(1 for d in deg.values() if d == 3) > 1:
        return False
    n4, n5 = labs.count(4), labs.count(5)
    if 3 in deg.values():                             # D_n or E_n: simply laced
        if n4 or n5:
            return False
        c = next(u for u in comp if deg[u] == 3)
        arms = []
        for v in adj[c]:
            length, prev, cur = 1, c, v
            while deg[cur] == 2:
                nxt = next(w for w in adj[cur] if w != prev)
                prev, cur, length = cur, nxt, length + 1
            arms.append(length)
        arms.sort()
        return arms[:2] == [1, 1] or arms in ([1, 2, 2], [1, 2, 3], [1, 2, 4])
    # a path
    ends = [u for u in comp if deg[u] == 1]
    order, prev, cur = [ends[0]], None, ends[0]
    while len(order) < k:
        nxt = next(w for w in adj[cur] if w != prev)
        prev, cur = cur, nxt
        order.append(cur)
    plabs = [lab(order[i], order[i + 1]) for i in range(k - 1)]
    if n4 == 0 and n5 == 0:
        return True                                   # A_n
    if n4 == 1 and n5 == 0:
        i = plabs.index(4)
        return i in (0, k - 2) or (k == 4 and i == 1)  # B_n, F_4
    if n5 == 1 and n4 == 0:
        i = plabs.index(5)
        return k in (3, 4) and i in (0, k - 2)       # H_3, H_4
    return False


# ------------------------------------------------ minimal non-elliptic labellings
def minimal_non_elliptic(k):
    """All labellings of the complete graph on nodes 0..k-1 (labels in LABELS) whose
    diagram is non-elliptic while every proper subdiagram is elliptic.  Enumerated
    by backtracking: an edge label is fixed only if every subset whose edges are then
    all labelled, other than the whole node set, is elliptic."""
    nodes = list(range(k))
    edges = list(itertools.combinations(nodes, 2))
    # the proper subsets completed by labelling edge number i
    done_at = {i: [] for i in range(len(edges))}
    for r in range(3, k):
        for S in itertools.combinations(nodes, r):
            last = max(edges.index(e) for e in itertools.combinations(S, 2))
            done_at[last].append(S)
    out, lab = [], {}

    def L(u, v):
        return lab[(min(u, v), max(u, v))]

    def rec(i):
        if i == len(edges):
            if not elliptic(nodes, L):
                out.append(dict(lab))
            return
        for m in LABELS:
            lab[edges[i]] = m
            if all(elliptic(S, L) for S in done_at[i]):
                rec(i + 1)
        del lab[edges[i]]
    rec(0)
    return out


COS = {2: sp.Integer(0), 3: sp.Rational(1, 2), 4: sp.sqrt(2) / 2,
       5: (1 + sp.sqrt(5)) / 4}


def det_sign(k, lab):
    """Exact sign of the Gram determinant, labels <= 5 (entries in Q(sqrt2, sqrt5))."""
    G = sp.eye(k)
    for (u, v), m in lab.items():
        G[u, v] = G[v, u] = -COS[m]
    d = sp.nsimplify(sp.expand(G.det(method="berkowitz")))
    d = sp.radsimp(d)
    if d == 0 or sp.expand(d) == 0:
        return 0
    # d is a nonzero element of a number field; 60 digits decide its sign
    v = sp.N(d, 60)
    assert abs(v) > sp.Float("1e-40"), f"near-zero determinant {d}"
    return 1 if v > 0 else -1


def lanner_labellings(k):
    """Labellings of nodes 0..k-1 that are Lannér: minimal non-elliptic and not
    parabolic.  Returns (list, number of isomorphism classes)."""
    out = []
    cache = {}
    for lab in minimal_non_elliptic(k):
        if k == 3:
            p, q, r = (lab[(0, 1)], lab[(0, 2)], lab[(1, 2)])
            s = Fraction(1, p) + Fraction(1, q) + Fraction(1, r)
            if s < 1:
                out.append(lab)                       # s == 1 is parabolic
            continue
        assert all(m <= 5 for m in lab.values()), "label > 5 in order >= 4"
        key = _canon(k, lab)
        if key not in cache:
            cache[key] = det_sign(k, lab)
        if cache[key] < 0:
            out.append(lab)
    classes = {_canon(k, lab) for lab in out}
    return out, len(classes)


def _canon(k, lab):
    return min(tuple(lab[(min(p[u], p[v]), max(p[u], p[v]))]
                     for u, v in itertools.combinations(range(k), 2))
               for p in itertools.permutations(range(k)))


# ------------------------------------------------------------------ the encoding
class Enc:
    def __init__(self, t, mne, lanner):
        self.mf = [tuple(sorted(m)) for m in t["missing_faces"]]
        self.dashed = {m for m in self.mf if len(m) == 2}
        big = [set(m) for m in self.mf if len(m) >= 3]
        assert all(len(m) <= 5 for m in self.mf), "a missing face of size > 5"
        self.pairs = [e for e in itertools.combinations(range(N), 2)
                      if e not in self.dashed]
        self.var = {}
        for e in self.pairs:
            for m in LABELS:
                self.var[(e, m)] = len(self.var) + 1
        self.nv = len(self.var)
        self.clauses = []
        # (1) exactly one label per ordinary pair
        for e in self.pairs:
            lits = [self.var[(e, m)] for m in LABELS]
            self.clauses.append(lits)
            self.clauses += [[-a, -b] for a, b in itertools.combinations(lits, 2)]
        # faces: every subset of [10] containing no missing face
        mfs = [set(m) for m in self.mf]
        self.faces = [S for r in range(3, D + 1)
                      for S in itertools.combinations(range(N), r)
                      if not any(M <= set(S) for M in mfs)]
        self.vertices = [S for S in self.faces if len(S) == D]
        # (2) forbid minimal non-elliptic labellings on every face
        for S in self.faces:
            for lab in mne[len(S)]:
                self.clauses.append([-self.var[((S[u], S[v]), m)]
                                     for (u, v), m in lab.items()])
        # (3) every missing face of size >= 3 carries a Lannér labelling
        for M in [m for m in self.mf if len(m) >= 3]:
            sel = []
            for lab in lanner[len(M)]:
                y = self.new()
                sel.append(y)
                for (u, v), m in lab.items():
                    self.clauses.append([-y, self.var[((M[u], M[v]), m)]])
            self.clauses.append(sel)
        # (5) Sigma is connected (Proposition 2.1): for every split of the facets
        # into two nonempty parts that no dashed edge crosses, some ordinary pair
        # across the split carries a label other than 2.
        self.cuts = 0
        for mask in range(1, 2 ** (N - 1)):
            U = {i for i in range(N) if mask >> i & 1}
            if any((a in U) != (b in U) for a, b in self.dashed):
                continue
            cross = [e for e in self.pairs if (e[0] in U) != (e[1] in U)]
            self.clauses.append([-self.var[(e, 2)] for e in cross])
            self.cuts += 1
        # (4) Burcroff Lemma 5.5(b): caps m <= 5
        self.capped = []
        for (u, v) in self.pairs:
            if any(u in M and v in M for M in big):
                continue
            for v1, v2 in ((u, v), (v, u)):
                if any(v2 in M and not any(v1 in d and (set(d) - {v1}) & M
                                           for d in self.dashed) for M in big):
                    self.capped.append((u, v))
                    self.clauses += [[-self.var[((u, v), 6)]], [-self.var[((u, v), 7)]]]
                    break

    def new(self):
        self.nv += 1
        return self.nv

    def decode(self, model):
        pos = {l for l in model if l > 0}
        return {e: next(m for m in LABELS if self.var[(e, m)] in pos)
                for e in self.pairs}

    def block(self, labelling):
        return [-self.var[(e, labelling[e])] for e in self.pairs]


def connected(labelling, dashed):
    g = nx.Graph()
    g.add_nodes_from(range(N))
    g.add_edges_from(e for e, m in labelling.items() if m != 2)
    g.add_edges_from(dashed)
    return nx.is_connected(g)


def automorphisms(mf):
    """Aut(T): facet permutations preserving the missing-face hypergraph."""
    g = nx.Graph()
    for i in range(N):
        g.add_node(("v", i), side=0)
    for k, m in enumerate(mf):
        g.add_node(("e", k), side=1)
        g.add_edges_from((("v", i), ("e", k)) for i in m)
    gm = iso.GraphMatcher(g, g, node_match=iso.categorical_node_match("side", None))
    return [[phi[("v", i)][1] for i in range(N)] for phi in gm.isomorphisms_iter()]


def image(labelling, perm):
    return {(min(perm[u], perm[v]), max(perm[u], perm[v])): m
            for (u, v), m in labelling.items()}


# ---------------------------------------------------------------- screened sets
_SCREENED = None


def screened_labellings(tid):
    """Every labelling of type `tid` that reached a screen in the run of record,
    as a map ordinary pair -> label, read from the instance dumps (the ones the exact
    certifiers refuted, plus P_{6,10})."""
    global _SCREENED
    if _SCREENED is None:
        _SCREENED = {}
        for lab_tid, lab in _read_dumps():
            _SCREENED.setdefault(lab_tid, []).append(lab)
    return _SCREENED.get(tid, [])


def _read_dumps():
    for f in ("plain_instances.jsonl", "wild_instances.jsonl"):
        for line in (ROOT / "runs/d6_n10/d6_valid" / f).read_text().splitlines():
            if not line.strip():
                continue
            r = json.loads(line)
            lab = {(min(i, j), max(i, j)): m for i, j, m in r["ordinary"]}
            for i, j in r.get("wild", []):
                lab[(min(i, j), max(i, j))] = 7
            yield int(r["tag"].split("|")[0]), lab


def refuted_ok():
    """The exact certificates must cover every screened labelling but one."""
    for f, key in (("plain_certificates.json", None), ("wild_certificates.json", None)):
        s = json.load(open(ROOT / "runs/d6_n10/d6_valid" / f))["summary"]
        assert s.get("undecided_among_rejected", 1) == 0, f"{f}: undecided rejections"
    return True


# --------------------------------------------------------------------- driving
def run_type(tid, t, mne, lanner, outdir, drat):
    enc = Enc(t, mne, lanner)
    blocked = 0
    screened = []
    state = {}
    for line in (ROOT / "runs/d6_n10/d6_valid/state.jsonl").read_text().splitlines():
        if line.strip():
            r = json.loads(line)
            state[r["key"]] = r
    n_screened = state[f"{tid}|"]["enum_count"]
    # Blocking clauses are STREAMED -- written to the CNF file and passed to the
    # solver one at a time, never held as Python lists -- because for the largest
    # types the Aut(T)-images of the screened labellings number in the millions.
    cnf = outdir / f"type{tid}.cnf"
    prf = outdir / f"type{tid}.drat"
    body = outdir / f"type{tid}.cnf.body"
    nclauses = 0
    # Glucose 4 emits text DRAT that drat-trim verifies.  (CaDiCaL 1.9.5 as bundled
    # with PySAT produced proofs drat-trim could not replay -- "no conflict" -- on
    # some of these formulas, so it is not used for the certificate.)
    s = Solver(name="glucose4", with_proof=True)
    fh = open(body, "w")

    def add(c):
        nonlocal nclauses
        s.add_clause(c)
        fh.write(" ".join(map(str, c)) + " 0\n")
        nclauses += 1

    for c in enc.clauses:
        add(c)
    blocked = 0
    if n_screened:
        screened = screened_labellings(tid)
        assert len(screened) == n_screened, (tid, len(screened), n_screened)
        # Non-vacuity: every screened labelling (which passed the pipeline's forward
        # checking) must satisfy this independent CNF.  A rejection here would mean
        # the encoding is stricter than the search -- a possible unsound constraint.
        chk = Solver(name="cadical195", bootstrap_with=enc.clauses)
        inconsistent = sum(
            1 for lab in screened
            if not chk.solve(assumptions=[enc.var[(e, lab[e])] for e in enc.pairs]))
        chk.delete()
        print(f"  [{tid}] non-vacuity: {len(screened)} screened labellings checked",
              file=sys.stderr, flush=True)
        assert inconsistent == 0, (
            f"type {tid}: {inconsistent} screened labellings violate the CNF")
        auts = automorphisms(enc.mf)
        seen = set()
        for lab in screened:
            for p in auts:
                im = image(lab, p)
                key = bytes(im[e] for e in enc.pairs)
                if key not in seen:
                    seen.add(key)
                    add(enc.block(im))
        blocked = len(seen)
        print(f"  [{tid}] {blocked} blocking clauses ({len(auts)} automorphisms)",
              file=sys.stderr, flush=True)
        del seen, screened
    t0 = time.time()
    disconnected = 0
    while s.solve():
        lab = enc.decode(s.get_model())
        if connected(lab, enc.dashed):
            s.delete()
            fh.close()
            return {"tid": tid, "verdict": "SAT", "solution": {f"{u}{v}": m for (u, v), m
                    in lab.items() if m != 2}}
        add(enc.block(lab))
        disconnected += 1
        if disconnected % 100 == 0:
            print(f"  [{tid}] {disconnected} disconnected solutions blocked",
                  file=sys.stderr, flush=True)
    fh.close()
    # copy the proof trace to disk without loading it
    import shutil
    s.solver.prfile.seek(0)
    with open(prf, "wb") as out:
        shutil.copyfileobj(s.solver.prfile, out)
    s.delete()
    sec = round(time.time() - t0, 1)
    with open(cnf, "w") as out:
        out.write(f"p cnf {enc.nv} {nclauses}\n")
        with open(body) as src:
            shutil.copyfileobj(src, out)
    body.unlink()
    proof_lines = sum(1 for _ in open(prf, "rb"))
    print(f"  [{tid}] UNSAT after {sec}s; checking {proof_lines} proof lines",
          file=sys.stderr, flush=True)
    chk = None
    if drat:
        r = subprocess.run([drat, str(cnf), str(prf)], capture_output=True, text=True)
        chk = "VERIFIED" if "s VERIFIED" in r.stdout else "FAILED"
    return {"tid": tid, "verdict": "UNSAT", "vars": enc.nv, "clauses": nclauses,
            "faces": len(enc.faces), "vertices": len(enc.vertices),
            "screened": n_screened, "screened_consistent_with_cnf": n_screened,
            "aut_order": len(automorphisms(enc.mf)), "blocking_clauses": blocked,
            "disconnected_blocked": disconnected, "capped": len(enc.capped),
            "proof_lines": proof_lines, "solve_seconds": sec, "drat_trim": chk}


def run_type_split(tid, t, mne, lanner, outdir, drat, tlim=200000):
    """The same certificate as run_type, for a type whose blocking clauses are too
    many to hold in one formula (types 57 and 81: 3.5e7 and 8.5e6 clauses).

    Two ordinary pairs e1, e2 are chosen and the formula is split into the 36 cases
    x[e1, l1] & x[e2, l2].  Each case gets the base clauses, those two unit clauses,
    and only the blocking clauses whose labelling has l1 at e1 and l2 at e2 (the
    others are satisfied by the units).  By the exactly-one clauses the 36 cases
    cover every assignment, so the formula is unsatisfiable iff every case is; each
    case's DRAT proof is checked separately.  Images are bucketed on disk and
    de-duplicated with sort -u, so memory stays flat."""
    import shutil
    enc = Enc(t, mne, lanner)
    state = {}
    for line in (ROOT / "runs/d6_n10/d6_valid/state.jsonl").read_text().splitlines():
        if line.strip():
            r = json.loads(line)
            state[r["key"]] = r
    n_screened = state[f"{tid}|"]["enum_count"]
    screened = screened_labellings(tid)
    assert len(screened) == n_screened
    chk = Solver(name="cadical195", bootstrap_with=enc.clauses)
    bad = sum(1 for lab in screened
              if not chk.solve(assumptions=[enc.var[(e, lab[e])] for e in enc.pairs]))
    chk.delete()
    assert bad == 0, f"type {tid}: {bad} screened labellings violate the CNF"
    print(f"  [{tid}] non-vacuity: {n_screened} screened labellings checked",
          file=sys.stderr, flush=True)
    P = enc.pairs
    idx = {e: i for i, e in enumerate(P)}
    auts = automorphisms(enc.mf)
    # pm[k][i] = index of the image of pair i under automorphism k
    pm = [[idx[(min(p[u], p[v]), max(p[u], p[v]))] for (u, v) in P] for p in auts]
    L = [bytes(lab[e] for e in P) for lab in screened]
    del screened
    # split pairs: the two whose labels vary most over a sample of images
    import random
    rnd = random.Random(tid)
    cnt = [[0] * 8 for _ in P]
    for _ in range(20000):
        lab, m = rnd.choice(L), rnd.choice(pm)
        img = bytearray(len(P))
        for i, j in enumerate(m):
            img[j] = lab[i]
        for j, x in enumerate(img):
            cnt[j][x] += 1
    spread = sorted(range(len(P)), key=lambda j: -sum(1 for x in cnt[j] if x))
    j1, j2 = spread[0], spread[1]
    bdir = outdir / f"type{tid}_buckets"
    bdir.mkdir(exist_ok=True)
    # Reuse buckets left by an interrupted run: if every case has a .raw or .uniq
    # file (or was already verified), skip the bucketing pass.
    done_before = {c for c in ARG.get("skipcases", "").split(",") if c}
    have = {f.stem for f in bdir.iterdir()}
    if not all(f"{a}{b}" in have or f"{a}{b}" in done_before
               for a in LABELS for b in LABELS):
        fhs = {(a, b): open(bdir / f"{a}{b}.raw", "wb") for a in LABELS for b in LABELS}
        for lab in L:
            for m in pm:
                img = bytearray(len(P))
                for i, j in enumerate(m):
                    img[j] = lab[i]
                fhs[(img[j1], img[j2])].write(bytes(img) + b"\n")
        for f in fhs.values():
            f.close()
        print(f"  [{tid}] images bucketed on pairs {P[j1]}, {P[j2]}", file=sys.stderr,
              flush=True)
    else:
        print(f"  [{tid}] reusing buckets in {bdir}", file=sys.stderr, flush=True)
    del L
    global _CASE
    _CASE = (tid, enc, P, j1, j2, bdir, outdir, drat, tlim)
    todo = [f"{a}{b}" for a in LABELS for b in LABELS if f"{a}{b}" not in done_before]
    jobs = int(ARG.get("jobs", 1))
    import multiprocessing as mp
    with mp.get_context("fork").Pool(processes=jobs) as pool:
        cases = list(pool.imap_unordered(_do_case, todo))
    cases += [{"case": c, "drat_trim": "VERIFIED", "blocking": None,
               "clauses": 0, "note": "verified in an earlier run"} for c in done_before]
    if any(c.get("verdict") == "SAT" for c in cases):
        return {"tid": tid, "verdict": "SAT",
                "case": [c["case"] for c in cases if c.get("verdict") == "SAT"]}
    total = sum(c["blocking"] or 0 for c in cases)
    allok = all(c["drat_trim"] == "VERIFIED" for c in cases)
    return {"tid": tid, "verdict": "UNSAT", "vars": enc.nv,
            "clauses": max(c["clauses"] for c in cases),
            "faces": len(enc.faces), "vertices": len(enc.vertices),
            "screened": n_screened, "screened_consistent_with_cnf": n_screened,
            "aut_order": len(auts), "blocking_clauses": total,
            "disconnected_blocked": 0, "capped": len(enc.capped),
            "split": [list(P[j1]), list(P[j2])], "cases": cases,
            "proof_lines": 0, "solve_seconds": 0,
            "drat_trim": "VERIFIED" if allok else "FAILED"}


def lex_leader_clauses(enc, perms, newvar):
    """CNF for 'L <=lex pi(L)' for every pi in perms, L the labelling as a vector
    over enc.pairs (in that order) with labels ordered 2 < 3 < ... < 6 < *.

    The standard chain encoding, one chain of auxiliary variables per permutation:
    e_0 holds; e_k implies L[k] <= I[k]; e_k and L[k] = I[k] imply e_{k+1}; where
    I = pi(L), I[k] = L[q[k]].  If L <=lex I the chain can stop at the first strict
    inequality; if L >lex I the equalities force e up to the first difference, where
    L[k] > I[k] contradicts it.  So, projected on L, the clauses say exactly
    'L is lexicographically least in its orbit'."""
    P = enc.pairs
    idx = {e: i for i, e in enumerate(P)}
    out = []
    for p in perms:
        # I[j] = L[q[j]] where pair q[j] is mapped by p onto pair j
        q = [None] * len(P)
        for i, (u, v) in enumerate(P):
            q[idx[(min(p[u], p[v]), max(p[u], p[v]))]] = i
        if q == list(range(len(P))):
            continue                                   # the identity
        e = newvar()
        out.append([e])
        for k in range(len(P)):
            if q[k] == k:
                continue
            a, b = P[k], P[q[k]]
            for la in LABELS:                          # L[k] <= I[k]
                for lb in LABELS:
                    if la > lb:
                        out.append([-e, -enc.var[(a, la)], -enc.var[(b, lb)]])
            nxt = newvar()
            for l in LABELS:                           # equal -> continue the chain
                out.append([-e, -enc.var[(a, l)], -enc.var[(b, l)], nxt])
            e = nxt
    return out


def run_type_sb(tid, t, mne, lanner, outdir, drat, tlim=200000):
    """The certificate for a type with a very large automorphism group, using a
    lex-leader symmetry-breaking predicate instead of blocking every image.

    Every constraint of the base CNF is determined by the type's missing-face
    hypergraph, so the formula F is invariant under Aut(T), acting on labellings by
    permuting pairs.  Let SB say that L is lexicographically least in its Aut(T)-
    orbit (lex_leader_clauses).  If F has a solution outside the orbits of the
    screened labellings, its orbit's least element solves F and SB and is not the
    least element of any screened orbit.  Hence: F and SB, with the least element of
    each screened orbit blocked, unsatisfiable  =>  every solution of F is an
    image of a screened labelling.  That formula's DRAT proof is checked by
    drat-trim.  What this adds to the trusted base is lex_leader_clauses, a
    standard encoding, and the canonical form below; as a check on both, every
    blocked representative is verified to satisfy F and SB."""
    import shutil
    enc = Enc(t, mne, lanner)
    state = {}
    for line in (ROOT / "runs/d6_n10/d6_valid/state.jsonl").read_text().splitlines():
        if line.strip():
            r = json.loads(line)
            state[r["key"]] = r
    n_screened = state[f"{tid}|"]["enum_count"]
    screened = screened_labellings(tid)
    assert len(screened) == n_screened
    chk = Solver(name="cadical195", bootstrap_with=enc.clauses)
    bad = sum(1 for lab in screened
              if not chk.solve(assumptions=[enc.var[(e, lab[e])] for e in enc.pairs]))
    chk.delete()
    assert bad == 0, f"type {tid}: {bad} screened labellings violate the CNF"
    print(f"  [{tid}] non-vacuity: {n_screened} screened labellings checked",
          file=sys.stderr, flush=True)
    P = enc.pairs
    idx = {e: i for i, e in enumerate(P)}
    auts = automorphisms(enc.mf)
    pm = [[idx[(min(p[u], p[v]), max(p[u], p[v]))] for (u, v) in P] for p in auts]
    # canonical representative of each screened orbit: the lexicographically least
    # image, labels compared as integers 2 < ... < 7 in the pair order of P
    reps = set()
    for lab in screened:
        L = [lab[e] for e in P]
        best = None
        for m in pm:
            img = [0] * len(P)
            for i, j in enumerate(m):
                img[j] = L[i]
            img = tuple(img)
            if best is None or img < best:
                best = img
        reps.add(best)
    del screened
    print(f"  [{tid}] {len(reps)} canonical representatives ({len(auts)} automorphisms)",
          file=sys.stderr, flush=True)
    sb = lex_leader_clauses(enc, auts, enc.new)
    # check: each representative satisfies F and SB (validates both the encoding of
    # SB and the canonical form)
    chk = Solver(name="cadical195", bootstrap_with=enc.clauses + sb)
    bad = sum(1 for r in reps
              if not chk.solve(assumptions=[enc.var[(P[j], r[j])] for j in range(len(P))]))
    chk.delete()
    assert bad == 0, f"type {tid}: {bad} canonical representatives violate F and SB"
    print(f"  [{tid}] {len(sb)} symmetry-breaking clauses; all representatives "
          f"satisfy F and SB", file=sys.stderr, flush=True)
    cnf, prf = outdir / f"type{tid}_sb.cnf", outdir / f"type{tid}_sb.drat"
    blocks = [[-enc.var[(P[j], r[j])] for j in range(len(P))] for r in sorted(reps)]
    allc = enc.clauses + sb + blocks
    s = Solver(name="glucose4", with_proof=True, bootstrap_with=allc)
    t0 = time.time()
    disconnected = 0
    while s.solve():
        lab = enc.decode(s.get_model())
        if connected(lab, enc.dashed):
            s.delete()
            return {"tid": tid, "verdict": "SAT", "mode": "lex-leader",
                    "solution": {f"{u}{v}": m for (u, v), m in lab.items() if m != 2}}
        c = enc.block(lab)
        s.add_clause(c)
        allc.append(c)
        disconnected += 1
    sec = round(time.time() - t0, 1)
    s.solver.prfile.seek(0)
    with open(prf, "wb") as out:
        shutil.copyfileobj(s.solver.prfile, out)
    s.delete()
    with open(cnf, "w") as fh:
        fh.write(f"p cnf {enc.nv} {len(allc)}\n")
        for c in allc:
            fh.write(" ".join(map(str, c)) + " 0\n")
    r = subprocess.run([drat, str(cnf), str(prf), "-t", str(tlim)],
                       capture_output=True, text=True) if drat else None
    ok = r is not None and "s VERIFIED" in r.stdout
    print(f"  [{tid}] UNSAT after {sec}s; drat-trim {'VERIFIED' if ok else 'FAILED'}",
          file=sys.stderr, flush=True)
    return {"tid": tid, "verdict": "UNSAT", "mode": "lex-leader", "vars": enc.nv,
            "clauses": len(allc), "faces": len(enc.faces),
            "vertices": len(enc.vertices), "screened": n_screened,
            "screened_consistent_with_cnf": n_screened, "aut_order": len(auts),
            "canonical_representatives": len(reps), "sb_clauses": len(sb),
            "blocking_clauses": len(blocks), "disconnected_blocked": disconnected,
            "capped": len(enc.capped), "proof_lines": sum(1 for _ in open(prf, "rb")),
            "solve_seconds": sec, "drat_trim": "VERIFIED" if ok else "FAILED"}


_CASE = None


def _do_case(case):
    """One of the 36 cases of run_type_split (run in a worker process)."""
    import shutil
    tid, enc, P, j1, j2, bdir, outdir, drat, tlim = _CASE
    a, b = int(case[0]), int(case[1])
    raw, uniq = bdir / f"{case}.raw", bdir / f"{case}.uniq"
    if raw.exists():
        subprocess.run(["sort", "-u", "-o", str(uniq), str(raw)], check=True,
                       env={"LC_ALL": "C", "PATH": "/usr/bin:/bin"})
        raw.unlink()
    cnf = outdir / f"type{tid}_{case}.cnf"
    prf = outdir / f"type{tid}_{case}.drat"
    s = Solver(name="glucose4", with_proof=True)
    n = nb = 0
    with open(cnf.with_suffix(".body"), "w") as fh:
        def add(c):
            nonlocal n
            s.add_clause(c)
            fh.write(" ".join(map(str, c)) + " 0\n")
            n += 1
        for c in enc.clauses:
            add(c)
        add([enc.var[(P[j1], a)]])
        add([enc.var[(P[j2], b)]])
        if uniq.exists():
            with open(uniq, "rb") as src:
                for line in src:
                    img = line.rstrip(b"\n")
                    add([-enc.var[(P[j], img[j])] for j in range(len(P))])
                    nb += 1
        while s.solve():
            lab = enc.decode(s.get_model())
            if connected(lab, enc.dashed):
                s.delete()
                return {"case": case, "verdict": "SAT"}
            add(enc.block(lab))
    s.solver.prfile.seek(0)
    with open(prf, "wb") as out:
        shutil.copyfileobj(s.solver.prfile, out)
    s.delete()
    with open(cnf, "w") as out:
        out.write(f"p cnf {enc.nv} {n}\n")
        with open(cnf.with_suffix(".body")) as src:
            shutil.copyfileobj(src, out)
    cnf.with_suffix(".body").unlink()
    if uniq.exists():
        uniq.unlink()
    r = subprocess.run([drat, str(cnf), str(prf), "-t", str(tlim)],
                       capture_output=True, text=True) if drat else None
    ok = r is not None and "s VERIFIED" in r.stdout
    print(f"  [{tid}] case {case}: {nb} blocking, {'VERIFIED' if ok else 'FAILED'}",
          file=sys.stderr, flush=True)
    if ok:
        cnf.unlink()
        prf.unlink()
    return {"case": case, "clauses": n, "blocking": nb,
            "drat_trim": "VERIFIED" if ok else "FAILED"}


def link_certificates():
    """Close the loop with the exact refutation certificates:
    every labelling the certificate blocks is read from the refutation dump
    (screened_labellings), and here each one is matched, labelling by labelling, to
    an exact certificate in the certifiers' artifacts -- or is the single labelling
    the cascade passed, P_{6,10}.  The counts are added to the artifact."""
    run = ROOT / "runs/d6_n10/d6_valid"

    def key(r, tid=None):
        lab = {(min(i, j), max(i, j)): m for i, j, m in r["ordinary"]}
        for i, j in r.get("wild", []):
            lab[(min(i, j), max(i, j))] = 7
        return (int(r["tag"].split("|")[0]), tuple(sorted(lab.items())))

    certified, accepted = set(), set()
    for r in json.load(open(run / "plain_certificates.json"))["instances"]:
        if r["certified_unrealizable"]:
            certified.add(key(r))
        elif r["cascade_decision"] == "candidate":
            accepted.add(key(r))
    shards = sorted(run.glob("wild_certificates.shard*of*.json"))
    assert shards, "wildcard certificate shards (full records) not found"
    for f in shards:
        for r in json.load(open(f))["instances"]:
            if r["certified_infeasible"]:
                certified.add(key(r))
    assert len(accepted) == 1, f"{len(accepted)} accepted labellings, expected 1"
    art = ROOT / "runs/d6_n10/tree_sat.json"
    ts = json.loads(art.read_text())
    per, missing = {}, 0
    for tid_s, rec in ts.items():
        tid = int(tid_s)
        if not rec.get("screened"):
            continue
        n_cert = n_acc = n_bad = 0
        for lab in screened_labellings(tid):
            k = (tid, tuple(sorted(lab.items())))
            if k in certified:
                n_cert += 1
            elif k in accepted:
                n_acc += 1
            else:
                n_bad += 1
        rec["screened_certified"] = n_cert
        rec["screened_accepted"] = n_acc
        rec["screened_uncertified"] = n_bad
        missing += n_bad
        per[tid] = (n_cert, n_acc, n_bad)
    art.write_text(json.dumps(ts, indent=1) + "\n")
    tot = [sum(v[i] for v in per.values()) for i in range(3)]
    print(f"{len(per)} types with screened labellings: {tot[0]} carry an exact "
          f"certificate, {tot[1]} is the accepted labelling, {tot[2]} have neither")
    return 1 if missing else 0


def merge_logs(paths):
    """Rebuild the artifact from the per-type JSON lines of several runs (the types
    were certified in parallel processes).  A type certified twice must agree."""
    got = {}
    for p in paths:
        for line in Path(p).read_text().splitlines():
            if line.startswith('{"tid"'):
                r = json.loads(line)
                old = got.get(r["tid"])
                if old and (old["verdict"], old.get("drat_trim"),
                            old.get("blocking_clauses")) != (
                        r["verdict"], r.get("drat_trim"), r.get("blocking_clauses")):
                    raise SystemExit(f"type {r['tid']}: runs disagree")
                got[r["tid"]] = r
    art = ROOT / "runs/d6_n10/tree_sat.json"
    art.write_text(json.dumps({str(k): got[k] for k in sorted(got)}, indent=1) + "\n")
    print(f"merged {len(got)} types from {len(paths)} logs -> {art.relative_to(ROOT)}")
    return 0


def main():
    if "merge" in ARG:
        return merge_logs(ARG["merge"].split(","))
    if "--link" in sys.argv[1:]:
        return link_certificates()
    types = {t["type_id"]: t for t in
             json.load(open(ROOT / "runs/d6_n10/stage2/types.json"))}
    flags = json.load(open(ROOT / "runs/d6_n10/type_flags.json"))["d6_n10"]
    valid = sorted(t for t in types if not flags[str(t)]["has_missing_face_of_size_d"])
    sel = ARG.get("tids", "all")
    tids = valid if sel == "all" else [int(x) for x in sel.split(",")]
    outdir = Path(ARG.get("out", ROOT / "runs/d6_n10/tree_sat"))
    outdir.mkdir(parents=True, exist_ok=True)
    drat = ARG.get("drat")
    refuted_ok()

    mne = {k: minimal_non_elliptic(k) for k in range(3, D + 1)}
    lanner = {}
    for k in (3, 4, 5):
        lanner[k], ncls = lanner_labellings(k)
        print(f"order {k}: {len(mne[k])} minimal non-elliptic labellings, "
              f"{len(lanner[k])} Lannér ({ncls} up to isomorphism)", flush=True)
        if k in (4, 5):
            want = {4: 9, 5: 5}[k]
            assert ncls == want, f"order {k}: {ncls} Lannér classes, Lannér lists {want}"
    print(f"order 6: {len(mne[6])} minimal non-elliptic labellings", flush=True)


    results = []
    split = {int(x) for x in ARG.get("split", "").split(",") if x}
    for tid in tids:
        sbset = {int(x) for x in ARG.get("sb", "").split(",") if x}
        r = (run_type_sb(tid, types[tid], mne, lanner, outdir, drat) if tid in sbset
             else run_type_split(tid, types[tid], mne, lanner, outdir, drat)
             if tid in split else run_type(tid, types[tid], mne, lanner, outdir, drat))
        results.append(r)
        print(json.dumps(r), flush=True)
    art = ROOT / "runs/d6_n10/tree_sat.json"
    prev = json.loads(art.read_text()) if art.exists() else {}
    prev.update({str(r["tid"]): r for r in results})
    art.write_text(json.dumps(prev, indent=1) + "\n")
    bad = [r["tid"] for r in results
           if r["verdict"] != "UNSAT" or (drat and r["drat_trim"] != "VERIFIED")]
    print(f"\n{len(results)} types: {len(results) - len(bad)} UNSAT"
          + (" with DRAT proof verified" if drat else "") + f"; problems: {bad}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
