#!/usr/bin/env python3
"""Parallel d=4 Stage-4 runner: per type, enumerate + structured-screen single-thread
(collect candidate starts), then solve the candidates (mpmath refine + minpoly) across a
worker pool. Fixes the budget-limitation that capped the solve-heavy types (e.g. type 5
recovered 17 -> 92 distinct when fully processed)."""
import sys, json, time
from pathlib import Path
import numpy as np
import multiprocessing as mp
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline.stage4_gram import (_build_vertex_groups, _build_face_tuple_groups,
    _build_lanner_groups, enumerate_labels_backtrack, _build_minor_index,
    _structured_screen, _refine_mpmath, _recognize_minpoly_and_verify, _GRAM_FLOAT,
    _numerical_screen, _build_gram_numpy, _check_signature_float)
from sympy import symbols
from run_d4 import canonical_key

N, D = 8, 4

def _solve_candidate(args):
    """Worker: given (label_assign dict, dotted, x_start list), run the exact solve."""
    la, dotted, x0 = args
    sl = [symbols(f"x_{p[0]}_{p[1]}", positive=True) for p in dotted]
    xh = _refine_mpmath(np.array(x0), la, dotted, N, D, dps=100)
    if xh is None:
        return []
    out = []
    for so in _recognize_minpoly_and_verify(xh, sl, la, dotted, N, D, dps=100):
        out.append({"label_assignment": {str(p): v for p, v in la.items()}, "dot_values": so})
    return out

def candidates_for_type(t, enum_timeout=1800):
    mf = [frozenset(m) for m in t["missing_faces"]]
    dotted = [tuple(sorted(m)) for m in mf if len(m) == 2]
    ordinary = frozenset((i, j) for i in range(N) for j in range(i+1, N)) - frozenset(dotted)
    vs = [frozenset(v) for v in t["vertex_sets"]]
    vg = _build_vertex_groups(ordinary, vs)
    ft = _build_face_tuple_groups(ordinary, mf, N, max_size=4)
    ex = {frozenset(v) for v, _ in vg}
    vg = vg + [(v, p) for v, p in ft if frozenset(v) not in ex]
    lg = _build_lanner_groups(ordinary, mf)
    mi = _build_minor_index(dotted, N, D)
    cands = []
    for la in enumerate_labels_backtrack(ordinary, vg, lg, max_count=10**9,
                                         timeout=enum_timeout, label_indices=None):
        ordf = {p: _GRAM_FLOAT[m] for p, m in la.items()}
        if not dotted:
            G = _build_gram_numpy(ordf, [], np.array([]), N)
            ev = np.linalg.eigvalsh(G)
            if int(np.sum(np.abs(ev) > 1e-8)) == D+1 and _check_signature_float(G, D):
                cands.append((dict(la), dotted, []))
            continue
        dec, xs = _structured_screen(ordf, dotted, mi, N, D)
        if dec is None:
            xa, res = _numerical_screen(ordf, dotted, N, D, residual_threshold=1e-6)
            if res <= 1e-6 and _check_signature_float(_build_gram_numpy(ordf, dotted, xa, N), D, tol=1e-3):
                cands.append((dict(la), dotted, list(xa)))
        elif dec:
            for s in xs:
                cands.append((dict(la), dotted, [s[p] for p in dotted]))
    return cands

def main():
    survivors = json.load(open("runs/d4_n8/stage3/surviving_types.json"))
    nproc = max(1, mp.cpu_count() - 2)
    print(f"d=4 PARALLEL run: {len(survivors)} types, {nproc} solve workers")
    t0 = time.time()
    all_distinct = set()
    per_type = {}
    with mp.get_context("fork").Pool(processes=nproc) as pool:
        for t in survivors:
            tid = t["type_id"]; te = time.time()
            cands = candidates_for_type(t)
            raw = []
            for r in pool.imap_unordered(_solve_candidate, cands, chunksize=4):
                raw.extend(r)
            keys = set(canonical_key(r, N) for r in raw)
            per_type[tid] = len(keys); all_distinct |= keys
            print(f"  type {tid}: {len(cands)} cands -> {len(raw)} raw -> {len(keys)} distinct "
                  f"({time.time()-te:.0f}s, cum {len(all_distinct)})", flush=True)
    print(f"\n{'='*56}\nd=4 PARALLEL RESULT: {len(all_distinct)} distinct polytopes "
          f"(target 348) in {(time.time()-t0)/60:.1f} min\n{'='*56}")
    Path("runs/d4_n8/stage4_par").mkdir(parents=True, exist_ok=True)
    (Path("runs/d4_n8/stage4_par")/"summary.json").write_text(json.dumps(
        {"distinct": len(all_distinct), "target": 348, "per_type": per_type}, indent=2))

if __name__ == "__main__":
    main()
