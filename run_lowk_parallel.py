#!/usr/bin/env python3
"""Partition-parallel Stage-4 for low-k types whose label-assignment space is too large
for single-thread enumerate-then-screen.  The label search is partitioned by fixing the
first P pairs of the deterministic ordering (prefix=); every prefix runs on a worker pool,
covering the whole space with no overlap.  Passing candidates are then solved and deduped.

Usage: python3 run_lowk_parallel.py <tid[,tid,...]> [branch_timeout_s] [P]
"""
import sys, json, time, itertools
from pathlib import Path
import numpy as np
import multiprocessing as mp
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline.stage4_gram import (_build_vertex_groups, _build_face_tuple_groups,
    _build_lanner_groups, enumerate_labels_backtrack, _build_minor_index,
    _structured_screen, _refine_mpmath, _recognize_minpoly_and_verify, _GRAM_FLOAT,
    _numerical_screen, _build_gram_numpy, _check_signature_float, VALID_LABELS)
from pipeline.utils.gale_exact import AffineGale
from run_d4 import canonical_key
from sympy import symbols
N, D = 8, 4

def _groups(t):
    mf=[frozenset(m) for m in t["missing_faces"]]
    dotted=[tuple(sorted(m)) for m in mf if len(m)==2]
    ordinary=frozenset((i,j) for i in range(N) for j in range(i+1,N))-frozenset(dotted)
    if t.get("vertex_sets"):
        vs=[frozenset(v) for v in t["vertex_sets"]]
    else:
        vs=AffineGale([tuple(p) for p in t["example_points"]],
                      frozenset(t["example_positive"]), D).vertex_sets()
    vg=_build_vertex_groups(ordinary, vs)
    ft=_build_face_tuple_groups(ordinary, mf, N, max_size=4)
    ex={frozenset(v) for v,_ in vg}; vg=vg+[(v,p) for v,p in ft if frozenset(v) not in ex]
    lg=_build_lanner_groups(ordinary, mf)
    mi=_build_minor_index(dotted, N, D)
    return ordinary, vg, lg, mi, dotted

def _enum_worker(args):
    """Enumerate one prefix-branch + screen; return (candidates, n_assigns, hit_timeout)."""
    t, prefix, branch_timeout = args
    ordinary, vg, lg, mi, dotted = _groups(t)
    cands=[]; n=0; t0=time.time()
    for la in enumerate_labels_backtrack(ordinary, vg, lg, max_count=10**9,
                                         timeout=branch_timeout, prefix=list(prefix)):
        n+=1
        ordf={p:_GRAM_FLOAT[m] for p,m in la.items()}
        dec,xs=_structured_screen(ordf, dotted, mi, N, D)
        if dec:
            for s in xs:
                cands.append((dict(la), dotted, [s[p] for p in dotted]))
        elif dec is None:
            xa,res=_numerical_screen(ordf, dotted, N, D, residual_threshold=1e-6)
            if res<=1e-6 and _check_signature_float(_build_gram_numpy(ordf,dotted,xa,N),D,tol=1e-3):
                cands.append((dict(la), dotted, list(xa)))
    el=time.time()-t0
    return cands, n, (el >= branch_timeout-2)

def _solve_worker(args):
    la, dotted, x0 = args
    if not dotted:
        return [{"label_assignment":{str(p):v for p,v in la.items()}, "dot_values":{}}]
    sl=[symbols(f"x_{p[0]}_{p[1]}", positive=True) for p in dotted]
    xh=_refine_mpmath(np.array(x0), la, dotted, N, D, dps=100)
    if xh is None: return []
    out=[]
    for so in _recognize_minpoly_and_verify(xh, sl, la, dotted, N, D, dps=100, recover_minpoly=False):
        out.append({"label_assignment":{str(p):v for p,v in la.items()}, "dot_values":so})
    return out

def main():
    tids=[int(x) for x in sys.argv[1].split(",")]
    branch_timeout=float(sys.argv[2]) if len(sys.argv)>2 else 2400
    P=int(sys.argv[3]) if len(sys.argv)>3 else 2
    survivors={t["type_id"]:t for t in json.load(open("runs/d4_n8/stage2/types.json"))}
    nproc=max(1, mp.cpu_count()-1)
    nlab=len(VALID_LABELS)
    print(f"low-k parallel: types={tids} P={P} ({nlab**P} prefixes) workers={nproc} "
          f"branch_timeout={branch_timeout}s", flush=True)
    for tid in tids:
        t=survivors[tid]; te=time.time()
        prefixes=[(t, pre, branch_timeout) for pre in itertools.product(range(nlab), repeat=P)]
        all_cands=[]; tot_assign=0; any_timeout=False; nonempty=0
        with mp.get_context("fork").Pool(processes=nproc) as pool:
            for cands,n,hit in pool.imap_unordered(_enum_worker, prefixes, chunksize=1):
                all_cands.extend(cands); tot_assign+=n
                if n>0: nonempty+=1
                if hit: any_timeout=True
            raw=[]
            for r in pool.imap_unordered(_solve_worker, all_cands, chunksize=4):
                raw.extend(r)
        keys=set(canonical_key(r,N) for r in raw)
        status="COMPLETE" if not any_timeout else "INCOMPLETE(timeout)"
        print(f"type {tid}: assigns={tot_assign} cands={len(all_cands)} raw={len(raw)} "
              f"-> {len(keys)} DISTINCT  [{status}, {nonempty} nonempty prefixes] "
              f"in {(time.time()-te)/60:.1f}min", flush=True)

if __name__=="__main__":
    main()
