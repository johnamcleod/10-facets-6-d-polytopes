#!/usr/bin/env python3
"""Resolve the d=6 block-paste-INCOMPLETE types with the backtracking Stage-4 solver.

Block-paste explodes on these types because it enumerates labels {2..7} + wildcard
{7,8,9,10,12}; but Burcroff bounds d>=5 dihedral angles to m in {2,3,4,5}, so the
correct label set is exactly VALID_LABELS[0:4].  process_type_stage4 with
label_indices=(0,1,2,3) uses bitmask forward-checking over that set -> sound, complete
for d=6, and tractable (type 2: 74s definitive 0 vs block-paste 15min INCOMPLETE).

Each distinct>0 result is a candidate SECOND compact hyperbolic Coxeter 6-polytope
with 10 facets.  Recovers P^B6 on type 379 as the soundness anchor.
"""
import json, sys, time
from pathlib import Path
import multiprocessing as mp
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline.stage4_gram import process_type_stage4
from pipeline.stage4_blockpaste import _setup
from run_d4 import canonical_key

D = 6
TIMEOUT = 900.0
MAXA = 5_000_000
LABELS = (0, 1, 2, 3)     # VALID_LABELS indices -> m in {2,3,4,5}, complete for d=6


def _with_vertex_sets(t):
    V, *_ = _setup(t)
    t = dict(t)
    t["vertex_sets"] = [sorted(v) for v in V]
    return t


def _worker(t):
    tid = t["type_id"]
    t0 = time.time()
    try:
        res = process_type_stage4(t, D, max_assignments=MAXA, enum_timeout=TIMEOUT,
                                  solve_timeout=300.0, label_indices=LABELS, verbose=False)
    except Exception as e:
        return tid, {"type_id": tid, "distinct": -1, "error": repr(e)[:200],
                     "elapsed_s": round(time.time() - t0, 1)}
    n = 10
    keys = set(canonical_key(r, n) for r in res)
    return tid, {"type_id": tid, "distinct": len(keys), "raw": len(res),
                 "elapsed_s": round(time.time() - t0, 1)}


def main():
    ids = json.load(open("runs/d6_n10/blockpaste/incomplete_ids.json"))
    types = {t["type_id"]: t for t in json.load(open("runs/d6_n10/stage2/types.json"))}
    todo = [_with_vertex_sets(types[i]) for i in ids]
    outdir = Path("runs/d6_n10/backtrack_incomplete")
    outdir.mkdir(parents=True, exist_ok=True)
    nproc = max(1, mp.cpu_count() - 1)
    print(f"backtracking resolve: {len(todo)} types, labels m in {{2,3,4,5}}, "
          f"timeout={TIMEOUT}s workers={nproc}", flush=True)
    summary = {}
    with mp.get_context("fork").Pool(processes=nproc) as pool:
        for tid, r in pool.imap_unordered(_worker, todo):
            tag = "FOUND" if r["distinct"] > 0 else ("ERR" if r["distinct"] < 0 else "empty")
            print(f"type {tid}: distinct={r['distinct']} raw={r.get('raw','-')} "
                  f"[{tag}] in {r['elapsed_s']}s", flush=True)
            summary[tid] = r
            (outdir / "summary.json").write_text(json.dumps(summary, indent=2))
    found = {k: v for k, v in summary.items() if v["distinct"] > 0}
    errs = {k: v for k, v in summary.items() if v["distinct"] < 0}
    print(f"\nDONE. realizing types={sorted(found)}  errors={sorted(errs)}", flush=True)
    print(f"Results -> {outdir}/summary.json", flush=True)


if __name__ == "__main__":
    main()
