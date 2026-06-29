#!/usr/bin/env python3
"""Block-pasting Stage-4 driver (Ma-Zheng method) for d=4 / n=8.

For each combinatorial type:
  1. Generate candidate ordinary-label vectors by block-pasting + library-join
     (``pipeline.stage4_blockpaste.paste_candidates``), partitioned over the 6 labels of a
     seed column so the work fans out across cores with bounded memory.
  2. Expand each candidate's ``7`` wildcards to concrete {7,8,9,10,12}, screen each via the
     shared ``stage4_gram.screen_candidate`` (sound structured cascade + numerical fallback),
     keeping only screen-passing solver starts.
  3. Solve+verify each start exactly (``_refine_mpmath`` + ``_recognize_minpoly_and_verify``:
     signature (4,1) + Jacobian-rank isolation), dedup via ``run_d4.canonical_key``.

This replaces the brute enumerate-then-screen path for candidate generation; the exact
solver is reused unchanged.  Output JSON matches ``run_lowk_parallel`` for downstream reuse.

Usage: python3 run_blockpaste_parallel.py <tid[,tid,...]|all> [nworkers]
"""
import sys, json, time
from pathlib import Path
import numpy as np
import multiprocessing as mp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline.stage4_blockpaste import (paste_candidates, expand_label_assignments,
                                        _build_constraints, _order_vertices)
from pipeline.stage4_gram import (screen_candidate, _build_minor_index, _refine_mpmath,
                                  _recognize_minpoly_and_verify)
from run_d4 import canonical_key
from sympy import symbols

N, D = 8, 4


def _seed_part_col(t):
    """First column of the greedily-ordered seed vertex -- the partition column."""
    _, _, ctx = _build_constraints(t)
    order = _order_vertices(ctx["V"], ctx["cols_of"])
    return ctx["cols_of"](ctx["V"][order[0]])[0]


def _paste_screen_worker(args):
    """Paste one partition (seed column pinned to a value), expand wildcards, screen.
    Returns (list of (la, dotted, x0) starts, n_label_assignments_screened)."""
    t, part_col, val = args
    dotted = [tuple(sorted(m)) for m in t["missing_faces"] if len(m) == 2]
    mi = _build_minor_index(dotted, N, D)
    cands, ordinary = paste_candidates(t, forced={part_col: val})
    starts = []
    n_la = 0
    for la in expand_label_assignments(cands, ordinary):
        n_la += 1
        starts.extend(screen_candidate(la, dotted, mi, N, D))
    return starts, n_la


def _solve_worker(args):
    la, dotted, x0 = args
    if not dotted:
        return [{"label_assignment": {str(p): v for p, v in la.items()}, "dot_values": {}}]
    sl = [symbols(f"x_{p[0]}_{p[1]}", positive=True) for p in dotted]
    xh = _refine_mpmath(np.array(x0), la, dotted, N, D, dps=100)
    if xh is None:
        return []
    out = []
    for so in _recognize_minpoly_and_verify(xh, sl, la, dotted, N, D, dps=100,
                                            recover_minpoly=False):
        out.append({"label_assignment": {str(p): v for p, v in la.items()}, "dot_values": so})
    return out


def process_type(tid, survivors, nproc, outdir):
    t = survivors[tid]
    te = time.time()
    part_col = _seed_part_col(t)
    part_args = [(t, part_col, v) for v in (2, 3, 4, 5, 6, 7)]

    all_starts = []
    tot_la = 0
    with mp.get_context("fork").Pool(processes=min(nproc, 6)) as pool:
        for starts, n_la in pool.imap_unordered(_paste_screen_worker, part_args):
            all_starts.extend(starts)
            tot_la += n_la
    raw = []
    with mp.get_context("fork").Pool(processes=nproc) as pool:
        for r in pool.imap_unordered(_solve_worker, all_starts, chunksize=8):
            raw.extend(r)
    keys = set(canonical_key(r, N) for r in raw)
    elapsed = (time.time() - te) / 60
    print(f"type {tid}: label_assigns={tot_la} screen_pass={len(all_starts)} "
          f"raw={len(raw)} -> {len(keys)} DISTINCT in {elapsed:.1f}min", flush=True)
    res = {"type_id": tid, "distinct": len(keys), "polytopes": raw,
           "label_assigns": tot_la, "screen_pass": len(all_starts),
           "elapsed_min": round(elapsed, 2)}
    (outdir / f"type_{tid}.json").write_text(json.dumps(res, indent=2))
    return res


def main():
    arg = sys.argv[1]
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 else max(1, mp.cpu_count() - 1)
    survivors = {t["type_id"]: t for t in json.load(open("runs/d4_n8/stage2/types.json"))}
    if arg == "all":
        tids = sorted(survivors)
    else:
        tids = [int(x) for x in arg.split(",")]
    outdir = Path("runs/d4_n8/blockpaste")
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"block-paste Stage 4: types={tids} workers={nproc}", flush=True)
    summary = {}
    for tid in tids:
        r = process_type(tid, survivors, nproc, outdir)
        summary[tid] = {k: r[k] for k in ("distinct", "label_assigns", "screen_pass",
                                          "elapsed_min")}
        (outdir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\nResults saved to {outdir}/", flush=True)


if __name__ == "__main__":
    main()
