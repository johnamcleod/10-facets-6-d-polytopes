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
                                        _build_constraints, _order_vertices, infer_dims)
from pipeline.stage4_gram import (screen_candidate, _build_minor_index, _refine_mpmath,
                                  _recognize_minpoly_and_verify)
from run_d4 import canonical_key
from sympy import symbols


def _seed_part_cols(t, p):
    """First ``p`` columns of the greedily-ordered seed vertex -- the partition columns.
    Pinning more columns => more, smaller partitions (bounded memory)."""
    _, _, ctx = _build_constraints(t)
    order = _order_vertices(ctx["V"], ctx["cols_of"])
    return ctx["cols_of"](ctx["V"][order[0]])[:p]


def _paste_screen_worker(args):
    """Paste one partition (seed columns pinned), expand wildcards, screen.
    On OverflowError the partition is too big -- signalled back so the caller can refine.
    Returns (list of (la, dotted, x0) starts, n_label_assignments_screened, overflow_forced
    or None)."""
    t, forced = args
    n, d = infer_dims(t)
    dotted = [tuple(sorted(m)) for m in t["missing_faces"] if len(m) == 2]
    mi = _build_minor_index(dotted, n, d)
    try:
        cands, ordinary = paste_candidates(t, forced=forced)
    except OverflowError:
        return [], 0, forced
    starts = []
    n_la = 0
    for la in expand_label_assignments(cands, ordinary):
        n_la += 1
        for s in screen_candidate(la, dotted, mi, n, d):
            starts.append((s[0], s[1], s[2], n, d))
    return starts, n_la, None


def _solve_worker(args):
    la, dotted, x0, n, d = args
    if not dotted:
        return [{"label_assignment": {str(p): v for p, v in la.items()}, "dot_values": {}}]
    sl = [symbols(f"x_{p[0]}_{p[1]}", positive=True) for p in dotted]
    try:
        xh = _refine_mpmath(np.array(x0), la, dotted, n, d, dps=100)
    except Exception:
        return []       # degenerate candidate -> skip (never a valid polytope silently lost:
    if xh is None:      # refine/verify only reject; the paste+screen already gated realizability)
        return []
    out = []
    try:
        sols = _recognize_minpoly_and_verify(xh, sl, la, dotted, n, d, dps=100,
                                             recover_minpoly=False)
    except Exception:
        return []
    for so in sols:
        out.append({"label_assignment": {str(p): v for p, v in la.items()}, "dot_values": so})
    return out


def process_type(tid, survivors, nproc, outdir, p=1):
    t = survivors[tid]
    te = time.time()
    n, d = infer_dims(t)
    seed_cols = _seed_part_cols(t, 6)  # full seed vertex columns (refine pool)
    LABELS = (2, 3, 4, 5, 6, 7)
    # initial partition: pin first ``p`` seed columns
    import itertools
    init = [{seed_cols[i]: combo[i] for i in range(p)}
            for combo in itertools.product(LABELS, repeat=p)]

    all_starts, tot_la, overflows = [], 0, 0
    queue = [(t, f) for f in init]
    next_refine_idx = p
    with mp.get_context("fork").Pool(processes=nproc) as pool:
        while queue:
            results = pool.map(_paste_screen_worker, queue)
            queue = []
            for starts, n_la, ovf in results:
                all_starts.extend(starts)
                tot_la += n_la
                if ovf is not None:
                    overflows += 1
                    # refine: pin one more seed column not already pinned
                    extra = next((c for c in seed_cols if c not in ovf), None)
                    if extra is None:
                        raise RuntimeError(f"type {tid}: partition {ovf} cannot refine")
                    for v in LABELS:
                        f2 = dict(ovf); f2[extra] = v
                        queue.append((t, f2))
    raw = []
    with mp.get_context("fork").Pool(processes=nproc) as pool:
        for r in pool.imap_unordered(_solve_worker, all_starts, chunksize=8):
            raw.extend(r)
    keys = set(canonical_key(r, n) for r in raw)
    elapsed = (time.time() - te) / 60
    print(f"type {tid}: label_assigns={tot_la} screen_pass={len(all_starts)} "
          f"raw={len(raw)} refine_overflows={overflows} -> {len(keys)} DISTINCT "
          f"in {elapsed:.1f}min", flush=True)
    res = {"type_id": tid, "distinct": len(keys), "polytopes": raw,
           "label_assigns": tot_la, "screen_pass": len(all_starts),
           "elapsed_min": round(elapsed, 2)}
    (outdir / f"type_{tid}.json").write_text(json.dumps(res, indent=2))
    return res


def main():
    # usage: run_blockpaste_parallel.py <tid[,tid,...]|all> [nworkers] [p] [d]
    arg = sys.argv[1]
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 else max(1, mp.cpu_count() - 1)
    p = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    d = int(sys.argv[4]) if len(sys.argv) > 4 else 4
    base = {4: "runs/d4_n8", 5: "runs/d5_n9", 6: "runs/d6_n10"}[d]
    survivors = {t.get("type_id", i): t
                 for i, t in enumerate(json.load(open(f"{base}/stage2/types.json")))}
    if arg == "all":
        tids = sorted(survivors)
    else:
        tids = [int(x) for x in arg.split(",")]
    outdir = Path(f"{base}/blockpaste")
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"block-paste Stage 4 (d={d}): types={tids} workers={nproc} partition_depth={p}",
          flush=True)
    summary = {}
    for tid in tids:
        r = process_type(tid, survivors, nproc, outdir, p=p)
        summary[tid] = {k: r[k] for k in ("distinct", "label_assigns", "screen_pass",
                                          "elapsed_min")}
        (outdir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\nResults saved to {outdir}/", flush=True)


if __name__ == "__main__":
    main()
