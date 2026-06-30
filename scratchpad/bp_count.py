"""Paste-only candidate counts for given types (no solve) -- feasibility gauge.

For each type: block-paste (partitioned 6-way over a seed column), report the candidate
count, wildcard distribution, and the estimated #label-assignments after wildcard
expansion (7 -> 5 concrete values).  Tells us whether the gap types are tractable for the
exact solve.
"""
import sys, json, time
import numpy as np
sys.path.insert(0, ".")
from pipeline.stage4_blockpaste import paste_candidates, _build_constraints, _order_vertices

T = {t["type_id"]: t for t in json.load(open("runs/d4_n8/stage2/types.json"))}
tids = [int(x) for x in sys.argv[1].split(",")]
CAP = int(sys.argv[2]) if len(sys.argv) > 2 else 6_000_000

for tid in tids:
    t = T[tid]
    dotted = [m for m in t["missing_faces"] if len(m) == 2]
    _, _, ctx = _build_constraints(t)
    order = _order_vertices(ctx["V"], ctx["cols_of"])
    pc = ctx["cols_of"](ctx["V"][order[0]])[0]
    t0 = time.time()
    parts, overflow = [], False
    for v in (2, 3, 4, 5, 6, 7):
        try:
            c, ordn = paste_candidates(t, max_candidates=CAP, forced={pc: v})
            if c.shape[0]:
                parts.append(c)
        except OverflowError as e:
            overflow = True
            print(f"type {tid}: PARTITION col{pc}={v} OVERFLOW ({e})", flush=True)
    dt = time.time() - t0
    if overflow and not parts:
        print(f"type {tid} (k={len(dotted)}): OVERFLOW, no count in {dt:.0f}s", flush=True)
        continue
    cands = np.unique(np.concatenate(parts, axis=0), axis=0) if parts else np.zeros((0, 1))
    nwild = (cands == 7).sum(axis=1)
    u, cc = np.unique(nwild, return_counts=True)
    dist = dict(zip(u.tolist(), cc.tolist()))
    est = int(sum(cnt * (5 ** w) for w, cnt in dist.items()))
    flag = " (PARTIAL-overflow)" if overflow else ""
    print(f"type {tid} (k={len(dotted)}): {cands.shape[0]} candidates, "
          f"wild-dist={dist}, est_expanded={est:,} in {dt:.0f}s{flag}", flush=True)
