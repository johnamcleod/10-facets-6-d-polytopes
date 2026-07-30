#!/usr/bin/env python3
"""Re-run the cheap d=6 types with orbit symmetry breaking DISABLED.

Symmetry breaking discards a partial labelling when it is not lexicographically
minimal among the images of its prefix under the automorphisms stabilising the
prefix's domain.  The soundness argument for that is the standard lex-leader
sketch, and a sketch is not a proof: a failure would be silent, discarding an orbit
whose canonical member is never scheduled.

The cheap way to settle it empirically is to turn the pruning off and check that the
verdicts do not change.  That is what this does, for every d=6 type that cost under
one worker-hour in the run of record.  The expensive types are excluded because
without pruning they do not terminate in reasonable time, which is why the device
exists.

Run:  python3 paper/checks/symmetry_breaking_off.py
"""
from __future__ import annotations

import collections
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pipeline.stage4_gram import process_type_stage4          # noqa: E402
from pipeline.stage4_blockpaste import _setup                 # noqa: E402
from pipeline.utils.automorphisms import compute_aut_group     # noqa: E402
from run_d4 import canonical_key                              # noqa: E402

D, NF, BUDGET = 6, 10, 900.0
TYPES = {t["type_id"]: t for t in
         json.load(open(ROOT / "runs/d6_n10/stage2/types.json"))}


def committed():
    st = {}
    for d in ("d6_tangency", "d6_rest"):
        for line in (ROOT / f"runs/d6_n10/{d}/state.jsonl").read_text().splitlines():
            line = line.strip()
            if line:
                r = json.loads(line)
                st[r.pop("key")] = r
    per, keys = collections.defaultdict(float), collections.defaultdict(set)
    for k, v in st.items():
        t = int(k.split("|")[0])
        per[t] += v.get("sec", 0)
        keys[t].update(v.get("keys", []))
    return per, keys


def work(task):
    tid, use_aut = task
    t = dict(TYPES[tid])
    V, *_ = _setup(t)
    t["vertex_sets"] = [sorted(v) for v in V]
    nn = 1 + max(max(v) for v in V)
    aut = compute_aut_group(V, nn) if use_aut else None
    so = {}
    t0 = time.time()
    res = process_type_stage4(t, D, max_assignments=50_000_000,
                              enum_timeout=BUDGET, solve_timeout=BUDGET / 2,
                              wildcard=True, use_burcroff_55b=True,
                              automorphisms=aut, stats_out=so, verbose=False) or []
    return (tid, use_aut, len({str(canonical_key(r, NF)) for r in res}),
            bool(so.get("exhausted")), so.get("enum_count", 0),
            round(time.time() - t0, 1))


def main():
    per, keys = committed()
    cheap = sorted(t for t, s in per.items() if s < 3600)
    print(f"d=6 types costing under one worker-hour in the run of record: {len(cheap)}")
    print(f"excluded as too expensive without pruning: "
          f"{sorted(t for t, s in per.items() if s >= 3600)}\n")
    tasks = [(t, False) for t in cheap]
    nproc = max(1, mp.cpu_count() - 1)
    out = {}
    with mp.get_context("fork").Pool(processes=nproc) as pool:
        for tid, ua, dist, exh, enum, sec in pool.imap_unordered(work, tasks):
            out[tid] = (dist, exh, enum, sec)
            ref = len(keys[tid])
            ok = dist == ref and exh
            print(f"  tid {tid:>4}  pruning OFF: {dist} polytopes, exhausted={exh}, "
                  f"enum={enum:,}, {sec}s   committed={ref}   "
                  f"{'OK' if ok else 'MISMATCH'}", flush=True)
    bad = [t for t, (d, e, _, _) in out.items() if d != len(keys[t]) or not e]
    print(f"\n{len(out)} types re-run with symmetry breaking disabled")
    print("RESULT:", "every verdict reproduced without the pruning"
          if not bad else f"MISMATCH on {bad}")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
