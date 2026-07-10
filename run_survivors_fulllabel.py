#!/usr/bin/env python3
"""Re-run the 54 facet-filter survivors of d=6 with FULL labels (label_indices=None).
Earlier verdicts used the WRONG restricted labels {2,3,4,5}; d=5/d=6 realizations use
labels up to 10/12 (e.g. d=5 type 302 has an m=10 edge).  Any survivor with distinct>0
is a candidate 2nd compact Coxeter 6-polytope (besides P^B6=379).
"""
import json, sys, time
sys.path.insert(0, '.')
import multiprocessing as mp
from pipeline.stage4_gram import process_type_stage4
from pipeline.stage4_blockpaste import _setup
from run_d4 import canonical_key

SURV = [8, 12, 17, 34, 36, 38, 40, 51, 55, 59, 60, 61, 69, 70, 92, 103, 120, 127, 132,
        140, 154, 159, 162, 168, 173, 206, 214, 218, 220, 229, 234, 239, 255, 265, 273,
        284, 286, 287, 295, 297, 308, 315, 317, 320, 329, 332, 344, 352, 354, 356, 360,
        378, 379, 382]
TYPES = {t['type_id']: t for t in json.load(open('runs/d6_n10/stage2/types.json'))}


ENUM_TO = 3600.0
SOLVE_TO = 600.0


def work(tid):
    try:
        t = dict(TYPES[tid]); V, *_ = _setup(t); t['vertex_sets'] = [sorted(v) for v in V]
        t0 = time.time()
        res = process_type_stage4(t, 6, max_assignments=50_000_000,
                                  enum_timeout=ENUM_TO, solve_timeout=SOLVE_TO,
                                  label_indices=None, verbose=False) or []
        el = time.time() - t0
        keys = set(canonical_key(r, 10) for r in res)
        # process_type_stage4 breaks its loop at enum_timeout+solve_timeout; if we ran
        # that long the enumeration did NOT exhaust -> a 0 here is NOT rigorous.
        complete = el < (ENUM_TO + SOLVE_TO) * 0.97
        return (tid, len(keys), len(res), round(el, 1), None, complete)
    except Exception as e:
        return (tid, -1, 0, 0.0, f"{type(e).__name__}: {e}", False)


def main():
    import os
    out = {}
    with mp.get_context("fork").Pool(processes=max(1, mp.cpu_count() - 1)) as pool:
        for tid, dist, raw, sec, err, complete in pool.imap_unordered(work, SURV):
            if err:
                tag = "ERROR: " + err
            elif dist > 0:
                tag = "**REALIZES**"
            elif complete:
                tag = "empty(COMPLETE)"
            else:
                tag = "empty(TIMEOUT-not-rigorous)"
            print(f"d6 tid {tid}: distinct={dist} raw={raw} [{tag}] in {sec}s", flush=True)
            out[tid] = {"distinct": dist, "raw": raw, "sec": sec, "err": err, "complete": complete}
            json.dump(out, open("runs/d6_n10/survivors_fulllabel.json", "w"), indent=2)
    realizing = sorted(k for k, v in out.items() if v["distinct"] > 0)
    errored = sorted(k for k, v in out.items() if v["distinct"] < 0)
    timeouts = sorted(k for k, v in out.items() if v["distinct"] == 0 and not v["complete"])
    print(f"\nSURVIVORS-DONE. realizing={realizing}  errored={errored}  "
          f"timeout-not-rigorous={timeouts}", flush=True)


if __name__ == "__main__":
    main()
