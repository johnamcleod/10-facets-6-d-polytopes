#!/usr/bin/env python3
"""Validate process_type_stage4 soundness at d=5 with FULL labels (label_indices=None).
Must reproduce the Ma-Zheng per-type census: line 319->22, 312->1, 322->18, 302->6,
313->3, 284->1  (total 51).  tid<->line: 0->319,5->312,6->322,19->302,29->313,63->284.
"""
import json, sys, time
sys.path.insert(0, '.')
import multiprocessing as mp
from pipeline.stage4_gram import process_type_stage4
from pipeline.stage4_blockpaste import _setup
from run_d4 import canonical_key

D5 = {t['type_id']: t for t in json.load(open('runs/d5_n9/stage2/types.json'))}
EXPECT = {0: 22, 5: 1, 6: 18, 19: 6, 29: 3, 63: 1}


def work(tid):
    try:
        t = dict(D5[tid])
        V, *_ = _setup(t)
        t['vertex_sets'] = [sorted(v) for v in V]
        t0 = time.time()
        res = process_type_stage4(t, 5, max_assignments=50_000_000,
                                  enum_timeout=1800.0, solve_timeout=600.0,
                                  label_indices=None, verbose=False)
        if res is None:
            res = []
        keys = set()
        for r in res:
            k = canonical_key(r, 9)
            keys.add(k)
        return (tid, len(keys), round(time.time() - t0, 1), None)
    except Exception as e:
        import traceback
        return (tid, -1, 0.0, f"{type(e).__name__}: {e}")


def main():
    tids = [63, 5, 29, 19, 6, 0]        # easy->hard
    with mp.get_context("fork").Pool(6) as pool:
        for tid, dist, sec, err in pool.imap_unordered(work, tids):
            if err:
                print(f"d5 tid {tid}: ERROR {err}", flush=True)
            else:
                exp = EXPECT[tid]
                tag = "OK" if dist == exp else f"MISMATCH(exp {exp})"
                print(f"d5 tid {tid}: distinct={dist} [{tag}] in {sec}s", flush=True)
    print("FULLLABEL-DONE", flush=True)


if __name__ == "__main__":
    main()
