#!/usr/bin/env python3
"""Definitive wildcard-mode validation: run ALL 109 d=5/n=9 types (k>=2) through
process_type_stage4 in wildcard mode ({2..6,7=any m>=7}, no label cap) and require
the census to come out EXACTLY as Ma-Zheng's: 6 realizing types
{0:22, 5:1, 6:18, 19:6, 29:3, 63:1} = 51 polytopes, every other type a rigorous
(exhausted, no wild_unbounded) zero.

This validates in one shot: the enumeration alphabet, the wildcard range analysis
(m=10 edges must be DISCOVERED, not assumed), the exact certification, and the
exhaustion accounting — the full soundness chain reused for the d=6 classification.

Usage: python3 validate_d5_wildcard.py [nproc]
"""
import json, sys, time
from pathlib import Path
import multiprocessing as mp
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline.stage4_gram import process_type_stage4
from pipeline.stage4_blockpaste import _setup
from run_d4 import canonical_key

EXPECT = {0: 22, 5: 1, 6: 18, 19: 6, 29: 3, 63: 1}
OUT = Path("runs/d5_n9/wildcard_validation.json")


def work(t):
    tid = t['type_id']
    try:
        t = dict(t); V, *_ = _setup(t); t['vertex_sets'] = [sorted(v) for v in V]
        so = {}
        t0 = time.time()
        res = process_type_stage4(t, 5, max_assignments=50_000_000,
                                  enum_timeout=7200.0, solve_timeout=1800.0,
                                  wildcard=True, stats_out=so, verbose=False) or []
        keys = set(canonical_key(r, 9) for r in res)
        high = sorted(set(int(m) for r in res
                          for m in r['label_assignment'].values() if int(m) >= 7))
        return (tid, len(keys), bool(so.get('exhausted')),
                bool(so.get('wild_unbounded')), so.get('wild_assignments', 0),
                high, round(time.time() - t0, 1), None)
    except Exception as e:
        return (tid, -1, False, False, 0, [], 0.0, f"{type(e).__name__}: {e}")


def main():
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else max(1, mp.cpu_count() - 1)
    types = json.load(open('runs/d5_n9/stage2/types.json'))
    done = json.loads(OUT.read_text()) if OUT.exists() else {}
    todo = [t for t in types if str(t['type_id']) not in done]
    # heaviest known types last-started = worst packing; put them first
    heavy = {0: 3, 6: 2, 19: 1}
    todo.sort(key=lambda t: -heavy.get(t['type_id'], 0))
    print(f"d5 wildcard validation: {len(todo)} types to run ({len(done)} cached), "
          f"nproc={nproc}", flush=True)
    with mp.get_context("fork").Pool(processes=nproc) as pool:
        for tid, dist, exh, unb, wilds, high, sec, err in pool.imap_unordered(work, todo):
            exp = EXPECT.get(tid, 0)
            ok = (dist == exp) and exh and not unb and not err
            tag = "OK" if ok else f"FAIL(exp {exp}, exh={exh}, unb={unb}, err={err})"
            print(f"d5 tid {tid}: distinct={dist} wild={wilds} high={high} "
                  f"[{tag}] in {sec}s", flush=True)
            done[str(tid)] = {"distinct": dist, "exhausted": exh, "unbounded": unb,
                              "high_labels": high, "sec": sec, "err": err}
            OUT.write_text(json.dumps(done, indent=1))
    total = sum(v["distinct"] for v in done.values() if v["distinct"] > 0)
    bad = {k: v for k, v in done.items()
           if v["distinct"] != EXPECT.get(int(k), 0) or not v["exhausted"]
           or v["unbounded"] or v["err"]}
    print(f"\nD5-WILDCARD-DONE. total={total} (expect 51)  "
          f"problem-types={sorted(bad) if bad else 'NONE'}", flush=True)


if __name__ == "__main__":
    main()
