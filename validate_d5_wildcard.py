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
import pipeline.stage4_gram as _sg
from pipeline.stage4_gram import process_type_stage4
from pipeline.stage4_blockpaste import _setup
from pipeline.utils.automorphisms import compute_aut_group
from run_d4 import canonical_key

EXPECT = {0: 22, 5: 1, 6: 18, 19: 6, 29: 3, 63: 1}

# Pass "55b" as the second argument to additionally enable Burcroff Lemma 5.5(b)
# low-weight caps -- the one solver flag the d=6 run sets that the baseline d=5
# validation did not exercise.  Separate output file so the two runs coexist.
_SMOKE = "smoke" in sys.argv[1:]
USE_55B = "55b" in sys.argv[1:]
# out=NAME writes to runs/d5_n9/NAME.json instead of the default.  This matters:
# main() RESUMES from OUT, so re-running after a code change against the committed
# output would silently skip all 109 cached types and report the OLD numbers as if
# freshly computed.  Always pass out= when validating a code change.
_OUT_ARG = next((a[4:] for a in sys.argv[1:] if a.startswith("out=")), None)
OUT = Path(f"runs/d5_n9/{_OUT_ARG}.json") if _OUT_ARG else Path(
    "runs/d5_n9/wildcard_validation_55b.json" if USE_55B
    else "runs/d5_n9/wildcard_validation.json")


def work(t):
    tid = t['type_id']
    try:
        for _k in _sg.REFINE_STATS: _sg.REFINE_STATS[_k] = 0
        for _k in _sg.WILD_STATS: _sg.WILD_STATS[_k] = 0
        t = dict(t); V, *_ = _setup(t); t['vertex_sets'] = [sorted(v) for v in V]
        n = 1 + max(max(v) for v in V)
        automorphisms = compute_aut_group(V, n)
        so = {}
        t0 = time.time()
        res = process_type_stage4(t, 5, max_assignments=50_000_000,
                                  enum_timeout=7200.0, solve_timeout=1800.0,
                                  wildcard=True, stats_out=so,
                                  use_burcroff_55b=USE_55B,
                                  automorphisms=automorphisms, verbose=False) or []
        keys = set(canonical_key(r, 9) for r in res)
        high = sorted(set(int(m) for r in res
                          for m in r['label_assignment'].values() if int(m) >= 7))
        return (tid, len(keys), bool(so.get('exhausted')),
                bool(so.get('wild_unbounded')), so.get('wild_assignments', 0),
                high, round(time.time() - t0, 1), None, dict(_sg.REFINE_STATS),
                so.get('wild_assignments', 0), dict(_sg.WILD_STATS))
    except Exception as e:
        return (tid, -1, False, False, 0, [], 0.0, f"{type(e).__name__}: {e}", {}, 0, {})


def main():
    nproc = next((int(a) for a in sys.argv[1:] if a.isdigit()),
                 max(1, mp.cpu_count() - 1))
    types = json.load(open("runs/d5_n9/stage2/types.json"))
    if _SMOKE: types = [t for t in types if t["type_id"] in (29, 63)]
    done = json.loads(OUT.read_text()) if OUT.exists() else {}
    if done:
        print(f"!! RESUMING from {OUT}: {len(done)} types are CACHED and will NOT "
              f"be recomputed.\n!! If you changed the solver, pass out=NEWNAME "
              f"instead -- otherwise this reports stale results.", flush=True)
    todo = [t for t in types if str(t['type_id']) not in done]
    # heaviest known types last-started = worst packing; put them first
    heavy = {0: 3, 6: 2, 19: 1}
    todo.sort(key=lambda t: -heavy.get(t['type_id'], 0))
    print(f"d5 wildcard validation{' [+Burcroff 5.5(b)]' if USE_55B else ''}: "
          f"{len(todo)} types to run ({len(done)} cached), nproc={nproc}",
          flush=True)
    with mp.get_context("fork").Pool(processes=nproc) as pool:
        agg = {}
        for (tid, dist, exh, unb, wilds, high, sec, err, rst, wassign,
             wstats) in pool.imap_unordered(work, todo):
            for _k, _v in rst.items(): agg[_k] = agg.get(_k, 0) + _v
            exp = EXPECT.get(tid, 0)
            ok = (dist == exp) and exh and not unb and not err
            tag = "OK" if ok else f"FAIL(exp {exp}, exh={exh}, unb={unb}, err={err})"
            print(f"d5 tid {tid}: distinct={dist} wild={wilds} high={high} "
                  f"[{tag}] in {sec}s", flush=True)
            done[str(tid)] = {"distinct": dist, "exhausted": exh, "unbounded": unb,
                              "high_labels": high, "sec": sec, "err": err,
                              "refine": rst, "wild_assignments": wassign,
                              "wild": {k: v for k, v in wstats.items() if v}}
            OUT.write_text(json.dumps(done, indent=1))
    total = sum(v["distinct"] for v in done.values() if v["distinct"] > 0)
    bad = {k: v for k, v in done.items()
           if v["distinct"] != EXPECT.get(int(k), 0) or not v["exhausted"]
           or v["unbounded"] or v["err"]}
    print(f"\n_refine_mpmath outcomes over the whole d=5 census: {agg}")
    resc = agg.get('rescued_exception',0)+agg.get('rescued_overshoot',0)
    print(f"solutions RESCUED by the refinement bugfix: {resc}  "
          f"(non-zero => the pre-fix code discarded genuine solutions here)")
    print(f"\nD5-WILDCARD-DONE. total={total} (expect 51)  "
          f"problem-types={sorted(bad) if bad else 'NONE'}", flush=True)


if __name__ == "__main__":
    main()
