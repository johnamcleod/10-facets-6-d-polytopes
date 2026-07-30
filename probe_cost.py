#!/usr/bin/env python3
"""Bounded cost probe for a wildcard-mode classification run.

The full driver (run_survivors_rigorous.py) answers "does this type realize
anything", and for a hard type that can take hours or split into thousands of
subtrees.  This script answers the cheaper question needed BEFORE committing that
compute: which types are cheap, which are the tail, how fast the enumerator runs,
and -- now that the screen carries branch counters -- whether the bounded-box
floating-point fallback is ever reached in this dimension.

Every type is run with the SAME flags as the real driver (wildcard labels,
Burcroff 5.5(b) caps, orbit symmetry breaking) but with a short fixed budget, so
the whole probe is bounded.

Phase A  every candidate type, budget seconds each.
Phase B  for each type that did NOT exhaust in phase A, its 6 depth-1 subtrees,
         budget seconds each.  The fraction of children that exhaust is a
         first-order signal for whether prefix splitting will terminate quickly
         (shallow-but-wide) or the tree is genuinely deep.

Usage:
    python3 probe_cost.py [--d D] [--budget SECONDS] [--nproc N] [--phase-a-only]

Output: a per-type table, the aggregate branch counters, and an explicit
statement of what the probe does and does not determine.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pipeline.stage4_gram as _sg                            # noqa: E402
from pipeline.stage4_gram import process_type_stage4          # noqa: E402
from pipeline.stage4_blockpaste import _setup                 # noqa: E402
from pipeline.utils.automorphisms import compute_aut_group     # noqa: E402
from run_d4 import canonical_key                              # noqa: E402

N_LAB = 6              # wildcard alphabet {2,3,4,5,6,*}
_CFG: dict = {}


def _types_path(d):
    return f"runs/d{d}_n{d + 4}/stage2/types.json"


def work(task):
    tid, prefix = task
    d, budget = _CFG["d"], _CFG["budget"]
    n = d + 4
    try:
        for k in _sg.REFINE_STATS: _sg.REFINE_STATS[k] = 0
        t = dict(_CFG["types"][tid])
        V, *_ = _setup(t)
        t["vertex_sets"] = [sorted(v) for v in V]
        nn = 1 + max(max(v) for v in V)
        so = {}
        t0 = time.time()
        res = process_type_stage4(
            t, d, max_assignments=50_000_000,
            enum_timeout=budget, solve_timeout=budget / 2.0,
            wildcard=True, use_burcroff_55b=True,
            prefix=tuple(prefix) or None,
            automorphisms=compute_aut_group(V, nn),
            stats_out=so) or []
        el = time.time() - t0
        return dict(tid=tid, prefix=prefix, sec=round(el, 1),
                    distinct=len(set(str(canonical_key(r, n)) for r in res)),
                    exhausted=bool(so.get("exhausted")),
                    enum=so.get("enum_count", 0),
                    screened=so.get("screened", 0),
                    passed=so.get("passed_screen", 0),
                    exact=so.get("exact_attempts", 0),
                    wild=so.get("wild_assignments", 0),
                    branches=so.get("screen_branches") or {},
                    refine=dict(_sg.REFINE_STATS),
                    err=None)
    except Exception as e:
        return dict(tid=tid, prefix=prefix, sec=0.0, distinct=0, exhausted=False,
                    enum=0, screened=0, passed=0, exact=0, wild=0, branches={},
                    refine={}, err=f"{type(e).__name__}: {e}")


def run_batch(tasks, nproc):
    out = []
    with mp.get_context("fork").Pool(processes=nproc) as pool:
        for r in pool.imap_unordered(work, tasks):
            out.append(r)
            tag = "exhausted" if r["exhausted"] else "TIMEOUT"
            pfx = "" if not r["prefix"] else f" pfx={r['prefix']}"
            print(f"  tid {r['tid']:>4}{pfx}: {tag:>9}  {r['sec']:>6.1f}s  "
                  f"enum={r['enum']:>9,}  screened={r['screened']:>9,}  "
                  f"exact={r['exact']:>6,}  distinct={r['distinct']}"
                  + (f"  ERROR {r['err']}" if r["err"] else ""), flush=True)
    return out


def agg_branches(rows):
    c = Counter()
    for r in rows:
        c.update(r["branches"])
    return dict(c)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--d", type=int, default=4)
    ap.add_argument("--budget", type=float, default=120.0)
    ap.add_argument("--nproc", type=int, default=max(1, mp.cpu_count() - 1))
    ap.add_argument("--phase-a-only", action="store_true")
    a = ap.parse_args()

    types = {t["type_id"]: t for t in json.load(open(_types_path(a.d)))}
    fp = Path(f"runs/d{a.d}_n{a.d + 4}/facet_profile_survivors.json")
    cands = (json.loads(fp.read_text())["survivors"] if fp.exists()
             else sorted(types))
    _CFG.update(d=a.d, budget=a.budget, types=types)

    print(f"cost probe: d={a.d}, {len(cands)} candidate types, "
          f"{a.budget:g}s budget each, nproc={a.nproc}")
    print(f"worst-case wall clock: phase A "
          f"{len(cands) * a.budget * 1.5 / a.nproc / 60:.0f} min\n")

    print("--- phase A: whole types ---")
    A = run_batch([(t, ()) for t in cands], a.nproc)
    done = [r for r in A if r["exhausted"] and not r["err"]]
    hard = [r for r in A if not r["exhausted"] and not r["err"]]
    errs = [r for r in A if r["err"]]

    print(f"\nphase A: {len(done)}/{len(cands)} types exhausted within budget, "
          f"{len(hard)} hit the budget, {len(errs)} errored")
    if done:
        print(f"   cheap types cost {sum(r['sec'] for r in done)/60:.1f} min total; "
              f"slowest {max(r['sec'] for r in done):.0f}s")
        realizing = {r['tid']: r['distinct'] for r in done if r['distinct']}
        print(f"   realizers among them: {realizing or 'none'}")
    if hard:
        print(f"   the tail: {sorted(r['tid'] for r in hard)}")
        thr = [r["enum"] / r["sec"] for r in hard if r["sec"] > 0]
        if thr:
            print(f"   enumeration throughput on the tail: "
                  f"{min(thr):,.0f}-{max(thr):,.0f} labellings/s")

    B = []
    if hard and not a.phase_a_only:
        print(f"\n--- phase B: depth-1 subtrees of the {len(hard)} tail types ---")
        B = run_batch([(r["tid"], (i,)) for r in hard for i in range(N_LAB)],
                      a.nproc)
        per = {}
        for r in B:
            per.setdefault(r["tid"], []).append(r["exhausted"])
        print("\nphase B: children exhausting within budget, per tail type")
        for tid, flags in sorted(per.items()):
            print(f"   tid {tid:>4}: {sum(flags)}/{len(flags)}")
        tot = sum(sum(v) for v in per.values()); n_all = sum(len(v) for v in per.values())
        print(f"   overall {tot}/{n_all} depth-1 subtrees exhausted at this budget")

    print("\n--- _refine_mpmath outcome counters ---")
    rs = Counter()
    for r in A + B: rs.update(r.get("refine") or {})
    print("   " + json.dumps(dict(rs)))
    resc = rs["rescued_exception"] + rs["rescued_overshoot"]
    print(f"   solutions RESCUED by the 2026-07-28 refinement bugfix: {resc}")
    print("   (non-zero means the pre-fix code would have discarded that many "
          "genuine solutions here)")

    print("\n--- screen branch counters (aggregate) ---")
    br = agg_branches(A + B)
    print("   " + (json.dumps(br) if br else "(none recorded)"))
    fb = br.get("numerical_fallback", 0)
    print(f"   bounded-box [1.001,1000] fallback reached: {fb} times")
    for k in ("stuck_no_pair", "stuck_zero_resultant", "stuck_minor_free",
              "stuck_leaves"):
        if br.get(k):
            print(f"   value-dependent stall {k}: {br[k]}")
    if br.get("cascade_pair"):
        print(f"   pair-resultant (np.polyfit) path used: {br['cascade_pair']} times")

    out = Path(f"runs/d{a.d}_n{a.d + 4}/cost_probe.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"budget": a.budget, "phase_a": A, "phase_b": B,
                               "branches": br}, indent=1, default=str))
    print(f"\nwritten: {out}")
    print("\nWHAT THIS DOES NOT DETERMINE: for the tail types the probe gives a")
    print("lower bound only -- it measures throughput and how readily subtrees")
    print("exhaust, not total tree size, so it cannot produce a total-hours")
    print("figure for them. It is a triage tool, not a cost oracle.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
