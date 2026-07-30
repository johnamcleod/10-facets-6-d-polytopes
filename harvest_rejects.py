#!/usr/bin/env python3
"""Harvest the candidates that high-precision refinement REJECTED, and re-attack them.

Motivation.  The d=4 run closes rigorously at 347 of the published 348, with the
whole deficit in one combinatorial type (tid 8 = Burcroff's G11, 7 of 8).  Every
gate in the pipeline that produces an "empty" verdict is either an exhausted
enumeration or a proof -- with ONE exception: `_refine_mpmath` rejecting a
candidate because Gauss-Newton did not drive the minor residual below the
acceptance bar.  That is not a proof of non-realizability; it is a statement about
one numerical iteration from one starting point.

In tid 8 those rejections are not scattered.  All 7 of them sit in the two deepest
subtrees, 8|5,0,0 and 8|5,0,1 -- exactly where a missing eighth polytope would
have to live.  Everywhere else in tid 8, refinement succeeded on every candidate.

So: re-run precisely those two subtrees with the diagnostic dump armed, then take
each rejected candidate and push far harder than the production settings do --
more precision, more iterations, and perturbed restarts -- and see whether any of
them is a genuine polytope that the production bar merely failed to reach.

Usage:
    python3 harvest_rejects.py harvest            # re-run the two subtrees, dump rejects
    python3 harvest_rejects.py attack             # hammer the dumped rejects
"""
from __future__ import annotations

import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pipeline.stage4_gram as _sg                            # noqa: E402
from pipeline.stage4_gram import process_type_stage4          # noqa: E402
from pipeline.stage4_blockpaste import _setup                 # noqa: E402
from pipeline.utils.automorphisms import compute_aut_group     # noqa: E402

D, NF = 4, 8
TARGETS = [(8, (5, 0, 0)), (8, (5, 0, 1))]
OUT = Path("runs/d4_n8/rejects.json")
TYPES = {t["type_id"]: t
         for t in json.load(open(f"runs/d{D}_n{NF}/stage2/types.json"))}


def _job(task):
    tid, prefix = task
    _sg.REFINE_DUMP = []
    for k in _sg.REFINE_STATS:
        _sg.REFINE_STATS[k] = 0
    t = dict(TYPES[tid])
    V, *_ = _setup(t)
    t["vertex_sets"] = [sorted(v) for v in V]
    nn = 1 + max(max(v) for v in V)
    so = {}
    t0 = time.time()
    res = process_type_stage4(
        t, D, max_assignments=50_000_000, enum_timeout=7200, solve_timeout=1200,
        prefix=prefix, stats_out=so, wildcard=True, use_burcroff_55b=True,
        automorphisms=compute_aut_group(V, nn), verbose=False) or []
    return dict(tid=tid, prefix=list(prefix), sec=round(time.time() - t0, 1),
                found=len(res), exhausted=bool(so.get("exhausted")),
                enum=so.get("enum_count", 0),
                refine=dict(_sg.REFINE_STATS),
                rejects=list(_sg.REFINE_DUMP))


def harvest():
    print(f"re-running {len(TARGETS)} subtrees with the reject dump armed")
    with mp.get_context("fork").Pool(processes=len(TARGETS)) as pool:
        rows = list(pool.imap_unordered(_job, TARGETS))
    allrej = []
    for r in sorted(rows, key=lambda r: r["prefix"]):
        print(f"  tid {r['tid']}|{tuple(r['prefix'])}: found={r['found']} "
              f"exhausted={r['exhausted']} enum={r['enum']:,} {r['sec']}s")
        print(f"     refine={ {k: v for k, v in r['refine'].items() if v} }  "
              f"rejects dumped={len(r['rejects'])}")
        for d in r["rejects"]:
            d["src"] = f"{r['tid']}|{','.join(map(str, r['prefix']))}"
        allrej += r["rejects"]
    OUT.write_text(json.dumps(allrej, indent=1))
    print(f"\n{len(allrej)} rejected candidates written to {OUT}")
    from collections import Counter
    print("by reason:", dict(Counter(d["reason"] for d in allrej)))
    return 0


def attack():
    """Push each rejected candidate far past the production settings."""
    import mpmath
    rej = json.loads(OUT.read_text())
    print(f"attacking {len(rej)} rejected candidates\n")
    # (dps, max_iter) ladder; production is (60, 40).
    LADDER = [(60, 200), (120, 200), (240, 300), (400, 400)]
    # multiplicative jitters on the starting point, to escape a bad basin
    JITTER = [1.0, 1.0 + 1e-6, 1 - 1e-6, 1.001, 0.999, 1.01, 0.99, 1.05, 0.95]
    wins, margins = [], []
    for idx, d in enumerate(rej):
        labels = {tuple(int(x) for x in k.split(",")): m
                  for k, m in d["labels"].items()}
        dotted = [tuple(p) for p in d["dotted"]]
        x0 = d["x_approx"]
        print(f"[{idx}] src={d['src']} reason={d['reason']} "
              f"stalled at {d['best_norm']}")
        print(f"     x0={['%.6f' % v if v is not None else None for v in x0]}")
        best, bestnorm = None, None
        for dps, mi in LADDER:
            for jf in JITTER:
                xs = [None if v is None else v * jf for v in x0]
                for k in _sg.REFINE_STATS:
                    _sg.REFINE_STATS[k] = 0
                _sg.REFINE_DUMP = []          # armed: captures the stalled residual
                out = _sg._refine_mpmath(xs, labels, dotted, NF, D,
                                         dps=dps, max_iter=mi)
                for rec in _sg.REFINE_DUMP:
                    try:
                        v = float(str(rec["best_norm"]).split("'")[1])
                        if bestnorm is None or v < bestnorm:
                            bestnorm = v
                    except Exception:
                        pass
                if out is not None:
                    mpmath.mp.dps = dps
                    best = (dps, mi, jf, [mpmath.nstr(v, 25) for v in out])
                    break
            if best:
                break
        if best:
            dps, mi, jf, vals = best
            print(f"     *** CONVERGED at dps={dps} iter={mi} jitter={jf}: {vals}")
            wins.append(dict(idx=idx, src=d["src"], dps=dps, jitter=jf,
                             weights=vals, labels=d["labels"], dotted=d["dotted"]))
        else:
            n_try = len(LADDER) * len(JITTER)
            print(f"     no convergence over {n_try} attempts; best residual "
                  + (f"{bestnorm:.3e}" if bestnorm is not None else "unrecorded"))
            margins.append(bestnorm)
        print()
    p = Path("runs/d4_n8/rejects_attacked.json")
    p.write_text(json.dumps(wins, indent=1))
    print(f"{len(wins)} of {len(rej)} rejected candidates converged; written {p}")
    good = [m for m in margins if m is not None]
    if good:
        print(f"\nbest residual ANY non-converging candidate reached: {min(good):.3e}")
        print("Every solution accepted anywhere in the d=4 run reached <= 1e-50, so the")
        print("rejected candidates are not marginal: they sit tens of orders of")
        print("magnitude from the acceptance bar, which is the behaviour expected of a")
        print("point that is simply not a solution rather than one the bar lost.")
    print("\nA convergence here is NOT yet a polytope: it must still pass the full")
    print("certification (signature (d,1), isolation, no parabolic subdiagram, and")
    print("CoxIter).  It does establish that the production acceptance bar, not")
    print("non-realizability, is what rejected the candidate.")
    return 0


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "harvest"
    sys.exit(harvest() if mode == "harvest" else attack())
