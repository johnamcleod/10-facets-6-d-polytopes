#!/usr/bin/env python3
"""Margin audit of the screen's floating-point decisions (Section 8, item 4).

Every "this labelling is not realizable" verdict is ultimately a float64 quantity
compared against a fixed tolerance.  The verdict is trustworthy exactly to the
extent that the quantity is not NEAR the tolerance.  This script instruments
every such comparison and reports how much room there was.

Comparisons instrumented, with the threshold the code uses:

  cascade, _quad_roots_gt1  (pin a dashed weight from a vanishing 8x8 minor,
                             which is a quadratic  a x^2 + b x + c  in that weight)
      |a| / scale            vs 1e-12    "is the minor really quadratic in x"
      |c| / scale            vs 1e-9     "is the minor identically zero"
      discriminant / scale   vs 0        "are the roots real"
      |root - 1|             vs 1e-7     "is the root a legal weight x > 1"
  cascade leaf, _check_signature_float
      smallest |eigenvalue| that must vanish   vs 1e-7 / 1e-6
  wildcard, _wild_feasible
      attained multistart L-BFGS-B floor       vs 1e-6

For each we report the minimum and low percentiles over all decisions made.  A
minimum many orders of magnitude above the threshold means an a-priori
backward-error bound on the float64 determinant/eigenvalue computation would
certify every one of these verdicts without redoing the search in exact
arithmetic; a minimum close to the threshold identifies exactly which decisions
need exact escalation.

Usage:
    python3 checks/screen_margins.py 308 59 320 ...        # d=6 types
    python3 checks/screen_margins.py --d5 5 63 29 19       # d=5 types

With no type arguments it runs the cheapest d=6 survivor types (those whose full
classification run cost under one CPU-hour), which is a few CPU-minutes.
"""
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import pipeline.stage4_gram as sg                       # noqa: E402
from pipeline.stage4_blockpaste import _setup           # noqa: E402
from pipeline.utils.automorphisms import compute_aut_group  # noqa: E402

# d=6 survivor types whose full run cost under 1 CPU-hour 
CHEAP_D6 = [8, 12, 17, 36, 51, 59, 60, 69, 127, 140, 154, 162, 168, 214, 220,
            239, 255, 273, 287, 308, 317, 320, 332, 344, 352, 354, 356, 382]

REC = {"quad_a": [], "quad_c": [], "quad_disc": [], "quad_root": [],
       "leaf_eig": [], "wild_resid": []}

THRESH = {"quad_a": 1e-12, "quad_c": 1e-9, "quad_disc": 0.0,
          "quad_root": 1e-7, "leaf_eig": 1e-7, "wild_resid": 1e-6}

_orig_quad = sg._quad_roots_gt1
_orig_sig = sg._check_signature_float
_orig_wild = sg._wild_feasible


def quad(Gsub, ii, jj, eps=1e-7):
    def det_at(x):
        M = Gsub.copy()
        M[ii, jj] = M[jj, ii] = -x
        return np.linalg.det(M)
    c = det_at(0.0)
    a = (det_at(2.0) - 2.0 * det_at(1.0) + c) / 2.0
    b = det_at(1.0) - c - a
    scale = max(1.0, abs(a), abs(b), abs(c))
    REC["quad_a"].append(abs(a) / scale)
    REC["quad_c"].append(abs(c) / scale)
    if abs(a) >= 1e-12:
        disc = b * b - 4 * a * c
        REC["quad_disc"].append(abs(disc) / max(1.0, b * b, abs(4 * a * c)))
        if disc >= 0:
            r = math.sqrt(disc)
            for x in ((-b + r) / (2 * a), (-b - r) / (2 * a)):
                REC["quad_root"].append(abs(x - 1.0))
    return _orig_quad(Gsub, ii, jj, eps)


def sigchk(G, d, tol=1e-8):
    ev = np.sort(np.linalg.eigvalsh(G))
    nz = G.shape[0] - d - 1                    # eigenvalues that must vanish
    if nz > 0:
        REC["leaf_eig"].append(min(abs(ev[i]) for i in range(1, 1 + nz)))
    return _orig_sig(G, d, tol)


def wild(*a, **k):
    out = _orig_wild(*a, **k)
    try:
        REC["wild_resid"].append(float(out[0]))
    except Exception:
        pass
    return out


sg._quad_roots_gt1 = quad
sg._check_signature_float = sigchk
sg._wild_feasible = wild


def run(tids, d):
    n = d + 4
    src = (ROOT / f"runs/d{d}_n{n}/stage2/types.json")
    TYPES = {t["type_id"]: t for t in json.load(open(src))}
    tot_screened = tot_wild = 0
    t0 = time.time()
    for tid in tids:
        t = dict(TYPES[tid])
        V, *_ = _setup(t)
        t["vertex_sets"] = [sorted(v) for v in V]
        nn = 1 + max(max(v) for v in V)
        so = {}
        res = sg.process_type_stage4(
            t, d, max_assignments=50_000_000, enum_timeout=3600.0,
            solve_timeout=600.0, wildcard=True, use_burcroff_55b=True,
            automorphisms=compute_aut_group(V, nn), stats_out=so) or []
        tot_screened += so.get("screened", 0)
        tot_wild += so.get("wild_assignments", 0)
        print(f"  tid {tid:>4}: screened={so.get('screened',0):>7} "
              f"wild={so.get('wild_assignments',0):>7} found={len(res)} "
              f"exhausted={so.get('exhausted')}", flush=True)
    return tot_screened, tot_wild, time.time() - t0


def main():
    args = [a for a in sys.argv[1:] if a != "--d5"]
    d = 5 if "--d5" in sys.argv[1:] else 6
    tids = [int(a) for a in args] if args else CHEAP_D6
    print(f"margin audit: d={d}, {len(tids)} types\n")
    screened, wild_n, secs = run(tids, d)
    print(f"\n{screened:,} labellings screened "
          f"({wild_n:,} wildcard-bearing) in {secs/60:.1f} min\n")
    print(f"{'comparison':<12} {'n':>10} {'min':>12} {'0.1 pct':>12} "
          f"{'median':>12} {'threshold':>11} {'min/thresh':>12}")
    worst = None
    for key, vals in REC.items():
        if not vals:
            print(f"{key:<12} {'(never evaluated)':>10}")
            continue
        v = np.array(vals)
        th = THRESH[key]
        ratio = (v.min() / th) if th > 0 else float("inf")
        print(f"{key:<12} {len(v):>10,} {v.min():>12.3e} "
              f"{np.percentile(v,0.1):>12.3e} {np.median(v):>12.3e} "
              f"{th:>11.0e} "
              f"{('inf' if th == 0 else f'{ratio:.1e}'):>12}")
        if th > 0 and (worst is None or ratio < worst[1]):
            worst = (key, ratio)
    print()
    if worst:
        print(f"Tightest margin: {worst[0]} clears its threshold by a factor "
              f"{worst[1]:.1e}.")
        print("An a-priori float64 error bound below that factor would certify "
              "every")
        print("decision recorded here without any exact re-computation.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
