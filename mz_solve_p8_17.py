#!/usr/bin/env python3
"""Solve Ma-Zheng's OWN candidate list for P8_17 (= our tid 8 = Burcroff's G11).

This is the decisive experiment for the one open discrepancy in our d=4 census.
Our run closes tid 8 rigorously at 7 polytopes; the published count is 8.  The two
halves of our pipeline that could lose a polytope are the ENUMERATOR (plus its
forward-checks, symmetry breaking and structured screen) and the SOLVER (feasibility
scan, Gauss-Newton refinement, exact certification).

Ma-Zheng's repository ships their intermediate candidate list for exactly this
combinatorial type: 325,957 admissible label vectors, produced by their independent
implementation of a different method.  Feeding that list to OUR solver, with our
enumerator and screen entirely bypassed, separates the two halves:

  * 8 polytopes  =>  their enumeration reaches a labelling ours does not, and the
                     defect is in our enumerator / screen / symmetry breaking.
  * 7 polytopes  =>  two independent enumerations and our solver agree on 7, and the
                     discrepancy lies in the published count or in a later stage of
                     their own pipeline -- not in a candidate we pruned.

Bypassing the screen is the whole point: `_solve_wild_assignment` is called on every
row, so no row is discarded on any heuristic ground before it is actually solved.

Usage:  python3 mz_solve_p8_17.py [nproc] [limit]
"""
from __future__ import annotations

import itertools
import json
import multiprocessing as mp
import re
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np                                              # noqa: E402
import pipeline.stage4_gram as _sg                              # noqa: E402
from pipeline.stage4_gram import _solve_wild_assignment    # noqa: E402
from run_d4 import canonical_key                                # noqa: E402
import sympy as sp                                              # noqa: E402

ROOT = Path(__file__).resolve().parent
MZ = ROOT / "scratchpad/HCPdm"
N, D, OUR_TID, MZ_LINE = 8, 4, 8, 17
KEY_ALLIN = [10 * a + b for a, b in itertools.combinations(range(1, N + 1), 2)]


def mz_vertex_sets(line_no):
    lines = [l for l in (MZ / "polytopeDATA/4d8m.txt").read_text().splitlines()
             if l.strip()]
    return [tuple(sorted(int(x) - 1 for x in g.split(",")))
            for g in re.findall(r"\[([\d,]+)\]", lines[line_no - 1])]


def dotted_of(vsets, n=N):
    on = [{k for k, v in enumerate(vsets) if i in v} for i in range(n)]
    return [(i, j) for i, j in itertools.combinations(range(n), 2)
            if not (on[i] & on[j])]


def build():
    """Return (sigma_inv, cols, our_dotted, rows)."""
    types = {t["type_id"]: t for t in
             json.load(open(ROOT / "runs/d4_n8/stage2/types.json"))}
    ours_v = [tuple(sorted(v)) for v in types[OUR_TID]["vertex_sets"]]
    theirs_v = mz_vertex_sets(MZ_LINE)
    tgt = {frozenset(v) for v in theirs_v}
    sigma = next(p for p in itertools.permutations(range(N))
                 if {frozenset(p[x] for x in v) for v in ours_v} == tgt)
    sinv = [0] * N
    for i, x in enumerate(sigma):
        sinv[x] = i
    our_dot = dotted_of(ours_v)
    their_dot = dotted_of(theirs_v)
    # the alignment must carry our dotted pairs onto theirs
    assert {tuple(sorted((sigma[i], sigma[j]))) for i, j in our_dot} == set(their_dot)
    their_inf = {10 * (i + 1) + (j + 1) for i, j in their_dot}
    cols = [k for k in KEY_ALLIN if k not in their_inf]
    f = MZ / f"output/P8_{MZ_LINE}/P8_{MZ_LINE}_LSIEr1_per.txt"
    rows = [tuple(int(x) for x in l.split())
            for l in f.read_text().splitlines() if len(l.split()) == len(cols)]
    return sinv, cols, our_dot, rows


SINV, COLS, OUR_DOT, ROWS = build()
# their column -> our (i,j) pair, precomputed
COLPAIR = [tuple(sorted((SINV[c // 10 - 1], SINV[c % 10 - 1]))) for c in COLS]
SYMS = {p: sp.Symbol(f"x_{p[0]}_{p[1]}", positive=True) for p in OUR_DOT}
SYM_LIST = [SYMS[p] for p in OUR_DOT]


def work(chunk):
    found, nsolved, nwild, err = {}, 0, 0, 0
    for ri in chunk:
        row = ROWS[ri]
        lab = {COLPAIR[c]: int(v) for c, v in enumerate(row)}
        wild = [p for p, m in lab.items() if m >= 7]
        if wild:
            nwild += 1
        try:
            flags = {}
            res = _solve_wild_assignment(
                lab, wild, OUR_DOT, SYM_LIST, N, D,
                m_scan_max=100, dps=100, deadline_s=600.0, flags=flags)
        except Exception:
            err += 1
            continue
        for labels_int, sol in res:
            rec = {"type_id": OUR_TID,
                   "label_assignment": {str(k): v for k, v in labels_int.items()},
                   "dot_values": sol}
            try:
                k = str(canonical_key(rec, N))
            except Exception:
                k = json.dumps(rec, sort_keys=True, default=str)
            found[k] = rec
            nsolved += 1
    return found, nsolved, nwild, err


def main():
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else max(1, mp.cpu_count() - 1)
    limit = int(sys.argv[2]) if len(sys.argv) > 2 else len(ROWS)
    idx = list(range(min(limit, len(ROWS))))
    print(f"Ma-Zheng P8_{MZ_LINE} candidate rows: {len(ROWS):,}  "
          f"(solving {len(idx):,} with nproc={nproc})")
    print(f"our dotted pairs: {OUR_DOT}")
    print(f"rows carrying a wildcard (m>=7): "
          f"{sum(1 for r in ROWS if max(r) >= 7):,}")
    print("NOTE: the structured screen is BYPASSED -- every row goes to the solver.\n")
    CH = 200
    chunks = [idx[i:i + CH] for i in range(0, len(idx), CH)]
    allf, tot, totw, tote, done = {}, 0, 0, 0, 0
    t0 = time.time()
    with mp.get_context("fork").Pool(processes=nproc) as pool:
        for f, ns, nw, er in pool.imap_unordered(work, chunks):
            allf.update(f); tot += ns; totw += nw; tote += er
            done += 1
            if done % 25 == 0 or done == len(chunks):
                el = time.time() - t0
                rate = done / el
                print(f"  {done}/{len(chunks)} chunks  {el/60:.1f} min  "
                      f"eta {(len(chunks)-done)/rate/60:.1f} min  "
                      f"distinct={len(allf)}  raw={tot}  errors={tote}", flush=True)
    out = ROOT / "runs/d4_n8/mz_p8_17_solved.json"
    out.write_text(json.dumps(list(allf.values()), indent=1, default=str))
    print(f"\nDISTINCT POLYTOPES from Ma-Zheng's own candidate list: {len(allf)}")
    print(f"  raw certified solutions {tot}, wildcard rows {totw}, errors {tote}")
    print(f"  written {out}")
    print("\nOur run's independent verdict for this type was 7; the published count "
          "is 8.")
    if len(allf) == 8:
        print("=> 8 here means OUR ENUMERATOR pruned the eighth polytope: their "
              "candidate\n   list contains a labelling our search never reached.")
    elif len(allf) == 7:
        print("=> 7 here means our solver and TWO independent enumerations agree. "
              "The\n   eighth polytope is not in a candidate we pruned; the "
              "discrepancy lies\n   in the published count or downstream in their "
              "own pipeline.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
