#!/usr/bin/env python3
"""Re-run types under random relabellings of the facets, with symmetry breaking ON.

The lex-leader symmetry breaking is argued sound, and on 280 of the 304 searched
types its verdicts are reproduced with the device switched off.  The other 23
types (and the realizing type 379, whose unpruned run is expensive) are too large
to enumerate unpruned.  This is a cheap consistency check
(second-round report, item 3.2): relabel the facets by a random permutation, which
changes the enumeration order, the prefix partition and the automorphism group
the pruner works with (it is conjugated), and require the same verdict.  A buggy
pruner that discarded orbits would be unlikely to discard the same ones under
every relabelling.

It is evidence, not a proof; the verdict is required to be identical (same number
of polytopes, exhausted) for every relabelling that finishes within its budget.

Run:  python3 checks/symmetry_relabel.py [seeds=3] [budget=7200] [nproc=7]
      [tids=5,16,...]    (default: the 23 budget-exceeding types, plus 379)
Artifact: runs/d6_n10/symmetry_relabel.json
"""
from __future__ import annotations

import json
import multiprocessing as mp
import random
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline.stage4_gram import process_type_stage4          # noqa: E402
from pipeline.stage4_blockpaste import _setup                 # noqa: E402
from pipeline.utils.automorphisms import compute_aut_group     # noqa: E402
from run_d4 import canonical_key                              # noqa: E402

D, NF = 6, 10
ARG = {a.split("=", 1)[0]: a.split("=", 1)[1] for a in sys.argv[1:] if "=" in a}
SEEDS = int(ARG.get("seeds", 3))
BUDGET = float(ARG.get("budget", 7200))
NPROC = int(ARG.get("nproc", 7))
DEFAULT = [5, 16, 42, 57, 66, 68, 72, 81, 102, 170, 173, 184, 187, 189, 198, 208,
           228, 253, 263, 268, 274, 358, 363, 379]
TIDS = [int(x) for x in ARG["tids"].split(",")] if "tids" in ARG else DEFAULT
TYPES = {t["type_id"]: t for t in json.load(open(ROOT / "runs/d6_n10/stage2/types.json"))}
REF = {}
for line in (ROOT / "runs/d6_n10/d6_valid/state.jsonl").read_text().splitlines():
    if line.strip():
        r = json.loads(line)
        REF[r["key"]] = r


def relabelled(tid, seed):
    """The type with facet i renamed perm[i]; vertex sets carried along exactly."""
    t = dict(TYPES[tid])
    V, *_ = _setup(t)
    perm = list(range(NF))
    random.Random(1000 * tid + seed).shuffle(perm)
    u = {"type_id": tid,
         "missing_faces": [sorted(perm[i] for i in m) for m in t["missing_faces"]],
         "vertex_sets": [sorted(perm[i] for i in v) for v in V],
         "p_count": t["p_count"]}
    return u, perm


def work(task):
    tid, seed = task
    u, perm = relabelled(tid, seed)
    V, *_ = _setup(u)
    aut = compute_aut_group(V, NF)
    so = {}
    t0 = time.time()
    res = process_type_stage4(u, D, max_assignments=50_000_000,
                              enum_timeout=BUDGET, solve_timeout=600.0,
                              wildcard=True, use_burcroff_55b=True,
                              automorphisms=aut, stats_out=so, verbose=False) or []
    n = len({str(canonical_key(r, NF)) for r in res})
    exh = bool(so.get("exhausted")) and not so.get("wild_unbounded")
    return (tid, seed, perm, n, exh, so.get("enum_count", 0),
            round(time.time() - t0, 1))


def main():
    tasks = [(t, s) for t in TIDS for s in range(1, SEEDS + 1)]
    out = {}
    with mp.get_context("fork").Pool(processes=NPROC) as pool:
        for tid, seed, perm, n, exh, enum, sec in pool.imap_unordered(work, tasks):
            ref = REF[f"{tid}|"]
            rn = len(ref.get("keys", []))
            tag = ("reproduced" if exh and n == rn else
                   "inconclusive (budget)" if not exh else "DISAGREES")
            out.setdefault(str(tid), []).append(
                {"seed": seed, "perm": perm, "polytopes": n, "exhausted": exh,
                 "labellings": enum, "seconds": sec, "reference_polytopes": rn,
                 "reference_labellings": ref.get("enum_count", 0), "verdict": tag})
            print(f"  tid {tid:>4} seed {seed}: {n} polytopes, exhausted={exh}, "
                  f"enum={enum:,}, {sec}s   reference={rn}   {tag}", flush=True)
    runs = [r for v in out.values() for r in v]
    bad = [r for r in runs if r["verdict"] == "DISAGREES"]
    inc = [r for r in runs if r["verdict"].startswith("inconclusive")]
    art = ROOT / "runs/d6_n10/symmetry_relabel.json"
    art.write_text(json.dumps({"seeds": SEEDS, "budget": BUDGET, "types": TIDS,
                               "runs": out, "disagreements": len(bad),
                               "inconclusive": len(inc),
                               "reproduced": len(runs) - len(bad) - len(inc)},
                              indent=1) + "\n")
    print(f"\n{len(runs)} relabelled runs over {len(TIDS)} types: "
          f"{len(runs) - len(bad) - len(inc)} reproduced, {len(inc)} over budget, "
          f"{len(bad)} disagreements")
    print(f"artifact -> {art.relative_to(ROOT)}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
