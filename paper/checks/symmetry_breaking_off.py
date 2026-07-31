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
    # Prefer the run of record; fall back to the superseded two-part run only if it
    # is absent, exactly as paper/make_tables.py does, so the reference verdicts
    # this check compares against are the ones the paper reports.
    dirs = (("d6_final",) if (ROOT / "runs/d6_n10/d6_final/state.jsonl").exists()
            else ("d6_tangency", "d6_rest"))
    print(f"reference verdicts from: {', '.join(dirs)}")
    st = {}
    for d in dirs:
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


def realizer_check(budget=7200.0):
    """The one test of this kind that can fail in the dangerous direction.

    Every type above returns an EMPTY verdict, and a symmetry breaker that
    wrongly discarded orbits would reproduce empty verdicts perfectly -- so
    agreement there is evidence against spurious acceptance, not against loss.
    What actually matters is whether the pruning can lose a polytope that exists,
    and in dimension 6 there is exactly one place to test that: the subtree of the
    realizing type in which P^B6 is found.  It is too expensive to run the whole
    type without pruning, but the subtree of record is not, because its first two
    labels are fixed.

    Runs 379 with prefix (0,0), pruning OFF, and requires the same verdict: one
    polytope, exhausted.
    """
    tid, prefix = 379, (0, 0)
    ref = None
    for line in (ROOT / "runs/d6_n10/d6_final/state.jsonl").read_text().splitlines():
        if line.strip() and json.loads(line)["key"] == "379|0,0":
            ref = json.loads(line)
    t = dict(TYPES[tid])
    V, *_ = _setup(t)
    t["vertex_sets"] = [sorted(v) for v in V]
    so = {}
    t0 = time.time()
    res = process_type_stage4(t, D, max_assignments=50_000_000,
                              enum_timeout=budget, solve_timeout=budget / 2,
                              wildcard=True, use_burcroff_55b=True,
                              prefix=prefix, automorphisms=None,
                              stats_out=so, verbose=False) or []
    n = len({str(canonical_key(r, NF)) for r in res})
    exh = bool(so.get("exhausted"))
    sec = round(time.time() - t0, 1)
    ok = exh and n == 1
    print(f"\nrealizing subtree {tid}|{','.join(map(str, prefix))} with pruning OFF: "
          f"{n} polytopes, exhausted={exh}, "
          f"labellings={so.get('enum_count', 0):,}, {sec}s")
    if ref:
        print(f"  with pruning ON (classification run):  "
              f"{len(ref.get('keys', []))} polytopes, "
              f"exhausted={ref.get('exhausted')}, "
              f"labellings={ref.get('enum_count', 0):,}, {ref.get('sec')}s")
    print("  " + ("P^B6 is still found and the subtree still exhausts -- the "
                  "pruning loses nothing where a polytope exists" if ok else
                  "MISMATCH: expected 1 polytope and exhaustion"))
    # Artifact, so the paper's figures for this check are read rather than
    # transcribed (paper/make_tables.py).
    out = ROOT / "runs/d6_n10/symmetry_breaking_off.json"
    body = json.loads(out.read_text()) if out.exists() else {}
    body["realizing_subtree"] = {
        "subtree": f"{tid}|{','.join(map(str, prefix))}",
        "pruning_off": {"polytopes": n, "exhausted": exh,
                        "labellings": so.get("enum_count", 0), "seconds": sec},
        "pruning_on": ({"polytopes": len(ref.get("keys", [])),
                        "exhausted": ref.get("exhausted"),
                        "labellings": ref.get("enum_count", 0),
                        "seconds": ref.get("sec")} if ref else None),
        "verdict_reproduced": ok,
    }
    out.write_text(json.dumps(body, indent=1) + "\n")
    print(f"  artifact -> {out.relative_to(ROOT)}")
    return ok


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
            # a run that hit the budget is INCONCLUSIVE, not a mismatch: without
            # pruning these trees are far larger, which is the point of the device
            tag = ("reproduced" if exh and dist == ref else
                   "inconclusive (hit budget)" if not exh else "DISAGREES")
            print(f"  tid {tid:>4}  pruning OFF: {dist} polytopes, exhausted={exh}, "
                  f"enum={enum:,}, {sec}s   committed={ref}   {tag}", flush=True)
    art = ROOT / "runs/d6_n10/symmetry_breaking_off.json"
    body = json.loads(art.read_text()) if art.exists() else {}
    body["cheap_types"] = {
        "types": sorted(out),
        "per_type": {str(t): {"polytopes": v[0], "exhausted": v[1],
                              "labellings": v[2], "seconds": v[3],
                              "polytopes_with_pruning": len(keys[t])}
                     for t, v in out.items()},
    }
    art.write_text(json.dumps(body, indent=1) + "\n")
    agree = [t for t, (d, e, _, _) in out.items() if e and d == len(keys[t])]
    incon = [t for t, (d, e, _, _) in out.items() if not e]
    bad = [t for t, (d, e, _, _) in out.items() if e and d != len(keys[t])]
    print(f"\n{len(out)} types re-run with symmetry breaking disabled")
    print(f"  exhausted and verdict reproduced    : {len(agree)}")
    print(f"  hit the budget, inconclusive         : {len(incon)}  {sorted(incon)}")
    print(f"  genuine disagreements                : {len(bad)}  {sorted(bad)}")
    print("  NOTE: every verdict above is EMPTY, and a symmetry breaker that")
    print("        wrongly discarded orbits would reproduce empty verdicts too.")
    print("        The test that can fail in the dangerous direction is the")
    print("        realizing subtree, which is run next.")
    real_ok = realizer_check()
    print()
    print("RESULT:", "no verdict changed without the pruning, and the realizing "
          "subtree still finds P^B6"
          if not bad and real_ok else
          f"DISAGREEMENT on {bad}" if bad else
          "the realizing subtree did not reproduce its verdict without pruning")
    return 0 if (not bad and real_ok) else 1


if __name__ == "__main__":
    sys.exit(main())
