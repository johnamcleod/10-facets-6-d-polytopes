#!/usr/bin/env python3
"""Rigorous full-label resolve of the 54 d=6 facet-filter survivors.

Replaces run_survivors_fulllabel.py, whose "empty(COMPLETE)" tag was a wall-clock
heuristic: 19 types hit the 3600s enum timeout and were mislabeled complete —
including tid 379 = P^B6 itself (distinct=0 there is a BUG per CLAUDE.md §0).

This driver uses the real exhaustion flag (process_type_stage4 stats_out) and
partitions the label search by enumeration prefix (enumerate_labels_backtrack
`prefix`), refining any timed-out subtree one level deeper.  A type's verdict is
rigorous iff EVERY subtree reports exhausted=True; the type's polytopes are the
union of subtree results deduped by canonical_key.

Resumable: completed subtrees are persisted to STATE and skipped on restart.

Usage:  python3 run_survivors_rigorous.py [fast|suspect|all|<tid,tid,...>] [nproc]
"""
import json, sys, time
from pathlib import Path
import multiprocessing as mp
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline.stage4_gram import process_type_stage4, VALID_LABELS
from pipeline.stage4_blockpaste import _setup
from run_d4 import canonical_key

SURV = [8, 12, 17, 34, 36, 38, 40, 51, 55, 59, 60, 61, 69, 70, 92, 103, 120, 127, 132,
        140, 154, 159, 162, 168, 173, 206, 214, 218, 220, 229, 234, 239, 255, 265, 273,
        284, 286, 287, 295, 297, 308, 315, 317, 320, 329, 332, 344, 352, 354, 356, 360,
        378, 379, 382]
# The 19 that hit the 3600s enum timeout in the fulllabel run (0 there = not rigorous).
SUSPECT = [34, 40, 55, 60, 61, 103, 120, 127, 132, 159, 173, 206, 284, 286, 295, 315,
           329, 344, 379]
FAST = [t for t in SURV if t not in SUSPECT]

TYPES = {t['type_id']: t for t in json.load(open('runs/d6_n10/stage2/types.json'))}
OUTDIR = Path("runs/d6_n10/survivors_rigorous")
STATE = OUTDIR / "state.json"
REALIZER_DIR = OUTDIR / "realizers"

N_LAB = len(VALID_LABELS)            # 10 — full labels, label_indices=None
MAX_DEPTH = 6
MAXA = 50_000_000
SOLVE_TO = 600.0

def enum_to(depth):
    return 3600.0 if depth == 0 else 1800.0


def _key(tid, prefix):
    return f"{tid}|{','.join(map(str, prefix))}"


def work(task):
    tid, prefix = task
    try:
        t = dict(TYPES[tid]); V, *_ = _setup(t); t['vertex_sets'] = [sorted(v) for v in V]
        so = {}
        t0 = time.time()
        res = process_type_stage4(t, 6, max_assignments=MAXA,
                                  enum_timeout=enum_to(len(prefix)),
                                  solve_timeout=SOLVE_TO,
                                  label_indices=None, prefix=tuple(prefix) or None,
                                  stats_out=so, verbose=False) or []
        el = round(time.time() - t0, 1)
        keys = sorted(set(str(canonical_key(r, 10)) for r in res))
        return (tid, prefix, keys, res, bool(so.get("exhausted")),
                so.get("enum_count", 0), el, None)
    except Exception as e:
        return (tid, prefix, [], [], False, 0, 0.0, f"{type(e).__name__}: {e}")


def main():
    sel = sys.argv[1] if len(sys.argv) > 1 else "all"
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 else max(1, mp.cpu_count() - 1)
    if sel == "fast":
        tids = FAST
    elif sel == "suspect":
        tids = SUSPECT
    elif sel == "all":
        tids = SUSPECT + FAST          # suspects first: they gate the result
    else:
        tids = [int(x) for x in sel.split(",")]

    OUTDIR.mkdir(parents=True, exist_ok=True)
    REALIZER_DIR.mkdir(parents=True, exist_ok=True)
    state = json.loads(STATE.read_text()) if STATE.exists() else {}

    # 379 (P^B6 anchor) always first if present.
    tids = sorted(set(tids), key=lambda t: (t != 379, t))
    queue = [(tid, ()) for tid in tids]

    wave = 0
    while queue:
        todo = [(tid, pfx) for tid, pfx in queue
                if _key(tid, pfx) not in state]
        skipped = len(queue) - len(todo)
        print(f"--- wave {wave}: {len(todo)} subtree tasks ({skipped} cached), "
              f"nproc={nproc} ---", flush=True)
        next_queue = []
        if todo:
            with mp.get_context("fork").Pool(processes=nproc) as pool:
                for tid, pfx, keys, res, exh, ecount, el, err in \
                        pool.imap_unordered(work, todo):
                    k = _key(tid, pfx)
                    if err:
                        print(f"tid {tid} pfx={pfx}: ERROR {err}", flush=True)
                        state[k] = {"error": err, "exhausted": False, "sec": el}
                    else:
                        tag = ("EXHAUSTED" if exh else
                               f"TIMEOUT->refine d{len(pfx)+1}")
                        if keys:
                            tag = "**REALIZES** " + tag
                            fn = REALIZER_DIR / f"tid{tid}_{'_'.join(map(str,pfx)) or 'root'}.json"
                            fn.write_text(json.dumps(res, indent=1, default=str))
                        print(f"tid {tid} pfx={pfx}: distinct={len(keys)} "
                              f"enum={ecount} [{tag}] in {el}s", flush=True)
                        state[k] = {"keys": keys, "exhausted": exh,
                                    "enum_count": ecount, "sec": el}
                    STATE.write_text(json.dumps(state, indent=1))
        # refinement pass over the whole queue (incl. cached entries)
        for tid, pfx in queue:
            st = state.get(_key(tid, pfx), {})
            if st.get("exhausted") or st.get("error"):
                continue
            if len(pfx) >= MAX_DEPTH:
                print(f"tid {tid} pfx={pfx}: UNRESOLVED at max depth", flush=True)
                continue
            next_queue.extend((tid, tuple(pfx) + (li,)) for li in range(N_LAB))
        queue = next_queue
        wave += 1

    # final per-type verdicts.  A subtree is covered if it was exhausted itself
    # or (having timed out) all 10 of its children are covered — the refinement
    # tree partitions the search exactly, so root coverage = rigorous verdict.
    def covered(tid, pfx):
        st = state.get(_key(tid, pfx))
        if st is None or st.get("error"):
            return False
        if st.get("exhausted"):
            return True
        if len(pfx) >= MAX_DEPTH:
            return False
        return all(covered(tid, tuple(pfx) + (li,)) for li in range(N_LAB))

    verdict = {}
    for tid in tids:
        ks = [s for key, s in state.items() if key.startswith(f"{tid}|")]
        keys = sorted(set(k for s in ks for k in s.get("keys", [])))
        rigorous = bool(ks) and covered(tid, ())
        verdict[tid] = {"distinct": len(keys), "rigorous": rigorous,
                        "subtrees": len(ks)}
        print(f"VERDICT tid {tid}: distinct={len(keys)} "
              f"{'RIGOROUS' if rigorous else 'NOT-RIGOROUS'}", flush=True)
    (OUTDIR / "verdicts.json").write_text(json.dumps(verdict, indent=1))
    realizing = sorted(t for t, v in verdict.items() if v["distinct"] > 0)
    nonrig = sorted(t for t, v in verdict.items() if not v["rigorous"])
    print(f"\nRIGOROUS-DONE. realizing={realizing} not-rigorous={nonrig}", flush=True)


if __name__ == "__main__":
    main()
