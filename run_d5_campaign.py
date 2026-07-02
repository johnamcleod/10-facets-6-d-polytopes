#!/usr/bin/env python3
"""Resumable nights/weekends campaign to finish the 11 intractable d=5 types.

The block-paste run left 11 d=5 types INCOMPLETE (aborted by a per-type budget); their
candidate sets are finite but huge.  This harness partitions each into small independent
work-units and grinds them over many short, resumable sessions (run when the machine is
idle).  A filesystem queue makes it crash/sleep-safe:

    runs/d5_campaign/todo/<uid>.json   pending unit  {type_id, forced}
    runs/d5_campaign/wip/<uid>.json    claimed (atomic rename); requeued on restart
    runs/d5_campaign/done/<uid>.json   result {type_id, forced, polytopes, n_la, ...}

A unit = one block-paste partition (seed columns pinned via ``forced``).  If a unit's paste
overflows the candidate cap it AUTO-SPLITS: it pins one more seed column and enqueues the 6
sub-units, so arbitrarily hard partitions are handled without losing work.  Each unit is
atomic and idempotent.

Usage:
    python3 run_d5_campaign.py init [p]                # build the todo queue (depth p, default 3)
    python3 run_d5_campaign.py run  [minutes] [nproc]  # grind pending units until the deadline
    python3 run_d5_campaign.py status                  # per-type + total progress vs census 51
"""
import sys, json, os, time, itertools, signal
from pathlib import Path
import numpy as np
import multiprocessing as mp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pipeline.stage4_blockpaste as bp
from pipeline.stage4_blockpaste import (paste_candidates, expand_label_assignments,
                                        _build_constraints, _order_vertices, infer_dims)
from pipeline.stage4_gram import (screen_candidate, _build_minor_index, _refine_mpmath,
                                  _recognize_minpoly_and_verify)
from run_d4 import canonical_key
from sympy import symbols

bp.USE_L4_BASIS = True                       # d=5 uses the prism-base convention
D5 = "runs/d5_n9"
CAMP = Path(f"{D5}/d5_campaign")
CAP = 3_000_000                              # per-partition candidate cap -> auto-split above
INCOMPLETE = [1, 4, 8, 15, 16, 20, 30, 31, 53, 55, 60]

_TYPES = None
def types():
    global _TYPES
    if _TYPES is None:
        _TYPES = {i: t for i, t in enumerate(json.load(open(f"{D5}/stage2/types.json")))}
    return _TYPES


def _seed_cols(t):
    _, _, ctx = _build_constraints(t)
    order = _order_vertices(ctx["V"], ctx["cols_of"])
    return ctx["cols_of"](ctx["V"][order[0]])          # 10 ordinary cols of the seed vertex


def _uid(tid, forced):
    return f"t{tid}_" + "_".join(f"{c}-{v}" for c, v in sorted(forced.items())) or f"t{tid}_"


# ---------------------------------------------------------------------------
# process one unit (one partition): paste -> expand -> screen -> solve
# ---------------------------------------------------------------------------
def _solve_start(la, dotted, x0, n, d):
    if not dotted:
        return [{"label_assignment": {str(p): v for p, v in la.items()}, "dot_values": {}}]
    sl = [symbols(f"x_{p[0]}_{p[1]}", positive=True) for p in dotted]
    try:
        xh = _refine_mpmath(np.array(x0), la, dotted, n, d, dps=100)
        if xh is None:
            return []
        out = []
        for so in _recognize_minpoly_and_verify(xh, sl, la, dotted, n, d, dps=100,
                                                recover_minpoly=False):
            out.append({"label_assignment": {str(p): v for p, v in la.items()},
                        "dot_values": so})
        return out
    except Exception:
        return []


def process_unit(tid, forced):
    """Return a result dict, OR {'split': col} if the paste overflowed and the unit should
    be split by pinning one more seed column."""
    t = types()[tid]
    n, d = infer_dims(t)
    dotted = [tuple(sorted(m)) for m in t["missing_faces"] if len(m) == 2]
    _, _, ctx = _build_constraints(t)
    pool = list(_seed_cols(t)) + [c for c in range(ctx["Nr"]) if c not in _seed_cols(t)]
    try:
        cands, ordinary = paste_candidates(t, forced=forced, max_candidates=CAP)
    except OverflowError:
        extra = next((c for c in pool if c not in forced), None)   # pin one more column
        return {"split": extra}
    mi = _build_minor_index(dotted, n, d)
    starts, n_la = [], 0
    for la in expand_label_assignments(cands, ordinary):
        n_la += 1
        starts.extend(screen_candidate(la, dotted, mi, n, d))
    raw = []
    for la, dot, x0 in starts:
        raw.extend(_solve_start(la, dot, x0, n, d))
    return {"type_id": tid, "forced": {str(k): v for k, v in forced.items()},
            "polytopes": raw, "n_la": n_la, "screen_pass": len(starts)}


# ---------------------------------------------------------------------------
# filesystem queue
# ---------------------------------------------------------------------------
def _dirs():
    for s in ("todo", "wip", "done"):
        (CAMP / s).mkdir(parents=True, exist_ok=True)
    return CAMP / "todo", CAMP / "wip", CAMP / "done"


def cmd_init(p=3):
    todo, wip, done = _dirs()
    LABELS = (2, 3, 4, 5, 6, 7)
    n_units = 0
    for tid in INCOMPLETE:
        sc = _seed_cols(types()[tid])[:p]
        for combo in itertools.product(LABELS, repeat=p):
            forced = {sc[i]: combo[i] for i in range(p)}
            uid = _uid(tid, forced)
            if (done / f"{uid}.json").exists() or (todo / f"{uid}.json").exists():
                continue
            (todo / f"{uid}.json").write_text(json.dumps({"type_id": tid, "forced": forced}))
            n_units += 1
    print(f"init: {n_units} new units queued (depth p={p}) for types {INCOMPLETE}", flush=True)
    print(f"  todo={len(list(todo.glob('*.json')))} done={len(list(done.glob('*.json')))}",
          flush=True)


def _claim(todo, wip):
    for f in todo.iterdir():
        if f.suffix != ".json":
            continue
        dest = wip / f.name
        try:
            os.rename(f, dest)              # atomic; only one worker wins
            return dest
        except (FileNotFoundError, OSError):
            continue
    return None


def _worker(deadline):
    todo, wip, done = CAMP / "todo", CAMP / "wip", CAMP / "done"
    while time.time() < deadline:
        claimed = _claim(todo, wip)
        if claimed is None:
            return
        u = json.loads(claimed.read_text())
        tid = u["type_id"]; forced = {int(k): v for k, v in u["forced"].items()}
        res = process_unit(tid, forced)
        if "split" in res and res["split"] is not None:
            # enqueue 6 sub-units pinning one more column
            for v in (2, 3, 4, 5, 6, 7):
                f2 = dict(forced); f2[res["split"]] = v
                uid2 = _uid(tid, f2)
                if not (done / f"{uid2}.json").exists():
                    (todo / f"{uid2}.json").write_text(
                        json.dumps({"type_id": tid, "forced": f2}))
            (done / claimed.name).write_text(json.dumps(
                {"type_id": tid, "forced": u["forced"], "split_into": 6, "polytopes": []}))
        else:
            (done / claimed.name).write_text(json.dumps(res))
        claimed.unlink(missing_ok=True)


def cmd_run(minutes=480, nproc=None):
    todo, wip, done = _dirs()
    for f in wip.glob("*.json"):              # requeue anything left mid-flight last time
        os.rename(f, todo / f.name)
    nproc = nproc or max(1, mp.cpu_count() - 1)
    deadline = time.time() + minutes * 60
    n0 = len(list(done.glob("*.json")))
    print(f"run: {len(list(todo.glob('*.json')))} pending, {nproc} workers, "
          f"deadline {minutes}min", flush=True)
    procs = [mp.get_context("fork").Process(target=_worker, args=(deadline,))
             for _ in range(nproc)]
    for pr in procs: pr.start()
    for pr in procs: pr.join()
    print(f"run done: {len(list(done.glob('*.json'))) - n0} units completed this session; "
          f"{len(list(todo.glob('*.json')))} still pending", flush=True)


def cmd_status():
    done = CAMP / "done"
    # combine campaign results with the already-COMPLETE d=5 types
    raw_by_type = {}
    for f in (done.glob("*.json") if done.exists() else []):
        r = json.loads(f.read_text())
        raw_by_type.setdefault(r["type_id"], []).extend(r.get("polytopes", []))
    # per-type pending
    todo = CAMP / "todo"
    pend = {}
    for f in (todo.glob("*.json") if todo.exists() else []):
        pend[json.loads(f.read_text())["type_id"]] = pend.get(
            json.loads(f.read_text())["type_id"], 0) + 1
    n = 9
    total = 0
    print("=== d=5 campaign status ===")
    # already-complete (non-incomplete) types from the original run
    import glob
    base = 0
    for tf in glob.glob(f"{D5}/blockpaste/type_*.json"):
        tid = int(tf.split("type_")[1].split(".")[0])
        if tid > 108 or tid in INCOMPLETE:
            continue
        base += json.load(open(tf))["distinct"]
    print(f"already-complete types total: {base}")
    for tid in INCOMPLETE:
        raw = raw_by_type.get(tid, [])
        dist = len(set(canonical_key(r, n) for r in raw))
        p = pend.get(tid, 0)
        total += dist
        print(f"  type {tid}: {dist} distinct so far, {p} units pending"
              + ("  [DONE]" if p == 0 else ""))
    print(f"\nCAMPAIGN total (incomplete types): {total}  + base {base} = {base + total}  "
          f"(census 51)")
    tp = len(list(todo.glob('*.json'))) if todo.exists() else 0
    dn = len(list(done.glob('*.json'))) if done.exists() else 0
    print(f"units: {dn} done, {tp} pending")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "status"
    if cmd == "init":
        cmd_init(int(sys.argv[2]) if len(sys.argv) > 2 else 3)
    elif cmd == "run":
        cmd_run(float(sys.argv[2]) if len(sys.argv) > 2 else 480,
                int(sys.argv[3]) if len(sys.argv) > 3 else None)
    elif cmd == "status":
        cmd_status()
    else:
        print(__doc__)
