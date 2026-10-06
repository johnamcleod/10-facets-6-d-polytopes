#!/usr/bin/env python3
"""Re-run selected subtrees of the classification with the wildcard-instance dump
enabled, producing the artifact the exact certifier consumes.

The numerical joint-feasibility screen of Section 5.4 rejects every
wildcard-bearing labelling, but the run of record recorded only the COUNT of
those rejections (state.jsonl `diag.wild.joint_infeasible`).  To replace that
screen by exact infeasibility certificates we need the instances themselves, so
this driver re-runs the subtrees that produce them with
`pipeline.stage4_gram.WILD_DUMP_PATH` set, writing one JSON record per
wildcard-bearing labelling: its ordinary labels, dashed pairs, wild pairs, and
the numerical residual the screen used to reject it.

The subtree parameters (timeouts, alphabet, symmetry pruning, Burcroff 5.5b) are
taken from run_survivors_rigorous, not restated, so the dumped instances are the
ones the run of record actually screened.

Usage:
    python3 run_wild_dump.py                       # d=6: the subtrees of record
    python3 run_wild_dump.py d=5                   # the whole d=5 census
    python3 run_wild_dump.py d=6 379|0,0           # one explicit subtree
    python3 run_wild_dump.py d=5 out=d5_wild.jsonl

Output (default):  runs/d{D}_n{D+4}/wild_instances.jsonl
"""
import json
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

D = next((int(a.split("=")[1]) for a in sys.argv[1:] if a.startswith("d=")), 6)
NF = D + 4

# run_survivors_rigorous reads sys.argv at import time (dimension, wildcard flag,
# output dir), so present it the command line it expects.  "wildcard" is not
# optional here: the wildcard alphabet is what creates wildcard instances at all.
_argv = sys.argv
sys.argv = [_argv[0], "all", "1", "wildcard", f"d={D}"]
import run_survivors_rigorous as R          # noqa: E402
import pipeline.stage4_gram as _sg          # noqa: E402
sys.argv = _argv

# The subtrees that produced complete labellings in the run of record.  For d=6
# the run of record (d6_exact, exact forward-checking gates) exhausts every type
# at the root, and only type 379's root reaches a complete labelling at all: all
# 952 labellings, 406 of them wildcard-bearing, live in it.  The dump of record
# was made from that root ("379|" in its .meta.json), so the root is the default
# here and a re-run reproduces the committed artifact exactly.  (In the earlier
# float-gate run, d6_final, the same 952 labellings all lay inside its subtree
# 379|0,0, the other 70 of its 71 subtrees being emptied by forward checking.)
DEFAULT_SUBTREES = {
    6: [(379, ())],
}


def subtrees_from_state(d):
    """Every subtree of the run of record with a nonzero labelling count."""
    for name in ("d6_exact", "d6_final", "survivors_wildcard", "survivors_rigorous"):
        p = ROOT / f"runs/d{d}_n{d + 4}/{name}/state.jsonl"
        if not p.exists():
            continue
        out = []
        for line in p.read_text().splitlines():
            if not line.strip():
                continue
            rec = json.loads(line)
            if (rec.get("diag") or {}).get("wild_assignments"):
                tid, _, pfx = rec["key"].partition("|")
                prefix = tuple(int(x) for x in pfx.split(",") if x != "")
                out.append((int(tid), prefix))
        if out:
            print(f"subtree list from runs/d{d}_n{d + 4}/{name}/state.jsonl")
            return out
    return []


def _reset(path):
    """Clear a sink AND its completion marker before writing.

    The marker must go first: leaving a previous run's marker beside a
    freshly-truncated dump would advertise a complete artifact that is being
    rewritten, and a consumer checking only for the marker's presence would read a
    prefix.  (Consumers also compare the record count, which is what actually
    catches this, but there is no reason to leave the trap lying around.)
    """
    meta = path.with_suffix(path.suffix + ".meta.json")
    if meta.exists():
        meta.unlink()
    if path.exists():
        path.unlink()


def parse_subtree(arg):
    tid, _, pfx = arg.partition("|")
    return int(tid), tuple(int(x) for x in pfx.split(",") if x != "")


def main():
    args = [a for a in sys.argv[1:] if not a.startswith(("d=", "out="))]
    out = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("out=")),
               f"runs/d{D}_n{NF}/wild_instances.jsonl")
    out_path = ROOT / out if not os.path.isabs(out) else Path(out)

    if args:
        subtrees = [parse_subtree(a) for a in args]
    else:
        subtrees = subtrees_from_state(D) or DEFAULT_SUBTREES.get(D)
        if not subtrees:
            # No dumped state to learn from: run every required type from the
            # root, which is what the census driver itself does.
            tids = sorted(set(R.SURV) - set(R.ESSELMANN_3FREE_KILLED)
                          - set(R.LEMMA_KILLED))
            subtrees = [(t, ()) for t in tids]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    _reset(out_path)
    _sg.WILD_DUMP_PATH = str(out_path)

    # The labellings without a wildcard go to a companion sink in the same run:
    # they are decided by a different code path (the double-precision cascade),
    # they are the other half of the screen caveat, and re-running the subtree
    # twice to collect them separately would be wasteful and would risk the two
    # dumps describing different runs.
    plain_path = out_path.with_name(out_path.name.replace("wild_", "plain_")
                                    if "wild_" in out_path.name
                                    else "plain_" + out_path.name)
    _reset(plain_path)
    _sg.PLAIN_DUMP_PATH = str(plain_path)

    print(f"d={D}, {len(subtrees)} subtree(s) -> {out_path}")
    t0 = time.time()
    total_wild = 0
    for tid, prefix in subtrees:
        _sg.WILD_DUMP_TAG = R._key(tid, prefix)
        t1 = time.time()
        res = R.work((tid, prefix))
        diag = res[8] or {}
        n_wild = diag.get("wild_assignments", 0)
        total_wild += n_wild
        print(f"  {R._key(tid, prefix):<16} labellings={res[5]:<8} "
              f"wild={n_wild:<6} exhausted={res[4]}  realizers={len(res[2])}  "
              f"{time.time() - t1:.0f}s"
              + (f"  ERROR {res[7]}" if res[7] else ""), flush=True)

    n_lines = sum(1 for _ in open(out_path)) if out_path.exists() else 0
    print(f"\n{total_wild} wildcard-bearing labellings reported by the counters, "
          f"{n_lines} records dumped, {time.time() - t0:.0f}s total")
    if n_lines != total_wild:
        print("MISMATCH: dump and counters disagree", file=sys.stderr)
        return 1
    # Completion marker.  The dump is appended to as the run proceeds, so a
    # consumer that merely finds the file cannot tell a finished dump from a
    # running one -- and certifying a partial dump would report "every instance
    # certified" about a prefix.  The certifier requires this marker and checks
    # the record count against it.
    n_plain = sum(1 for _ in open(plain_path)) if plain_path.exists() else 0
    print(f"{n_plain} non-wildcard labellings dumped to {plain_path}")
    info = {
        "complete": True,
        "d": D,
        "records": n_lines,
        "wild_assignments_from_counters": total_wild,
        "plain_records": n_plain,
        "subtrees": [f"{t}|{','.join(map(str, p))}" for t, p in subtrees],
        "seconds": round(time.time() - t0, 1),
    }
    for path, n in ((out_path, n_lines), (plain_path, n_plain)):
        meta = path.with_suffix(path.suffix + ".meta.json")
        meta.write_text(json.dumps(dict(info, records=n), indent=1) + "\n")
        print(f"wrote {meta}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
