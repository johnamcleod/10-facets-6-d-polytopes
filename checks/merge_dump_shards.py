#!/usr/bin/env python3
"""Merge the sharded instance dumps of the lemmas=valid run into one sink.

run_wild_dump.py writes one wildcard sink and one plain sink per invocation, each
with a completion marker recording its record count.  The 304-type run was dumped
in several shards, run in parallel over disjoint subtree lists.  This script:

  * refuses to merge unless every shard's marker is present, says complete, and
    matches its file's line count;
  * checks that the shards' subtree lists are disjoint and together cover every
    subtree of the run of record with a nonzero labelling count, and nothing else;
  * checks the merged record counts against the run of record's own counters
    (labellings, and wildcard-bearing labellings, summed over those subtrees);
  * writes runs/d6_n10/d6_valid/{wild,plain}_instances.jsonl with markers.

Run:  python3 checks/merge_dump_shards.py
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "runs/d6_n10/d6_valid"
DUMP = RUN / "dump"


def records(path):
    return [l for l in path.read_text().splitlines() if l.strip()]


def meta(path):
    return json.loads(path.with_suffix(path.suffix + ".meta.json").read_text())


def main():
    state = {}
    for l in (RUN / "state.jsonl").read_text().splitlines():
        if l.strip():
            r = json.loads(l)
            state[r.pop("key")] = r
    want = {k for k, v in state.items() if v.get("enum_count", 0)}
    want_all = sum(state[k]["enum_count"] for k in want)
    want_wild = sum((state[k].get("diag") or {}).get("wild_assignments", 0) for k in want)

    shards = sorted(DUMP.glob("wild_shard*.jsonl"))
    seen, wild, plain = set(), [], []
    for w in shards:
        p = w.with_name(w.name.replace("wild_", "plain_"))
        for f in (w, p):
            m = meta(f)
            n = len(records(f))
            if not (m.get("complete") and m.get("records") == n):
                print(f"INCOMPLETE shard {f.name}: marker {m.get('records')}, "
                      f"file {n}", file=sys.stderr)
                return 2
        subs = set(meta(w)["subtrees"])
        if subs & seen:
            print(f"overlapping subtrees in {w.name}: {sorted(subs & seen)}",
                  file=sys.stderr)
            return 2
        seen |= subs
        wild += records(w)
        plain += records(p)

    if seen != want:
        print(f"shards cover {len(seen)} subtrees, the run has {len(want)} with "
              f"labellings; missing {sorted(want - seen)}, extra {sorted(seen - want)}",
              file=sys.stderr)
        return 2
    if len(wild) != want_wild or len(wild) + len(plain) != want_all:
        print(f"merged {len(wild)} wild + {len(plain)} plain, run of record "
              f"counts {want_wild} wild of {want_all}", file=sys.stderr)
        return 2

    for name, recs in (("wild_instances.jsonl", wild), ("plain_instances.jsonl", plain)):
        out = RUN / name
        out.write_text("\n".join(recs) + "\n")
        out.with_suffix(out.suffix + ".meta.json").write_text(json.dumps({
            "complete": True, "d": 6, "records": len(recs),
            "subtrees": sorted(seen), "merged_from": [s.name for s in shards],
        }, indent=1) + "\n")
        print(f"{name}: {len(recs)} records")
    print(f"RESULT: {len(shards)} shards merged; {len(seen)} subtrees, "
          f"{len(wild) + len(plain)} labellings, matching the run of record")
    return 0


if __name__ == "__main__":
    sys.exit(main())
