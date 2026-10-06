#!/usr/bin/env python3
"""Generate the combinatorial types of simple d-polytopes with d+4 facets and p >= 2
from the order-type database, with the exact affine-Gale criterion, sharded and
resumable.

This is the generator of `pipeline/stage2_gale.py` (`process_order_type`: the exact
face criterion of `pipeline/utils/gale_exact.py`, the polytopality gate, missing
faces up to size 5, p >= 2), driven over record ranges of the `.chi` cache so that
the 14,309,547 order types on 10 points can be processed in pieces on whatever
cores are free.  The same function reproduces the published d=4 (30/30) and d=5
(109/109) censuses (`tests/test_generator_coverage.py`).

    python3 generate_types.py shard I [--d 6] [--size 50000]   # one record range
    python3 generate_types.py run [--d 6] [--max-workers 6] [--idle 1.5]
    python3 generate_types.py merge [--d 6] [--compare runs/d6_n10/stage2/types.json]

`run` starts shards at nice 19, one at a time, whenever at least --idle cores are
measured idle, up to --max-workers; finished shards are skipped, so it can be
stopped and restarted.  `merge` deduplicates the raw missing-face systems of all
shards exactly (colour refinement, then search over colour-preserving
permutations), writes `types.json`, and, with --compare, checks that the result is
the given type list up to relabelling of the facets.

Input: `data/aak/otypes{n}.chi`, the chirotope cache that `pipeline/c/aak_parse`
builds from the database file `otypes{n}.b16`.
Output: `runs/d{d}_n{n}/stage2_exact/`.
"""
from __future__ import annotations

import argparse
import json
import os
import struct
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

MAGIC = 0x41414B01


def chi_path(n):
    return ROOT / f"data/aak/otypes{n:02d}.chi"


def out_dir(d):
    return ROOT / f"runs/d{d}_n{d + 4}/stage2_exact"


def open_chi(n):
    """Memory-map the .chi cache: (num_records, chirotopes int8, points uint16)."""
    p = chi_path(n)
    with open(p, "rb") as f:
        magic, nn, num = struct.unpack("<IIQ", f.read(16))
    if magic != MAGIC or nn != n:
        raise ValueError(f"{p}: bad header (magic {magic:#x}, n {nn})")
    tri = n * (n - 1) * (n - 2) // 6
    row = tri + 4 * n
    body = np.memmap(p, dtype=np.uint8, mode="r", offset=16, shape=(num, row))
    return num, tri, body


def shard_file(d, i):
    return out_dir(d) / "shards" / f"shard_{i:04d}.json"


def do_shard(d, i, size):
    from pipeline.stage2_gale import (_batch_prefilter, _build_prefilter_tables,
                                      process_order_type)
    n = d + 4
    num, tri, body = open_chi(n)
    lo, hi = i * size, min((i + 1) * size, num)
    if lo >= num:
        raise SystemExit(f"shard {i} is past the end ({num} records)")
    chunk = np.array(body[lo:hi])
    chi = chunk[:, :tri].view(np.int8)
    pts = chunk[:, tri:].view(np.uint16).reshape(hi - lo, n, 2)
    mask = _batch_prefilter(chi, _build_prefilter_tables(n))
    t0 = time.time()
    raw = {}                       # labelled missing-face system -> entry
    for k in np.where(mask)[0]:
        p = [tuple(map(int, q)) for q in pts[k]]
        for _key, mf, ag in process_order_type(None, p, d, exact_dedup=False):
            key = tuple(sorted(tuple(sorted(m)) for m in mf))
            rec = lo + int(k)
            e = raw.get(key)
            if e is None:
                raw[key] = {"missing_faces": [list(m) for m in key],
                            "example_record": rec,
                            "example_points": [list(q) for q in p],
                            "example_positive": sorted(ag.positive),
                            "records": [rec]}
            elif e["records"][-1] != rec:
                e["records"].append(rec)
    f = shard_file(d, i)
    f.parent.mkdir(parents=True, exist_ok=True)
    tmp = f.with_suffix(".tmp")
    tmp.write_text(json.dumps({"d": d, "lo": lo, "hi": hi, "passed_prefilter":
                               int(mask.sum()), "seconds": round(time.time() - t0, 1),
                               "raw": list(raw.values())}, separators=(",", ":")))
    tmp.rename(f)                  # a shard file exists only once it is complete
    print(f"shard {i}: records {lo}..{hi - 1}, {len(raw)} raw systems, "
          f"{time.time() - t0:.0f}s", flush=True)


def idle_cores(sample=5):
    """Idle CPU, in cores, measured over `sample` seconds (macOS `top`, else load)."""
    try:
        out = subprocess.run(["top", "-l", "2", "-n", "0", "-s", str(sample)],
                             capture_output=True, text=True, timeout=60).stdout
        line = [l for l in out.splitlines() if l.startswith("CPU usage")][-1]
        pct = float(line.split(",")[-1].split("%")[0])
        return pct / 100 * os.cpu_count()
    except Exception:
        return os.cpu_count() - os.getloadavg()[0]


def run(d, size, max_workers, idle):
    n = d + 4
    num, _, _ = open_chi(n)
    nshards = (num + size - 1) // size
    todo = [i for i in range(nshards) if not shard_file(d, i).exists()]
    log = out_dir(d) / "logs"
    log.mkdir(parents=True, exist_ok=True)
    print(f"{nshards} shards of {size}, {len(todo)} to do", flush=True)
    procs = {}
    while todo or procs:
        for i, pr in list(procs.items()):
            if pr.poll() is not None:
                del procs[i]
                if pr.returncode != 0:
                    print(f"shard {i} failed (exit {pr.returncode}); see logs", flush=True)
        if todo and len(procs) < max_workers and idle_cores() >= idle:
            i = todo.pop(0)
            env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                       VECLIB_MAXIMUM_THREADS="1")
            with open(log / f"shard_{i:04d}.log", "w") as fh:
                procs[i] = subprocess.Popen(
                    ["nice", "-n", "19", sys.executable, __file__, "shard", str(i),
                     "--d", str(d), "--size", str(size)],
                    stdout=fh, stderr=subprocess.STDOUT, env=env)
            print(f"{time.strftime('%H:%M')} started shard {i} "
                  f"({len(procs)} running, {len(todo)} queued)", flush=True)
            time.sleep(90)         # let the load average register the new worker
        else:
            time.sleep(30)
    print("all shards done", flush=True)


def merge(d, size, compare):
    import networkx as nx
    from networkx.algorithms import isomorphism as iso
    from pipeline.utils.canonical import canonical_missing_face_hypergraph
    n = d + 4
    num, _, _ = open_chi(n)
    nshards = (num + size - 1) // size
    missing = [i for i in range(nshards) if not shard_file(d, i).exists()]
    if missing:
        raise SystemExit(f"{len(missing)} shards not done, e.g. {missing[:5]}")
    types, seconds = {}, 0.0
    for i in range(nshards):
        s = json.loads(shard_file(d, i).read_text())
        seconds += s["seconds"]
        for e in s["raw"]:
            key = canonical_missing_face_hypergraph(e["missing_faces"], n, exact=True)
            t = types.get(key)
            if t is None:
                types[key] = {"canonical_key": [list(m) for m in key],
                              "missing_faces": e["missing_faces"],
                              "p_count": sum(1 for m in e["missing_faces"] if len(m) == 2),
                              "example_points": e["example_points"],
                              "example_positive": e["example_positive"],
                              "source_order_type_ids": set(e["records"])}
            else:
                t["source_order_type_ids"].update(e["records"])
    out = []
    for k, t in enumerate(sorted(types.values(),
                                 key=lambda t: min(t["source_order_type_ids"]))):
        t["type_id"] = k
        t["source_order_type_ids"] = sorted(t["source_order_type_ids"])
        out.append(t)
    (out_dir(d) / "types.json").write_text(json.dumps(out, separators=(",", ":")))
    print(f"{len(out)} types from {num} order types, {seconds / 3600:.1f} CPU-hours "
          f"-> {(out_dir(d) / 'types.json').relative_to(ROOT)}")
    if not compare:
        return 0

    def graph(mf):
        g = nx.Graph()
        g.add_nodes_from((("v", i) for i in range(n)), side=0)
        for j, m in enumerate(mf):
            g.add_node(("e", j), side=1)
            g.add_edges_from((("v", i), ("e", j)) for i in m)
        return g

    def inv(mf):
        return tuple(sorted(len(m) for m in mf))

    ref = json.loads(Path(compare).read_text())
    nm = iso.categorical_node_match("side", None)
    pool = {}
    for t in out:
        pool.setdefault(inv(t["missing_faces"]), []).append(t)
    matched, unmatched_ref = set(), []
    for r in ref:
        gr = graph(r["missing_faces"])
        hit = next((t for t in pool.get(inv(r["missing_faces"]), [])
                    if t["type_id"] not in matched
                    and nx.is_isomorphic(gr, graph(t["missing_faces"]), node_match=nm)),
                   None)
        if hit is None:
            unmatched_ref.append(r["type_id"])
        else:
            matched.add(hit["type_id"])
    new = [t["type_id"] for t in out if t["type_id"] not in matched]
    res = {"generated": len(out), "reference": len(ref),
           "reference_not_generated": unmatched_ref, "generated_not_in_reference": new}
    (out_dir(d) / "comparison.json").write_text(json.dumps(res, indent=1) + "\n")
    print(f"compared with {compare}: {len(ref) - len(unmatched_ref)}/{len(ref)} matched; "
          f"missing from ours {unmatched_ref or 'none'}; new in ours {new or 'none'}")
    return 0 if not unmatched_ref and not new else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["shard", "run", "merge"])
    ap.add_argument("index", nargs="?", type=int)
    ap.add_argument("--d", type=int, default=6)
    ap.add_argument("--size", type=int, default=50_000)
    ap.add_argument("--max-workers", type=int, default=6)
    ap.add_argument("--idle", type=float, default=1.5,
                    help="start a shard only while at least this many cores are idle")
    ap.add_argument("--compare", default=None)
    a = ap.parse_args()
    if a.cmd == "shard":
        do_shard(a.d, a.index, a.size)
    elif a.cmd == "run":
        run(a.d, a.size, a.max_workers, a.idle)
    else:
        return merge(a.d, a.size, a.compare)
    return 0


if __name__ == "__main__":
    sys.exit(main())
