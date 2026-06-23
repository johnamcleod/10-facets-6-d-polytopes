#!/usr/bin/env python3
"""Full d=4 Stage-4 production run: enumerate + structured-screen + exact solve
over all surviving combinatorial types, then deduplicate raw configs to distinct
compact hyperbolic Coxeter 4-polytopes (target: 348).

Usage:
    python run_d4.py                      # all types, exhaustive
    python run_d4.py --types 5,6,18       # subset (dry run)
    python run_d4.py --workers 8 --timeout 7200
"""
import argparse, json, time, itertools
from pathlib import Path
from collections import defaultdict
import numpy as np

from pipeline.stage4_gram import run_stage4

N, D = 8, 4


def config_to_gram(r, n):
    G = np.eye(n)
    for es, lab in r["label_assignment"].items():
        i, j = eval(es)
        G[i, j] = G[j, i] = -np.cos(np.pi / lab)
    for es, dv in r.get("dot_values", {}).items():
        p = es.split("_"); i, j = int(p[1]), int(p[2])
        v = float(dv["value"]) if isinstance(dv, dict) else float(dv)
        G[i, j] = G[j, i] = -v
    return G


def is_valid(G, d, tol=1e-4):
    evals = np.linalg.eigvalsh(G)
    return int(np.sum(evals < -tol)) == 1 and int(np.sum(np.abs(evals) < tol)) == n_kernel(G.shape[0], d)


def n_kernel(n, d):
    return n - d - 1


def canonical_key(r, n, dot_eps=2e-4):
    L = np.zeros((n, n))
    for es, lab in r["label_assignment"].items():
        i, j = eval(es); L[i, j] = L[j, i] = lab
    G = config_to_gram(r, n)
    certs = [tuple(sorted(int(L[i, j]) for j in range(n) if j != i)) for i in range(n)]
    classes = defaultdict(list)
    for i, c in enumerate(certs):
        classes[c].append(i)
    groups = [classes[k] for k in sorted(classes)]
    total = 1
    for g in groups:
        for x in range(2, len(g) + 1):
            total *= x
    if total > 2_000_000:
        # coarse fallback: sorted multiset of entries (iso-invariant, not exact)
        return ("coarse", tuple(sorted(round(G[i, j], 6) for i in range(n)
                                       for j in range(i + 1, n))))

    def ek(v):
        return round(v, 6) if abs(v) <= 1 + 1e-10 else round(round(v / dot_eps) * dot_eps, 8)
    best = None
    for combo in itertools.product(*[itertools.permutations(g) for g in groups]):
        perm = [x for grp in combo for x in grp]
        flat = tuple(ek(G[perm[i], perm[j]]) for i in range(n) for j in range(i + 1, n))
        if best is None or flat < best:
            best = flat
    return ("exact", best)


def deduplicate(results, n, d):
    seen, invalid = set(), 0
    for r in results:
        G = config_to_gram(r, n)
        if not is_valid(G, d):
            invalid += 1
            continue
        seen.add(canonical_key(r, n))
    return len(seen), invalid


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--types", default=None, help="comma list of type_ids (default all)")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--timeout", type=float, default=7200.0, help="per-type budget (s)")
    ap.add_argument("--stage3", default="runs/d4_n8/stage3")
    ap.add_argument("--out", default="runs/d4_n8/stage4_prod")
    args = ap.parse_args()

    survivors = json.load(open(Path(args.stage3) / "surviving_types.json"))
    if args.types:
        want = set(int(x) for x in args.types.split(","))
        survivors = [t for t in survivors if t["type_id"] in want]
    print(f"d=4 production run: {len(survivors)} types, {args.workers} workers, "
          f"per-type budget {args.timeout}s")

    t0 = time.time()
    results = run_stage4(D, survivors, output_dir=args.out, verbose=True,
                         per_type_timeout=args.timeout, n_workers=args.workers,
                         label_indices=None, max_assignments=10**9,
                         face_tuples_max_size=4)
    elapsed = time.time() - t0
    distinct, invalid = deduplicate(results, N, D)
    print(f"\n{'='*56}")
    print(f"d=4 RESULT: {len(results)} raw configs -> {distinct} distinct polytopes "
          f"({invalid} invalid) in {elapsed/60:.1f} min")
    print(f"  target (Ma-Zheng / Burcroff): 348")
    print(f"{'='*56}")
    Path(args.out).mkdir(parents=True, exist_ok=True)
    (Path(args.out) / "summary.json").write_text(json.dumps(
        {"raw_configs": len(results), "distinct": distinct, "invalid": invalid,
         "target": 348, "elapsed_s": round(elapsed, 1),
         "n_types": len(survivors)}, indent=2))


if __name__ == "__main__":
    main()
