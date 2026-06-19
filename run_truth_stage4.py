#!/usr/bin/env python3
"""
Bypass runner: feed PUBLISHED combinatorial types (Ma-Zheng vertex flags)
straight into Stage 4, skipping the broken order-type -> Gale-diagram generator.

This isolates Stage 4 (label enumeration + exact Gram realizability) and tells
us how many of the known 348 (d=4) / 51 (d=5) polytopes the current solver can
recover GIVEN correct combinatorial types.

Usage:
    python run_truth_stage4.py --d 4 --timeout 120 --workers 8
"""

import argparse
import json
import time
from pathlib import Path

from pipeline.validate_coverage import load_truth_as_stage3
from pipeline.stage4_gram import run_stage4
from run_regression import deduplicate_gram_results

EXPECTED = {4: 348, 5: 51}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--d", type=int, required=True, choices=[4, 5])
    ap.add_argument("--timeout", type=float, default=120.0, help="per-type timeout (s)")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--max-assignments", type=int, default=2_000_000)
    ap.add_argument("--truth", default=None)
    ap.add_argument("--limit", type=int, default=None, help="only first N types (debug)")
    args = ap.parse_args()

    n = args.d + 4
    truth = args.truth or f"data/ground_truth/{args.d}d{n}m.txt"
    types = load_truth_as_stage3(truth, args.d)
    if args.limit:
        types = types[: args.limit]

    print(f"{'='*60}")
    print(f"Bypass Stage 4 from ground truth: d={args.d}, n={n}")
    print(f"  truth types: {len(types)}  (expect {EXPECTED[args.d]} polytopes)")
    print(f"  labels: ALL (incl. pi/7..pi/12);  per-type timeout {args.timeout}s")
    print(f"{'='*60}")

    out_dir = f"runs/d{args.d}_n{n}/stage4_truth"
    t0 = time.time()
    gram = run_stage4(
        args.d, types,
        output_dir=out_dir,
        verbose=True,
        per_type_timeout=args.timeout,
        n_workers=args.workers,
        label_indices=None,            # d=4 genuinely uses labels >= 7
        max_assignments=args.max_assignments,
        face_tuples_max_size=4,
    )
    elapsed = time.time() - t0

    unique, invalid = deduplicate_gram_results(gram, n, args.d)
    print(f"\n{'='*60}")
    print(f"  truth types:        {len(types)}")
    print(f"  raw Gram configs:   {len(gram)}")
    print(f"  invalid signature:  {invalid}")
    print(f"  DISTINCT polytopes: {unique}   (expected {EXPECTED[args.d]})")
    print(f"  elapsed:            {elapsed/60:.1f} min")
    print(f"{'='*60}")

    Path(out_dir).mkdir(parents=True, exist_ok=True)
    (Path(out_dir) / "truth_summary.json").write_text(json.dumps({
        "d": args.d, "n": n, "truth_types": len(types),
        "raw_configs": len(gram), "invalid": invalid,
        "unique_polytopes": unique, "expected": EXPECTED[args.d],
        "elapsed_s": round(elapsed, 1),
    }, indent=2))


if __name__ == "__main__":
    main()
