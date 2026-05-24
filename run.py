#!/usr/bin/env python3
"""
Main pipeline driver for compact hyperbolic Coxeter d-polytopes with d+4 facets.

Usage:
    python run.py --d 4          # Regression: d=4, expect 348 polytopes
    python run.py --d 5          # Regression: d=5, expect 51 polytopes
    python run.py --d 6          # The target: d=6, n=10
    python run.py --d 4 --stages 1,2,3   # Run specific stages only
"""

import argparse
import json
import sys
from pathlib import Path

from pipeline.stage1_order_types import run_stage1, load_stage1
from pipeline.stage2_gale import run_stage2, load_stage2
from pipeline.stage3_filters import run_stage3, load_stage3
from pipeline.stage4_gram import run_stage4


def stage_dirs(d):
    n = d + 4
    base = Path(f"runs/d{d}_n{n}")
    return {
        1: base / "stage1",
        2: base / "stage2",
        3: base / "stage3",
        4: base / "stage4_gram",
        5: base / "stage5_verified",
        6: base / "stage6_crosscheck",
        "result": base / "RESULT",
    }


def run_pipeline(d, stages=None, force=False, data_dir="data/aak"):
    if stages is None:
        stages = [1, 2, 3, 4]

    n = d + 4
    dirs = stage_dirs(d)
    print(f"\n{'='*60}")
    print(f"Pipeline: d={d}, n={n} (compact hyperbolic Coxeter polytopes)")
    print(f"{'='*60}\n")

    # Stage 1
    s1_dir = dirs[1]
    if 1 in stages:
        ot = run_stage1(d, output_dir=str(s1_dir), data_dir=data_dir)
    elif s1_dir.exists():
        print(f"Stage 1: loading from {s1_dir}")
        ot = load_stage1(str(s1_dir))
    else:
        print("Stage 1: not run and no cached results. Run with --stages 1")
        sys.exit(1)

    print(f"  Stage 1 done: {len(ot)} order types\n")

    # Stage 2
    s2_dir = dirs[2]
    if 2 in stages:
        types = run_stage2(d, ot, output_dir=str(s2_dir))
    elif s2_dir.exists():
        print(f"Stage 2: loading from {s2_dir}")
        types = load_stage2(str(s2_dir))
    else:
        print("Stage 2: not run. Run with --stages 1,2")
        sys.exit(1)

    print(f"  Stage 2 done: {len(types)} combinatorial types\n")

    # Stage 3
    s3_dir = dirs[3]
    if 3 in stages:
        survivors, elim_log = run_stage3(d, types, output_dir=str(s3_dir))
    elif s3_dir.exists():
        print(f"Stage 3: loading from {s3_dir}")
        survivors = load_stage3(str(s3_dir))
    else:
        print("Stage 3: not run. Run with --stages 1,2,3")
        sys.exit(1)

    print(f"  Stage 3 done: {len(survivors)} surviving types\n")

    # Stage 4
    if 4 in stages:
        gram_results = run_stage4(d, survivors, output_dir=str(dirs[4]))
        print(f"  Stage 4 done: {len(gram_results)} valid Gram configurations\n")
    else:
        gram_results = []

    # Summary
    print(f"\n{'='*60}")
    print(f"Summary for d={d}, n={n}:")
    print(f"  Order types: {len(ot)}")
    print(f"  Combinatorial types (after p>=2): {len(types)}")
    print(f"  Surviving types (after Stage 3): {len(survivors)}")
    print(f"  Valid Gram configurations: {len(gram_results)}")
    print(f"{'='*60}\n")

    return {
        "d": d, "n": n,
        "num_order_types": len(ot),
        "num_comb_types": len(types),
        "num_survivors": len(survivors),
        "num_gram": len(gram_results),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Compact hyperbolic Coxeter polytope classifier"
    )
    parser.add_argument("--d", type=int, default=6, help="Polytope dimension (default 6)")
    parser.add_argument(
        "--stages", type=str, default="1,2,3,4",
        help="Comma-separated list of stages to run (default: 1,2,3,4)"
    )
    parser.add_argument("--data-dir", type=str, default="data/aak",
                        help="Directory containing AAK order-type database files")
    parser.add_argument("--force", action="store_true",
                        help="Recompute even if cached results exist")
    args = parser.parse_args()

    stages = [int(s) for s in args.stages.split(",")]
    result = run_pipeline(args.d, stages=stages, force=args.force,
                          data_dir=args.data_dir)

    # Write summary
    dirs = stage_dirs(args.d)
    dirs["result"].mkdir(parents=True, exist_ok=True)
    (dirs["result"] / "summary.json").write_text(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
