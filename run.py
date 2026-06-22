#!/usr/bin/env python3
"""
Main pipeline driver for compact hyperbolic Coxeter d-polytopes with d+4 facets.

Usage:
    python run.py --d 4          # Regression: d=4, expect 348 polytopes
    python run.py --d 5          # Regression: d=5, expect 51 polytopes
    python run.py --d 6          # The target: d=6, n=10
    python run.py --d 4 --stages 2,3,4   # Skip stage 1 (default)

Stage 1 is skipped by default: Stage 2 streams directly from the AAK .chi
file via the C filter, avoiding materialising all order types as Python objects.
"""

import argparse
import json
import sys
import time
from pathlib import Path

from pipeline.stage2_gale import run_stage2_from_chi, load_stage2
from pipeline.stage3_filters import run_stage3, load_stage3
from pipeline.stage4_gram import run_stage4

PER_TYPE_TIMEOUT = 120.0   # seconds per type for Stage 4 (d=4/5); d=6 uses 600s

CHI_PATHS = {
    4: Path("data/aak/otypes08.chi"),
    5: Path("data/aak/otypes09.chi"),
    6: Path("data/aak/otypes10.chi"),
}


def stage_dirs(d):
    n = d + 4
    base = Path(f"runs/d{d}_n{n}")
    return {
        2: base / "stage2",
        3: base / "stage3",
        4: base / "stage4",
        "result": base / "RESULT",
    }


def run_pipeline(d, force=False):
    n = d + 4
    dirs = stage_dirs(d)
    t_total = time.time()

    chi_path = CHI_PATHS.get(d)
    if chi_path is None:
        print(f"No CHI_PATH configured for d={d}")
        sys.exit(1)
    if not chi_path.exists():
        print(f"Chi file not found: {chi_path}")
        print(f"Run the C parser first: pipeline/c/aak_parse <b16_file> {n} {chi_path}")
        sys.exit(1)

    print(f"\n{'='*60}")
    print(f"Pipeline: d={d}, n={n} (compact hyperbolic Coxeter polytopes)")
    print(f"{'='*60}\n")

    # Stage 2
    s2_dir = dirs[2]
    if not force and s2_dir.exists():
        print("Stage 2: loading cached results...")
        types = load_stage2(str(s2_dir))
        print(f"  Loaded {len(types):,} combinatorial types")
    else:
        if force and s2_dir.exists():
            import shutil; shutil.rmtree(s2_dir)
        # Python streaming path with the exact affine-Gale criterion
        # (gale_exact.AffineGale).  The C/gd2 path is disabled until the C
        # criterion is ported to match (the old C filter is buggy; see
        # memory stage2-generator-bugs).  d=6 (14.3M order types) will be slow
        # here and needs the C port before it is practical.
        types = run_stage2_from_chi(d, chi_path, output_dir=str(s2_dir),
                                    exact_dedup=True)
    print(f"  Stage 2 done: {len(types):,} combinatorial types\n")

    # Stage 3
    s3_dir = dirs[3]
    if not force and s3_dir.exists():
        print("Stage 3: loading cached results...")
        survivors = load_stage3(str(s3_dir))
        print(f"  Loaded {len(survivors):,} surviving types")
    else:
        if force and s3_dir.exists():
            import shutil; shutil.rmtree(s3_dir)
        survivors, _ = run_stage3(d, types, output_dir=str(s3_dir))
    print(f"  Stage 3 done: {len(survivors):,} surviving types\n")

    # Stage 4
    # For d>=5, Burcroff's theorem bounds dihedral angles to pi/m with m in {2,3,4,5}.
    # Restricting to these 4 labels enables bitmask forward-checking for size-5
    # sub-vertex groups (1,386 valid PD combos vs 4,096 limit), dramatically
    # reducing the search space.
    label_indices = (0, 1, 2, 3) if d >= 5 else None  # labels {2,3,4,5} for d>=5
    # d=6 backtracking is slow (hard type takes ~200s); use 600s timeout.
    per_type_timeout = 600.0 if d == 6 else PER_TYPE_TIMEOUT
    gram_results = run_stage4(
        d, survivors,
        output_dir=str(dirs[4]),
        verbose=True,
        per_type_timeout=per_type_timeout,
        n_workers=None,
        label_indices=label_indices,
    )
    print(f"  Stage 4 done: {len(gram_results)} valid Gram configurations\n")

    elapsed = time.time() - t_total
    print(f"\n{'='*60}")
    print(f"Summary for d={d}, n={n}:")
    print(f"  Combinatorial types:  {len(types):>8,}")
    print(f"  Stage 3 survivors:    {len(survivors):>8,}")
    print(f"  Valid Gram configs:   {len(gram_results):>8,}")
    print(f"  Total time:           {elapsed/60:.1f} min")
    print(f"{'='*60}\n")

    result = {
        "d": d, "n": n,
        "num_comb_types": len(types),
        "num_survivors": len(survivors),
        "num_gram": len(gram_results),
        "elapsed_s": round(elapsed, 1),
    }
    dirs["result"].mkdir(parents=True, exist_ok=True)
    (dirs["result"] / "summary.json").write_text(json.dumps(result, indent=2))
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Compact hyperbolic Coxeter polytope classifier"
    )
    parser.add_argument("--d", type=int, default=6, help="Polytope dimension (default 6)")
    parser.add_argument("--force", action="store_true",
                        help="Recompute even if cached results exist")
    args = parser.parse_args()
    run_pipeline(args.d, force=args.force)


if __name__ == "__main__":
    main()
