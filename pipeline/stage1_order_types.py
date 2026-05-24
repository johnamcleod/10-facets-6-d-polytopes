"""
Stage 1: Enumerate realizable planar order types on n points.

For d=4: n = d+4 = 8 points
For d=5: n = d+4 = 9 points
For d=6: n = d+4 = 10 points

Data source priority:
1. AAK order-type database (data/aak/)
2. Exhaustive generation via random sampling + deduplication (for small n)
3. Future: NumPSLA integration

Each order type is stored as:
  - A canonical chirotope tuple (the lexicographically minimal relabeling)
  - An explicit integer-coordinate realization
"""

import json
import struct
from pathlib import Path
from itertools import combinations, permutations
import random
import math

from pipeline.utils.chirotope import (
    compute_chirotope, chirotope_to_tuple, canonical_chirotope,
    verify_chirotope_from_points, relabel_chirotope
)
from pipeline.utils.manifest import write_manifest

# Known counts from the literature (Aichholzer-Aurenhammer-Krasser)
# These are the numbers of realizable order types on n points in general position.
KNOWN_COUNTS = {
    3: 1,
    4: 2,
    5: 3,
    6: 16,
    7: 135,
    8: 3315,
    9: 158817,
    10: 14309547,
    11: 2334512907,
}


def _load_aak_database(n, data_dir="data/aak"):
    """Attempt to load the AAK order-type database for n points.

    The database files from Aichholzer et al. are binary files with
    sequences of n (x,y) coordinate pairs (typically 2-byte integers).

    Returns list of (canonical_chi_tuple, points) or None if not found.
    """
    data_dir = Path(data_dir)
    # Common AAK filename patterns
    candidates = [
        data_dir / f"otypes{n:02d}.b16",
        data_dir / f"otypes{n}.b16",
        data_dir / f"order{n}.dat",
        data_dir / f"n{n}.ot",
    ]
    for path in candidates:
        if path.exists():
            return _parse_aak_file(path, n)
    return None


def _parse_aak_file(path, n):
    """Parse an AAK binary database file.

    Format: each entry is n pairs of 16-bit unsigned integers (x, y).
    """
    results = []
    data = path.read_bytes()
    record_size = n * 2 * 2  # n points, 2 coords, 2 bytes each
    num_records = len(data) // record_size

    for i in range(num_records):
        offset = i * record_size
        coords = struct.unpack_from(f'{n * 2}H', data, offset)
        points = [(coords[2*j], coords[2*j+1]) for j in range(n)]

        chi = compute_chirotope(points)
        if chi is None:
            continue  # skip degenerate
        chi_tuple = chirotope_to_tuple(chi, n)
        canonical, _ = canonical_chirotope(chi_tuple, n)
        results.append((canonical, points))

    return results


def _generate_random_order_types(n, num_samples=100000, seed=42, coord_range=1000):
    """Generate realizable order types by random sampling.

    This is NOT exhaustive for large n, but sufficient for validation
    and small cases. Returns dict: canonical_chi_tuple -> points.
    """
    rng = random.Random(seed)
    seen = {}

    for _ in range(num_samples):
        # Generate n random integer points
        while True:
            pts = [(rng.randint(-coord_range, coord_range),
                    rng.randint(-coord_range, coord_range)) for _ in range(n)]
            chi = compute_chirotope(pts)
            if chi is not None:
                break

        chi_tuple = chirotope_to_tuple(chi, n)
        canonical, _ = canonical_chirotope(chi_tuple, n)

        if canonical not in seen:
            seen[canonical] = pts

    return seen


def _generate_small_exhaustive(n, coord_range=100, num_attempts=5000, seed=0):
    """Generate order types by random sampling (not exhaustive for large n).

    Uses raw chirotope tuples as keys — fast, no n! canonicalization.
    Different labelings of the same order type may appear as separate entries;
    Stage 2 deduplication (on combinatorial types) handles this correctly.
    """
    rng = random.Random(seed)
    seen = {}

    for _ in range(num_attempts):
        pts = [(rng.randint(-coord_range, coord_range),
                rng.randint(-coord_range, coord_range)) for _ in range(n)]
        chi = compute_chirotope(pts)
        if chi is None:
            continue
        chi_tuple = chirotope_to_tuple(chi, n)
        if chi_tuple not in seen:
            seen[chi_tuple] = pts

    return seen


def run_stage1(d, output_dir=None, data_dir="data/aak", seed=42):
    """Run Stage 1: enumerate realizable order types on n = d+4 points.

    Returns:
        list of dicts with keys: 'chirotope', 'points', 'id'

    Writes output to output_dir/order_types/ if provided.
    """
    n = d + 4
    print(f"Stage 1: d={d}, n={n} (Gale dimension = {n-d-2})")

    # Try AAK database first
    print(f"  Trying AAK database for n={n}...")
    aak_results = _load_aak_database(n, data_dir)

    if aak_results is not None:
        print(f"  Loaded {len(aak_results)} order types from AAK database.")
        order_types = [
            {"id": i, "chirotope": list(chi), "points": list(pts)}
            for i, (chi, pts) in enumerate(aak_results)
        ]
        source = "aak_database"
    else:
        print(f"  AAK database not found. Using random sampling.")
        print(f"  NOTE: Sampling is NOT exhaustive. For a complete run,")
        print(f"        download the AAK database to data/aak/.")

        seen = _generate_small_exhaustive(n, seed=seed)

        print(f"  Found {len(seen)} distinct raw chirotopes by sampling.")
        order_types = [
            {"id": i, "chirotope": list(chi), "points": list(pts)}
            for i, (chi, pts) in enumerate(seen.items())
        ]
        source = "random_sampling"

    # Acceptance test: verify a sample of C(n,3) = 120 triple orientations
    print(f"  Running acceptance test (spot-checking up to 20 entries)...")
    sample_ids = random.Random(seed).sample(
        range(len(order_types)), min(20, len(order_types))
    )
    for idx in sample_ids:
        ot = order_types[idx]
        pts = [tuple(p) for p in ot["points"]]
        chi_stored = tuple(ot["chirotope"])
        assert verify_chirotope_from_points(chi_stored, pts), (
            f"Acceptance test FAILED: chirotope mismatch for entry {idx}"
        )
    print(f"  Acceptance test passed.")

    # Check count against known
    if n in KNOWN_COUNTS and source == "aak_database":
        expected = KNOWN_COUNTS[n]
        actual = len(order_types)
        if actual != expected:
            print(f"  WARNING: expected {expected} order types, got {actual}.")
        else:
            print(f"  Count matches known value: {actual} order types.")

    # Write output
    if output_dir is not None:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        (out / "order_types.json").write_text(
            json.dumps(order_types, indent=None, separators=(',', ':'))
        )
        write_manifest(out, "stage1", {"d": d, "n": n, "source": source},
                       {"num_order_types": len(order_types)})
        print(f"  Written to {out}/")

    return order_types


def load_stage1(output_dir):
    """Load Stage 1 results from disk."""
    path = Path(output_dir) / "order_types.json"
    return json.loads(path.read_text())
