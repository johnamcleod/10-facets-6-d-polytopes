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
    sequences of n (x,y) coordinate pairs (1-byte for .b08, 2-byte for .b16).

    Returns list of (canonical_chi_tuple, points) or None if not found.
    """
    data_dir = Path(data_dir)
    # (path, coord_bytes) — try b08 first for n<=8, b16 otherwise
    candidates = [
        (data_dir / f"otypes{n:02d}.b08", 1),
        (data_dir / f"otypes{n:02d}.b16", 2),
        (data_dir / f"otypes{n}.b16",     2),
        (data_dir / f"order{n}.dat",      2),
        (data_dir / f"n{n}.ot",           2),
    ]
    for path, coord_bytes in candidates:
        if path.exists():
            return _parse_aak_file(path, n, coord_bytes)
    return None


def _chi_path(b16_path):
    """Return the preprocessed .chi path for a given .b16 path."""
    return b16_path.with_suffix(".chi")


def _build_chi_file(b_path, n, coord_bytes=2):
    """Run the C parser to convert .b08/.b16 → chi.  Raises RuntimeError on failure."""
    import subprocess
    c_dir = Path(__file__).parent / "c"
    binary = c_dir / "aak_parse"
    if not binary.exists():
        raise RuntimeError(
            f"C parser not built.  Run:  cd {c_dir} && make"
        )
    chi = _chi_path(b_path)
    print(f"    Building {chi.name} via C parser (one-time)...", flush=True)
    result = subprocess.run(
        [str(binary), str(b_path), str(n), str(chi), str(coord_bytes)],
        capture_output=False,
    )
    if result.returncode != 0:
        raise RuntimeError(f"aak_parse failed with exit code {result.returncode}")
    return chi


def _parse_chi_file(chi_path, n):
    """Read a preprocessed .chi file produced by aak_parse into Python lists.

    File layout (little-endian):
      uint32  magic = 0x41414B01
      uint32  n
      uint64  num_valid_records
      then for each record:
        int8  [C(n,3)]  chirotope signs
        uint16[n*2]     (x,y) coordinates
    """
    import numpy as np
    import struct as _struct

    MAGIC = 0x41414B01
    num_triples = n * (n - 1) * (n - 2) // 6

    data = chi_path.read_bytes()
    magic, nn, num_records = _struct.unpack_from("<IIQ", data, 0)
    if magic != MAGIC:
        raise ValueError(f"Bad magic in {chi_path}: {magic:#010x}")
    if nn != n:
        raise ValueError(f".chi file has n={nn}, expected {n}")

    header_size = 4 + 4 + 8   # 16 bytes
    chi_bytes = num_triples
    pts_bytes = n * 2 * 2     # uint16 * n * 2
    row_bytes = chi_bytes + pts_bytes

    body = np.frombuffer(data, dtype=np.uint8, offset=header_size)
    expected = num_records * row_bytes
    if len(body) < expected:
        raise ValueError(f".chi file truncated: got {len(body)} body bytes, need {expected}")

    body = body[:expected].reshape(num_records, row_bytes)

    chi_raw = body[:, :chi_bytes].view(np.int8)                          # (N, C(n,3))
    pts_raw = body[:, chi_bytes:].view(np.uint16).reshape(num_records, n, 2)  # (N, n, 2)

    print(f"    Loaded {num_records:,} records from {chi_path.name}", flush=True)
    return list(zip(chi_raw.tolist(), pts_raw.tolist()))


def _parse_aak_file(path, n, coord_bytes=2):
    """Parse an AAK binary database file.

    Uses a pre-built .chi cache produced by the C parser (pipeline/c/aak_parse).
    If the cache does not exist, builds it first (one-time cost, ~seconds).
    """
    chi = _chi_path(path)
    if not chi.exists():
        _build_chi_file(path, n, coord_bytes)
    return _parse_chi_file(chi, n)


def fetch_chi_records(chi_path, n, record_ids):
    """Return point coordinates for specific record IDs from a .chi file.

    Uses random access (fixed-width rows) — O(1) per record.
    Returns dict: record_id -> list of (x, y) tuples.
    """
    import numpy as np
    import struct as _struct

    MAGIC = 0x41414B01
    num_triples = n * (n - 1) * (n - 2) // 6
    chi_bytes = num_triples
    pts_bytes = n * 2 * 2
    row_bytes = chi_bytes + pts_bytes
    header_size = 16

    data = Path(chi_path).read_bytes()
    magic, nn, num_records = _struct.unpack_from("<IIQ", data, 0)
    if magic != MAGIC:
        raise ValueError(f"Bad magic: {magic:#010x}")
    if nn != n:
        raise ValueError(f".chi has n={nn}, expected {n}")

    body = np.frombuffer(data, dtype=np.uint8, offset=header_size)
    body = body[:num_records * row_bytes].reshape(num_records, row_bytes)
    pts_raw = body[:, chi_bytes:].view(np.uint16).reshape(num_records, n, 2)

    return {rid: [tuple(map(int, p)) for p in pts_raw[rid]] for rid in record_ids}


def stream_chi_file(chi_path, n, chunk_size=50_000):
    """Yield (chi_chunk, pts_chunk, start, end, total) from a .chi file.

    chi_chunk : int8  array  (M, C(n,3))
    pts_chunk : uint16 array (M, n, 2)
    start/end : record indices for this chunk
    total     : total records in file

    Avoids loading all 14M records into Python objects at once.
    """
    import numpy as np
    import struct as _struct

    MAGIC = 0x41414B01
    num_triples = n * (n - 1) * (n - 2) // 6
    chi_bytes = num_triples
    pts_bytes = n * 2 * 2
    row_bytes = chi_bytes + pts_bytes

    data = chi_path.read_bytes()
    magic, nn, num_records = _struct.unpack_from("<IIQ", data, 0)
    if magic != MAGIC:
        raise ValueError(f"Bad magic in {chi_path}: {magic:#010x}")
    if nn != n:
        raise ValueError(f".chi file has n={nn}, expected {n}")

    body = np.frombuffer(data, dtype=np.uint8, offset=16)
    body = body[:num_records * row_bytes].reshape(num_records, row_bytes)

    for start in range(0, num_records, chunk_size):
        end = min(start + chunk_size, num_records)
        chunk = body[start:end]
        chi_chunk = chunk[:, :chi_bytes].view(np.int8).copy()
        pts_chunk = chunk[:, chi_bytes:].view(np.uint16).reshape(end - start, n, 2).copy()
        yield chi_chunk, pts_chunk, start, end, num_records


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


def _generate_small_exhaustive(n, coord_range=100, num_attempts=None, seed=0):
    # Default sample count scales down for larger n (Stage 2 cost grows with n)
    if num_attempts is None:
        num_attempts = max(500, 5000 - (n - 8) * 1500)
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
