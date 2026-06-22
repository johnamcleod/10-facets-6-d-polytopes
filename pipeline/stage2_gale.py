"""
Stage 2: Affine Gale diagrams → combinatorial types of simple d-polytopes.

For each order type from Stage 1:
  - Choose pairs {u,v} ⊆ {0,...,n-1} to be the "positive" points
  - Check: both positive points lie strictly inside conv(negative points)
  - Compute the missing-face list (combinatorial type)
  - Deduplicate by canonical missing-face hypergraph
  - Filter to simple polytopes (or bucket separately)

The deduplication key is the canonical form of the missing-face hypergraph,
computed via a graph-automorphism based canonicalization (or by trying all
n! relabelings for small n — feasible for n ≤ 10).
"""

import json
import time
from pathlib import Path
from itertools import combinations
from functools import lru_cache

import numpy as np

from pipeline.utils.gale_exact import AffineGale, _strictly_inside
from pipeline.utils.canonical import canonical_missing_face_hypergraph
from pipeline.utils.manifest import write_manifest


# ---------------------------------------------------------------------------
# Vectorised pre-filter for streaming Stage 2
# ---------------------------------------------------------------------------

def _build_prefilter_tables(n):
    """Precompute chirotope index/sign lookup for the batch interior test.

    For each C(n,2) positive pair (u,v) and candidate interior point p ∈ {u,v},
    we need chi(p, s_i, s_j) for all s_i, s_j ∈ S = {0..n-1}\\{u,v}, i ≠ j.

    Interior test:  p is NOT strictly inside conv(S)
                    iff ∃ s_i ∈ S : ∀ s_j ∈ S\\{s_i} : chi(p, s_i, s_j) = +1.

    Returns list of dicts, one per (u,v) pair, with numpy arrays:
      idx_u  / sign_u : shape (k, k-1)  for checking u inside conv(S)
      idx_v  / sign_v : shape (k, k-1)  for checking v inside conv(S)
    where k = n - 2.
    """
    triple_to_idx = {abc: t for t, abc in enumerate(combinations(range(n), 3))}

    def chi_lookup(p, q, r):
        a, b, c = sorted((p, q, r))
        cidx = triple_to_idx[(a, b, c)]
        perm = [p, q, r]
        target = [a, b, c]
        swaps = 0
        for i in range(3):
            if perm[i] != target[i]:
                j = perm.index(target[i], i)
                perm[i], perm[j] = perm[j], perm[i]
                swaps += 1
        return cidx, (1 if swaps % 2 == 0 else -1)

    k = n - 2
    tables = []
    for u, v in combinations(range(n), 2):
        S = [s for s in range(n) if s != u and s != v]

        def build(p):
            idx_arr  = np.empty((k, k - 1), dtype=np.int32)
            sign_arr = np.empty((k, k - 1), dtype=np.int8)
            for i, si in enumerate(S):
                jj = 0
                for j, sj in enumerate(S):
                    if j == i:
                        continue
                    cidx, csign = chi_lookup(p, si, sj)
                    idx_arr[i, jj]  = cidx
                    sign_arr[i, jj] = csign
                    jj += 1
            return idx_arr, sign_arr

        tables.append({
            'u': u, 'v': v,
            'idx_u': build(u)[0], 'sign_u': build(u)[1],
            'idx_v': build(v)[0], 'sign_v': build(v)[1],
        })
    return tables


def _batch_prefilter(chi_chunk, tables):
    """Boolean mask: True if the record might yield a valid (u,v) Gale pair.

    A record passes iff at least one (u,v) pair has both u and v strictly
    inside conv(S), tested via the chirotope interior criterion.
    """
    M = len(chi_chunk)
    any_valid = np.zeros(M, dtype=bool)

    for t in tables:
        # vals_u[m, i, j] = chi(u, S[i], S[j])  for record m
        vals_u = chi_chunk[:, t['idx_u']] * t['sign_u']   # (M, k, k-1)
        u_outside = np.any(np.all(vals_u == 1, axis=2), axis=1)

        vals_v = chi_chunk[:, t['idx_v']] * t['sign_v']
        v_outside = np.any(np.all(vals_v == 1, axis=2), axis=1)

        any_valid |= (~u_outside) & (~v_outside)
        if any_valid.all():
            break

    return any_valid


def _degree_sequence(mf_list, n):
    """Compute the degree sequence of each node in the missing-face hypergraph.

    Returns a tuple of sorted (degree, node) pairs used as a fast invariant.
    """
    deg = [0] * n
    for m in mf_list:
        for i in m:
            deg[i] += 1
    return tuple(sorted(deg, reverse=True))


def _canonicalize_missing_faces(mf_list, n):
    """Fast canonical form of the missing-face hypergraph.

    Uses degree-sequence-guided relabeling: assign new labels in order of
    decreasing degree (ties broken by local neighborhood structure).
    This is a heuristic — good enough for a fast dedup pass.
    Exact dedup (for final runs) uses _canonical_form_exact.
    """
    mf_sorted = tuple(sorted(tuple(sorted(m)) for m in mf_list))
    # Augment with degree sequence for better discrimination
    deg_seq = _degree_sequence(mf_list, n)
    return (deg_seq, mf_sorted)


def _canonical_form_exact(mf_list, n):
    """Exact canonical form via exhaustive relabeling (O(n!), n <= 10).

    Returns the lexicographically minimal representation over all
    permutations of {0,...,n-1}.

    For n=10: 10! = 3.6M permutations. With |mf| typically small (~5-15),
    this costs ~50M comparisons — feasible for a few hundred types but
    not for millions. Called only when exact_dedup=True.
    """
    from itertools import permutations as perms
    best = None
    for perm in perms(range(n)):
        relabeled = tuple(sorted(
            tuple(sorted(perm[i] for i in m))
            for m in mf_list
        ))
        if best is None or relabeled < best:
            best = relabeled
    return best


def process_order_type(chi_key, points, d, exact_dedup=False):
    """Process one order type: try all positive-pair choices.

    Returns list of (canonical_type_key, missing_faces, AffineGale) for valid
    choices.  Uses the exact integer-arithmetic affine-Gale criterion
    (`pipeline.utils.gale_exact.AffineGale`), validated against the d=4 (30/30)
    and d=5 (109/109 of the k>=2 types) Ma-Zheng ground truth on 2026-06-22.
    The earlier float (`gale.GaleDiagram`) / C (`c/stage2_filter`) criteria were
    both wrong (see memory `stage2-generator-bugs`).
    """
    n = len(points)
    pts = [tuple(p) for p in points]
    results = []

    for u, v in combinations(range(n), 2):
        neg_pts = [pts[i] for i in range(n) if i != u and i != v]

        # Fast exact pre-filter: both positive points strictly inside conv(neg).
        if not (_strictly_inside(pts[u], neg_pts) and _strictly_inside(pts[v], neg_pts)):
            continue

        ag = AffineGale(pts, (u, v), d)

        # Polytopality gate: every singleton a face, simple (no (d+1)-face),
        # and f0 >= Lower-Bound-Theorem bound.  Removes non-polytopal configs
        # that otherwise produce spurious high-k "types".
        if not ag.is_polytope():
            continue

        # Missing faces up to size 5 (Lannér bound, filter F2).
        mf = ag.missing_faces(max_size=5)
        if sum(1 for m in mf if len(m) == 2) < 2:   # F1 pre-check: p >= 2
            continue

        key = canonical_missing_face_hypergraph(mf, n, exact=exact_dedup)
        results.append((key, mf, ag))

    return results


def _gd2_path(chi_path):
    return Path(chi_path).with_suffix(".gd2")


def _build_gd2_file(chi_path, d):
    """Run the C stage2_filter to produce a .gd2 file. Raises on failure."""
    import subprocess
    c_dir = Path(__file__).parent / "c"
    binary = c_dir / "stage2_filter"
    if not binary.exists():
        raise RuntimeError(f"C filter not built.  Run:  cd {c_dir} && make")
    gd2 = _gd2_path(chi_path)
    print(f"  Running C stage2_filter → {gd2.name} ...", flush=True)
    result = subprocess.run(
        [str(binary), str(chi_path), str(d), str(gd2)],
        capture_output=False,
    )
    if result.returncode != 0:
        raise RuntimeError(f"stage2_filter failed (exit {result.returncode})")
    return gd2


def _read_gd2_file(gd2_path, n):
    """Read a .gd2 file into a list of (record_id, u, v, missing_faces).

    missing_faces: list of frozensets of facet indices.
    """
    import struct as _struct
    MAGIC = 0x47443200
    data = gd2_path.read_bytes()
    magic, nn, num_entries = _struct.unpack_from("<IIQ", data, 0)
    if magic != MAGIC:
        raise ValueError(f"Bad .gd2 magic: {magic:#010x}")
    if nn != n:
        raise ValueError(f".gd2 has n={nn}, expected {n}")

    results = []
    pos = 16
    for _ in range(num_entries):
        rec_id, u, v, nm = _struct.unpack_from("<IBBB", data, pos)
        pos += 7
        masks = _struct.unpack_from(f"<{nm}H", data, pos)
        pos += nm * 2
        mf = [frozenset(i for i in range(n) if (m >> i) & 1) for m in masks]
        results.append((rec_id, u, v, mf))
    print(f"  Loaded {num_entries:,} unique raw patterns from {gd2_path.name}",
          flush=True)
    return results


def run_stage2_from_gd2(d, chi_path, output_dir=None, exact_dedup=False,
                        verbose=True, rebuild=False):
    """Run Stage 2 using the precomputed C .gd2 file.

    Builds the .gd2 file via the C stage2_filter if it does not already exist
    (or if rebuild=True).  Python then deduplicates using WL canonical form.
    """
    from pipeline.stage1_order_types import fetch_chi_records

    chi_path = Path(chi_path)
    n = d + 4
    gd2 = _gd2_path(chi_path)

    print(f"Stage 2: d={d}, n={n}  (C-accelerated from {chi_path.name})")
    if not exact_dedup and n >= 9:
        print(f"  WARNING: inexact dedup for n={n}. Use exact_dedup=True for final runs.")

    if rebuild or not gd2.exists():
        _build_gd2_file(chi_path, d)

    t0 = time.time()
    entries = _read_gd2_file(gd2, n)

    seen_types = {}
    for rec_id, u, v, mf in entries:
        key = canonical_missing_face_hypergraph(mf, n, exact=exact_dedup)
        p_count = sum(1 for m in mf if len(m) == 2)
        if key not in seen_types:
            seen_types[key] = {
                "type_id": len(seen_types),
                "canonical_key": list(key) if isinstance(key, tuple) else key,
                "missing_faces": [sorted(m) for m in mf],
                "p_count": p_count,
                "source_order_type_ids": [int(rec_id)],
                "example_positive": sorted([u, v]),
                "_example_rec_id": int(rec_id),
            }
        else:
            seen_types[key]["source_order_type_ids"].append(int(rec_id))

    # Fetch coordinates for one representative record per type
    rep_ids = [t["_example_rec_id"] for t in seen_types.values()]
    print(f"  Fetching example coordinates for {len(rep_ids)} types...", flush=True)
    coords = fetch_chi_records(chi_path, n, rep_ids)
    for t in seen_types.values():
        t["example_points"] = coords[t.pop("_example_rec_id")]

    types_list = list(seen_types.values())
    for i, t in enumerate(types_list):
        t["type_id"] = i

    elapsed = time.time() - t0
    print(f"  {len(entries):,} raw patterns → {len(types_list)} types  ({elapsed:.1f}s)")

    if verbose:
        p_counts = {}
        for t in types_list:
            p = t["p_count"]
            p_counts[p] = p_counts.get(p, 0) + 1
        for p in sorted(p_counts):
            print(f"    p={p}: {p_counts[p]} types")

    if output_dir is not None:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        (out / "types.json").write_text(
            json.dumps(types_list, indent=None, separators=(',', ':'))
        )
        write_manifest(out, "stage2",
                       {"d": d, "n": n, "exact_dedup": exact_dedup,
                        "source": "c_stage2_filter"},
                       {"num_types": len(types_list)})
        print(f"  Written to {out}/")

    return types_list


def run_stage2_from_chi(d, chi_path, output_dir=None, exact_dedup=False,
                        verbose=True, chunk_size=50_000):
    """Run Stage 2 streaming directly from a .chi file.

    Avoids materialising all 14M order types as Python objects.  A vectorised
    numpy pre-filter (chirotope interior test) rejects records in bulk; only
    survivors reach the Python per-record logic.
    """
    from pipeline.stage1_order_types import stream_chi_file

    chi_path = Path(chi_path)
    n = d + 4
    print(f"Stage 2: d={d}, n={n}  (streaming from {chi_path.name})")
    print(f"  chunk_size={chunk_size:,}  exact_dedup={exact_dedup}")
    if not exact_dedup and n >= 9:
        print(f"  WARNING: using fast (inexact) dedup for n={n}. "
              f"Use exact_dedup=True for final runs.")

    t0 = time.time()
    print("  Building pre-filter tables...", flush=True)
    tables = _build_prefilter_tables(n)
    print(f"  Pre-filter tables ready ({len(tables)} pairs).  Starting stream...",
          flush=True)

    seen_types = {}
    total_in = 0
    total_passed = 0

    for chi_chunk, pts_chunk, start, end, total_records in \
            stream_chi_file(chi_path, n, chunk_size):

        mask = _batch_prefilter(chi_chunk, tables)
        passed_idx = np.where(mask)[0]
        total_in     += end - start
        total_passed += len(passed_idx)

        for i in passed_idx:
            pts      = [tuple(map(int, p)) for p in pts_chunk[i]]
            chi_key  = tuple(int(x) for x in chi_chunk[i])
            ot_id    = start + int(i)

            for key, mf, gd in process_order_type(chi_key, pts, d, exact_dedup):
                if key not in seen_types:
                    seen_types[key] = {
                        "type_id": len(seen_types),
                        "canonical_key": list(key) if isinstance(key, tuple) else key,
                        "missing_faces": [sorted(m) for m in mf],
                        "p_count": sum(1 for m in mf if len(m) == 2),
                        "source_order_type_ids": [ot_id],
                        "example_points": [list(p) for p in pts],
                        "example_positive": sorted(gd.positive),
                        "vertex_sets": [sorted(v) for v in gd.vertex_sets()],
                    }
                else:
                    seen_types[key]["source_order_type_ids"].append(ot_id)

        if verbose and (end % (chunk_size * 20) < chunk_size or end == total_records):
            elapsed = time.time() - t0
            rate    = total_in / elapsed
            remain  = (total_records - total_in) / rate if rate > 0 else 0
            pct_pass = 100.0 * total_passed / total_in if total_in else 0
            print(f"  {total_in:>10,}/{total_records:,}  ({100*total_in//total_records}%)  "
                  f"pass-rate={pct_pass:.3f}%  types={len(seen_types)}  "
                  f"~{remain/60:.1f} min remaining", flush=True)

    types_list = list(seen_types.values())
    for i, t in enumerate(types_list):
        t["type_id"] = i

    elapsed = time.time() - t0
    print(f"  Done in {elapsed:.1f}s.  "
          f"filter kept {total_passed:,}/{total_in:,} ({100*total_passed/total_in:.3f}%)  "
          f"→ {len(types_list)} distinct types.")

    if verbose:
        p_counts = {}
        for t in types_list:
            p = t["p_count"]
            p_counts[p] = p_counts.get(p, 0) + 1
        for p in sorted(p_counts):
            print(f"    p={p}: {p_counts[p]} types")

    if output_dir is not None:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        (out / "types.json").write_text(
            json.dumps(types_list, indent=None, separators=(',', ':'))
        )
        write_manifest(out, "stage2", {"d": d, "n": n, "exact_dedup": exact_dedup},
                       {"num_types": len(types_list)})
        print(f"  Written to {out}/")

    return types_list


def run_stage2(d, stage1_results, output_dir=None, exact_dedup=False, verbose=True):
    """Run Stage 2: build Gale diagrams and extract combinatorial types.

    Args:
        d: polytope dimension
        stage1_results: list of order-type dicts from Stage 1
        exact_dedup: use full n!-relabeling canonical form (slow but exact)
        output_dir: where to write output

    Returns:
        list of combinatorial type dicts with keys:
          'type_id', 'canonical_key', 'missing_faces', 'p_count',
          'source_order_type_ids', 'example_points', 'example_positive'
    """
    n = d + 4
    print(f"Stage 2: d={d}, n={n}")
    print(f"  Processing {len(stage1_results)} order types...")
    print(f"  Exact dedup: {exact_dedup}")
    if not exact_dedup and n >= 9:
        print(f"  WARNING: using fast (inexact) dedup for n={n}. "
              f"Use exact_dedup=True for final runs.")

    seen_types = {}  # canonical_key -> type dict

    for ot in stage1_results:
        pts = [tuple(p) for p in ot["points"]]
        chi_key = tuple(ot["chirotope"])
        ot_id = ot["id"]

        for key, mf, gd in process_order_type(chi_key, pts, d, exact_dedup):
            if key not in seen_types:
                seen_types[key] = {
                    "type_id": len(seen_types),
                    "canonical_key": list(key) if isinstance(key, tuple) else key,
                    "missing_faces": [sorted(m) for m in mf],
                    "p_count": sum(1 for m in mf if len(m) == 2),
                    "source_order_type_ids": [ot_id],
                    "example_points": [list(p) for p in pts],
                    "example_positive": sorted(gd.positive),
                    "vertex_sets": [sorted(v) for v in gd.vertex_sets()],
                }
            else:
                seen_types[key]["source_order_type_ids"].append(ot_id)

    types_list = list(seen_types.values())
    # Reassign sequential IDs
    for i, t in enumerate(types_list):
        t["type_id"] = i

    if verbose:
        print(f"  Found {len(types_list)} distinct combinatorial types "
              f"(after p>=2 filter).")
        p_counts = {}
        for t in types_list:
            p = t["p_count"]
            p_counts[p] = p_counts.get(p, 0) + 1
        for p in sorted(p_counts):
            print(f"    p={p}: {p_counts[p]} types")

    if output_dir is not None:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        (out / "types.json").write_text(
            json.dumps(types_list, indent=None, separators=(',', ':'))
        )
        write_manifest(out, "stage2", {"d": d, "n": n, "exact_dedup": exact_dedup},
                       {"num_types": len(types_list)})
        print(f"  Written to {out}/")

    return types_list


def load_stage2(output_dir):
    """Load Stage 2 results from disk."""
    path = Path(output_dir) / "types.json"
    return json.loads(path.read_text())
