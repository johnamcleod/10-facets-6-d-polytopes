"""
Coverage validator: compare the pipeline's Stage-3 combinatorial types against
published ground truth (Ma-Zheng HCPdm) up to facet relabeling.

Ground truth files (data/ground_truth/{4d8m,5d9m}.txt) list, one polytope per
line, the VERTEX FLAGS: each [...] is the set of facets incident to one vertex,
facets 1-indexed.  From these we derive the MINIMAL NON-FACE system (a.k.a.
missing faces): the minimal sets of facets with no common vertex.  That system
is exactly what Stage 2/3 must reproduce (Stage 4 consumes it directly), so it
is the right invariant to compare on.

Two combinatorial types are equal iff there is a permutation of the n facets
carrying one minimal-non-face system to the other.  We test this exactly:
color-refine the facets from the non-face hypergraph, then search only
color-preserving permutations (a tiny fraction of n! once colors split).

Usage:
    python -m pipeline.validate_coverage --d 4
    python -m pipeline.validate_coverage --d 4 --stage3 runs/d4_n8/stage3 \
        --truth data/ground_truth/4d8m.txt
"""

import argparse
import itertools
import json
from collections import Counter, defaultdict
from pathlib import Path


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def _detect_base(text):
    """Detect facet-label base (0- or 1-indexed) from the min label in the file.

    The HCPdm ground-truth files are inconsistent: 4d8m.txt is 1-indexed
    (labels 1..8) while 5d9m.txt is 0-indexed (labels 0..8).  We normalise to
    0-indexed by subtracting the detected base.
    """
    mn = None
    for tok in text.replace("[", " ").replace("]", " ").replace(",", " ").split():
        try:
            v = int(tok)
        except ValueError:
            continue
        mn = v if mn is None else min(mn, v)
    return mn if mn in (0, 1) else 0


def parse_vertex_flags(line, base=1):
    """Parse one ground-truth line into a list of frozensets (0-indexed facets)."""
    verts = []
    for tok in line.replace("[", "|").replace("]", "").split("|"):
        tok = tok.strip()
        if not tok:
            continue
        verts.append(frozenset(int(x) - base for x in tok.split(",")))
    return verts


def minimal_non_faces(vertices, n):
    """Minimal subsets S of {0..n-1} not contained in any vertex's facet set.

    A subset S is a 'face' iff S is contained in some vertex (the facets in S
    share that vertex).  We enumerate candidate non-faces by increasing size
    and keep only the minimal ones.
    """
    vertices = [set(v) for v in vertices]

    def is_face(S):
        return any(S <= v for v in vertices)

    minimal = []
    minimal_sets = []
    # Sizes 2..n: a singleton {i} is always a face (facet i exists).
    for size in range(2, n + 1):
        for combo in itertools.combinations(range(n), size):
            cs = set(combo)
            if is_face(cs):
                continue
            # minimal iff no already-found minimal non-face is a subset
            if any(m <= cs for m in minimal_sets):
                continue
            minimal.append(frozenset(combo))
            minimal_sets.append(cs)
    return minimal


def load_ground_truth(path, n):
    types = []
    text = Path(path).read_text()
    base = _detect_base(text)
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        verts = parse_vertex_flags(line, base)
        types.append({
            "nverts": len(verts),
            "mnf": minimal_non_faces(verts, n),
            "vertices": [tuple(sorted(v)) for v in verts],
        })
    return types


def load_truth_as_stage3(path, d):
    """Load ground-truth combinatorial types as Stage-4-ready type dicts.

    Bypasses the order-type -> Gale-diagram generator: the vertex flags ARE the
    vertex sets, and the minimal-non-face system IS the missing-face list.
    Each returned dict has the keys process_type_stage4 needs:
      type_id, missing_faces, vertex_sets, p_count.
    """
    n = d + 4
    types = []
    text = Path(path).read_text()
    base = _detect_base(text)
    for i, line in enumerate(text.splitlines()):
        line = line.strip()
        if not line:
            continue
        verts = parse_vertex_flags(line, base)
        # Sanity: simple d-polytope -> every vertex on exactly d facets.
        if any(len(v) != d for v in verts):
            raise ValueError(f"line {i}: vertex not incident to exactly d={d} facets")
        mnf = minimal_non_faces(verts, n)
        types.append({
            "type_id": i,
            "missing_faces": [sorted(m) for m in mnf],
            "vertex_sets": [sorted(v) for v in verts],
            "p_count": sum(1 for m in mnf if len(m) == 2),
        })
    return types


def load_our_types(stage_dir):
    """Load our combinatorial types from a Stage-2 (types.json) or Stage-3
    (surviving_types.json) directory.

    NOTE: the ground-truth files are the full set of simple (d,d+4)-polytope
    combinatorial candidates (4d8m is pre-filtered to p>=2; 5d9m is not), NOT
    the realized-compact subset.  So the GENERATOR gate is Stage-2 coverage;
    Stage 3 then legitimately prunes candidates toward realizability and is
    expected to drop some.  Point this at the Stage-2 dir to test the generator.
    """
    p2 = Path(stage_dir) / "types.json"
    p3 = Path(stage_dir) / "surviving_types.json"
    raw = json.loads((p3 if p3.exists() else p2).read_text())
    return [{"type_id": t["type_id"],
             "mnf": [frozenset(m) for m in t["missing_faces"]]} for t in raw]


# ---------------------------------------------------------------------------
# Exact isomorphism of minimal-non-face systems up to facet relabeling
# ---------------------------------------------------------------------------

def _refine_colors(mnf, n):
    """Color-refine facets from the minimal-non-face hypergraph.

    Initial color = multiset of sizes of non-faces a facet belongs to.
    Refine by appending the sorted multiset of co-members' colors.
    Returns a list `color[i]` of small ints, stable under relabeling.
    """
    incidence = [[] for _ in range(n)]   # incidence[i] = list of frozensets containing i
    for f in mnf:
        for i in f:
            incidence[i].append(f)

    colors = [tuple(sorted(len(f) for f in incidence[i])) for i in range(n)]
    colors = _renumber(colors)
    for _ in range(n):
        new = []
        for i in range(n):
            sig = (colors[i],)
            nbr = []
            for f in incidence[i]:
                nbr.append((len(f), tuple(sorted(colors[j] for j in f if j != i))))
            new.append((sig, tuple(sorted(nbr))))
        new = _renumber(new)
        if new == colors:
            break
        colors = new
    return colors


def _renumber(seq):
    order = {c: idx for idx, c in enumerate(sorted(set(seq)))}
    return [order[c] for c in seq]


def canonical_form(mnf, n, perm_cap=2_000_000):
    """Lexicographically minimal representation of `mnf` over color-preserving
    permutations of the n facets.  Exact canonical form (two systems are
    isomorphic iff their canonical forms are equal), provided the color-class
    product stays under `perm_cap`; otherwise raises to signal the caller to
    fall back to the coarse cert.
    """
    colors = _refine_colors(mnf, n)
    classes = defaultdict(list)
    for i, c in enumerate(colors):
        classes[c].append(i)
    class_keys = sorted(classes.keys())
    groups = [classes[k] for k in class_keys]

    total = 1
    for g in groups:
        total *= _factorial(len(g))
    if total > perm_cap:
        raise OverflowError(f"color-preserving perms {total} exceed cap {perm_cap}")

    best = None
    for perm_combo in itertools.product(*[itertools.permutations(g) for g in groups]):
        # Build relabeling: original facet -> new index.
        # Concatenate the permuted classes in fixed class order to assign 0..n-1.
        mapping = {}
        idx = 0
        for grp in perm_combo:
            for orig in grp:
                mapping[orig] = idx
                idx += 1
        relabeled = tuple(sorted(
            tuple(sorted(mapping[x] for x in f)) for f in mnf
        ))
        if best is None or relabeled < best:
            best = relabeled
    return best


def coarse_cert(mnf, n):
    """Iso-invariant (necessary, not sufficient) key for bucketing / fallback."""
    colors = _refine_colors(mnf, n)
    by_size = Counter(len(f) for f in mnf)
    edge_colors = tuple(sorted(
        tuple(sorted(colors[x] for x in f)) for f in mnf
    ))
    return (tuple(sorted(by_size.items())), tuple(sorted(colors)), edge_colors)


def _factorial(k):
    r = 1
    for i in range(2, k + 1):
        r *= i
    return r


def type_key(mnf, n):
    """Best available isomorphism key: exact canonical form, else coarse cert."""
    try:
        return ("exact", canonical_form(mnf, n))
    except OverflowError:
        return ("coarse", coarse_cert(mnf, n))


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------

def k_of(mnf):
    """k = number of size-2 minimal non-faces (= disjoint facet pairs)."""
    return sum(1 for f in mnf if len(f) == 2)


def compare(d, stage3_dir, truth_path, verbose=True):
    n = d + 4
    truth = load_ground_truth(truth_path, n)
    ours = load_our_types(stage3_dir)

    truth_keys = defaultdict(list)
    for t in truth:
        truth_keys[type_key(t["mnf"], n)].append(t)
    our_keys = defaultdict(list)
    for t in ours:
        our_keys[type_key(t["mnf"], n)].append(t)

    truth_set = set(truth_keys)
    our_set = set(our_keys)

    covered = truth_set & our_set
    missing = truth_set - our_set      # ground-truth types we fail to produce
    spurious = our_set - truth_set     # our types with no ground-truth match

    if verbose:
        print(f"\n{'='*64}")
        print(f"Coverage validation  d={d}, n={n}")
        print(f"  ground truth: {truth_path}")
        print(f"  our stage3:   {stage3_dir}")
        print(f"{'='*64}")
        print(f"Ground-truth combinatorial types:  {len(truth)} lines, "
              f"{len(truth_set)} distinct (up to relabeling)")
        print(f"Our Stage-3 survivor types:        {len(ours)} lines, "
              f"{len(our_set)} distinct")
        print(f"\n  COVERED  (truth types we produce):   {len(covered)}/{len(truth_set)}")
        print(f"  MISSING  (truth types we DROP):      {len(missing)}")
        print(f"  SPURIOUS (our types not in truth):   {len(spurious)}")

        def kdist(keys, bucket):
            c = Counter()
            for key in keys:
                reps = bucket[key]
                c[k_of(reps[0]["mnf"])] += 1
            return sorted(c.items())

        print(f"\n  k-distribution (k = # disjoint facet pairs):")
        print(f"    truth  all:      {kdist(truth_set, truth_keys)}")
        print(f"    truth  MISSING:  {kdist(missing, truth_keys)}")
        print(f"    ours   all:      {kdist(our_set, our_keys)}")
        print(f"    ours   SPURIOUS: {kdist(spurious, our_keys)}")

        approx = any(key[0] == "coarse" for key in truth_set | our_set)
        if approx:
            print("\n  NOTE: some types used the coarse (necessary-only) cert; "
                  "a 'covered' on those is not a guaranteed exact match.")

        if missing:
            print(f"\n  Example MISSING ground-truth types (first 5):")
            for key in list(missing)[:5]:
                rep = truth_keys[key][0]
                size2 = sorted(tuple(sorted(f)) for f in rep["mnf"] if len(f) == 2)
                print(f"    k={k_of(rep['mnf'])}  nverts={rep['nverts']}  "
                      f"disjoint_pairs={size2}")

    return {
        "d": d, "n": n,
        "truth_lines": len(truth), "truth_distinct": len(truth_set),
        "ours_lines": len(ours), "ours_distinct": len(our_set),
        "covered": len(covered), "missing": len(missing), "spurious": len(spurious),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--d", type=int, required=True, choices=[4, 5])
    ap.add_argument("--stage3", default=None, help="Stage-3 dir (default runs/d{d}_n{d+4}/stage3)")
    ap.add_argument("--truth", default=None, help="ground-truth file (default data/ground_truth/{d}d{d+4}m.txt)")
    args = ap.parse_args()
    n = args.d + 4
    stage3 = args.stage3 or f"runs/d{args.d}_n{n}/stage3"
    truth = args.truth or f"data/ground_truth/{args.d}d{n}m.txt"
    compare(args.d, stage3, truth)


if __name__ == "__main__":
    main()
