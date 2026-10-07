#!/usr/bin/env python3
"""Does the truncated (size <= 5) missing-face hypergraph determine the
full one?  Exact, combinatorial, per type.

For a simple 6-polytope with 10 facets, a set of <= 5 facets is a face iff it
contains no minimal non-face of size <= 5, so the truncated hypergraph K fixes all
faces of size <= 5 (labelled).  Vertices are the 6-sets that are faces; every one
of them has all 5-subsets faces, so vertices lie in C(K) = {6-sets all of whose
5-subsets are faces}.  Every edge (5-set face) of a simple polytope lies in exactly
two vertices.  We enumerate ALL subsets F of C(K) satisfying that, i.e. all
0/1 solutions of  sum_{S in C(K), S > R} x_S = 2  for every 5-face R.  The size-6
minimal non-faces are then exactly C(K) \\ F.  If the solution is unique for K, then
EVERY candidate (any order type, any positive pair) with truncated hypergraph equal
to K up to relabelling has full hypergraph equal to K + (C(K)\\F) up to the same
relabelling: the truncated class is a full class.

Used in paper section 3 (Proposition 3.2): deduplication on the size <= 5
hypergraph merges no two distinct types.

Run:  python3 checks/size6_completion.py [out.json]
"""
import itertools, json, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from pipeline.utils.gale_exact import AffineGale  # noqa: E402

N, D = 10, 6


def faces_upto5(K):
    Kf = [frozenset(m) for m in K]
    def is_face(S):
        return not any(m <= S for m in Kf)
    return is_face


def solve(K, cap=1000):
    is_face = faces_upto5(K)
    ridges = [frozenset(R) for R in itertools.combinations(range(N), 5)
              if is_face(frozenset(R))]
    C = [frozenset(S) for S in itertools.combinations(range(N), 6)
         if all(is_face(frozenset(R)) for R in itertools.combinations(S, 5))]
    inc = {R: [i for i, S in enumerate(C) if R <= S] for R in ridges}
    sols = []
    x = [None] * len(C)

    def dfs():
        if len(sols) >= cap:
            return
        # propagate
        changed = True
        trail = []
        ok = True
        while changed and ok:
            changed = False
            for R, lst in inc.items():
                ones = sum(1 for i in lst if x[i] == 1)
                und = [i for i in lst if x[i] is None]
                if ones > 2 or ones + len(und) < 2:
                    ok = False
                    break
                if und and ones == 2:
                    for i in und:
                        x[i] = 0; trail.append(i)
                    changed = True
                elif und and ones + len(und) == 2:
                    for i in und:
                        x[i] = 1; trail.append(i)
                    changed = True
        if ok:
            und = [i for i in range(len(C)) if x[i] is None]
            if not und:
                sols.append(frozenset(C[i] for i in range(len(C)) if x[i] == 1))
            else:
                i = und[0]
                for val in (1, 0):
                    x[i] = val
                    dfs()
                    x[i] = None
        for i in trail:
            x[i] = None

    dfs()
    return C, sols


def main():
    types = json.load(open(ROOT / "runs/d6_n10/stage2/types.json"))
    flags = json.load(open(ROOT / "runs/d6_n10/type_flags.json"))["d6_n10"]
    t0 = time.time()
    out = {}
    bad = []
    for t in types:
        tid = t["type_id"]
        K = t["missing_faces"]
        assert max(len(m) for m in K) <= 5
        C, sols = solve(K)
        if len(sols) != 1:
            bad.append((tid, len(sols)))
            out[tid] = {"solutions": len(sols)}
            continue
        F = sols[0]
        mf6 = sorted(sorted(S) for S in C if S not in F)
        # cross-check against the exact Gale criterion on the stored realisation
        ag = AffineGale([tuple(p) for p in t["example_points"]],
                        t["example_positive"], D)
        full = ag.missing_faces(max_size=N)
        g6 = sorted(sorted(m) for m in full if len(m) >= 6)
        gtr = sorted(sorted(m) for m in full if len(m) <= 5)
        agree_trunc = gtr == sorted(sorted(m) for m in K)
        agree6 = g6 == mf6
        verts = sorted(sorted(S) for S in F)
        agree_v = verts == sorted(sorted(v) for v in ag.vertex_sets())
        flag6 = flags[str(tid)]["has_missing_face_of_size_d"]
        out[tid] = {"solutions": 1, "C": len(C), "vertices": len(F),
                    "size6_missing": len(mf6), "flag_size6": flag6,
                    "gale_agrees_trunc": agree_trunc, "gale_agrees_size6": agree6,
                    "gale_agrees_vertices": agree_v,
                    "max_mf_size_gale": max(len(m) for m in full)}
        if not (agree_trunc and agree6 and agree_v and (len(mf6) > 0) == flag6):
            bad.append((tid, "mismatch", out[tid]))
    n83 = sum(1 for v in out.values() if v.get("size6_missing"))
    res = {"types": len(types), "unique_completion": sum(1 for v in out.values()
           if v["solutions"] == 1), "with_size6": n83, "problems": bad,
           "seconds": round(time.time() - t0, 1), "per_type": out}
    (Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "runs/d6_n10/size6_completion.json").write_text(
        json.dumps(res, indent=1, default=str))
    print({k: v for k, v in res.items() if k != "per_type"})


if __name__ == "__main__":
    main()
