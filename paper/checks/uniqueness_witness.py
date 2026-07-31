#!/usr/bin/env python3
"""The data behind ``exactly one polytope up to isometry''.

The refutation certificates establish that exactly one LABELLING of one type is
realizable.  That is not by itself uniqueness of the polytope, for two reasons, and
this script measures both.

1. A labelling need not determine its ultraparallel weights.  Each pinning step of
   the cascade solves a quadratic and admits up to two roots > 1, so the weight tree
   of a single labelling can have several leaves, each a candidate weight vector,
   and several could in principle pass.  For the accepted labelling the tree turns
   out to have exactly one leaf.

2. Several labellings of the same type can describe the same polytope, namely the
   orbit of one labelling under the combinatorial automorphism group Aut(T) of the
   type.  The orbit size is |Aut(T)| / |Stab|, and the search returns exactly those
   labellings when symmetry breaking is off.  Since a compact hyperbolic Coxeter
   polytope is determined up to isometry of H^d by its Gram matrix up to
   simultaneous row and column permutation, one orbit means one polytope.

Both are computed here from the realizer record and the type's missing-face
hypergraph, and written to runs/d6_n10/uniqueness_witness.json for the paper to
read.  Aut(T) is computed by brute force over all 10! permutations, which takes a
few seconds and needs no graph-isomorphism library to be trusted.

Run:  python3 paper/checks/uniqueness_witness.py
"""
import ast
import itertools
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

import pipeline.stage4_gram as sg   # noqa: E402

D, N, TID = 6, 10, 379


def realizer():
    for run in ("d6_exact", "d6_final", "d6_tangency"):
        for name in (f"tid{TID}_root.json", f"tid{TID}_0_0.json"):
            p = ROOT / f"runs/d{D}_n{N}/{run}/realizers/{name}"
            if p.exists():
                return json.loads(p.read_text())[0], f"{run}/realizers/{name}"
    raise SystemExit("no realizer record found")


def main():
    rec, src = realizer()
    labels = {ast.literal_eval(k): v for k, v in rec["label_assignment"].items()}
    dashed = tuple(sorted(tuple(int(x) for x in k.split("_")[1:])
                          for k in rec["dot_values"]))

    # (1) the weight tree of the accepted labelling
    ordinary = {p: sg._GRAM_FLOAT[m] for p, m in labels.items()}
    minor_index = sg._build_minor_index(dashed, N, D)
    stats = {}
    decision, sols = sg._structured_screen(ordinary, dashed, minor_index, N, D,
                                           stats=stats)

    # (2) Aut(T) and the stabiliser of the labelling in it
    types = {t["type_id"]: t for t in
             json.loads((ROOT / f"runs/d{D}_n{N}/stage2/types.json").read_text())}
    mf = {tuple(sorted(m)) for m in types[TID]["missing_faces"]}

    def is_aut(perm):
        return {tuple(sorted(perm[i] for i in m)) for m in mf} == mf

    aut = [p for p in itertools.permutations(range(N)) if is_aut(p)]

    dashed_set = set(dashed)

    def lab(i, j):
        p = (i, j) if i < j else (j, i)
        return ("dashed",) if p in dashed_set else ("ordinary", labels[p])

    stab = [p for p in aut
            if all(lab(p[i], p[j]) == lab(i, j)
                   for i in range(N) for j in range(i + 1, N))]

    out = {
        "source": src,
        "type": TID,
        "cascade_decision": decision,
        "weight_tree_leaves": stats.get("leaves"),
        "weight_tree_leaves_passing": stats.get("leaf_solutions"),
        "candidate_weight_vectors": len(sols or []),
        "aut_type_order": len(aut),
        "labelling_stabiliser_order": len(stab),
        "orbit_of_the_labelling": len(aut) // len(stab),
        "weights": {f"x{k}": float(v) for k, v in (sols or [{}])[0].items()},
    }
    (ROOT / f"runs/d{D}_n{N}/uniqueness_witness.json").write_text(
        json.dumps(out, indent=1) + "\n")

    print(f"realizer record: {src}")
    print(f"weight tree of the accepted labelling: "
          f"{out['weight_tree_leaves']} leaf/leaves, "
          f"{out['weight_tree_leaves_passing']} passing the leaf test")
    print(f"|Aut(T)| = {out['aut_type_order']}, "
          f"|Stab(labelling)| = {out['labelling_stabiliser_order']}, "
          f"orbit = {out['orbit_of_the_labelling']} label assignments")
    print(f"artifact -> runs/d{D}_n{N}/uniqueness_witness.json")

    problems = []
    if out["weight_tree_leaves"] != 1:
        problems.append(f"the weight tree has {out['weight_tree_leaves']} leaves; "
                        f"the paper's uniqueness paragraph must say so")
    if out["weight_tree_leaves_passing"] != 1:
        problems.append(f"{out['weight_tree_leaves_passing']} leaves pass the leaf "
                        f"test, so one labelling presents several candidate "
                        f"polytopes")
    if out["aut_type_order"] % out["labelling_stabiliser_order"]:
        problems.append("the stabiliser order does not divide |Aut(T)|")
    print()
    print("RESULT:", "the accepted labelling determines its weights uniquely, and "
          "its Aut(T)-orbit is one polytope presented "
          f"{out['orbit_of_the_labelling']} ways"
          if not problems else "PROBLEMS: " + "; ".join(problems))
    return 0 if not problems else 1


if __name__ == "__main__":
    sys.exit(main())
