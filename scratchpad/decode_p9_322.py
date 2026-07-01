"""Decode Ma-Zheng's 3 P9_322 label-vectors (change7) and check whether they are among
our 18 solutions' ordinary-label assignments, up to the automorphism group."""
import json, ast, itertools, sys
import numpy as np
sys.path.insert(0, ".")
from pipeline.validate_coverage import parse_vertex_flags, minimal_non_faces, _detect_base

n = 9
txt = open("scratchpad/HCPdm/polytopeDATA/5d9m.txt").read()
verts = [frozenset(v) for v in parse_vertex_flags(txt.splitlines()[321], _detect_base(txt))]
mf = [tuple(sorted(m)) for m in minimal_non_faces(verts, n)]
dotted0 = [m for m in mf if len(m) == 2]                    # 0-indexed dotted pairs
mf_set = frozenset(frozenset(m) for m in mf)
# automorphism group (facet perms preserving missing faces)
auts = [p for p in itertools.permutations(range(n))
        if frozenset(frozenset(p[i] for i in m) for m in mf) == mf_set]

# their column order: 1-indexed pairs 10*i+j sorted, minus dotted(infty)
infty = {10 * (a + 1) + (b + 1) for (a, b) in dotted0}
key_all = sorted(10 * i + j for i in range(1, n + 1) for j in range(i + 1, n + 1))
key_sel = [k for k in key_all if k not in infty]            # 30 ordinary pair-keys
assert len(key_sel) == 30, len(key_sel)
cols = [((k // 10) - 1, (k % 10) - 1) for k in key_sel]      # 0-indexed ordinary pairs

def label_matrix(labels):
    L = np.full((n, n), 2, dtype=int)
    for (i, j), m in zip(cols, labels):
        L[i, j] = L[j, i] = int(m)
    for (i, j) in dotted0:
        L[i, j] = L[j, i] = -1                              # dotted marker
    return L

def canon(L):
    best = None
    for s in auts:
        f = tuple(int(L[s[i], s[j]]) for i in range(n) for j in range(i + 1, n))
        best = f if best is None or f < best else best
    return best

their = []
for line in open("scratchpad/HCPdm/output/P9_322/P9_322_LSIE_per_change7.txt"):
    labs = [int(x) for x in line.split()]
    if len(labs) == 30:
        their.append(canon(label_matrix(labs)))
print("their canonical label-matrices:", len(their))

# our 18 solutions' ordinary label matrices
r = json.load(open("runs/d5_n9/blockpaste/type_322.json"))
from run_d4 import canonical_key
seen = {}
for p in r["polytopes"]:
    seen.setdefault(canonical_key(p, n), p)
ours = set()
for p in seen.values():
    L = np.full((n, n), 2, dtype=int)
    for e, m in p["label_assignment"].items():
        i, j = ast.literal_eval(e); L[i, j] = L[j, i] = int(m)
    for key in p.get("dot_values", {}):
        i, j = (int(x) for x in key.split("_")[1:]); L[i, j] = L[j, i] = -1
    ours.add(canon(L))
print("our distinct label-matrices:", len(ours))
print("their 3 that appear among ours:", sum(1 for t in their if t in ours), "/", len(their))
