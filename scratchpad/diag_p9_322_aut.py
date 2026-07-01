"""Is P9_322's 18 an under-dedup? Compute Aut(type) (facet perms preserving the missing-face
set) and re-canonicalize the raw solutions under the FULL automorphism group; compare to
canonical_key's 18 and the truth 3."""
import json, itertools
import numpy as np
import sys
sys.path.insert(0, ".")
from run_d4 import canonical_key, config_to_gram

r = json.load(open("runs/d5_n9/blockpaste/type_322.json"))
polys = r["polytopes"]
n = 9
mf = None
# rebuild missing faces from the type (same construction as the anchor)
from pipeline.validate_coverage import parse_vertex_flags, minimal_non_faces, _detect_base
txt = open("scratchpad/HCPdm/polytopeDATA/5d9m.txt").read()
verts = parse_vertex_flags(txt.splitlines()[321], _detect_base(txt))
mf = [frozenset(m) for m in minimal_non_faces(verts, n)]
mf_set = frozenset(mf)

# automorphism group: permutations of facets preserving the missing-face set
auts = []
for perm in itertools.permutations(range(n)):
    if frozenset(frozenset(perm[i] for i in m) for m in mf) == mf_set:
        auts.append(perm)
print(f"P9_322: {len(polys)} raw, |Aut(type)| = {len(auts)}")
print(f"canonical_key distinct = {len(set(canonical_key(p, n) for p in polys))}")

# re-dedup under the full Aut group: canonical form = min over sigma in Aut of the
# permuted, rounded Gram flattening
def ek(v):
    return round(v, 5) if abs(v) <= 1 + 1e-9 else round(v, 4)

def aut_key(p):
    G = config_to_gram(p, n)
    best = None
    for s in auts:
        flat = tuple(ek(G[s[i], s[j]]) for i in range(n) for j in range(i + 1, n))
        if best is None or flat < best:
            best = flat
    return best

print(f"Aut-group distinct  = {len(set(aut_key(p) for p in polys))}  (truth = 3)")
