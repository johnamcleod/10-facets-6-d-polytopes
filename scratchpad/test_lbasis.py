"""Empirical test: does enabling the l_basis prism saver reconcile the counts?
  d=5 P9_322 : expect 3 (Ma-Zheng)   [we get 18 without l_basis]
  d=4 type 6 : expect 49 (census P2)  [we get 49 without l_basis; 12 with it, earlier]
CoxIter-check the survivors."""
import json, sys, ast
from pathlib import Path
sys.path.insert(0, ".")
import pipeline.stage4_blockpaste as bp
bp.USE_L4_BASIS = True            # enable the prism-basis saver
from run_blockpaste_parallel import process_type
from pipeline.utils import coxiter
from pipeline.validate_coverage import parse_vertex_flags, minimal_non_faces, _detect_base

def coxiter_compact(res, n, d):
    ok = 0
    from run_d4 import canonical_key
    seen = {}
    for p in res["polytopes"]:
        seen.setdefault(canonical_key(p, n), p)
    for p in seen.values():
        la = {tuple(ast.literal_eval(e)): m for e, m in p["label_assignment"].items()}
        dotted = [tuple(int(x) for x in k.split("_")[1:]) for k in p.get("dot_values", {})]
        cc, dim = coxiter.check(la, dotted, n, d)
        if cc and dim == d:
            ok += 1
    return len(seen), ok

out = Path("runs/tmp_lbasis"); out.mkdir(parents=True, exist_ok=True)

# d=5 P9_322
txt = open("scratchpad/HCPdm/polytopeDATA/5d9m.txt").read()
verts = parse_vertex_flags(txt.splitlines()[321], _detect_base(txt))
n = 1 + max(f for v in verts for f in v)
t = {"type_id": 322, "vertex_sets": [sorted(v) for v in verts],
     "missing_faces": [sorted(m) for m in minimal_non_faces(verts, n)]}
r5 = process_type(322, {322: t}, nproc=7, outdir=out, p=2)
d5tot, d5cc = coxiter_compact(r5, 9, 5)
print(f"\nP9_322 with l_basis: {d5tot} distinct, {d5cc} CoxIter-compact  (Ma-Zheng 3)", flush=True)

# d=4 type 6
surv = {tt["type_id"]: tt for tt in json.load(open("runs/d4_n8/stage2/types.json"))}
r4 = process_type(6, surv, nproc=7, outdir=out, p=1)
d4tot, d4cc = coxiter_compact(r4, 8, 4)
print(f"type 6 with l_basis: {d4tot} distinct, {d4cc} CoxIter-compact  (census P2 49)", flush=True)
