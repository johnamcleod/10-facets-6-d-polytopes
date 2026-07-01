"""d=5 anchor: block-paste on Ma-Zheng P9_322 must give 3 (their 180 candidates -> 3
after signature).  If it gives >3, the block-paste/solver over-accepts for d=5."""
import sys, json, time
from pathlib import Path
sys.path.insert(0, ".")
from pipeline.validate_coverage import parse_vertex_flags, minimal_non_faces, _detect_base
from run_blockpaste_parallel import process_type

txt = open("scratchpad/HCPdm/polytopeDATA/5d9m.txt").read()
base = _detect_base(txt)
line = txt.splitlines()[321]          # P9_322 = line 322 (1-indexed) = index 321
verts = parse_vertex_flags(line, base)
n = 1 + max(f for v in verts for f in v)
mf = minimal_non_faces(verts, n)
t = {"type_id": 322,
     "vertex_sets": [sorted(v) for v in verts],
     "missing_faces": [sorted(m) for m in mf]}
k = sum(1 for m in mf if len(m) == 2)
print(f"P9_322: n={n} nverts={len(verts)} k(dotted)={k} "
      f"mf-sizes={sorted(len(m) for m in mf)}", flush=True)

t0 = time.time()
outdir = Path("runs/d5_n9/blockpaste"); outdir.mkdir(parents=True, exist_ok=True)
res = process_type(322, {322: t}, nproc=7, outdir=outdir, p=2)
print(f"\nP9_322 -> {res['distinct']} DISTINCT  (EXPECT 3)  "
      f"[screen_pass={res['screen_pass']} in {time.time()-t0:.0f}s]", flush=True)
print("ANCHOR", "PASS" if res["distinct"] == 3 else f"MISMATCH ({res['distinct']} != 3)")
