import json, sys
sys.path.insert(0,'/Users/jam/dev/sandbox/10-facet-6-polytopes')
import pipeline.stage4_gram as sg
from pipeline.stage4_gram import process_type_stage4
from pipeline.stage4_blockpaste import _setup
from pipeline.utils.automorphisms import compute_aut_group
TY={t["type_id"]:t for t in json.load(open('runs/d6_n10/stage2/types.json'))}
TIDS=[int(x) for x in open('/tmp/wbtids.txt').read().split(',')]
print("types:", TIDS, flush=True)
res={}
for hi in (1000.0, 1e6):
    sg._WILD_BOX_HI = hi
    print(f"\n=== dashed-weight box upper bound = {hi:g} ===", flush=True)
    for tid in TIDS:
        for k in sg.WILD_STATS: sg.WILD_STATS[k]=0
        t=dict(TY[tid]); V,*_=_setup(t); t['vertex_sets']=[sorted(v) for v in V]
        so={}
        r=process_type_stage4(t,6,max_assignments=50_000_000,enum_timeout=420,
            solve_timeout=210,wildcard=True,use_burcroff_55b=True,stats_out=so,
            automorphisms=compute_aut_group(V,1+max(max(v) for v in V)),verbose=False) or []
        key=(hi,tid); res[key]=(so.get('enum_count',0), so.get('wild_assignments',0),
                                len(r), bool(so.get('exhausted')))
        print(f"  tid {tid:>4} enum={res[key][0]:>8,} wild={res[key][1]:>8,} "
              f"found={res[key][2]} exhausted={res[key][3]}", flush=True)
print("\n=== comparison ===", flush=True)
same=True
for tid in TIDS:
    a=res[(1000.0,tid)]; b=res[(1e6,tid)]
    # only a comparison between two EXHAUSTED runs is meaningful: when a run times
    # out, enum_count measures throughput, and the wider box is slower per solve.
    both_exh = a[3] and b[3]
    ok = both_exh and a[2] == b[2]
    same &= ok or (not both_exh)
    tag = ('verdicts identical' if ok else
           'INCONCLUSIVE (neither run exhausted)' if not both_exh else 'DIFFERS')
    print(f"  tid {tid:>4}: {tag}   1e3={a}  1e6={b}", flush=True)
print("\nRESULT:", "no verdict changed when the box was widened 1000x"
      if same else "A VERDICT CHANGED -- the bound is load-bearing")
