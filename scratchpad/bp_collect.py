"""Collect block-paste per-type results and tabulate vs brute / census.

Reads runs/d4_n8/blockpaste/type_*.json, groups by k (#dotted = size-2 missing faces),
prints per-type block-paste distinct counts and the total, alongside the brute 334-run
per-type counts and the census per-k totals.
"""
import json, glob, os

T = {t["type_id"]: t for t in json.load(open("runs/d4_n8/stage2/types.json"))}

def k_of(tid):
    return sum(1 for m in T[tid]["missing_faces"] if len(m) == 2)

# Brute 334-run + lowk-completion per-type distinct (type 15 was 0 = incomplete).
BRUTE = {0:115,5:130,6:49, 1:3,2:0,3:2,7:15,18:1, 4:4,8:8,17:2,15:0,
         10:2,14:2,16:1,19:0, 20:0,24:0,22:0,26:0}
CENSUS_PER_K = {6:294, 5:23, 4:24, 3:7, 2:0}  # sums to 348

bp = {}
for f in sorted(glob.glob("runs/d4_n8/blockpaste/type_*.json")):
    r = json.load(open(f))
    bp[r["type_id"]] = r["distinct"]

print(f"{'tid':>3} {'k':>2} {'brute':>6} {'blockpaste':>11}")
done = sorted(bp)
for tid in done:
    b = BRUTE.get(tid, "-")
    flag = "" if (b == bp[tid] or b == "-") else "  <-- DIFF"
    print(f"{tid:>3} {k_of(tid):>2} {str(b):>6} {bp[tid]:>11}{flag}")

print(f"\nblock-paste types done: {len(bp)}/30   total distinct: {sum(bp.values())}")
# per-k totals (block-paste, done types only)
fromk = {}
for tid, d in bp.items():
    fromk.setdefault(k_of(tid), 0)
    fromk[k_of(tid)] += d
print("per-k (block-paste, done):", dict(sorted(fromk.items(), reverse=True)))
print("per-k census:            ", CENSUS_PER_K, " sum", sum(CENSUS_PER_K.values()))
missing = [tid for tid in range(30) if tid not in bp]
if missing:
    print("NOT YET DONE:", missing)
