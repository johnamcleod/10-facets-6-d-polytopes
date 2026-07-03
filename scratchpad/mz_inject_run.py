"""Run Ma-Zheng chcp59 on OUR d=5 type (by injecting vertex data), flag 2 and 1.
Bypasses the num->line file lookup by substituting our 0-indexed vertex_sets for the
file-parsed Vert (chcp59 then IncreaseOne's them to 1-indexed).
Compares chcp59 round-3 count to our pipeline basis-ON/OFF distinct paste-candidate count.
Usage: python3 scratchpad/mz_inject_run.py <tid>
"""
import subprocess, re, json, sys, tempfile, os, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HCP = ROOT + "/scratchpad/HCPdm"
SRC = HCP + "/pyFile/chcp59.py"

tid = int(sys.argv[1])
types = json.load(open(f"{ROOT}/runs/d5_n9/stage2/types.json"))
T = next(x for i, x in enumerate(types) if x.get("type_id", i) == tid)
V0 = T["vertex_sets"]                        # 0-indexed facets 0..8
inject = repr([sorted(v) for v in V0])       # chcp59 does IncreaseOne after
print(f"tid={tid}  nverts={len(V0)}  vertex_sets(0-idx)={V0}", flush=True)

src = open(SRC).read()
# substitute our Vert in place of the file-parsed one; keep the IncreaseOne+sort that follow
assert "Vert = list([convert_data5(_) for _ in Vert])" in src
src = src.replace("Vert = list([convert_data5(_) for _ in Vert])",
                  f"Vert = {inject}")

def run(flag):
    t = re.sub(r'^flag = 2.*$', f'flag = {flag}', src, flags=re.M)
    tf = tempfile.NamedTemporaryFile('w', suffix='.py', delete=False, dir='/tmp')
    tf.write(t); tf.close()
    t0 = time.time()
    try:
        p = subprocess.run([sys.executable, tf.name], cwd=HCP,
                           capture_output=True, text=True, timeout=900)
        m = re.search(r'data round 3:\s*\((\d+),', p.stdout)
        # also grab the printed l5_basis line for inspection
        lb = re.search(r'l5_basis \(.*', p.stdout)
        cnt = int(m.group(1)) if m else None
        info = lb.group(0) if lb else ""
        err = None if cnt is not None else (p.stderr.strip()[-400:] or p.stdout.strip()[-400:])
    except subprocess.TimeoutExpired:
        cnt, info, err = None, "", "TIMEOUT_900"
    os.unlink(tf.name)
    return cnt, info, err, time.time()-t0

for flag in (2, 1):
    cnt, info, err, dt = run(flag)
    print(f"  chcp59 flag={flag}: round3={cnt}  ({dt:.0f}s)  {info}", flush=True)
    if err: print(f"    ERR: {err}", flush=True)
