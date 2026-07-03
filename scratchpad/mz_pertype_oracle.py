"""Authoritative per-type d=5 census counts from Ma-Zheng chcp59.
For each num, run flag=2 (basis ON = census canonical) and flag=1 (basis OFF),
capture the 'data round 3' shape[0] = per-type distinct-count.
Run from HCPdm root.  Results -> scratchpad/mz_oracle_results.json
"""
import subprocess, re, json, sys, tempfile, os, time

HCP = os.path.dirname(os.path.abspath(__file__)) + "/HCPdm"
SRC = HCP + "/pyFile/chcp59.py"

# my_idx -> Ma-Zheng num  (from d5-gap-diagnosis memory)
IDX2NUM = {1:254,4:255,8:249,15:10,16:2,20:280,30:248,31:268,53:242,55:4,60:9}
# small/fast nums first; 322 anchor last (known 3/18)
NUMS = [2, 4, 9, 10, 242, 248, 249, 254, 255, 268, 280, 322]

src = open(SRC).read()
res = {}
for num in NUMS:
    res[num] = {}
    for flag in (2, 1):
        t = src
        t = re.sub(r'^num = 322.*$', f'num = {num}', t, flags=re.M)
        t = re.sub(r'^flag = 2.*$', f'flag = {flag}', t, flags=re.M)
        tf = tempfile.NamedTemporaryFile('w', suffix='.py', delete=False, dir='/tmp')
        tf.write(t); tf.close()
        t0 = time.time()
        try:
            p = subprocess.run([sys.executable, tf.name], cwd=HCP,
                               capture_output=True, text=True, timeout=1800)
            out = p.stdout
            m = re.search(r'data round 3:\s*\((\d+),', out)
            cnt = int(m.group(1)) if m else None
            err = None if cnt is not None else (p.stderr.strip()[-300:] or out.strip()[-300:])
        except subprocess.TimeoutExpired:
            cnt, err = None, "TIMEOUT_1800s"
        res[num][f"flag{flag}"] = cnt
        res[num].setdefault("err", {})[f"flag{flag}"] = err
        os.unlink(tf.name)
        print(f"num={num} flag={flag} -> {cnt}  ({time.time()-t0:.0f}s)  err={str(err)[:80]}", flush=True)
        json.dump(res, open(os.path.dirname(SRC)+"/../../mz_oracle_results.json","w"), indent=2)
print("DONE")
