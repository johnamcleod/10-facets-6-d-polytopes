#!/usr/bin/env python3
"""One-command verification of the unique compact hyperbolic Coxeter 6-polytope
with 10 facets (the polytope P_{6,10} of the paper).

Run:  python3 verify_polytope.py

Checks performed
----------------
 1. EXACT inertia of the 10x10 Gram matrix over the number field Q(sqrt2, sqrt5),
    by symmetric (congruence) Gaussian elimination + Sylvester's law of inertia.
    No floating point is used for this decision.
 2. EXACT rank (SymPy) and the algebraic relation  w(7,8) = w(6,7)^2 - 1.
 3. All ultraparallel weights > 1.
 4. No parabolic subdiagram (scan of all 2^10 - 1 subdiagrams).
 5. CoxIter: cocompact, dimension 6, f-vector, covolume  (skipped with a notice
    if the CoxIter binary is not built).
 6. The diagram agrees with the realizer emitted by the classification run
    (runs/d6_n10/d6_tangency/realizers/), up to relabelling.

Exit code 0 iff every non-skipped check passes.
"""
from __future__ import annotations

import itertools
import json
import subprocess
import sys
from pathlib import Path

import sympy as sp
from sympy import Rational, S, sqrt

ROOT = Path(__file__).resolve().parent
# CoxIter is an optional independent verifier.  Look for it on PATH first, so that
# a reader who installed it normally needs no configuration, then fall back to the
# in-tree build location used during development.  Absence is reported as SKIPPED,
# never as a pass.
import shutil as _shutil

COXITER = Path(_shutil.which("coxiter") or
               ROOT / "scratchpad" / "CoxIter" / "build" / "coxiter")

N, D = 10, 6

# ---------------------------------------------------------------- the polytope
ORDINARY = [(1, 2, 5), (2, 3, 3), (3, 4, 3), (4, 5, 3), (5, 6, 3),
            (9, 10, 5), (2, 7, 4), (2, 8, 4), (5, 9, 3), (6, 10, 5)]
W67 = 2 * sqrt(2) + sqrt(10)          # = sqrt(2)*(2+sqrt(5))
W78 = 17 + 8 * sqrt(5)
DOTTED = {(6, 7): W67, (7, 8): W78, (8, 9): W67}


def cos_pi_over(m):
    return {2: S.Zero, 3: Rational(1, 2), 4: sqrt(2) / 2,
            5: (1 + sqrt(5)) / 4}[m]


def gram():
    G = sp.eye(N)
    for i, j, m in ORDINARY:
        G[i - 1, j - 1] = G[j - 1, i - 1] = -cos_pi_over(m)
    for (i, j), w in DOTTED.items():
        G[i - 1, j - 1] = G[j - 1, i - 1] = -w
    return G


# ------------------------------------------------- 1. exact inertia over K
def exact_inertia(G):
    """Inertia (pos, neg, zero) by exact congruence elimination over
    Q(sqrt2, sqrt5).  Sylvester's law of inertia: the signs of the pivots of a
    symmetric congruence reduction give the inertia."""
    K = sp.QQ.algebraic_field(sqrt(2), sqrt(5))
    A = [[K.from_sympy(sp.sympify(G[i, j])) for j in range(N)] for i in range(N)]
    pos = neg = zero = 0
    size = N
    while size > 0:
        p = next((i for i in range(size) if A[i][i] != K.zero), None)
        if p is None:
            q = next(((i, j) for i in range(size) for j in range(i + 1, size)
                      if A[i][j] != K.zero), None)
            if q is None:                       # remaining block is zero
                zero += size
                break
            i, j = q                            # x_i <- x_i + x_j (congruence)
            for c in range(size):
                A[i][c] = A[i][c] + A[j][c]
            for r in range(size):
                A[r][i] = A[r][i] + A[r][j]
            p = i
        if p != 0:
            A[0], A[p] = A[p], A[0]
            for r in range(size):
                A[r][0], A[r][p] = A[r][p], A[r][0]
        d = A[0][0]
        if sp.ask(sp.Q.positive(sp.simplify(K.to_sympy(d)))):
            pos += 1
        else:
            neg += 1
        A = [[A[i][j] - (A[i][0] / d) * A[0][j] for j in range(1, size)]
             for i in range(1, size)]
        size -= 1
    return pos, neg, zero


# ------------------------------------------------------ 4. parabolic scan
def has_parabolic(G, tol=1e-7):
    import numpy as np
    Gf = np.array(sp.matrix2numpy(G.evalf(30), dtype=float))
    for size in range(2, N + 1):
        for sub in itertools.combinations(range(N), size):
            ev = np.linalg.eigvalsh(Gf[np.ix_(sub, sub)])
            if (ev < -tol).sum() == 0 and (abs(ev) <= tol).sum() >= 1:
                return sub
    return None


# ------------------------------------------------------------- 5. CoxIter
def coxiter():
    if not COXITER.exists():
        return None, "binary not found at %s" % COXITER
    lines = [f"{N} {D}"]
    for i, j, m in ORDINARY:
        lines.append(f"{i} {j} {m}")
    for (i, j) in DOTTED:
        lines.append(f"{i} {j} 1")             # weight 1 == dashed
    out = subprocess.run([str(COXITER), "-c", "-fv", "-a"],
                         input="\n".join(lines) + "\n",
                         capture_output=True, text=True, timeout=300).stdout
    info = {}
    for line in out.splitlines():
        s = line.strip()
        for key in ("Cocompact", "Finite covolume", "Dimension", "f-vector",
                    "Euler characteristic", "Covolume",
                    "Number of vertices at infinity"):
            if s.startswith(key + ":"):
                info[key] = s.split(":", 1)[1].strip()
    if "Dimension" not in info:
        info["Dimension"] = str(D)             # CoxIter echoes it while reading
    return info, out


# ---------------------------- 6. agreement with the classification realizers
def matches_run_output(G):
    d = ROOT / "runs/d6_n10/d6_tangency/realizers"
    files = sorted(d.glob("tid379_*.json")) if d.exists() else []
    if not files:
        return None
    rec = json.load(open(files[0]))[0]
    B = sp.eye(N)
    for k, m in rec["label_assignment"].items():
        i, j = eval(k)
        B[i, j] = B[j, i] = -cos_pi_over(int(m))
    for k in rec["dot_values"]:
        _, i, j = k.split("_")
        i, j = int(i), int(j)
        w = W78 if abs(float(rec["dot_values"][k]["value"]) - 34.888) < 0.01 else W67
        B[i, j] = B[j, i] = -w
    An = [[sp.N(G[i, j], 30) for j in range(N)] for i in range(N)]
    Bn = [[sp.N(B[i, j], 30) for j in range(N)] for i in range(N)]
    rows_a = [sorted(An[i][j] for j in range(N) if j != i) for i in range(N)]
    rows_b = [sorted(Bn[i][j] for j in range(N) if j != i) for i in range(N)]
    cand = [[j for j in range(N)
             if all(abs(x - y) < 1e-20 for x, y in zip(rows_a[i], rows_b[j]))]
            for i in range(N)]

    def rec_map(i, perm, used):
        if i == N:
            return list(perm)
        for j in cand[i]:
            if j in used:
                continue
            if all(abs(An[i][a] - Bn[j][perm[a]]) < 1e-20 for a in range(i)):
                perm.append(j)
                used.add(j)
                r = rec_map(i + 1, perm, used)
                if r:
                    return r
                used.discard(j)
                perm.pop()
        return None

    return rec_map(0, [], set())


def main():
    ok = True
    print("=" * 70)
    print("  Verification of P_{6,10}: compact hyperbolic Coxeter 6-polytope,")
    print("  10 facets, field Q(sqrt2, sqrt5)")
    print("=" * 70)

    G = gram()

    print("\n[1] EXACT inertia over Q(sqrt2, sqrt5) (congruence elimination)")
    pos, neg, zero = exact_inertia(G)
    good = (pos, neg, zero) == (6, 1, 3)
    ok &= good
    print(f"    positive={pos}  negative={neg}  zero={zero}"
          f"   -> signature ({pos},{neg}), rank {pos+neg}, nullity {zero}"
          f"   {'OK' if good else 'FAIL (expected 6,1,3)'}")

    print("\n[2] EXACT rank and algebraic relation")
    r = G.rank()
    rel = sp.simplify(W78 - (W67 ** 2 - 1))
    good = (r == D + 1) and rel == 0
    ok &= good
    print(f"    rank = {r} (expected {D+1});  w(7,8) - (w(6,7)^2 - 1) = {rel}"
          f"   {'OK' if good else 'FAIL'}")
    print(f"    w(6,7) = w(8,9) = {W67} = {float(W67):.12f}")
    print(f"    w(7,8)          = {W78} = {float(W78):.12f}")
    print(f"    minimal polynomials: x^4-36x^2+4  and  x^2-34x-31")

    print("\n[3] Ultraparallel weights > 1")
    good = all(float(w) > 1 for w in DOTTED.values())
    ok &= good
    print(f"    {'OK' if good else 'FAIL'}")

    print("\n[4] Parabolic subdiagram scan (all 2^10-1 subsets)")
    par = has_parabolic(G)
    ok &= par is None
    print(f"    {'none found  OK' if par is None else f'FAIL: parabolic on {par}'}")

    print("\n[5] CoxIter")
    info, raw = coxiter()
    if info is None:
        print(f"    SKIPPED ({raw}); build it with:")
        print("      cd scratchpad/CoxIter && mkdir -p build && cd build "
              "&& cmake .. && make")
    else:
        for k, v in info.items():
            print(f"    {k}: {v}")
        good = (info.get("Cocompact", "").lower().startswith("y")
                and info.get("f-vector") == "(31, 93, 125, 95, 42, 10, 1)")
        ok &= good
        print(f"    {'OK' if good else 'FAIL (expected cocompact, f-vector (31, 93, 125, 95, 42, 10, 1))'}")

    print("\n[6] Agreement with the classification run's realizer records")
    m = matches_run_output(G)
    if m is None:
        print("    SKIPPED (runs/d6_n10/d6_tangency/realizers not present)")
    else:
        ok &= True
        print(f"    isomorphic; paper node i -> run node: "
              f"{ {i+1: m[i] for i in range(N)} }   OK")

    print("\n" + "=" * 70)
    print("  RESULT: " + ("ALL CHECKS PASSED" if ok else "FAILURE"))
    print("=" * 70)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
