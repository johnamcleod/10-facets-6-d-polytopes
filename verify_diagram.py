"""
verify_diagram.py  –  Verify a Coxeter diagram is a compact hyperbolic d-polytope.

Handles:
  - Ordinary Coxeter edges  (i, j, m)        →  G_ij = −cos(π/m)
  - Dotted (ultraparallel) edges with exact weight  (i, j, "dotted", weight)
  - Dotted edges with unknown weight  (i, j, "dotted")  →  solved from rank condition

Usage: edit _RAW_EDGES, D, N below, then:
    python verify_diagram.py
"""

# ── INPUT ────────────────────────────────────────────────────────────────────
# 1-indexed nodes; m=2 (right angle) edges omitted; "dotted" for ultraparallel.
# For dotted edges with known exact weight, add the weight as a 4th element:
#   (i, j, "dotted", weight_expr)  where weight_expr is a sympy expression.
# For dotted edges with unknown weight, use (i, j, "dotted") — the solver will
# find the weight from the rank condition.

from sympy import sqrt, Rational
s2, s5, s10 = sqrt(2), sqrt(5), sqrt(10)

_RAW_EDGES = [
    # Compact hyperbolic Coxeter 6-polytope with 10 facets, d=6.
    # Field Q(sqrt2, sqrt5).
    # Verified: rank=7, signature=(6,1), no parabolic subdiagrams.
    # Dotted weights: w(6,7)=w(8,9)=sqrt(10)+2*sqrt(2)~5.99,  w(7,8)=17+8*sqrt(5)~34.89
    # Algebraic relation: w(7,8) = w(6,7)^2 - 1.
    (1,  2,  5),
    (2,  3,  3),
    (3,  4,  3),
    (4,  5,  3),
    (5,  6,  3),
    (6,  7,  "dotted", s10 + 2*s2),   # w ~ 5.99
    (7,  8,  "dotted", 17 + 8*s5),    # w ~ 34.89
    (8,  9,  "dotted", s10 + 2*s2),   # w ~ 5.99
    (9,  10, 5),
    (2,  7,  4),
    (2,  8,  4),
    (5,  9,  3),
    (6,  10, 5),
]

D = 6    # hyperbolic dimension
N = 10   # number of facets
# ─────────────────────────────────────────────────────────────────────────────

import sys, itertools

try:
    import numpy as np
    from numpy.linalg import eigvalsh
except ImportError:
    sys.exit("numpy is required.")

try:
    import sympy
    from sympy import (Matrix, Rational as Rat, sqrt as ssqrt, S, cos, pi,
                       symbols, zeros as sym_zeros, simplify, nsimplify)
    SYMPY = True
except ImportError:
    SYMPY = False
    print("WARNING: SymPy not available; exact verification skipped.")

try:
    import scipy.optimize as _opt
    SCIPY = True
except ImportError:
    SCIPY = False


# ---------------------------------------------------------------------------
# Exact Gram entry G_ij = −cos(π/m)
# ---------------------------------------------------------------------------

def gram_entry_exact(m):
    if m == 2:  return S.Zero
    if m == 3:  return Rat(-1, 2)
    if m == 4:  return -ssqrt(2) / 2
    if m == 5:  return -(1 + ssqrt(5)) / 4
    if m == 6:  return -ssqrt(3) / 2
    if m == 10: return -(ssqrt(5) - 1) / 4
    return -cos(pi / m)

def gram_entry_float(m):
    import math
    return -math.cos(math.pi / m)


# ---------------------------------------------------------------------------
# Parse edges (1-indexed → 0-indexed)
# ---------------------------------------------------------------------------

def parse_edges(raw):
    ordinary_float  = {}   # (i,j) → float
    ordinary_exact  = {}   # (i,j) → sympy
    dotted_known    = {}   # (i,j) → (float, sympy) — dotted with given weight
    dotted_unknown  = []   # [(i,j)] — dotted with weight to be solved

    for entry in raw:
        a, b = entry[0]-1, entry[1]-1
        m_or_d = entry[2]
        if a > b: a, b = b, a
        key = (a, b)

        if m_or_d == "dotted":
            if len(entry) == 4:
                weight_sym = entry[3]
                weight_fl  = float(weight_sym.evalf())
                dotted_known[key] = (weight_fl, weight_sym)
            else:
                dotted_unknown.append(key)
        else:
            m = m_or_d
            ordinary_float[key] = gram_entry_float(m)
            if SYMPY:
                ordinary_exact[key] = gram_entry_exact(m)

    return ordinary_float, ordinary_exact, dotted_known, dotted_unknown


# ---------------------------------------------------------------------------
# Build numpy Gram matrix
# ---------------------------------------------------------------------------

def build_gram_np(n, ordinary_float, dotted_known, dotted_unknown, x_vals):
    G = np.eye(n)
    for (i,j), v in ordinary_float.items():
        G[i,j] = G[j,i] = v
    for (i,j), (wf, _) in dotted_known.items():
        G[i,j] = G[j,i] = -wf
    for (i,j), x in zip(dotted_unknown, x_vals):
        G[i,j] = G[j,i] = -x
    return G


# ---------------------------------------------------------------------------
# Build exact SymPy Gram matrix
# ---------------------------------------------------------------------------

def build_gram_sym(n, ordinary_exact, dotted_known, dotted_unknown, x_vals_float):
    G = sympy.zeros(n)
    for i in range(n): G[i,i] = S.One
    for (i,j), v in ordinary_exact.items():
        G[i,j] = G[j,i] = v
    for (i,j), (_, ws) in dotted_known.items():
        G[i,j] = G[j,i] = -ws
    for (i,j), xf in zip(dotted_unknown, x_vals_float):
        xe = nsimplify(xf, [ssqrt(2), ssqrt(3), ssqrt(5)], tolerance=1e-9)
        G[i,j] = G[j,i] = -xe
        print(f"    dotted ({i+1},{j+1}): numeric={xf:.8f}  →  exact≈{xe}")
    return G


# ---------------------------------------------------------------------------
# Solve unknown dotted weights from rank = d+1 condition
# ---------------------------------------------------------------------------

def solve_dotted_unknown(n, d, ordinary_float, dotted_known, dotted_unknown):
    k = len(dotted_unknown)
    if k == 0:
        return np.array([]), 0.0

    def residual_vec(xs):
        G = build_gram_np(n, ordinary_float, dotted_known, dotted_unknown, xs)
        ev = np.sort(eigvalsh(G))
        return ev[1:3]   # 2nd and 3rd smallest → 0 for sig (d,1) rank d+1

    best_x, best_cost = None, 1e18
    starts = [(1.2,)*k, (1.5,)*k, (2.0,)*k, (phi_v,)*k, (5.0,)*k, (10.0,)*k]
    import math; phi_v = (1+math.sqrt(5))/2
    for x0 in starts:
        if len(x0) != k: x0 = (x0[0],)*k
        try:
            res = _opt.least_squares(
                residual_vec, x0,
                bounds=([1+1e-6]*k, [1e6]*k),
                method='trf', ftol=1e-14, xtol=1e-14, gtol=1e-14, max_nfev=50000)
            if res.cost < best_cost:
                best_cost = res.cost; best_x = res.x
        except Exception as e:
            print(f"    (optimizer error: {e})")
    return best_x, best_cost


# ---------------------------------------------------------------------------
# Parabolic subdiagram check
# ---------------------------------------------------------------------------

def check_parabolic(G_np, n, tol=1e-6):
    parabolic = []
    for size in range(2, n+1):
        for subset in itertools.combinations(range(n), size):
            sub = G_np[np.ix_(subset, subset)]
            ev = np.sort(eigvalsh(sub))
            if (ev < -tol).sum() == 0 and (np.abs(ev) <= tol).sum() >= 1:
                parabolic.append((subset, ev.tolist()))
    return parabolic


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 62)
    print(f"  Coxeter diagram verifier   N={N} facets, D={D}")
    print("=" * 62)

    ordinary_float, ordinary_exact, dotted_known, dotted_unknown = parse_edges(_RAW_EDGES)

    print(f"\nOrdinary edges ({len(ordinary_float)}):")
    for (i,j),v in sorted(ordinary_float.items()):
        m = next((e[2] for e in _RAW_EDGES
                  if sorted([e[0]-1,e[1]-1])==[i,j] and e[2]!="dotted"), "?")
        print(f"    {i+1} —[{m}]— {j+1}   G={v:.6f}")

    print(f"\nDotted edges with known weight ({len(dotted_known)}):")
    for (i,j),(wf,ws) in sorted(dotted_known.items()):
        print(f"    {i+1} ···· {j+1}   w={ws} ≈ {wf:.6f}  →  G={-wf:.6f}")

    if dotted_unknown:
        print(f"\nDotted edges with unknown weight ({len(dotted_unknown)}) — will solve:")
        for (i,j) in dotted_unknown:
            print(f"    {i+1} ···· {j+1}")
        print(f"\n[0] Solving for {len(dotted_unknown)} unknown dotted weight(s) …")
        x_vals, residual = solve_dotted_unknown(N, D, ordinary_float, dotted_known, dotted_unknown)
        if x_vals is None or residual > 1e-8:
            print(f"  FAIL: solver did not converge (residual={residual:.3e}).")
            return
        print(f"  Converged (residual={residual:.3e})")
        for (i,j),x in zip(dotted_unknown, x_vals):
            print(f"    w({i+1},{j+1}) = {x:.8f}")
    else:
        x_vals = np.array([])

    # Build float matrix
    G_np = build_gram_np(N, ordinary_float, dotted_known, dotted_unknown, x_vals)

    # ── Signature check ─────────────────────────────────────────────────────
    print(f"\n[1] Signature check")
    ev = np.sort(eigvalsh(G_np))
    neg = int((ev < -1e-6).sum())
    zer = int((np.abs(ev) <= 1e-6).sum())
    pos = int((ev > 1e-6).sum())
    rank = N - zer
    print(f"  Eigenvalues: {[f'{e:.5g}' for e in ev]}")
    print(f"  neg={neg}, zero={zer}, pos={pos}  →  rank={rank}, nullity={zer}")
    sig_ok = (neg == 1 and rank == D+1)
    print(f"  {'OK: signature ('+str(D)+',1) ✓' if sig_ok else 'FAIL: expected signature ('+str(D)+',1) with rank '+str(D+1)}'")

    # ── x > 1 check ─────────────────────────────────────────────────────────
    print(f"\n[2] Ultraparallel weight check (all dotted w > 1)")
    all_x_ok = True
    for (i,j),(wf,_) in sorted(dotted_known.items()):
        ok = wf > 1.0
        if not ok: all_x_ok = False
        print(f"    w({i+1},{j+1}) = {wf:.6f}  {'✓' if ok else 'FAIL'}")
    for (i,j),x in zip(dotted_unknown, x_vals):
        ok = x > 1.0
        if not ok: all_x_ok = False
        print(f"    w({i+1},{j+1}) = {x:.6f}  {'✓' if ok else 'FAIL'}")
    if all_x_ok:
        print("  All dotted weights > 1. ✓")

    # ── Parabolic subdiagram check ───────────────────────────────────────────
    print(f"\n[3] Parabolic subdiagram check (compactness)")
    print("  Scanning all 2^10 − 1 non-empty subsets … ", end="", flush=True)
    par = check_parabolic(G_np, N)
    print("done.")
    if par:
        print(f"  FAIL: {len(par)} parabolic subdiagram(s) found — NOT compact.")
        for sub, ev_sub in par[:5]:
            print(f"    nodes {tuple(s+1 for s in sub)}: min evals {[f'{e:.4g}' for e in ev_sub[:3]]}")
    else:
        print("  No parabolic subdiagrams. ✓")

    # ── Exact rank (SymPy) ───────────────────────────────────────────────────
    if SYMPY:
        print(f"\n[4] Exact rank (SymPy)")
        G_sym = build_gram_sym(N, ordinary_exact, dotted_known, dotted_unknown, x_vals)
        rk = G_sym.rank()
        rank_ok = (rk == D+1)
        print(f"  Exact rank = {rk}  {'✓' if rank_ok else 'FAIL (expected '+str(D+1)+')'}")
        # Check relation w(7,8) = w(6,7)^2 - 1 if applicable
        if len(dotted_known) >= 2:
            entries = sorted(dotted_known.items())
            ws_vals = [ws for _,(_, ws) in entries]
            if len(ws_vals) >= 2:
                from sympy import simplify
                rel = simplify(ws_vals[1] - (ws_vals[0]**2 - 1))
                print(f"  Algebraic check w(7,8) = w(6,7)²−1: residual = {rel}")
    else:
        rank_ok = sig_ok

    # ── Summary ─────────────────────────────────────────────────────────────
    print("\n" + "=" * 62)
    all_ok = sig_ok and all_x_ok and (len(par) == 0) and rank_ok
    if all_ok:
        print(f"  PASS: compact hyperbolic {D}-polytope with {N} facets ✓")
    else:
        print("  FAIL: one or more checks did not pass.")
    print("=" * 62)


if __name__ == "__main__":
    main()
