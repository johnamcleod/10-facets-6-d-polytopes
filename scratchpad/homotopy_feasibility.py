"""Feasibility probe for a per-type homotopy solver.

(1) A minimal total-degree homotopy tracker (numpy predictor-corrector) -- validates the
    method on a tiny polynomial system with known solutions, and shows the path count is the
    product of the equation degrees (Bezout).
(2) The per-type rank system's path count, computed for a real combinatorial type, to see
    whether homotopy's path count is tractable or reproduces the enumeration blowup.

The per-type system (variables = ordinary Gram entries c_e + dotted y_f):
  - cyclotomic constraint per ordinary edge:  prod_{m in VALID}(c_e + cos(pi/m)) = 0  (deg 10)
  - rank condition: every (d+2)-minor of the Gram vanishes
The total-degree (Bezout) path count is dominated by 10^(#ordinary edges).
"""
import json, itertools, cmath
import numpy as np


# ---------------------------------------------------------------------------
# (1) minimal total-degree homotopy on a tiny system (validation)
# ---------------------------------------------------------------------------
def total_degree_homotopy(F, dF, degs, steps=2000):
    """Track all product-of-degrees paths of the total-degree homotopy
    H(x,t) = (1-t) * gamma * G(x) + t * F(x),  G_i(x) = x_i^{d_i} - 1  (roots of unity start).
    Returns the finite solutions found (deduped)."""
    n = len(degs)
    gamma = complex(0.6, 0.8)                      # random unit gamma (avoid real singularities)
    # start solutions: all tuples of d_i-th roots of unity
    starts = itertools.product(*[[cmath.exp(2j * cmath.pi * k / degs[i]) for k in range(degs[i])]
                                 for i in range(n)])
    paths = 0
    sols = []
    for s in starts:
        paths += 1
        x = np.array(s, dtype=complex)
        ok = True
        for st in range(1, steps + 1):
            t = st / steps
            # H = (1-t)*gamma*(x^d - 1) + t*F  ; corrector: Newton on H(.,t)=0
            for _ in range(8):
                G = np.array([x[i] ** degs[i] - 1 for i in range(n)], dtype=complex)
                Hv = (1 - t) * gamma * G + t * F(x)
                Jg = np.diag([degs[i] * x[i] ** (degs[i] - 1) for i in range(n)])
                J = (1 - t) * gamma * Jg + t * dF(x)
                try:
                    dx = np.linalg.solve(J, -Hv)
                except np.linalg.LinAlgError:
                    ok = False; break
                x = x + dx
                if np.linalg.norm(dx) < 1e-12:
                    break
            if not ok or not np.all(np.isfinite(x)):
                ok = False; break
        if ok and np.linalg.norm(F(x)) < 1e-6:
            sols.append(x)
    # dedupe
    uniq = []
    for x in sols:
        if not any(np.linalg.norm(x - u) < 1e-6 for u in uniq):
            uniq.append(x)
    return uniq, paths


def _tiny_demo():
    # F = { x^2 - 2 , x*y - 1 }  -> solutions (sqrt2, 1/sqrt2), (-sqrt2, -1/sqrt2); Bezout 4
    F = lambda v: np.array([v[0] ** 2 - 2, v[0] * v[1] - 1], dtype=complex)
    dF = lambda v: np.array([[2 * v[0], 0], [v[1], v[0]]], dtype=complex)
    sols, paths = total_degree_homotopy(F, dF, [2, 2])
    print(f"tiny demo  x^2-2, xy-1:  tracked {paths} paths (Bezout 2*2=4), "
          f"found {len(sols)} finite solutions")
    for x in sorted(sols, key=lambda z: z[0].real):
        print(f"    x={x[0].real:+.4f}{x[0].imag:+.4f}i  y={x[1].real:+.4f}{x[1].imag:+.4f}i")


# ---------------------------------------------------------------------------
# (2) per-type path count
# ---------------------------------------------------------------------------
def per_type_path_count(t):
    n = 1 + max(max(m) for m in t["missing_faces"])
    dotted = [m for m in t["missing_faces"] if len(m) == 2]
    n_ord = n * (n - 1) // 2 - len(dotted)
    VALID = 10                                     # |VALID_LABELS|
    bez = VALID ** n_ord                           # dominant term (degree-10 per ordinary edge)
    return n, n_ord, len(dotted), bez


if __name__ == "__main__":
    _tiny_demo()
    print()
    T = {t["type_id"]: t for t in json.load(open("runs/d4_n8/stage2/types.json"))}
    for tid in (8, 20, 6):
        n, no, k, bez = per_type_path_count(T[tid])
        print(f"d=4 type {tid}: {no} ordinary edges, k={k} dotted -> "
              f"total-degree path count ~ 10^{no}  (= {bez:.1e})")
    print("d=6 low-k: ~40 ordinary edges -> ~10^40 paths")
