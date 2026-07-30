"""ctypes wrapper for pipeline/c/wild_kernel.c — a fast small-matrix eigenvalue
kernel for the wildcard rank-deficiency objective (see stage4_gram._wild_feasible).

Computes IDENTICAL math to the pure-numpy objective_grad (same Jacobi-derived
formula: f = sum of squares of the num_zero smallest-|eigenvalue| eigenvalues,
grad_e = -4 * sum_i lambda_i * v[i0,i] * v[j0,i]) — this module exists purely to
remove numpy/scipy Python-dispatch overhead from a call made ~10^5-10^6 times
per subtree solve. If the shared library is unavailable (not built, wrong
platform), `available()` returns False and callers must fall back to the numpy
path — never silently change behavior.

Use Prepared(...) to pin the invariant arguments (G_base, ui, uj) to ctypes
buffers ONCE per L-BFGS-B run — re-marshaling those on every call (as a naive
per-call wrapper would) turned out to dominate the wall time as much as the
numpy path it replaced (measured: only 84->100 assignments/60s). Only `u`
(the actual optimization variable) is cast per call.
"""
from __future__ import annotations

import ctypes
from pathlib import Path

import numpy as np

_LIB_PATH = Path(__file__).resolve().parents[2] / "pipeline" / "c" / "libwildkernel.dylib"

_lib = None
_load_error = None
_c_double_p = ctypes.POINTER(ctypes.c_double)
_c_int_p = ctypes.POINTER(ctypes.c_int)
try:
    _lib = ctypes.CDLL(str(_LIB_PATH))
    _lib.wild_objective_grad.argtypes = [
        _c_double_p, ctypes.c_int,
        _c_int_p, _c_int_p, ctypes.c_int,
        _c_double_p, ctypes.c_int,
        _c_double_p, _c_double_p,
    ]
    _lib.wild_objective_grad.restype = ctypes.c_int
    _lib.wild_objective.argtypes = [
        _c_double_p, ctypes.c_int,
        _c_int_p, _c_int_p, ctypes.c_int,
        _c_double_p, ctypes.c_int,
    ]
    _lib.wild_objective.restype = ctypes.c_double
except OSError as e:
    _lib = None
    _load_error = e


def available() -> bool:
    return _lib is not None


class Prepared:
    """Pins G_base/ui/uj/output buffers to ctypes pointers ONCE; .objective and
    .objective_grad are closures over locals (not bound-method attribute
    lookups) that marshal only the k-length `u` vector per call and re-use a
    single pre-allocated ctypes.c_double output cell (no per-call ctypes
    object construction) — this matters because L-BFGS-B calls these
    thousands of times per subtree and measured profiling showed even cheap
    per-call Python/ctypes bookkeeping adds up to real wall time here."""

    def __init__(self, G_base: np.ndarray, ui: np.ndarray, uj: np.ndarray,
                num_zero: int):
        if _lib is None:
            raise RuntimeError(f"wild_kernel library not available: {_load_error}")
        n = G_base.shape[0]
        k = len(ui)
        self.n, self.k, self.num_zero = n, k, num_zero
        G = np.ascontiguousarray(G_base, dtype=np.float64)
        uia = np.ascontiguousarray(ui, dtype=np.int32)
        uja = np.ascontiguousarray(uj, dtype=np.int32)
        Gp = G.ctypes.data_as(_c_double_p)
        uip = uia.ctypes.data_as(_c_int_p)
        ujp = uja.ctypes.data_as(_c_int_p)
        grad_buf = np.empty(k, dtype=np.float64)
        gradp = grad_buf.ctypes.data_as(_c_double_p)
        u_buf = np.empty(k, dtype=np.float64)
        up = u_buf.ctypes.data_as(_c_double_p)
        f_cell = ctypes.c_double(0.0)
        fp = ctypes.pointer(f_cell)          # valid pointer to f_cell for the
                                              # lifetime of this Prepared instance
        c_obj = _lib.wild_objective
        c_grad = _lib.wild_objective_grad

        def objective(u):
            u_buf[:] = u
            return c_obj(Gp, n, uip, ujp, k, up, num_zero)

        def objective_grad(u):
            u_buf[:] = u
            rc = c_grad(Gp, n, uip, ujp, k, up, num_zero, fp, gradp)
            if rc != 0:
                raise RuntimeError(f"wild_objective_grad failed (n={n} > MAXN?)")
            return f_cell.value, grad_buf

        self.objective = objective
        self.objective_grad = objective_grad
        # keep references alive (ctypes pointers don't keep the numpy arrays
        # they point into alive on their own)
        self._keepalive = (G, uia, uja, grad_buf, u_buf, f_cell)


def objective_grad(G_base: np.ndarray, ui: np.ndarray, uj: np.ndarray,
                   u: np.ndarray, num_zero: int):
    """One-shot convenience wrapper (re-marshals every argument) — prefer
    Prepared for repeated calls with the same G_base/ui/uj (e.g. inside one
    L-BFGS-B run). Returns (f: float, grad: np.ndarray[len(u)])."""
    p = Prepared(G_base, ui, uj, num_zero)
    return p.objective_grad(u)


def objective(G_base: np.ndarray, ui: np.ndarray, uj: np.ndarray,
              u: np.ndarray, num_zero: int) -> float:
    p = Prepared(G_base, ui, uj, num_zero)
    return p.objective(u)
