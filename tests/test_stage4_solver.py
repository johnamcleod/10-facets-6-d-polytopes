"""Fast regression tests for the Stage-4 field-agnostic exact solver.

Guards the 2026-06-22 fix that replaced the fixed-√-basis recognition (which
could not represent cos(π/m) for m>=7 and mis-identified even 17+8√5) with
minimal-polynomial recovery + high-precision (d,1)-signature certification.
"""
import math
import pytest

mpmath = pytest.importorskip("mpmath")
from pipeline.stage4_gram import (
    _recover_minpoly, _numerical_screen, _refine_mpmath,
    _recognize_minpoly_and_verify, _GRAM_FLOAT,
)
from sympy import symbols


def test_minpoly_field_agnostic():
    """Minimal polynomials recovered for cubic, Q(√5), and Q(√2,√5) values.
    Coeffs [c0..cD], sum ci x^i = 0, leading coeff normalised > 0."""
    mpmath.mp.dps = 100
    assert _recover_minpoly(2 * mpmath.cos(mpmath.pi / 7)) == [1, -2, -1, 1]   # x³−x²−2x+1
    assert _recover_minpoly(17 + 8 * mpmath.sqrt(5)) == [-31, -34, 1]          # x²−34x−31
    assert _recover_minpoly(2 * mpmath.sqrt(2) + mpmath.sqrt(10)) == [4, 0, -36, 0, 1]  # x⁴−36x²+4


# P^B6 ordinary labels (0-indexed); the solver must recover its exact dotted weights.
PB6_ORDINARY = {(0, 1): 5, (1, 2): 3, (2, 3): 3, (3, 4): 3, (4, 5): 3,
                (8, 9): 5, (1, 6): 4, (1, 7): 4, (4, 8): 3, (5, 9): 5}
PB6_DOTTED = [(5, 6), (6, 7), (7, 8)]


def test_pb6_solver_recovers_exact_weights():
    """Given P^B6's labels, the refine+recognize path recovers the exact
    minimal polynomials and certifies signature (6,1).  ~2s."""
    n, d = 10, 6
    sym_list = [symbols(f"x_{p[0]}_{p[1]}", positive=True) for p in PB6_DOTTED]
    ordf = {p: _GRAM_FLOAT[m] for p, m in PB6_ORDINARY.items()}

    x_approx, resid = _numerical_screen(ordf, PB6_DOTTED, n, d)
    assert resid < 1e-10
    x_hp = _refine_mpmath(x_approx, PB6_ORDINARY, PB6_DOTTED, n, d, dps=100)
    assert x_hp is not None, "high-precision Gauss-Newton failed to converge"

    sols = _recognize_minpoly_and_verify(x_hp, sym_list, PB6_ORDINARY,
                                         PB6_DOTTED, n, d, dps=100)
    assert len(sols) == 1
    sol = sols[0]
    assert sol["x_5_6"]["minpoly"] == [4, 0, -36, 0, 1]      # 2√2+√10
    assert sol["x_7_8"]["minpoly"] == [4, 0, -36, 0, 1]
    assert sol["x_6_7"]["minpoly"] == [-31, -34, 1]          # 17+8√5
    # sanity on the decimal values
    assert abs(float(sol["x_5_6"]["value"]) - (2 * math.sqrt(2) + math.sqrt(10))) < 1e-9
    assert abs(float(sol["x_6_7"]["value"]) - (17 + 8 * math.sqrt(5))) < 1e-9
