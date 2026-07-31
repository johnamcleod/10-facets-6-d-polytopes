"""Exact arithmetic in the multiquadratic field K = Q(sqrt2, sqrt3, sqrt5), plus
rational interval arithmetic with outward rounding.

Every Gram entry of a compact Coxeter polytope with dihedral angles pi/m,
m in {2,...,6}, lies in K:

    -cos(pi/2) = 0      -cos(pi/3) = -1/2       -cos(pi/4) = -sqrt2/2
    -cos(pi/5) = -(1+sqrt5)/4                   -cos(pi/6) = -sqrt3/2

so the coefficients of every minor of such a Gram matrix lie in K as well, and
can be computed with no rounding at all.  K has degree 8 over Q with the basis
{ sqrt(prod S) : S subset of {2,3,5} }, indexed here by the bitmask of S
(bit 0 -> 2, bit 1 -> 3, bit 2 -> 5), so the basis order is

    1, sqrt2, sqrt3, sqrt6, sqrt5, sqrt10, sqrt15, sqrt30

-- note sqrt6 before sqrt5; prefer `sqrt_basis(radicand)` to positional
literals.  The bitmask indexing makes multiplication a one-liner:

    e_S * e_T = (prod of the primes in S cap T) * e_{S xor T}.

This module deliberately does NOT use sympy: the certifier evaluates thousands
of 8x8 determinants over K, and Fraction tuples are about two orders of
magnitude faster than Expr trees for that.

The interval layer (`IV`) exists so that an exact element of K, or an exact
polynomial over K, can be enclosed rigorously by rationals: all arithmetic
rounds OUTWARD onto a fixed denominator, so a computed enclosure is always a
superset of the true range.  That one-sidedness is what makes an emptiness
verdict from interval evaluation a proof rather than an estimate.
"""
from __future__ import annotations

from fractions import Fraction
from math import isqrt

PRIMES = (2, 3, 5)
DIM = 8                                    # 2^3 basis elements

# FACT[mask] = product of the primes selected by `mask`; also the radicand of the
# basis element e_mask = sqrt(FACT[mask]).
FACT = tuple(
    (2 if m & 1 else 1) * (3 if m & 2 else 1) * (5 if m & 4 else 1)
    for m in range(DIM)
)

ZERO = (Fraction(0),) * DIM
ONE = (Fraction(1),) + (Fraction(0),) * (DIM - 1)


def alg(*coeffs) -> tuple:
    """Element of K from up to 8 rational coefficients in basis order."""
    c = list(coeffs) + [0] * (DIM - len(coeffs))
    return tuple(Fraction(x) for x in c)


def rat(q) -> tuple:
    return (Fraction(q),) + (Fraction(0),) * (DIM - 1)


def sqrt_basis(radicand: int) -> tuple:
    """The basis element sqrt(radicand) for radicand in {1,2,3,5,6,10,15,30}."""
    mask = FACT.index(radicand)
    return tuple(Fraction(1) if i == mask else Fraction(0) for i in range(DIM))


def add(u, v):
    return tuple(a + b for a, b in zip(u, v))


def sub(u, v):
    return tuple(a - b for a, b in zip(u, v))


def neg(u):
    return tuple(-a for a in u)


def mul(u, v):
    res = [Fraction(0)] * DIM
    for i, ui in enumerate(u):
        if not ui:
            continue
        for j, vj in enumerate(v):
            if not vj:
                continue
            f = FACT[i & j]
            res[i ^ j] += ui * vj * f if f != 1 else ui * vj
    return tuple(res)


def scale(u, q):
    q = Fraction(q)
    return tuple(a * q for a in u)


def is_zero(u) -> bool:
    return not any(u)


# SIGN[e][mask] = the sign the automorphism e of K/Q puts on the basis element
# e_mask.  The Galois group is (Z/2)^3: each of sqrt2, sqrt3, sqrt5 may flip, and
# the sign on e_mask is the product of the flips of the primes it contains.
SIGN = tuple(
    tuple(-1 if bin(e & mask).count("1") % 2 else 1 for mask in range(DIM))
    for e in range(DIM)
)


def conj(u, e):
    """Apply the automorphism indexed by the flip-mask `e`."""
    if e == 0:
        return u
    s = SIGN[e]
    return tuple(c * s[m] if c else c for m, c in enumerate(u))


def inv(u):
    """Inverse in K, via the conjugates: 1/u = (prod_{sigma != 1} sigma(u)) / N(u),
    where the norm N(u) = prod_sigma sigma(u) is rational.  Cheaper than solving
    the 8x8 multiplication system, and the conjugates are sign flips, so the only
    real work is six field multiplications."""
    if is_zero(u):
        raise ZeroDivisionError("inverse of 0 in Q(sqrt2,sqrt3,sqrt5)")
    v = ONE
    for e in range(1, DIM):
        v = mul(v, conj(u, e))
    nrm = mul(u, v)
    if any(nrm[1:]) or not nrm[0]:
        raise ArithmeticError(f"norm is not a nonzero rational: {nrm}")
    return scale(v, Fraction(1, 1) / nrm[0])


def det(rows) -> tuple:
    """Determinant of a square matrix over K by Gaussian elimination.

    Exact: K is a field, so elimination with division is legitimate and no
    tolerance enters.  Used on 8x8 principal minors, where the cost is
    negligible next to the interpolation that calls it.
    """
    n = len(rows)
    M = [list(r) for r in rows]
    sign = 1
    acc = ONE
    for c in range(n):
        piv = next((r for r in range(c, n) if not is_zero(M[r][c])), None)
        if piv is None:
            return ZERO
        if piv != c:
            M[c], M[piv] = M[piv], M[c]
            sign = -sign
        p = M[c][c]
        acc = mul(acc, p)
        pinv = inv(p)
        for r in range(c + 1, n):
            if is_zero(M[r][c]):
                continue
            f = mul(M[r][c], pinv)
            M[r] = [sub(a, mul(f, b)) for a, b in zip(M[r], M[c])]
    return acc if sign > 0 else neg(acc)


# ---------------------------------------------------------------------------
# Rational interval arithmetic with outward rounding.
# ---------------------------------------------------------------------------

# Every operation rounds its endpoints outward onto denominator DEN.  Without
# this the denominators of a Horner evaluation grow without bound and the
# certifier slows to a crawl; with it, enclosures stay valid (only ever wider)
# at a fixed cost per operation.
DEN = 1 << 64


def _floor_den(q: Fraction) -> Fraction:
    return Fraction((q.numerator * DEN) // q.denominator, DEN)


def _ceil_den(q: Fraction) -> Fraction:
    return Fraction(-((-q.numerator * DEN) // q.denominator), DEN)


class IV:
    """Closed rational interval [lo, hi]; all operations round outward."""

    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None, _raw=False):
        if hi is None:
            hi = lo
        lo = Fraction(lo)
        hi = Fraction(hi)
        if not _raw:
            lo = _floor_den(lo)
            hi = _ceil_den(hi)
        if lo > hi:
            raise ValueError(f"empty interval [{lo}, {hi}]")
        self.lo, self.hi = lo, hi

    def __repr__(self):
        return f"IV({float(self.lo):.6g}, {float(self.hi):.6g})"

    def __add__(self, o):
        o = _as_iv(o)
        return IV(self.lo + o.lo, self.hi + o.hi)

    __radd__ = __add__

    def __neg__(self):
        return IV(-self.hi, -self.lo, _raw=True)

    def __sub__(self, o):
        return self + (-_as_iv(o))

    def __rsub__(self, o):
        return _as_iv(o) + (-self)

    def __mul__(self, o):
        o = _as_iv(o)
        a, b, c, d = self.lo, self.hi, o.lo, o.hi
        # sign-case-free: four products
        p = (a * c, a * d, b * c, b * d)
        return IV(min(p), max(p))

    __rmul__ = __mul__

    def __pow__(self, k: int):
        r = IV(1)
        for _ in range(k):
            r = r * self
        return r

    def contains_zero(self) -> bool:
        return self.lo <= 0 <= self.hi

    def width(self) -> Fraction:
        return self.hi - self.lo

    def mid(self) -> Fraction:
        return (self.lo + self.hi) / 2


def _as_iv(o):
    return o if isinstance(o, IV) else IV(o)


# Rational enclosures of the eight basis elements sqrt(FACT[mask]), tight to
# ~10^-30.  isqrt gives the floor of an integer square root, so
# [isqrt(k*S^2)/S, (isqrt(k*S^2)+1)/S] encloses sqrt(k) rigorously.
_SCALE = 10 ** 30


def _sqrt_iv(k: int) -> IV:
    if k == 1:
        return IV(1)
    r = isqrt(k * _SCALE * _SCALE)
    return IV(Fraction(r, _SCALE), Fraction(r + 1, _SCALE))


SQRT_IV = tuple(_sqrt_iv(FACT[m]) for m in range(DIM))


def alg_iv(u) -> IV:
    """Rigorous rational enclosure of an element of K."""
    acc = IV(0)
    for m, coeff in enumerate(u):
        if coeff:
            acc = acc + SQRT_IV[m] * coeff
    return acc


def alg_float(u) -> float:
    """Float value of an element of K, for diagnostics only -- never for a verdict."""
    from math import sqrt
    return float(sum(float(c) * sqrt(FACT[m]) for m, c in enumerate(u) if c))
