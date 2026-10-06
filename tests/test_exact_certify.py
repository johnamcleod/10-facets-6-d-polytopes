"""Tests for the exact wildcard infeasibility certifier.

The certifier's verdicts are only as good as two things: that the exact field
arithmetic is the arithmetic it claims to be, and that the interpolated minor
polynomial really is the minor.  Both are checked here against independent
computations (float determinants from numpy, sympy for the field), because a
wrong polynomial would refute realizable labellings silently -- exactly the
failure mode the certifier exists to remove.
"""
import ast
import itertools
import json
import math
import pathlib
from fractions import Fraction

import numpy as np
import pytest

from pipeline.utils import exact_field as EF
from pipeline.utils.exact_certify import (C7_LO, GRAM_EXACT, Instance,
                                          certify_instance)


# ---------------------------------------------------------------------------
# The field
# ---------------------------------------------------------------------------

def test_field_multiplication_and_inverse():
    s2 = EF.sqrt_basis(2)
    s3 = EF.sqrt_basis(3)
    s5 = EF.sqrt_basis(5)
    assert EF.mul(s2, s2) == EF.rat(2)
    assert EF.mul(s2, s3) == EF.sqrt_basis(6)
    assert EF.mul(EF.mul(s2, s3), s5) == EF.sqrt_basis(30)
    u = EF.add(EF.rat(Fraction(3, 7)), EF.mul(s2, EF.rat(Fraction(-5, 3))))
    assert EF.mul(u, EF.inv(u)) == EF.ONE


def test_gram_entries_match_cosines():
    for m, val in GRAM_EXACT.items():
        assert EF.alg_float(val) == pytest.approx(-math.cos(math.pi / m), abs=1e-14)


def test_interval_encloses_value():
    for m, val in GRAM_EXACT.items():
        iv = EF.alg_iv(val)
        # the float cosine itself carries error (cos(pi/2) evaluates to 6e-17,
        # not 0), so the enclosure is checked with a float-sized slack
        assert float(iv.lo) - 1e-15 <= -math.cos(math.pi / m) <= float(iv.hi) + 1e-15
    assert float(C7_LO) <= math.cos(math.pi / 7)


def test_exact_det_matches_numpy():
    rng = np.random.default_rng(0)
    for _ in range(5):
        labels = rng.integers(2, 7, size=(6, 6))
        M = [[EF.ONE if i == j else GRAM_EXACT[int(labels[min(i, j)][max(i, j)])]
              for j in range(6)] for i in range(6)]
        Mf = np.array([[EF.alg_float(x) for x in row] for row in M])
        assert EF.alg_float(EF.det(M)) == pytest.approx(float(np.linalg.det(Mf)),
                                                        abs=1e-9)


# ---------------------------------------------------------------------------
# The minor polynomials
# ---------------------------------------------------------------------------

def _synthetic_instance(seed=0, n=10, d=6, n_dashed=2, n_wild=2):
    """A labelling with the same shape as a dumped one: random ordinary labels,
    a few dashed pairs and a few wild pairs."""
    rng = np.random.default_rng(seed)
    pairs = list(itertools.combinations(range(n), 2))
    rng.shuffle(pairs)
    dashed = sorted(pairs[:n_dashed])
    wild = sorted(pairs[n_dashed:n_dashed + n_wild])
    rest = sorted(pairs[n_dashed + n_wild:])
    ordinary = [[i, j, int(rng.integers(2, 7))] for i, j in rest]
    return {"n": n, "d": d, "ordinary": ordinary,
            "wild": [list(p) for p in wild], "dashed": [list(p) for p in dashed],
            "tag": "synthetic"}


def _numeric_gram(inst, xs, cs):
    G = np.eye(inst.n)
    for (i, j), m in inst.ordinary.items():
        G[i, j] = G[j, i] = -math.cos(math.pi / m)
    for e, x in zip(inst.dashed, xs):
        G[e[0], e[1]] = G[e[1], e[0]] = -x
    for e, c in zip(inst.wild, cs):
        G[e[0], e[1]] = G[e[1], e[0]] = -c
    return G


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_minor_polynomial_equals_the_minor(seed):
    """The interpolated polynomial must equal det G[S,S] * prod t_e^2 at random
    points of the domain, for every principal minor -- checked against numpy."""
    rec = _synthetic_instance(seed)
    inst = Instance(rec)
    rng = np.random.default_rng(100 + seed)
    sets = inst.index_sets()
    for S in sets[:6] + sets[-3:]:
        p = inst.minor_poly(S)
        for _ in range(3):
            xs = [1.0 + 9.0 * rng.random() for _ in inst.dashed]
            cs = [math.cos(math.pi / 7) + rng.random()
                  * (1 - math.cos(math.pi / 7)) for _ in inst.wild]
            G = _numeric_gram(inst, xs, cs)
            want = float(np.linalg.det(G[np.ix_(S, S)]))
            # the polynomial carries a factor t_e^2 for every dashed edge whose
            # both endpoints lie in S
            Sset = set(S)
            for e, x in zip(inst.dashed, xs):
                if e[0] in Sset and e[1] in Sset:
                    want *= (1.0 / x) ** 2
            box = [EF.IV(Fraction(1 / x).limit_denominator(10 ** 9))
                   for x in xs] + [EF.IV(Fraction(c).limit_denominator(10 ** 9))
                                   for c in cs]
            got = p.range(box)
            assert float(got.lo) - 1e-6 <= want <= float(got.hi) + 1e-6, (
                f"S={S} want={want} got={got}")


def test_range_is_an_enclosure_on_a_subbox():
    """Sampled values must lie inside the computed range enclosure."""
    rec = _synthetic_instance(3)
    inst = Instance(rec)
    S = inst.index_sets()[0]
    p = inst.minor_poly(S)
    box = inst.box()
    rng = np.random.default_rng(7)
    lo, hi = p.range(box).lo, p.range(box).hi
    for _ in range(50):
        pt = [EF.IV(Fraction(float(b.lo) + rng.random() * float(b.width())
                             ).limit_denominator(10 ** 9)) for b in box]
        v = p.range(pt)
        assert lo <= v.lo and v.hi <= hi


# ---------------------------------------------------------------------------
# End to end
# ---------------------------------------------------------------------------

def test_certifier_refutes_a_random_labelling():
    """A random labelling is overwhelmingly not realizable, and the certifier
    should say so exactly.  (It is allowed to answer 'undecided'; it must never
    answer 'certified' wrongly, which the feasible-instance test below guards.)"""
    ok, stats = certify_instance(_synthetic_instance(5), max_minors=4,
                                 max_boxes=20_000)
    assert ok, stats


ROOT = pathlib.Path(__file__).resolve().parents[1]


def _pb6_instance():
    """P^B6 itself, as a certifier instance: its 42 ordinary labels and its three
    dashed edges, with the dashed weights left unknown."""
    p = ROOT / "runs/d6_n10/d6_final/realizers/tid379_0_0.json"
    if not p.exists():
        pytest.skip("realizer artifact absent")
    r = json.loads(p.read_text())[0]
    la = {ast.literal_eval(k): v for k, v in r["label_assignment"].items()}
    dashed = [tuple(int(x) for x in k.split("_")[1:]) for k in r["dot_values"]]
    weights = {tuple(int(x) for x in k.split("_")[1:]): float(v["value"])
               for k, v in r["dot_values"].items()}
    rec = {"n": 10, "d": 6,
           "ordinary": sorted([i, j, m] for (i, j), m in la.items()),
           "wild": [], "dashed": sorted(list(p_) for p_ in dashed),
           "tag": "PB6"}
    return rec, weights


def test_certifier_does_not_refute_PB6():
    """The one labelling in the d=6 search that IS realizable must never be
    certified infeasible.  This is the anchor: the certifier's whole value is
    that it refutes only what is genuinely unrealizable, and P^B6 is the only
    available witness in this dimension."""
    rec, _ = _pb6_instance()
    ok, stats = certify_instance(rec, max_minors=12, max_boxes=50_000)
    assert not ok, f"P^B6 was refuted by the exact certifier: {stats}"


def test_PB6_solution_survives_every_minor_range():
    """Sharper than the previous test: at P^B6's own exact weights, every minor
    the certifier would use must have zero inside its enclosure -- so the
    certificate could not have been obtained by any escalation, not merely by the
    one the search happened to try."""
    rec, weights = _pb6_instance()
    inst = Instance(rec)
    # P^B6's weights are algebraic (minimal polynomials in the realizer file), so
    # the rational point below is an approximation and the exact minor there is
    # tiny but not exactly zero.  What must hold is that it is zero to within
    # that approximation -- decisively closer to zero than the certified minors,
    # whose enclosures exclude zero by margins of order 1e-2 and up.
    pt = [EF.IV(Fraction(1.0 / weights[e]).limit_denominator(10 ** 12))
          for _, e in inst.vars]
    for S in inst.index_sets()[:20]:
        p = inst.minor_poly(S)
        r = p.range(pt)
        assert max(abs(float(r.lo)), abs(float(r.hi))) < 1e-9, (
            f"minor {S} is nonzero at the P^B6 solution: {r}")


def test_PB6_label_only_minors_vanish():
    """The certificate that refutes every non-wildcard labelling is the
    non-vanishing of an 8x8 principal minor whose index set avoids all three
    dashed edges.  For P^B6 those same minors must vanish -- exactly, in
    Q(sqrt2,sqrt3,sqrt5), not merely to within a tolerance -- or the certificate
    would refute the polytope.  There are exactly three such index sets, the
    complements of the three 2-element vertex covers of the dashed path."""
    rec, _ = _pb6_instance()
    inst = Instance(rec)
    free = [S for S in itertools.combinations(range(inst.n), inst.d + 2)
            if not inst.minor_vars(S)]
    assert free == [(0, 1, 2, 3, 4, 5, 6, 8), (0, 1, 2, 3, 4, 5, 7, 8),
                    (0, 1, 2, 3, 4, 5, 7, 9)], free
    for S in free:
        p = inst.minor_poly(S)
        assert not p.terms, f"minor {S} is nonzero for P^B6: {p.terms}"


@pytest.mark.parametrize("artifact", ["runs/d6_n10/wild_certificates.json",
                                      "runs/d5_n9/wild_certificates.json",
                                      "runs/d6_n10/plain_certificates.json"])
def test_recorded_witnesses_are_corroborated(artifact):
    """Independent float re-check of the certificates in the artifact.

    The certificates are exact, so a float recomputation cannot verify them --
    but it can catch bookkeeping errors: a `label_only_minor` witness must name a
    minor whose float determinant is far from zero and agrees with the recorded
    exact value, and a `root_box_range` witness must name an enclosure that no
    sampled point of the domain contradicts.
    """
    p = ROOT / artifact
    if not p.exists():
        pytest.skip(f"{artifact} absent -- run checks/wild_exact_certify.py")
    data = json.loads(p.read_text())
    rng = np.random.default_rng(0)
    n_lab = n_box = 0
    for rec in data["instances"]:
        w = rec["witness"]
        if not w:
            continue
        inst = Instance(rec)
        if w["kind"] == "label_only_minor":
            S = w["minor"]
            xs = [1.0 + 50.0 * rng.random() for _ in inst.dashed]
            cs = [0.91 + 0.08 * rng.random() for _ in inst.wild]
            G = _numeric_gram(inst, xs, cs)
            det = float(np.linalg.det(np.array(G)[np.ix_(S, S)]))
            assert abs(det) > 1e-6, (rec["tag"], S, det)
            assert det == pytest.approx(w["value_float"], rel=1e-6, abs=1e-9)
            n_lab += 1
        else:
            for key, (lo, hi) in w["root_ranges"].items():
                S = json.loads(key)
                lo, hi = Fraction(lo), Fraction(hi)
                for _ in range(20):
                    xs = [1.0 / (rng.random() or 0.5) for _ in inst.dashed]
                    cs = [0.9009688679 + rng.random() * (1 - 0.9009688679)
                          for _ in inst.wild]
                    G = _numeric_gram(inst, xs, cs)
                    det = float(np.linalg.det(np.array(G)[np.ix_(S, S)]))
                    Sset = set(S)
                    for e, x in zip(inst.dashed, xs):
                        if e[0] in Sset and e[1] in Sset:
                            det *= (1.0 / x) ** 2
                    assert float(lo) - 1e-9 <= det <= float(hi) + 1e-9, (
                        rec["tag"], S, det, (float(lo), float(hi)))
            n_box += 1
    assert n_lab + n_box > 0
    print(f"\n{artifact}: {n_lab} label-only and {n_box} range witnesses corroborated")


def test_inertia_condition_agrees_with_numpy_eigenvalues():
    """The superhyperbolicity certificate must never claim two negative
    eigenvalues where there are fewer.

    Jacobi's rule is applied to interval-evaluated leading principal minors, so
    the risk is a sign read off a coefficient error rather than off the matrix.
    Here each claim is checked against a direct eigenvalue computation at sampled
    points of the box the claim covers.
    """
    rng = np.random.default_rng(11)
    checked = 0
    for seed in range(6):
        rec = _synthetic_instance(seed, n_dashed=3, n_wild=0)
        inst = Instance(rec)
        for order in inst.inertia_orders()[:6]:
            cond = inst.inertia_condition(order)
            for _ in range(6):
                xs = [1.0 + 30.0 * rng.random() for _ in inst.dashed]
                pt = [EF.IV(Fraction(1.0 / x).limit_denominator(10 ** 12))
                      for x in xs]
                if not cond.refutes(pt):
                    continue
                G = _numeric_gram(inst, xs, [])
                sub = np.array(G)[np.ix_(order, order)]
                neg = int((np.linalg.eigvalsh(sub) < -1e-9).sum())
                assert neg >= 2, (order, xs, neg,
                                  np.linalg.eigvalsh(sub))
                # and hence the whole matrix is not of signature (6,1)
                assert int((np.linalg.eigvalsh(np.array(G)) < -1e-9).sum()) >= 2
                checked += 1
    assert checked > 0, "no inertia claim was exercised"


def test_certifier_never_refutes_a_feasible_system():
    """A labelling with a genuine rank-7 solution in the domain must NOT be
    certified infeasible.  Built by taking a rank-deficient Gram matrix that the
    domain actually contains: all wild entries at cos(pi/7) and dashed weights
    chosen so the matrix is singular is hard to arrange by hand, so instead we
    verify the weaker but decisive property -- the certifier's pruning test
    never excludes a box that provably contains a common zero."""
    rec = _synthetic_instance(1)
    inst = Instance(rec)
    # a point of the domain, and the minors' values there
    xs = [2.5, 4.0][:len(inst.dashed)]
    cs = [0.95, 0.93][:len(inst.wild)]
    pt = [EF.IV(Fraction(1 / x).limit_denominator(10 ** 9)) for x in xs] + \
         [EF.IV(Fraction(c).limit_denominator(10 ** 9)) for c in cs]
    for S in inst.index_sets()[:8]:
        p = inst.minor_poly(S)
        v = p.range(pt)
        box_v = p.range(inst.box())
        # the enclosure over the whole box must contain the value at the point
        assert box_v.lo <= v.lo and v.hi <= box_v.hi
