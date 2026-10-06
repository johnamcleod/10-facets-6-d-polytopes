"""Exact infeasibility certificates for wildcard-bearing labellings.

WHAT IS BEING PROVED.  A wildcard-bearing labelling fixes the ordinary Gram
entries -cos(pi/m), m in {2,...,6}, leaves the dashed (ultraparallel) weights
x_e > 1 free, and leaves each wild entry free as c_e = cos(pi/m_e) with
m_e >= 7, i.e. c_e in [cos(pi/7), 1).  The labelling can only be realized if the
Gram matrix G has rank <= d+1, which forces every (d+2)x(d+2) PRINCIPAL minor of
G to vanish.  So exhibiting a set of principal minors with no common zero over
the whole domain refutes the labelling outright -- no local search, no bounded
weight range, no tolerance.

WHY THE DOMAIN IS COMPACT.  The wild entries already range over the compact
interval [cos(pi/7), 1); the integer scan happens only after feasibility, so an
infeasibility certificate for the continuous relaxation proves that no integer
m_e >= 7 works, and no unbounded quantifier over labels arises.  The dashed
weights are the unbounded direction, and the substitution

    x_e = 1 / t_e,      x_e in (1, infinity)  <->  t_e in (0, 1),

removes it: a minor is a polynomial of degree <= 2 in each x_e (each unknown
occupies exactly the two symmetric positions (i,j), (j,i) of one edge), so
multiplying by t_e^2 turns it into a polynomial in t_e.  Every minor thus
becomes a polynomial on the CLOSED rational box

    t_e in [0, 1],      c_e in [c7lo, 1],      c7lo <= cos(pi/7) rational,

which contains the true domain.  Certifying emptiness on the closed box
therefore certifies it on the domain; the enlargement is in the safe direction.

HOW EMPTINESS IS CERTIFIED.  Interval branch-and-bound.  A box is discarded as
soon as some minor's interval range excludes 0.  All coefficient arithmetic is
exact in Q(sqrt2,sqrt3,sqrt5) (pipeline/utils/exact_field.py) and every interval
operation rounds outward, so a discarded box provably contains no common zero.
If every box is discarded the certificate is complete.  If the box budget is
exhausted the answer is "not certified" -- never a silent success.

This module is the machinery; checks/wild_exact_certify.py is the driver
that runs it over the dumped instances.
"""
from __future__ import annotations

import itertools
from fractions import Fraction

from pipeline.utils import exact_field as EF
from pipeline.utils.exact_field import IV

# -cos(pi/m) as exact elements of K = Q(sqrt2,sqrt3,sqrt5).  Built by radicand
# via sqrt_basis, not by basis position: the basis is indexed by the BITMASK of
# the prime subset (mask 3 = sqrt6, not sqrt5), so positional literals are easy
# to get wrong and this is the load-bearing table of the whole certifier.
GRAM_EXACT = {
    2: EF.rat(0),
    3: EF.rat(Fraction(-1, 2)),
    4: EF.scale(EF.sqrt_basis(2), Fraction(-1, 2)),        # -sqrt2/2
    5: EF.add(EF.rat(Fraction(-1, 4)),
              EF.scale(EF.sqrt_basis(5), Fraction(-1, 4))),  # -(1+sqrt5)/4
    6: EF.scale(EF.sqrt_basis(3), Fraction(-1, 2)),        # -sqrt3/2
}

# A rational lower bound for cos(pi/7) = 0.90096886790241912...  cos(pi/7) is
# cubic over Q, hence not in K, so the box uses a rational bound BELOW it; that
# only enlarges the domain, which is the safe direction for a refutation.
C7_LO = Fraction(9009688679, 10 ** 10)


# ---------------------------------------------------------------------------
# Sparse multivariate polynomials with coefficients in K.
# Represented as {exponent tuple: K element}; degree <= 2 in each variable.
# ---------------------------------------------------------------------------

def poly_from_values(values, nvars, nodes=(0, 1, 2)):
    """Tensor-product interpolation of a degree-<=2-per-variable polynomial.

    `values` maps each multi-index in nodes^nvars (as a tuple of node indices)
    to the K-value of the polynomial there.  Exact: the Vandermonde inverse for
    three distinct rational nodes is rational, so no rounding occurs.
    """
    # Univariate inverse Vandermonde for nodes (v0,v1,v2): coefficients of
    # p(v) = a0 + a1 v + a2 v^2 from (f(v0), f(v1), f(v2)).
    v0, v1, v2 = (Fraction(v) for v in nodes)
    W = _vandermonde_inverse(v0, v1, v2)

    cur = dict(values)
    for axis in range(nvars):
        nxt = {}
        keys = sorted({k[:axis] + k[axis + 1:] for k in cur})
        for k in keys:
            # zero coefficients must be KEPT here: a later axis looks up all
            # three nodes of this one, and dropping a vanishing intermediate
            # coefficient makes that lookup fail (or, worse, silently wrong).
            fs = [cur.get(k[:axis] + (j,) + k[axis:], EF.ZERO) for j in range(3)]
            for deg in range(3):
                coeff = EF.ZERO
                for j in range(3):
                    if W[deg][j]:
                        coeff = EF.add(coeff, EF.scale(fs[j], W[deg][j]))
                nxt[k[:axis] + (deg,) + k[axis:]] = coeff
        cur = nxt
    return {k: v for k, v in cur.items() if not EF.is_zero(v)}


def _vandermonde_inverse(v0, v1, v2):
    """Rows deg=0,1,2 of the inverse of [[1,v,v^2]] at the three nodes."""
    V = [[Fraction(1), v, v * v] for v in (v0, v1, v2)]
    # invert the 3x3 exactly
    A = [row[:] + [Fraction(1) if i == j else Fraction(0) for j in range(3)]
         for i, row in enumerate(V)]
    for c in range(3):
        piv = next(r for r in range(c, 3) if A[r][c])
        A[c], A[piv] = A[piv], A[c]
        p = A[c][c]
        A[c] = [x / p for x in A[c]]
        for r in range(3):
            if r != c and A[r][c]:
                f = A[r][c]
                A[r] = [a - f * b for a, b in zip(A[r], A[c])]
    inv = [[A[r][3 + c] for c in range(3)] for r in range(3)]
    # inv maps values -> coefficients: coeff[deg] = sum_j inv[deg][j] f_j
    return inv


def poly_reciprocal_vars(poly, recip_axes, nvars):
    """Substitute v -> 1/v on the given axes and clear denominators.

    Multiplies by v^2 on each reciprocal axis, i.e. maps the exponent a to 2-a
    there.  Since the polynomial has degree <= 2 in every variable this is exact
    and stays polynomial; the resulting polynomial has the same zero set on
    v != 0 as the original had on the reciprocal point.
    """
    out = {}
    for exps, coeff in poly.items():
        e = list(exps)
        for ax in recip_axes:
            e[ax] = 2 - e[ax]
        out[tuple(e)] = coeff
    return out


class IPoly:
    """A polynomial with interval coefficients, ready for range evaluation."""

    __slots__ = ("terms", "nvars", "name", "exact", "index_set")

    def __init__(self, poly, nvars, name="", index_set=None):
        self.terms = [(exps, EF.alg_iv(c)) for exps, c in poly.items()
                      if not EF.is_zero(c)]
        self.nvars = nvars
        self.name = name
        self.index_set = index_set
        self.exact = poly

    def __bool__(self):
        return bool(self.terms)

    def is_nonzero_constant(self):
        """True if this polynomial is a nonzero constant.

        Decided on the EXACT field element, not on its interval enclosure: a
        nonzero element of K is nonzero, whereas an enclosure of a very small
        value could still straddle zero and lose a valid certificate.  The
        interval is kept only to report the witness.
        """
        return (len(self.terms) == 1
                and all(e == 0 for e in self.terms[0][0])
                and not EF.is_zero(self.exact[self.terms[0][0]]))

    def range(self, box):
        """Rigorous enclosure of the polynomial's range over `box`.

        Naive monomial form, with the powers of each variable cached per box.
        Wider than a centered form, but every wildcard instance here dies at or
        very near the root box, so the extra width costs nothing in practice and
        the form is simple enough to be obviously sound.
        """
        pw = [[IV(1), v, v * v] for v in box]
        acc = IV(0)
        for exps, coeff in self.terms:
            term = coeff
            for ax, e in enumerate(exps):
                if e:
                    term = term * pw[ax][e]
            acc = acc + term
        return acc


# ---------------------------------------------------------------------------
# Instance -> minors
# ---------------------------------------------------------------------------

class Instance:
    """One wildcard-bearing labelling, as dumped by the enumerator.

    Variables are ordered: the dashed weights first (as t_e = 1/x_e), then the
    wild entries (as c_e).  That is the order used by every box and polynomial
    in this module.
    """

    def __init__(self, rec):
        self.n = rec["n"]
        self.d = rec["d"]
        self.ordinary = {(i, j): m for i, j, m in rec["ordinary"]}
        self.wild = [tuple(p) for p in rec["wild"]]
        self.dashed = [tuple(p) for p in rec["dashed"]]
        self.tag = rec.get("tag", "")
        self.residual = rec.get("joint_residual")
        self.vars = [("t", e) for e in self.dashed] + [("c", e) for e in self.wild]
        self.var_index = {e: i for i, (_, e) in enumerate(self.vars)}
        self.nvars = len(self.vars)
        self.minor_size = self.d + 2

    def box(self):
        """The closed rational box containing the domain."""
        return [IV(0, 1) if kind == "t" else IV(C7_LO, 1)
                for kind, _ in self.vars]

    def entry_exact(self, i, j):
        """The (i,j) Gram entry as a K element, or None if it is a variable."""
        if i == j:
            return EF.ONE
        p = (i, j) if i < j else (j, i)
        if p in self.var_index:
            return None
        m = self.ordinary.get(p)
        if m is None:
            raise KeyError(f"pair {p} is neither labelled, dashed nor wild")
        if m not in GRAM_EXACT:
            raise ValueError(f"label {m} on {p} is not in the exact table "
                             f"(labels 2..6 only; wild edges must be listed as wild)")
        return GRAM_EXACT[m]

    def minor_vars(self, S):
        """Indices (into self.vars) of the unknowns occurring in minor S."""
        Sset = set(S)
        return [self.var_index[e] for _, e in self.vars
                if e[0] in Sset and e[1] in Sset]

    def index_sets(self):
        """Principal (d+2)-subsets, cheapest first: fewest unknowns, and among
        those, dashed-free ones first (a dashed unknown costs a reciprocal axis
        whose box reaches to x = infinity)."""
        sets = list(itertools.combinations(range(self.n), self.minor_size))
        def key(S):
            vs = self.minor_vars(S)
            n_t = sum(1 for v in vs if self.vars[v][0] == "t")
            return (len(vs), n_t, S)
        return sorted(sets, key=key)

    def inertia_orders(self, max_size=5):
        """Orderings of principal submatrices worth testing for two negative
        eigenvalues, best first.

        A single dashed edge $\\{i,j\\}$ already contributes one negative
        eigenvalue: the $2\\times2$ block $[[1,-x],[-x,1]]$ has determinant
        $1-x^2<0$ for $x>1$.  Two DISJOINT dashed edges therefore give two
        negatives unless the coupling between them is strong enough to pull one
        back, which is the classical source of superhyperbolicity, so those
        orderings come first; extensions by a further facet come next.
        """
        orders = []
        D = self.dashed
        for a in range(len(D)):
            for b in range(a + 1, len(D)):
                e, f = D[a], D[b]
                if set(e) & set(f):
                    continue                       # not disjoint: only one block
                orders.append((e[0], e[1], f[0], f[1]))
        # extensions: a disjoint pair plus one more facet
        base = list(orders)
        if max_size >= 5:
            for o in base:
                for v in range(self.n):
                    if v not in o:
                        orders.append(o + (v,))
        # a single dashed edge plus two further facets, as a fallback for
        # instances with no two disjoint dashed edges
        for e in D:
            for v in range(self.n):
                if v in e:
                    continue
                for w in range(v + 1, self.n):
                    if w in e:
                        continue
                    orders.append((e[0], e[1], v, w))
        return orders

    def inertia_condition(self, order, need=2):
        """Leading principal minors of the submatrix on `order`, as polynomials."""
        polys = [self.minor_poly(order[:m]) for m in range(2, len(order) + 1)]
        return InertiaCondition(polys, order, need=need)

    def minor_poly(self, S):
        """The minor det G[S,S], as a polynomial in the box variables.

        Computed by exact tensor-product interpolation over K: the determinant
        is evaluated exactly at 3^r integer nodes in the r unknowns of this
        minor (degree <= 2 in each), then the dashed axes are inverted by
        x = 1/t and denominators cleared.  Returns (IPoly, list of global var
        indices) -- the polynomial's axes are the FULL variable list, so all
        minors of an instance share one box.
        """
        vs = self.minor_vars(S)
        r = len(vs)
        local = {self.vars[v][1]: k for k, v in enumerate(vs)}

        base = [[self.entry_exact(i, j) for j in S] for i in S]
        values = {}
        for node in itertools.product(range(3), repeat=r):
            M = [row[:] for row in base]
            for a, i in enumerate(S):
                for b, j in enumerate(S):
                    if M[a][b] is None:
                        p = (i, j) if i < j else (j, i)
                        # entry is -value; nodes 0,1,2 are interpolation nodes
                        M[a][b] = EF.rat(-node[local[p]])
            values[node] = EF.det(M)

        poly = poly_from_values(values, r)
        # lift the local axes to global axes
        lifted = {}
        for exps, coeff in poly.items():
            g = [0] * self.nvars
            for k, e in enumerate(exps):
                g[vs[k]] = e
            lifted[tuple(g)] = coeff
        recip = [v for v in vs if self.vars[v][0] == "t"]
        lifted = poly_reciprocal_vars(lifted, recip, self.nvars)
        return IPoly(lifted, self.nvars, name=f"minor{S}", index_set=tuple(S))


# ---------------------------------------------------------------------------
# Conditions that refute a box
# ---------------------------------------------------------------------------
#
# A labelling is realizable only at a point where BOTH rank Gr <= d+1 and the
# signature is (d,1).  Accordingly there are two ways to refute a box, and the
# branch and bound prunes on either:
#
#   RankCondition     some (d+2)-principal minor is nonzero throughout the box,
#                     so rank Gr >= d+2 there;
#   InertiaCondition  some PRINCIPAL SUBMATRIX has at least two negative
#                     eigenvalues throughout the box, so Gr is superhyperbolic
#                     there.  The number of negative eigenvalues of a principal
#                     submatrix never exceeds that of the whole matrix, so two
#                     negatives anywhere below rule out signature (d,1).
#
# The wildcard labellings of Section 5.4 all fall to the first kind alone.  The
# second is what the non-wildcard labellings need: their weights are DETERMINED
# by the rank condition, so rank-deficient points of the domain generally do
# exist and the rejection is about the signature at those points, not about rank.
#
# Inertia is read off by Jacobi's rule: if the leading principal minors
# D_1,...,D_k of a symmetric matrix are all nonzero, the number of negative
# eigenvalues equals the number of sign changes in 1, D_1, ..., D_k.  Each D_i is
# one of our exact minor polynomials, so "all nonzero on this box" is decided by
# the same outward-rounded interval evaluation, and the sign sequence with it.
# Because the polynomials carry a factor t_e^2 >= 0, their signs agree with the
# determinants' on t_e > 0; at t_e = 0 the factor does not vanish the polynomial
# but evaluates it at the leading coefficient in x_e, which is the sign as
# x_e -> infinity, so the closed box is handled without a separate argument.


def _sign(iv):
    if iv.lo > 0:
        return 1
    if iv.hi < 0:
        return -1
    return 0


class RankCondition:
    __slots__ = ("poly", "name")

    def __init__(self, poly):
        self.poly = poly
        self.name = f"rank{list(poly.index_set)}"

    def refutes(self, box):
        return not self.poly.range(box).contains_zero()

    def report(self, box):
        r = self.poly.range(box)
        return {"kind": "rank", "minor": list(self.poly.index_set),
                "range": [str(r.lo), str(r.hi)]}


class InertiaCondition:
    """At least `need` negative eigenvalues in a principal submatrix."""

    __slots__ = ("polys", "order", "need", "name")

    def __init__(self, polys, order, need=2):
        self.polys = polys              # leading minors, sizes 1..k, in order
        self.order = tuple(order)
        self.need = need
        self.name = f"inertia{list(order)}"

    def _negatives(self, box):
        signs = [1]
        for p in self.polys:
            s = _sign(p.range(box))
            if s == 0:
                return None             # a vanishing leading minor: Jacobi's rule
            signs.append(s)             # does not apply, so claim nothing
        return sum(1 for a, b in zip(signs, signs[1:]) if a != b)

    def refutes(self, box):
        neg = self._negatives(box)
        return neg is not None and neg >= self.need

    def report(self, box):
        return {"kind": "inertia", "submatrix": list(self.order),
                "negative_eigenvalues": self._negatives(box),
                "leading_minor_ranges": [[str(p.range(box).lo),
                                          str(p.range(box).hi)]
                                         for p in self.polys]}


# ---------------------------------------------------------------------------
# Branch and bound
# ---------------------------------------------------------------------------

def certify_empty(conditions, box, max_boxes=200_000,
                  min_width=Fraction(1, 2 ** 24)):
    """Try to prove that no point of `box` is realizable.

    `conditions` are RankCondition / InertiaCondition objects (or bare IPolys,
    treated as rank conditions); a box is discarded as soon as ONE of them
    refutes it.  Returns (certified: bool, stats: dict).  `certified` is True
    only when every box in the subdivision was discarded by such a refutation,
    which is a proof.  False means undecided (box budget or minimum width
    reached) and never asserts feasibility.
    """
    conditions = [c if hasattr(c, "refutes") else RankCondition(c)
                  for c in conditions]
    stack = [box]
    boxes = 0
    pruned_by = {}
    max_depth_width = None
    while stack:
        B = stack.pop()
        boxes += 1
        if boxes > max_boxes:
            return False, {"boxes": boxes, "reason": "box budget",
                           "pruned_by": pruned_by}
        killed = None
        for c in conditions:
            if c.refutes(B):
                killed = c.name
                break
        if killed is not None:
            pruned_by[killed] = pruned_by.get(killed, 0) + 1
            continue
        # cannot exclude: split the widest axis
        widths = [b.width() for b in B]
        ax = max(range(len(B)), key=lambda i: widths[i])
        if widths[ax] < min_width:
            return False, {"boxes": boxes, "reason": "min width",
                           "box": [(str(b.lo), str(b.hi)) for b in B],
                           "pruned_by": pruned_by}
        mid = B[ax].mid()
        left = list(B)
        right = list(B)
        left[ax] = IV(B[ax].lo, mid)
        right[ax] = IV(mid, B[ax].hi)
        stack.append(left)
        stack.append(right)
        max_depth_width = widths[ax]
    return True, {"boxes": boxes, "pruned_by": pruned_by,
                  "last_width": str(max_depth_width) if max_depth_width else None}


def certify_instance(rec, max_minors=12, max_boxes=200_000, escalate=(1, 2, 4, 8, 12)):
    """Certify one dumped wildcard instance infeasible, or report undecided.

    Minors are added to the constraint set cheapest-first (fewest unknowns) and
    the branch-and-bound is retried after each escalation, so an instance that a
    single minor refutes costs a single determinant interpolation.

    Returns (certified, stats).  On success stats["witness"] records enough to
    re-check the certificate by hand or by an independent implementation: which
    principal minors were used, and either the exact nonzero value of a
    label-only minor or the rational range enclosure that excludes zero.
    """
    inst = Instance(rec)
    box = inst.box()
    sets = inst.index_sets()
    polys = []
    stats = {"nvars": inst.nvars, "n_dashed": len(inst.dashed),
             "n_wild": len(inst.wild), "residual": inst.residual,
             "minors_used": 0, "boxes": 0}
    seen = 0
    for target in escalate:
        if target > max_minors:
            break
        while len(polys) < target and seen < len(sets):
            S = sets[seen]
            seen += 1
            p = inst.minor_poly(S)
            if not p:
                continue                      # identically zero: no information
            if p.is_nonzero_constant():
                # The strongest certificate available: this 8x8 principal minor
                # involves no unknown at all, so the ordinary labels alone force
                # rank >= 8 and neither the wild entries nor the dashed weights
                # can rescue the labelling.
                exps, iv = p.terms[0]
                val = p.exact[exps]
                stats.update(minors_used=len(polys) + 1, boxes=1,
                             certificate=f"constant nonzero minor {S}",
                             witness={"kind": "label_only_minor",
                                      "minor": list(S),
                                      "value_exact": [str(c) for c in val],
                                      "value_enclosure": [str(iv.lo), str(iv.hi)],
                                      "value_float": EF.alg_float(val)})
                return True, stats
            polys.append(p)
        if not polys:
            break
        ok, bb = certify_empty(polys, box, max_boxes=max_boxes)
        stats["minors_used"] = len(polys)
        stats["boxes"] = bb.get("boxes", 0)
        stats["pruned_by"] = bb.get("pruned_by")
        if ok:
            stats["certificate"] = f"{len(polys)} minors, {bb['boxes']} boxes"
            root = [(p.index_set, p.range(box)) for p in polys]
            stats["witness"] = {
                "kind": ("root_box_range" if bb["boxes"] == 1
                         else "interval_branch_and_bound"),
                "minors": [list(p.index_set) for p in polys],
                "boxes": bb["boxes"],
                "root_ranges": {str(list(S)): [str(r.lo), str(r.hi)]
                                for S, r in root},
                "box": [[str(b.lo), str(b.hi)] for b in box],
            }
            return True, stats
        stats["reason"] = bb.get("reason")
        if seen >= len(sets):
            break
    return False, stats


def certify_instance_full(rec, max_minors=8, max_inertia=40, max_boxes=200_000):
    """Certify a labelling unrealizable using BOTH refutation kinds.

    Used for the labellings that carry no wildcard.  There the dashed weights are
    determined by the rank condition, so rank-deficient points of the domain
    usually do exist and no rank minor can refute the labelling; what fails at
    those points is the signature.  Superhyperbolicity is therefore tested
    alongside: a principal submatrix with two negative eigenvalues throughout a
    box refutes it just as well as a nonvanishing $(d+2)$-minor does.

    Escalation order is cheapest-first and stops at the first success:
      1. rank minors alone (this settles a labelling refuted on rank);
      2. each candidate inertia condition ALONE, root box only -- the common
         case, and much cheaper than subdividing;
      3. rank minors together with all inertia conditions, with subdivision.
    """
    inst = Instance(rec)
    box = inst.box()
    stats = {"nvars": inst.nvars, "n_dashed": len(inst.dashed),
             "n_wild": len(inst.wild), "boxes": 0, "minors_used": 0}

    # 1. rank alone
    ok, s = certify_instance(rec, max_minors=max_minors, max_boxes=max_boxes)
    if ok:
        s["route"] = "rank"
        return True, s

    # 2. one inertia condition at a time, root box only
    conds = []
    for order in inst.inertia_orders()[:max_inertia]:
        c = inst.inertia_condition(order)
        conds.append(c)
        if c.refutes(box):
            stats.update(boxes=1, route="inertia",
                         certificate=f"superhyperbolic on {list(order)}",
                         witness=dict(c.report(box),
                                      kind="superhyperbolic_root_box",
                                      box=[[str(b.lo), str(b.hi)] for b in box]))
            return True, stats

    # 3. everything together, with subdivision
    sets = inst.index_sets()
    polys = [p for p in (inst.minor_poly(S) for S in sets[:max_minors]) if p]
    all_conds = [RankCondition(p) for p in polys] + conds
    if not all_conds:
        return False, stats
    ok, bb = certify_empty(all_conds, box, max_boxes=max_boxes)
    stats["boxes"] = bb.get("boxes", 0)
    stats["pruned_by"] = bb.get("pruned_by")
    stats["minors_used"] = len(polys)
    if ok:
        stats.update(route="mixed",
                     certificate=f"{len(all_conds)} conditions, "
                                 f"{bb['boxes']} boxes",
                     witness={"kind": "mixed_branch_and_bound",
                              "conditions": [c.name for c in all_conds],
                              "boxes": bb["boxes"],
                              "pruned_by": bb.get("pruned_by"),
                              "box": [[str(b.lo), str(b.hi)] for b in box]})
        return True, stats
    stats["reason"] = bb.get("reason")
    return False, stats
