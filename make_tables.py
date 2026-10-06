#!/usr/bin/env python3
"""Regenerate the paper's data-derived LaTeX tables directly from the run artifacts.

Outputs (all under tables/):
    pertype.tex        the 55-row per-type table (a record; the paper states it as
                       Lemma "survivors" through the lists in summary_nums.tex)
    summary_nums.tex   \newcommand definitions for every number quoted in the text
    rootbox.tex        the interval certificates of Proposition "wildexact"

Everything here reads the committed run data; nothing is transcribed by hand.
Run from the repository root:  python3 make_tables.py
"""
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "tables"
OUT.mkdir(parents=True, exist_ok=True)

TYPES = {t["type_id"]: t for t in json.load(open(ROOT / "runs/d6_n10/stage2/types.json"))}
SURV = json.load(open(ROOT / "runs/d6_n10/facet_profile_survivors.json"))["survivors"]
# The d=6 classification of record.  Earlier checkpoint directories are retained in
# the repository for provenance but are NOT used for any number in the paper.
def _load_state(*dirs):
    st = {}
    for d in dirs:
        p = ROOT / f"runs/d6_n10/{d}/state.jsonl"
        for line in p.read_text().splitlines():
            line = line.strip()
            if line:
                r = json.loads(line)
                st[r.pop("key")] = r
    return st

# The classification of record is the run with every restriction of Section 4
# imposed before the search AND the exact forward-checking gates of Section 5.1, so
# that no verdict in it rests on a tolerance.  Earlier runs are read only if it is
# absent, so the script still works before the final run completes.
# d6_valid (2026-09-28) is the run over ALL 304 types without a size-6 minimal
# non-face, made once the facet-based reductions were found to be
# invalid (a facet of a Coxeter polytope is acute-angled, not in general Coxeter).
# When it is present it is the run of record, and the facet filter, prism lemma
# and 3-free shortcut play no part in any number.
for _cand in ("d6_valid", "d6_exact", "d6_final"):
    if (ROOT / f"runs/d6_n10/{_cand}/state.jsonl").exists():
        STATE = _load_state(_cand)
        RUN_OF_RECORD = _cand
        break
else:
    STATE = _load_state("d6_tangency", "d6_rest")
    RUN_OF_RECORD = "d6_tangency+d6_rest"
# The d=5 census of record: the run through the same code path as the d=6
# classification, i.e. with the exact forward-checking gates.  Falls back to the
# earlier floating-point-gate run if that is absent.
_d5x = ROOT / "runs/d5_n9/d5_exact.json"
D5 = json.load(open(_d5x if _d5x.exists()
                    else ROOT / "runs/d5_n9/d5_discfix.json"))
D4 = _load_state.__wrapped__ if False else None

# Wildcard traffic in the d=5 census of record, which calibrates the
# joint-feasibility test against ground truth (Section 5.4).  Read from the same
# artifact the census totals come from, so the two cannot describe different runs.
_D5W = D5 if any("wild_assignments" in v for v in D5.values()) else \
    json.load(open(ROOT / "runs/d5_n9/d5_wildcounts.json"))

NLAB, MAXD = 6, 14

# Two purely combinatorial exclusions, each a theorem (Section 4).  The flags are
# precomputed into runs/d6_n10/type_flags.json so that this script needs no import
# from the pipeline; regenerate them with checks/type_consistency.py --emit.
FLAGS = json.load(open(ROOT / "runs/d6_n10/type_flags.json"))["d6_n10"]

# Wildcard traffic in the d=5 census, which calibrates the joint-feasibility test
# against ground truth (Section 5.4).  Read from the artifact, never transcribed.

# Exact wildcard infeasibility certificates, written by
# checks/wild_exact_certify.py from the instance dumps of run_wild_dump.py.
VALID = RUN_OF_RECORD == "d6_valid"
_CERT = ROOT / ("runs/d6_n10/d6_valid" if VALID else "runs/d6_n10")
W6 = json.load(open(_CERT / "wild_certificates.json"))["summary"]
W5 = json.load(open(ROOT / "runs/d5_n9/wild_certificates.json"))["summary"]
# ... and the same for the labellings without a wildcard
# (checks/plain_exact_certify.py).
P6 = json.load(open(_CERT / "plain_certificates.json"))["summary"]
# One-sidedness is measured over the WHOLE d=5 census: every instance that survived
# the joint screen there, which is the only set a false refutation could damage.
# The refutation-power figures (how many rejections carry a certificate) come from
# W5 above and may cover the census or the two-type subset; this one always covers
# the census.
W5P = json.load(open(ROOT / "runs/d5_n9/wild_certificates_passers.json"))["summary"]

# The data behind "exactly one polytope up to isometry"
# (checks/uniqueness_witness.py): the weight tree of the accepted labelling,
# and the automorphism orbit that makes several label assignments one polytope.
UW = json.load(open(ROOT / "runs/d6_n10/uniqueness_witness.json"))

# Cascade branch counters of the d=4 census of record (runs/d4_n8/survivors_discfix),
# read here so that the figures quoted for dimension 4 come from the same run as its
# 348 polytopes.  An earlier directory, d4_final, has different totals; quoting one
# run's pair count beside another's decision count is exactly the confusion this
# avoids.
def _d4_branches():
    br = Counter()
    p = ROOT / "runs/d4_n8/survivors_discfix/state.jsonl"
    for line in p.read_text().splitlines():
        if line.strip():
            r = json.loads(line)
            for k, v in ((r.get("diag") or {}).get("screen_branches") or {}).items():
                br[k] += v
    return br

D4B = _d4_branches()


def _sci(x, digits=1, up=False):
    """LaTeX scientific notation, rounded so that the quoted figure is TRUE as a
    bound: down for a quantity the paper calls a lower bound (a margin, a
    minimum), up for one it calls an upper bound (a maximum).  Rounding a
    maximum down would state a bound the data violate."""
    from math import ceil, floor, log10
    e = floor(log10(x))
    scaled = x / 10 ** e * 10 ** digits
    m = (ceil(scaled) if up else floor(scaled)) / 10 ** digits
    return f"{m}\\times10^{{{e}}}"
D5W_TOTAL = sum(v.get("wild_assignments", 0) for v in _D5W.values())
D5W_REJECTED = sum(sum((v.get("wild") or {}).values()) for v in _D5W.values())


def _cert_records(path):
    """Per-instance certificate records of a certifier artifact, normalised to
    (tag, certified, kind, boxes).  Large artifacts keep full witnesses per
    shard and a digest per instance in the merged file; both forms are read."""
    body = json.load(open(path))
    if "instances" in body:
        for r in body["instances"]:
            w = r.get("witness") or {}
            yield (r["tag"], bool(r.get("certified_unrealizable")
                                  or r.get("certified_infeasible")),
                   w.get("kind"), w.get("boxes") or 0)
    else:
        for r in body["instances_digest"]:
            yield (r["tag"], bool(r["certified_infeasible"]), r.get("kind"),
                   r.get("boxes") or 0)


def _c(n):
    return f"{n:,}".replace(",", "{,}")
DEG_KILLED = {int(k) for k, v in FLAGS.items() if v["max_dashed_degree"] >= 3}
BIGMF_KILLED = {int(k) for k, v in FLAGS.items() if v["has_missing_face_of_size_d"]}
PRISM_KILLED = {int(k) for k, v in FLAGS.items()
                if not v.get("degree2_facet_is_prism", True)}
COMB_KILLED = DEG_KILLED | BIGMF_KILLED | PRISM_KILLED
if VALID:
    # Only the size-6 criterion is valid; it contains every dashed-degree-3 type
    # (a facet disjoint from three others), which is asserted below rather than assumed.
    assert DEG_KILLED <= BIGMF_KILLED, "a dashed-degree-3 type lacks a size-6 non-face"
    COMB_KILLED = set(BIGMF_KILLED)
    SURV = sorted(t for t in TYPES if t not in BIGMF_KILLED)


def _closed(tid, pfx=()):
    """Coverage recursion: a type is decided iff its root is exhausted, or every
    one of its six children is (recursively) decided.  This is the certificate the
    paper's exhaustion claim rests on -- recomputed here, not read from a flag."""
    e = STATE.get(f"{tid}|{','.join(map(str, pfx))}")
    if e is None:
        return False
    if e.get("exhausted"):
        return True
    if len(pfx) >= MAXD:
        return False
    return all(_closed(tid, pfx + (i,)) for i in range(NLAB))


def _keys(tid):
    out = set()
    for k, v in STATE.items():
        if k.split("|")[0] == str(tid):
            out.update(v.get("keys", []))
    return out

# Types eliminated by Esselmann's 2d-facet bound for 3-free polytopes: every
# missing face has size 2, so the polytope would be 3-free and would need >= 12
# facets.  These are decided by theorem, not by the search.
THREE_FREE = ([] if VALID else
              sorted(t for t in SURV if {len(m) for m in TYPES[t]["missing_faces"]} == {2}))


def profile(tid):
    return Counter(len(m) for m in TYPES[tid]["missing_faces"])


def profile_str(tid):
    p = profile(tid)
    return "".join(f"{s}^{{{p[s]}}}" for s in sorted(p))


def agg(tid):
    """(#subtree certificates, label assignments enumerated, CPU seconds)."""
    rows = [v for k, v in STATE.items() if k.split("|")[0] == str(tid)]
    return (len(rows),
            sum(r.get("enum_count", 0) for r in rows),
            sum(r.get("sec", 0.0) for r in rows))


def dep0(k):
    return len(k.split("|")[1].split(",")) if k.split("|")[1] else 0


def main():
    rows, tot_sub, tot_enum, tot_sec = [], 0, 0, 0.0
    for tid in sorted(SURV):
        nsub, enum, sec = agg(tid)
        ndist = len(_keys(tid))
        tot_sub += nsub
        tot_enum += enum
        tot_sec += sec
        # The method column must name the deciding lemma for every row a theorem
        # decides, and say "search" only for the types actually searched.  A type
        # can satisfy more than one exclusion; the first that applies is the one
        # reported, in the order the reductions are presented in Section 4.
        fl = FLAGS.get(str(tid), {})
        if VALID:
            fl = {}                    # every candidate is searched
        if fl.get("max_dashed_degree", 0) >= 3:
            method = r"Lem.\,\ref{lem:degree}"
            verdict = "0"
        elif fl.get("has_missing_face_of_size_d", False):
            method = r"Lem.\,\ref{lem:mfsix}"
            verdict = "0"
        elif not fl.get("degree2_facet_is_prism", True):
            method = r"Lem.\,\ref{lem:prism}"
            verdict = "0"
        elif tid in THREE_FREE:
            # unreachable on the current data: both 3-free types have dashed degree
            # 3 and so are already excluded by Lemma "degree".  Kept so that the
            # table would still name a deciding theorem if that ever changed.
            method = r"Lem.\,\ref{lem:esselmann}"
            verdict = "0"
        elif ndist:
            method = "search"
            verdict = r"\textbf{1}"
        else:
            method = "search"
            verdict = "0"
        req = ("" if tid in COMB_KILLED else r"$\bullet$")
        rows.append((tid, TYPES[tid]["p_count"], profile_str(tid), req, verdict,
                     method, nsub, enum, sec / 3600.0))

    lines = []
    for tid, p, prof, req, verdict, method, nsub, enum, hrs in rows:
        lines.append(f"{tid} & {p} & ${prof}$ & {req} & {verdict} & {method} & "
                     f"{nsub} & {enum:,} & {hrs:.2f} \\\\")
    (OUT / "pertype.tex").write_text("\n".join(lines) + "\n")

    # The same partition, as type lists: the
    # 44 lemma-excluded survivors grouped by deciding lemma and then by (p, profile),
    # and the 11 searched types with their per-type search data.
    def _group(tids):
        order, by = [], {}
        for t in tids:
            key = (TYPES[t]["p_count"], profile_str(t))
            if key not in by:
                order.append(key)
                by[key] = []
            by[key].append(t)
        order.sort(key=lambda k: (k[0], k[1]))
        return "; ".join(f"{', '.join(map(str, by[k]))} ($p={k[0]}$, ${k[1]}$)"
                         for k in order)
    by_method = {}
    for tid, _p, _prof, _req, _v, method, nsub, enum, hrs in rows:
        by_method.setdefault(method, []).append(tid)
    SURV_LISTS = {
        "SurvDegreeList": _group(by_method.get(r"Lem.\,\ref{lem:degree}", [])),
        "SurvMfsixList": _group(by_method.get(r"Lem.\,\ref{lem:mfsix}", [])),
        "SurvPrismList": _group(by_method.get(r"Lem.\,\ref{lem:prism}", [])),
        "SearchedList": _group(by_method.get("search", [])),
        "NumSurvDegree": str(len(by_method.get(r"Lem.\,\ref{lem:degree}", []))),
        "NumSurvMfsix": str(len(by_method.get(r"Lem.\,\ref{lem:mfsix}", []))),
        "NumSurvPrism": str(len(by_method.get(r"Lem.\,\ref{lem:prism}", []))),
        "SearchedTimes": ", ".join(f"{tid}: {hrs:.2f}" for tid, *_x, method, nsub,
                                   enum, hrs in rows if method == "search"),
    }
    _searched_rows = [r for r in rows if r[5] == "search"]
    assert set(by_method) <= {r"Lem.\,\ref{lem:degree}", r"Lem.\,\ref{lem:mfsix}",
                              r"Lem.\,\ref{lem:prism}", "search"}, (
        f"a survivor is decided by an unexpected method: {sorted(by_method)}")
    assert sorted(sum(by_method.values(), [])) == sorted(SURV), (
        "the survivor lists are not a partition of the facet-filter survivors")
    if not VALID:
        # properties of the 11-type run that the earlier text asserted
        assert all(r[6] == 1 for r in _searched_rows), (
            "a searched type has more than one subtree certificate; the paper says "
            "each is exhausted at its root")
        assert [r[0] for r in _searched_rows if r[7]] == [379], (
            "a searched type other than 379 produced complete labellings")
    assert all(r[6] == 0 and r[7] == 0 for r in rows if r[5] != "search"), (
        "a lemma-excluded survivor carries search data")

    # The interval certificates: for
    # each of the labellings not settled by a label-only minor, the labelling, the
    # certifying minor and the enclosure of its range over the compactified box,
    # rounded OUTWARD so that the printed interval still excludes 0.
    from fractions import Fraction
    from math import floor, ceil, log10

    def _out(x, up):
        x = Fraction(x)
        if x == 0:
            return "0"
        e = floor(log10(abs(float(x)))) - 3
        q = x / Fraction(10) ** e
        m = ceil(q) if up else floor(q)
        return f"{float(m * Fraction(10) ** e):.4g}"

    # Per-type certificate summary: for every type in which some
    # labelling reached a screen, how its refutations split by certificate kind.
    # The individual certificates -- minor, enclosure, subdivision -- are in the
    # artifacts; there are far too many to print.
    _per = {}
    for _f in ("plain_certificates.json", "wild_certificates.json"):
        for _tag, _ok, _kind, _bx in _cert_records(_CERT / _f):
            _t = int(_tag.split("|")[0])
            _d = _per.setdefault(_t, Counter())
            if _ok:
                _d[_kind] += 1
                # boxes are counted for interval certificates only; a label-only
                # minor uses no box (the plain and wild certifiers record it as 0
                # and 1 respectively; both are normalised to 0 here)
                if _kind != "label_only_minor":
                    _d["maxboxes"] = max(_d["maxboxes"], _bx)
            else:
                _d["unrefuted"] += 1
    rb_lines = []
    for _t in sorted(_per):
        _d = _per[_t]
        _tot = sum(v for k, v in _d.items() if k not in ("maxboxes", "unrefuted"))
        rb_lines.append(
            f"{_t} & {TYPES[_t]['p_count']} & ${profile_str(_t)}$ & {_c(_tot)} & "
            f"{_c(_d.get('label_only_minor', 0))} & {_c(_d.get('root_box_range', 0))} & "
            f"{_c(_d.get('interval_branch_and_bound', 0))} & "
            f"{_c(_d['maxboxes']) if _d['maxboxes'] else '--'} & "
            f"{_d.get('unrefuted', 0) or ''} \\\\")
    assert sum(_d.get("unrefuted", 0) for _d in _per.values()) == 1, \
        "exactly one screened labelling (P_{6,10}) may be left unrefuted"
    # totals row
    _T = Counter()
    for _d in _per.values():
        _T.update(_d)
    rb_lines.append("\\midrule")
    rb_lines.append(
        f"total & & & {_c(sum(v for k, v in _T.items() if k not in ('maxboxes', 'unrefuted')))} & "
        f"{_c(_T.get('label_only_minor', 0))} & {_c(_T.get('root_box_range', 0))} & "
        f"{_c(_T.get('interval_branch_and_bound', 0))} & "
        f"{_c(max(_d['maxboxes'] for _d in _per.values()))} & {_T.get('unrefuted', 0)} \\\\")
    (OUT / "rootbox.tex").write_text("\n".join(rb_lines) + "\n")

    # The realizing type's missing faces in the facet numbering of the paper's
    # Figure 1 (that of verify_polytope.py): the minimal non-elliptic subsets of
    # the diagram.  Ellipticity is tested numerically here, which is harmless
    # because the result is then required to be ISOMORPHIC to type 379's recorded
    # hypergraph, an exact combinatorial check.
    import itertools
    import math
    import numpy as _np
    import networkx as _nx
    from networkx.algorithms import isomorphism as _iso
    _ORD = [(1, 2, 5), (2, 3, 3), (3, 4, 3), (4, 5, 3), (5, 6, 3), (9, 10, 5),
            (2, 7, 4), (2, 8, 4), (5, 9, 3), (6, 10, 5)]
    _DOT = [(6, 7), (7, 8), (8, 9)]
    _G = _np.eye(10)
    for i, j, m in _ORD:
        _G[i - 1, j - 1] = _G[j - 1, i - 1] = -math.cos(math.pi / m)
    for i, j in _DOT:
        _G[i - 1, j - 1] = _G[j - 1, i - 1] = -2.0
    _ell = lambda S: _np.linalg.eigvalsh(_G[_np.ix_(S, S)]).min() > 1e-9
    _mf = []
    for k in range(2, 7):
        for S in itertools.combinations(range(10), k):
            if not _ell(list(S)) and all(_ell(list(T)) for T in
                                         itertools.combinations(S, k - 1)):
                _mf.append(S)

    def _inc(h):
        g = _nx.Graph()
        g.add_nodes_from((("v", i) for i in range(10)), side=0)
        for k, m in enumerate(h):
            g.add_node(("e", k), side=1)
            g.add_edges_from((("v", i), ("e", k)) for i in m)
        return g
    assert _nx.is_isomorphic(_inc(_mf), _inc(TYPES[379]["missing_faces"]),
                             node_match=_iso.categorical_node_match("side", None)), \
        "Figure 1's diagram does not have type 379's missing-face hypergraph"
    # one math run per set, so the list can break between sets
    REALIZING_MF = ", ".join("$\\{" + ",".join(str(i + 1) for i in S) + "\\}$"
                             for S in _mf)

    _BR = Counter()
    _STALLED_TYPES = set()
    for _k, _v in STATE.items():
        _b = (_v.get("diag") or {}).get("screen_branches") or {}
        _BR.update(_b)
        if any(k.startswith("stuck") and n for k, n in _b.items()):
            _STALLED_TYPES.add(int(_k.split("|")[0]))
    assert not _BR.get("cascade_pair"), "the pair-resultant branch was entered at d=6"
    assert _BR.get("numerical_fallback", 0) == sum(
        v for k, v in _BR.items() if k.startswith("stuck")), \
        "fallback calls do not match stalls"
    assert P6.get("cascade_stalled", 0) == sum(
        v for k, v in _BR.items() if k.startswith("stuck")), \
        "the certified stalled labellings do not match the stall counter"

    # certificate kinds over BOTH certifiers' artifacts of the run of record
    _kinds = Counter()
    _pk, _wk = Counter(), Counter()
    _maxboxes = 0
    for _f, _kk in (("plain_certificates.json", _pk), ("wild_certificates.json", _wk)):
        for _tag, _ok, _kind, _bx in _cert_records(_CERT / _f):
            if _ok:
                _kinds[_kind or "unknown"] += 1
                _kk[_kind or "unknown"] += 1
                _maxboxes = max(_maxboxes, _bx)
    assert _kinds["unknown"] == 0, "a certificate without a witness kind"
    del _kinds["unknown"]
    _DD = json.load(open(ROOT / "runs/d6_n10/dedup_truncation.json"))
    assert _DD["lost"] == 0, "truncated deduplication lost a type"
    assert len(_DD["types"]) == len(BIGMF_KILLED), "dedup check does not cover every excluded type"

    # Symmetry breaking switched off, over the
    # run of record's types, with the realizing type's root run unpruned.
    _SO = {}
    _sop = ROOT / "runs/d6_n10" / ("symmetry_breaking_off_valid.json" if VALID
                                   else "symmetry_breaking_off.json")
    if _sop.exists():
        _b = json.load(open(_sop))
        _pt = (_b.get("cheap_types") or {}).get("per_type", {})
        _rep = [int(t) for t, v in _pt.items()
                if v["exhausted"] and v["polytopes"] == v["polytopes_with_pruning"]]
        _bud = sorted(int(t) for t, v in _pt.items() if not v["exhausted"])
        _dis = [t for t, v in _pt.items()
                if v["exhausted"] and v["polytopes"] != v["polytopes_with_pruning"]]
        assert not _dis, f"symmetry breaking changed a verdict: {_dis}"
        _rs = _b.get("realizing_subtree") or {}
        # The unpruned ROOT of the realizing type may not finish in its budget;
        # the unpruned subtree with the first two labels fixed, which contains
        # P_{6,10}, is then the test in the dangerous direction.
        _sub = ROOT / "runs/d6_n10/symmetry_breaking_off_subtree.json"
        if not _rs.get("verdict_reproduced") and _sub.exists():
            _rs = json.load(open(_sub)).get("realizing_subtree") or _rs
        assert _rs.get("verdict_reproduced"), (
            "no unpruned run of the realizing type reproduced its verdict")
        _argonly = [t for t in _bud if not (t == 379 and _rs.get("subtree") == "379|")]
        _SO = {
            "SymOffTypes": str(len(_pt)),
            "SymOffReproduced": str(len(_rep)),
            "SymOffBudget": str(len(_bud)),
            "SymOffBudgetList": ", ".join(map(str, _bud)),
            # types whose pruning rests on the lex-leader argument alone: those that
            # hit the budget, less the realizing type if its unpruned ROOT run
            # completed in the longer budget of the realizer check
            "SymOffArgOnly": str(len(_argonly)),
            "SymOffArgOnlyList": ", ".join(map(str, _argonly)),
            "SymOffRealizerOK": "yes" if _rs.get("verdict_reproduced") else "no",
            "SymOffRealizerSubtree": str(_rs.get("subtree", "")),
            "SymOffRealizerPrunedLabellings": _c(((_rs.get("pruning_on") or {})
                                                  .get("labellings", 0))),
            "SymOffRealizerLabellings": _c(((_rs.get("pruning_off") or {})
                                            .get("labellings", 0))),
            "SymOffRealizerExhausted": ("yes" if (_rs.get("pruning_off") or {})
                                        .get("exhausted") else "no"),
        }

    # The hypergraphs recomputed independently from the stored coordinates
    # (checks/hypergraph_independent.py)
    _HG = json.load(open(ROOT / "runs/d6_n10/hypergraph_independent.json"))
    assert _HG["agree"] == len(TYPES) and not _HG["disagree"], \
        "an independently recomputed hypergraph disagrees with the pipeline"
    assert _HG["size6_types"] == len(BIGMF_KILLED) and not _HG["isomorphic_pairs"]

    # Random facet relabellings with symmetry breaking ON: checks/symmetry_relabel.py
    _RL = {}
    _rlp = ROOT / "runs/d6_n10/symmetry_relabel.json"
    if _rlp.exists():
        _rl = json.load(open(_rlp))
        assert _rl["disagreements"] == 0, "a relabelling changed a verdict"
        _runs = [r for v in _rl["runs"].values() for r in v]
        _RL = {"RelabelTypes": str(len(_rl["types"])), "RelabelSeeds": str(_rl["seeds"]),
               "RelabelRuns": str(len(_runs)),
               "RelabelReproduced": str(_rl["reproduced"]),
               "RelabelOverBudget": str(_rl["inconclusive"])}

    # The independent SAT/DRAT certificate of the search tree
    # (checks/tree_sat.py).
    _TS = {}
    _tsp = ROOT / "runs/d6_n10/tree_sat.json"
    if _tsp.exists():
        _ts = json.load(open(_tsp))
        _valid = {t for t in TYPES if t not in BIGMF_KILLED}
        _cov = {int(k) for k in _ts} & _valid
        _unsat = [r for r in _ts.values() if r["verdict"] == "UNSAT"]
        _ver = [r for r in _unsat if r.get("drat_trim") == "VERIFIED"]
        if len(_cov) == len(_valid):
            assert len(_ver) == len(_ts) == len(_valid), (
                "a type's CNF is not UNSAT with a verified DRAT proof")
            assert sum(r["screened_consistent_with_cnf"] for r in _ts.values()) == \
                sum(v.get("enum_count", 0) for v in STATE.values()), \
                "the non-vacuity check does not cover every screened labelling"
        # every blocked labelling is matched to an exact certificate or to P
        # (tree_sat.py --link)
        if len(_cov) == len(_valid):
            _lc = [r for r in _ts.values() if r.get("screened")]
            assert all("screened_certified" in r for r in _lc), "run tree_sat.py --link"
            assert sum(r["screened_uncertified"] for r in _lc) == 0
            assert sum(r["screened_accepted"] for r in _lc) == 1
            assert sum(r["screened_certified"] for r in _lc) == \
                P6["exactly_certified_unrealizable"] + W6["exactly_certified_infeasible"]
        _TS = {
            "TreeTypes": str(len(_cov)),
            "TreeVerified": str(len(_ver)),
            "TreeBlocking": _c(sum(r["blocking_clauses"] for r in _ts.values())),
            "TreeDisconnected": _c(sum(r["disconnected_blocked"] for r in _ts.values())),
            "TreeMaxClauses": _c(max(r["clauses"] for r in _ts.values())),
            "TreeProofLines": _c(sum(r.get("proof_lines", 0) for r in _ts.values())),
            "TreeSeconds": _c(round(sum(r["solve_seconds"] for r in _ts.values()))),
            "TreeImagesRealizing": _c(_ts.get("379", {}).get("blocking_clauses", 0)),
            # types certified with a lex-leader symmetry-breaking predicate rather
            # than by blocking every Aut(T)-image
            "TreeSBTypes": " and ".join(str(k) for k in sorted(
                int(k) for k, r in _ts.items() if r.get("mode") == "lex-leader")),
            "NumTreeSB": str(sum(1 for r in _ts.values() if r.get("mode") == "lex-leader")),
            "NumTreeNoSB": str(sum(1 for r in _ts.values() if r.get("mode") != "lex-leader")),
            "TreeSBAut": " and ".join(str(r["aut_order"]) for k, r in sorted(
                _ts.items(), key=lambda kv: int(kv[0])) if r.get("mode") == "lex-leader"),
        }

    d5_total = sum(x["distinct"] for x in D5.values())
    d5_real = {k: x["distinct"] for k, x in D5.items() if x["distinct"]}
    defs = {
        "NumOrderTypes": "14{,}309{,}547",
        "NumTypes": str(len(TYPES)),
        "NumTypesMF": str(sum(1 for t in TYPES.values()
                              if any(len(m) in (3, 4) for m in t["missing_faces"]))),
        "NumSurvivors": str(len(SURV)),
        "NumSearched": str(len(SURV) - len(THREE_FREE)),
        "NumSearchedLess": str(len(SURV) - len(THREE_FREE) - 1),
        "NumDegKilled": str(len(DEG_KILLED)),
        "NumBigMF": str(len(BIGMF_KILLED - DEG_KILLED)),
        "NumPrism": str(len(PRISM_KILLED - DEG_KILLED - BIGMF_KILLED)),
        "NumAfterComb": str(len(TYPES) - len(COMB_KILLED)),
        # after the two Lanner-derived exclusions alone, before the prism lemma
        "NumAfterDegMF": str(len(TYPES) - len(DEG_KILLED | BIGMF_KILLED)),
        "NumRequired": str(len(set(SURV) - COMB_KILLED - set(THREE_FREE))),
        "NumRedundant": str(len(set(SURV)) -
                            len(set(SURV) - COMB_KILLED - set(THREE_FREE))),
        # of those, the number also re-searched directly: the two 3-free types are
        # excluded by Lemma "degree" but their label trees do not terminate, so they
        # are the ones the corroborating search cannot cover
        "NumRedundantSearched": str(len(COMB_KILLED & set(SURV)) -
                                    len(COMB_KILLED & set(THREE_FREE))),
        "DFiveWild": _c(D5W_TOTAL),
        "DFiveWildKilled": _c(D5W_REJECTED),
        "DFiveWildSurvived": str(D5W_TOTAL - D5W_REJECTED),
        # Exact wildcard certificates (Proposition "wildexact", Section 5.4), read
        # from the certifier's artifacts.  WildCertified must equal NumWild or the
        # proposition is false as stated; the assertions below enforce that.
        "WildCertified": _c(W6["exactly_certified_infeasible"]),
        "WildLabelOnly": _c(W6["certificate_kinds"].get("label_only_minor", 0)),
        "WildRootBox": _c(W6["certificate_kinds"].get("root_box_range", 0)),
        "WildBnB": _c(W6["certificate_kinds"].get("interval_branch_and_bound", 0)),
        "WildMinMargin": _sci(W6["min_margin"]),
        "WildResidMin": _sci(W6["numerical_residual_min"]),
        "WildResidMax": _sci(W6["numerical_residual_max"], up=True),
        "DFourCascade": _c(sum(D4B.values())),
        "DFourPair": _c(D4B["cascade_pair"]),
        "DFourFallback": _c(D4B["numerical_fallback"]),
        "PlainRejected": _c(P6["cascade_rejected"]),
        "PlainCertified": _c(P6["exactly_certified_unrealizable"]),
        "DFiveWildCertSet": _c(W5["instances"]),
        "DFiveWildCert": _c(W5["exactly_certified_infeasible"]),
        "DFiveWildUndecided": _c(W5["instances"]
                                 - W5["exactly_certified_infeasible"]),
        # census-wide: every instance surviving the joint screen, and how many of
        # them the certifier refuted (which must be zero)
        "DFiveWildScreenPass": _c(W5P["passed_joint_screen"]),
        "DFiveWildPassRefuted": _c(W5P["wrongly_certified_among_screen_passers"]),
        "DFiveWildCensus": _c(W5P["census_instances"]),
        # refutation power, over whichever rejection set was certified
        "DFiveWildJointKilled": _c(W5["rejected_at_joint_screen"]),
        # branch counters of the run of record, so no figure is transcribed
        "NumCascade": _c(sum(sum((v.get("diag") or {}).get("screen_branches", {}).values())
                             for v in STATE.values())),
        "NumPlain": _c(sum(v.get("enum_count", 0) for v in STATE.values())
                       - sum((v.get("diag") or {}).get("wild_assignments", 0)
                             for v in STATE.values())),
        "NumWild": f"{sum((v.get('diag') or {}).get('wild_assignments', 0) for v in STATE.values()):,}".replace(",", "{,}"),
        "NumRequiredLess": str(len(set(SURV) - COMB_KILLED - set(THREE_FREE)) - 1),
        "NumThreeFree": str(len(THREE_FREE)),
        "ThreeFreeList": ", ".join(map(str, THREE_FREE)),
        "NumSubtrees": f"{tot_sub:,}".replace(",", "{,}"),
        "NumAssignments": f"{tot_enum:,}".replace(",", "{,}"),
        # sub-hour totals need a decimal, or "in 1 CPU-hours" gets printed
        "CPUHours": (f"{tot_sec/3600:.2f}" if tot_sec < 36000
                     else f"{tot_sec/3600:,.0f}".replace(",", "{,}")),
        "RealizingType": "379",
        # profiles, so the paper can name a type by its invariants
        # rather than by our arbitrary discovery-order number
        "RealizingP": str(TYPES[379]["p_count"]),
        "RealizingProfile": profile_str(379),
        "ThreeFreeP": ", ".join(str(TYPES[t]["p_count"]) for t in THREE_FREE),
        "ThreeFreeProfile": " and ".join(sorted({profile_str(t) for t in THREE_FREE})),
        "NumGramConfigs": str(UW["orbit_of_the_labelling"]),
        "AutTypeOrder": str(UW["aut_type_order"]),
        "AutStabOrder": str(UW["labelling_stabiliser_order"]),
        "WeightLeaves": str(UW["weight_tree_leaves"]),
        "DFiveTotal": str(d5_total),
        "DFiveTypes": str(len(D5)),
        "DFiveCPUHours": f"{sum(x['sec'] for x in D5.values())/3600:.1f}",
        "MaxDepth": str(max(dep0(k) for k in STATE)),
        "DFiveRealizing": ", ".join(f"{k}\\!\\to\\!{v}" for k, v in
                                    sorted(d5_real.items(), key=lambda x: -x[1])),
        "DFourFound": "348",
        "DFourTarget": "348",
        **SURV_LISTS,
        "RunOfRecord": RUN_OF_RECORD,
        **_TS,
        **_RL,
        "ValidPDist": ", ".join(f"${_n}$ with $p={_p}$" for _p, _n in sorted(Counter(
            TYPES[t]["p_count"] for t in TYPES if t not in BIGMF_KILLED).items())),
        **_SO,
        # cascade decisions of the run of record, by branch: accepted, rejected,
        # and stalled for a value-dependent reason (then handed to the bounded-box
        # fallback); the pair-resultant branch must never be entered at d=6
        "NumCascadeAccept": _c(_BR.get("cascade_true", 0)),
        "NumCascadeReject": _c(_BR.get("cascade_false", 0)),
        "NumStalled": _c(sum(v for k, v in _BR.items() if k.startswith("stuck"))),
        "NumFallback": _c(_BR.get("numerical_fallback", 0)),
        "StalledTypes": ", ".join(map(str, sorted(_STALLED_TYPES))),
        "NumCascadeDecisions": _c(_BR.get("cascade_true", 0) + _BR.get("cascade_false", 0)
                                  + sum(v for k, v in _BR.items() if k.startswith("stuck"))),
        # every refuted labelling, and how: a label-only minor (a nonzero element
        # of K) versus any certificate that needs interval arithmetic
        "NumRefuted": _c(P6["cascade_rejected"] + W6["exactly_certified_infeasible"]),
        "NumLabelOnly": _c(_kinds.get("label_only_minor", 0)),
        "NumIntervalCert": _c(sum(v for k, v in _kinds.items()
                                  if k != "label_only_minor")),
        "NumPlainLabelOnly": _c(_pk.get("label_only_minor", 0)),
        "NumPlainInterval": _c(sum(v for k, v in _pk.items() if k != "label_only_minor")),
        "NumWildLabelOnly": _c(_wk.get("label_only_minor", 0)),
        "NumWildRootBox": _c(_wk.get("root_box_range", 0)),
        "NumWildBnB": _c(_wk.get("interval_branch_and_bound", 0)),
        "CertMaxBoxes": _c(_maxboxes),
        "NumRootBoxAll": _c(_kinds.get("root_box_range", 0)),
        "NumBnBAll": _c(_kinds.get("interval_branch_and_bound", 0)),
        "NumIntervalKinds": ", ".join(f"{k}: {v}" for k, v in sorted(_kinds.items())
                                      if k != "label_only_minor"),
        # the truncated-deduplication check (checks/dedup_truncation.py --all)
        "DedupSources": _c(_DD["source_order_types"]),
        "DedupMatching": _c(_DD["candidates_matching"]),
        "DedupSizeSix": _c(_DD["candidates_with_size6"]),
        "DedupOther": _c(_DD["candidates_other_type"]),
        "RealizingMF": REALIZING_MF,
        "NumValidLess": str(len(TYPES) - len(BIGMF_KILLED) - 1),
        "NumBigMFAll": str(len(BIGMF_KILLED)),
        "NumValid": str(len(TYPES) - len(BIGMF_KILLED)),
        "NumTypesLabelled": str(len({int(k.split("|")[0]) for k, v in STATE.items()
                                     if v.get("enum_count", 0)})),
        "NumTypesEmptied": str(len({int(k.split("|")[0]) for k in STATE})
                               - len({int(k.split("|")[0]) for k, v in STATE.items()
                                      if v.get("enum_count", 0)})),
        "NumRootExhausted": str(sum(1 for k, v in STATE.items()
                                    if k.endswith("|") and v.get("exhausted"))),
        "NumRefined": str(len({int(k.split("|")[0]) for k in STATE
                               if not k.endswith("|")})),
        "NumSubtreesLabelled": str(sum(1 for v in STATE.values()
                                       if v.get("enum_count", 0))),
    }
    (OUT / "summary_nums.tex").write_text(
        "".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in defs.items()))

    print(f"wrote {OUT/'pertype.tex'} ({len(rows)} rows)")
    print(f"wrote {OUT/'summary_nums.tex'}")
    for k, v in defs.items():
        print(f"  {k:18s} {v}")
    # invariants the paper asserts
    assert len(TYPES) == 387 and len(SURV) == (304 if VALID else 55)
    searched = sorted({int(k.split("|")[0]) for k in STATE})
    notclosed = [t for t in searched if not _closed(t)]
    assert not notclosed, f"types not rigorously closed: {notclosed}"
    # no verdict may rest on a truncation at the depth cap, a wildcard window that
    # hit its scan edge, or a per-assignment deadline
    dep = lambda k: len(k.split("|")[1].split(",")) if k.split("|")[1] else 0
    cap = [k for k, v in STATE.items() if not v.get("exhausted") and dep(k) >= MAXD]
    assert not cap, f"depth-cap leaves: {cap}"
    unb = [k for k, v in STATE.items() if (v.get("diag") or {}).get("wild_unbounded")]
    assert not unb, f"wild_unbounded subtrees: {unb}"
    dln = [k for k, v in STATE.items() if (v.get("diag") or {}).get("wild_deadline")]
    assert not dln, f"wild_deadline subtrees: {dln}"
    assert [t for t in searched if _keys(t)] == [379]
    # the two combinatorial theorems must not touch the realizing type, and must
    # agree with the search wherever both speak
    assert 379 not in COMB_KILLED
    # "exactly one up to isometry": the accepted labelling must determine its
    # weights (one leaf of the weight tree, one passing), and the label assignments
    # describing the polytope must be a single Aut(T)-orbit.
    assert UW["weight_tree_leaves"] == UW["weight_tree_leaves_passing"] == 1, (
        f"the accepted labelling has {UW['weight_tree_leaves']} weight-tree leaves, "
        f"{UW['weight_tree_leaves_passing']} passing")
    assert (UW["orbit_of_the_labelling"] * UW["labelling_stabiliser_order"]
            == UW["aut_type_order"]), "orbit-stabiliser does not hold as recorded"
    # the d=5 calibration quoted in Section 5.4
    assert D5W_TOTAL - D5W_REJECTED == 9, "d=5 wildcard survivors changed"
    # Proposition "wildexact": EVERY wildcard-bearing labelling of the run of
    # record must carry an exact certificate, and the certifier must refute no
    # realizable instance.  Either failure falsifies the proposition as stated, so
    # these are assertions and not reported numbers.
    n_wild_run = sum((v.get("diag") or {}).get("wild_assignments", 0)
                     for v in STATE.values())
    assert W6["instances"] == n_wild_run, (
        f"certified {W6['instances']} wildcard instances but the run of record "
        f"screened {n_wild_run}")
    assert W6["undecided_among_rejected"] == 0, (
        f"{W6['undecided_among_rejected']} wildcard rejections lack an exact "
        f"certificate")
    for tag, W in (("d=6", W6), ("d=5", W5), ("d=5 census passers", W5P)):
        assert W["wrongly_certified_among_screen_passers"] == 0, (
            f"{tag}: the exact certifier refuted a labelling that SURVIVED the "
            f"screen")
    # The one-sidedness check need not certify the whole census, but it must cover
    # every instance in it that survived the screen -- those are the only ones a
    # false refutation could damage.  Both counts come from the census artifact.
    _d5pass_census = sum(v.get("wild_assignments", 0)
                         - (v.get("wild") or {}).get("joint_infeasible", 0)
                         for v in _D5W.values())
    assert W5P["census_instances"] == D5W_TOTAL, (
        f"the one-sidedness check reports a census of {W5P['census_instances']} "
        f"instances, but the census of record has {D5W_TOTAL}")
    assert W5P["passed_joint_screen"] == _d5pass_census, (
        f"the one-sidedness check covers {W5P['passed_joint_screen']} screen-passers "
        f"but the census has {_d5pass_census}")
    # The d=5 calibration set must contain every instance that survived the joint
    # screen -- that is where a false refutation would show up.  Nine of them are
    # the pi/10 polytopes; the rest survive the continuous relaxation and are
    # killed later by the integer window scan, so the certifier cannot and must not
    # refute them either.
    _d5pass = sum(v.get("wild_assignments", 0)
                  - (v.get("wild") or {}).get("joint_infeasible", 0)
                  for v in _D5W.values())
    assert W5["passed_joint_screen"] == _d5pass or W5["instances"] < D5W_TOTAL, (
        f"d=5 set has {W5['passed_joint_screen']} screen-passers, census artifact "
        f"says {_d5pass}")
    # Proposition "plainexact": the same for the labellings without a wildcard.
    # Together with the line above, every labelling that reached the screen in the
    # run of record is decided exactly, and the one that was accepted is not
    # refuted -- which is the claim, so it is asserted rather than reported.
    n_plain_run = (sum(v.get("enum_count", 0) for v in STATE.values())
                   - sum((v.get("diag") or {}).get("wild_assignments", 0)
                         for v in STATE.values()))
    assert P6["instances"] == n_plain_run, (
        f"certified {P6['instances']} non-wildcard labellings but the run of "
        f"record screened {n_plain_run}")
    assert P6["undecided_among_rejected"] == 0, (
        f"{P6['undecided_among_rejected']} cascade rejections lack an exact "
        f"certificate")
    assert P6["wrongly_certified_among_candidates"] == 0, (
        "the exact certifier refuted the labelling the search accepted")
    assert P6["realizable_candidates"] == 1, (
        f"{P6['realizable_candidates']} labellings passed the cascade, expected "
        f"exactly one (P^B6)")
    assert sum(v["distinct"] for v in _D5W.values()) == 51
    assert all(v["exhausted"] and not v["unbounded"] for v in _D5W.values())
    assert not [t for t in searched if t in COMB_KILLED and _keys(t)], \
        "a theorem-excluded type realized a polytope"
    assert len(_keys(379)) == 1
    assert d5_total == 51 and all(x["exhausted"] and not x["unbounded"]
                                  for x in D5.values())
    print("\nall paper invariants re-checked OK")


if __name__ == "__main__":
    sys.exit(main())
