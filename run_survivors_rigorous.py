#!/usr/bin/env python3
"""Rigorous full-label resolve of the 54 d=6 facet-filter survivors.

Replaces run_survivors_fulllabel.py, whose "empty(COMPLETE)" tag was a wall-clock
heuristic: 19 types hit the 3600s enum timeout and were mislabeled complete —
including tid 379 = P^B6 itself (distinct=0 there is a BUG per CLAUDE.md §0).

This driver uses the real exhaustion flag (process_type_stage4 stats_out) and
partitions the label search by enumeration prefix (enumerate_labels_backtrack
`prefix`), refining any timed-out subtree one level deeper.  A type's verdict is
rigorous iff EVERY subtree reports exhausted=True; the type's polytopes are the
union of subtree results deduped by canonical_key.

Resumable: completed subtrees are persisted to STATE and skipped on restart.

Usage:  python3 run_survivors_rigorous.py [fast|suspect|all|<tid,tid,...>] [nproc] [wildcard] [d=D] [out=DIR]

With the "wildcard" flag the solve is the RIGOROUS label treatment (Ma-Zheng
Prop 3.5): labels enumerate over {2,...,6,7} with 7 = "any m >= 7" resolved by
continuous-c range analysis + integer instantiation — no a-priori label cap.
Without it, labels are the {2,...,10,12} alphabet (a-posteriori d=4/d=5 set;
verdicts are then relative to that assumption).  Separate state dirs.

`d=D` selects the dimension (default 6, so every existing invocation is
unaffected).  For D != 6 the driver reads runs/dD_n{D+4}/stage2/types.json, keeps
its own state directory, and takes the type list from the facet-profile survivor
file if one exists, else from all types in the file.  The d=6 fast/suspect
partition and the Esselmann 3-free elimination are dimension-aware; see below.
"""
import json, sys, time
from pathlib import Path
import multiprocessing as mp
sys.path.insert(0, str(Path(__file__).resolve().parent))
import pipeline.stage4_gram as _sg
from pipeline.stage4_gram import process_type_stage4, VALID_LABELS, WILDCARD_INDICES
from pipeline.stage4_blockpaste import _setup
from pipeline.utils.automorphisms import compute_aut_group
from run_d4 import canonical_key

WILDCARD = "wildcard" in sys.argv[1:]
# Dimension: default 6 so all existing command lines behave exactly as before.
D = next((int(a.split("=")[1]) for a in sys.argv[1:] if a.startswith("d=")), 6)
NF = D + 4                       # number of facets

SURV = [8, 12, 17, 34, 36, 38, 40, 51, 55, 59, 60, 61, 69, 70, 92, 103, 120, 127, 132,
        140, 154, 159, 162, 168, 173, 206, 214, 218, 220, 229, 234, 239, 255, 265, 273,
        284, 286, 287, 295, 297, 308, 315, 317, 320, 329, 332, 344, 352, 354, 356, 360,
        378, 379, 382]
# The 19 that hit the 3600s enum timeout in the fulllabel run (0 there = not rigorous).
SUSPECT = [34, 40, 55, 60, 61, 103, 120, 127, 132, 159, 173, 206, 284, 286, 295, 315,
           329, 344, 379]
FAST = [t for t in SURV if t not in SUSPECT]

# Killed by THEOREM, not enumeration.  Types 159 and 329 have missing-face
# profile (2,2,2,2,2,2): every missing face has size 2, i.e. they are 3-FREE.
#
# CORRECTED JUSTIFICATION (2026-07-27).  This used to cite Burcroff
# arXiv:2201.03437 Theorem 8.1 ("no compact Coxeter 6-polytope with 10 facets
# having missing faces of orders only 2 and 5").  That citation is NOT adequate:
# her proof of Thm 8.1 analyses only the two types I1 and I2 she lists, and both
# contain a missing face of size 5, so the argument does not reach the 3-free
# case.  The correct authority is ESSELMANN (Über kompakte hyperbolische
# Coxeter-Polytope mit wenigen Facetten, Bielefeld 1994, Lemma 6.7; quoted as
# Burcroff Lemma 10.6): a 3-free compact Coxeter d-polytope has at least 2d
# facets, with equality only for the d-cube.  For d = 6 that is 12 > 10, so no
# 3-free type can realize.  Unconditional, and it makes the whole d=6 chain
# independent of Burcroff Thm 8.1 (among the 54 survivors the only types lacking
# a size-3/4 missing face are exactly these two 3-free ones).
#
# These are also the two types with no Lannér subdiagram constraints, hence the
# deepest enumeration grinders; their partial enumeration (139 subtrees each,
# ~8.5M and ~11.9M labellings, zero candidates) is retained as independent
# consistency evidence only -- neither enumeration is complete.
ESSELMANN_3FREE_KILLED = [159, 329]

_RUNDIR = f"runs/d{D}_n{NF}"
TYPES = {t['type_id']: t for t in json.load(open(f'{_RUNDIR}/stage2/types.json'))}

if D != 6:
    # Take the candidate list from the facet-profile survivors if that filter has
    # been run for this dimension, else every generated type.
    _fp = Path(f"{_RUNDIR}/facet_profile_survivors.json")
    SURV = (json.loads(_fp.read_text())["survivors"] if _fp.exists()
            else sorted(TYPES))
    SUSPECT, FAST = [], SURV      # no pre-classification for other dimensions

# Esselmann (1994, Lemma 6.7): a 3-free compact Coxeter d-polytope has at least
# 2d facets, with equality only for the d-cube.  So the elimination applies only
# when n < 2d.  For d=6, n=10 < 12 -> 3-free types are killed.  For d=4, n=8 = 2d
# -> a 3-free type could be the 4-cube, which does admit compact Coxeter
# structures, so those types must be SEARCHED, not killed.
if NF < 2 * D:
    ESSELMANN_3FREE_KILLED = [t for t in SURV
                              if {len(m) for m in TYPES[t]['missing_faces']} == {2}]
else:
    ESSELMANN_3FREE_KILLED = []

# Two further combinatorial exclusions, each a theorem, both consequences of
# Lannér's classification (paper Lemmas 4.2 and 4.3):
#
#   LEMMA 4.2 (bounded dashed degree).  A facet disjoint from t others is a compact
#   Coxeter (d-1)-polytope with n-1-t facets.  For d-1 >= 5 such a polytope is not
#   a simplex, so it has at least d+1 facets and t <= n-d-2.  For d=6, n=10 this
#   caps the dashed-edge degree at 2.
#
#   LEMMA 4.3 (no missing face of size d).  It would need a Lannér subdiagram of
#   order d, and Lannér diagrams exist only in orders 2..5.  Detected by comparing
#   the vertex list with the d-subsets containing no RECORDED missing face, since
#   the generator records minimal non-faces only up to size 5.
#
# The flags are precomputed by paper/checks/type_consistency.py --emit into
# runs/d{D}_n{NF}/type_flags.json; if that artifact is absent the exclusions are
# simply not applied, so the driver still runs (conservatively) without it.
LEMMA_KILLED = []
_flags_path = Path(f"{_RUNDIR}/type_flags.json")
if _flags_path.exists():
    _fl = json.loads(_flags_path.read_text()).get(f"d{D}_n{NF}", {})
    LEMMA_KILLED = sorted(
        t for t in SURV
        if (_fl.get(str(t), {}).get("max_dashed_degree", 0) >= 3
            or _fl.get(str(t), {}).get("has_missing_face_of_size_d", False)))

# Output directory.  `out=NAME` selects a different one, which is how a FRESH run
# is started: without it the driver resumes the existing state and every cached
# subtree is skipped, so a "re-run" would silently re-emit the old verdicts
# without recomputing anything.  A loud warning is printed if a populated state is
# resumed (see main()).
_OUT = next((a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("out=")),
            "survivors_wildcard" if WILDCARD else "survivors_rigorous")
OUTDIR = Path(f"{_RUNDIR}/{_OUT}")
STATE = OUTDIR / "state.json"
STATE_JSONL = OUTDIR / "state.jsonl"
REALIZER_DIR = OUTDIR / "realizers"

N_LAB = len(WILDCARD_INDICES) if WILDCARD else len(VALID_LABELS)   # 6 / 10
# 6 -> 8 (2026-07-17): the deep (2^5,3,5)-profile cores reach depth-6 timeouts.
# 8 -> 10 (2026-07-25): tid34 pfx=(0,1,0,0,0,1,0,3) hit the depth-8 cap itself
# (TIMEOUT->refine d9 tag, 1800s, still not exhausted) -> would surface as
# UNRESOLVED under the old cap. Extending preemptively rather than waiting for
# the rest of wave 8 (5202 tasks) to finish and discover more such cases.
# 10 -> 14 (2026-07-26): with orbit-symmetry pruning deployed, wave 10 (the
# depth-10 leaves) already shows 3 tid34 prefixes hitting the full 1800s
# timeout AT the cap itself (TIMEOUT->refine d11) -- these are genuinely
# canonical branches (not symmetric duplicates, or they'd have resolved in
# ~0.1s), so they need real further depth, not just a bigger timeout.
# Extending well past the observed failure point rather than incrementally,
# to avoid a third restart. Refinement continues losslessly from cached
# state on relaunch.
MAX_DEPTH = 14
MAXA = 50_000_000
SOLVE_TO = 600.0

def enum_to(depth):
    return 3600.0 if depth == 0 else 1800.0


def _key(tid, prefix):
    return f"{tid}|{','.join(map(str, prefix))}"


# Orbit-based symmetry-breaking (pipeline/utils/automorphisms.py): validated
# exactly against all 4 d=5 anchors (identical distinct/exhausted, real
# speedup) and against a direct 180s-window benchmark on the two most
# resistant d=6 types (tid34/tid103: full 180s timeout -> 0.1s exhausted).
# Sound because it does a full-vector lex comparison against every
# domain-stabilizing automorphism image of the CURRENT fixed prefix, never an
# independent per-pair inequality -- a non-canonical prefix's entire subtree
# is provably redundant with the (separately-scheduled) canonical prefix's
# subtree, so reporting it empty+exhausted immediately loses no solutions.
# Computed once per type before the pool forks (cheap: VF2 on <=42 ordinary
# pairs), so per-task overhead is just a lookup.
_AUT_CACHE: dict = {}


def _aut_for(tid):
    if tid not in _AUT_CACHE:
        t = dict(TYPES[tid]); V, *_ = _setup(t)
        n = 1 + max(max(v) for v in V)
        _AUT_CACHE[tid] = compute_aut_group(V, n)
    return _AUT_CACHE[tid]


def work(task):
    tid, prefix = task
    try:
        # Reset the per-subtree diagnostic counters.  These are PERSISTED with the
        # verdict (2026-07-28): the previous run recorded only keys/exhausted/
        # enum_count/sec, which made it impossible to answer afterwards which gate
        # rejected what -- the reason the refinement-bug exposure had to be
        # estimated by re-probing instead of read off the artifacts.
        for _D in (_sg.CERT_STATS, _sg.WILD_STATS, _sg.REFINE_STATS):
            for _k in _D:
                _D[_k] = 0
        t = dict(TYPES[tid]); V, *_ = _setup(t); t['vertex_sets'] = [sorted(v) for v in V]
        so = {}
        t0 = time.time()
        res = process_type_stage4(t, D, max_assignments=MAXA,
                                  enum_timeout=enum_to(len(prefix)),
                                  solve_timeout=SOLVE_TO,
                                  label_indices=None, prefix=tuple(prefix) or None,
                                  stats_out=so, wildcard=WILDCARD,
                                  use_burcroff_55b=WILDCARD,
                                  automorphisms=_aut_for(tid), verbose=False) or []
        el = round(time.time() - t0, 1)
        keys = sorted(set(str(canonical_key(r, NF)) for r in res))
        # a wild window hitting the scan edge means the subtree verdict cannot be
        # closed by instantiation — treat as not exhausted (never silently rigorous)
        exhausted = bool(so.get("exhausted")) and not so.get("wild_unbounded")
        diag = {
            "screened": so.get("screened", 0),
            "passed_screen": so.get("passed_screen", 0),
            "exact_attempts": so.get("exact_attempts", 0),
            "wild_assignments": so.get("wild_assignments", 0),
            "wild_unbounded": bool(so.get("wild_unbounded")),
            "wild_deadline": bool(so.get("wild_deadline")),
            "screen_branches": so.get("screen_branches") or {},
            "cert": {k: v for k, v in _sg.CERT_STATS.items() if v},
            "wild": {k: v for k, v in _sg.WILD_STATS.items() if v},
            "refine": {k: v for k, v in _sg.REFINE_STATS.items() if v},
        }
        return (tid, prefix, keys, res, exhausted,
                so.get("enum_count", 0), el, None, diag)
    except Exception as e:
        return (tid, prefix, [], [], False, 0, 0.0, f"{type(e).__name__}: {e}", {})


def main():
    sel = sys.argv[1] if len(sys.argv) > 1 else "all"
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 else max(1, mp.cpu_count() - 1)
    if sel == "fast":
        tids = FAST
    elif sel == "suspect":
        tids = SUSPECT
    elif sel == "all":
        tids = SUSPECT + FAST          # suspects first: they gate the result
    else:
        tids = [int(x) for x in sel.split(",")]

    OUTDIR.mkdir(parents=True, exist_ok=True)
    REALIZER_DIR.mkdir(parents=True, exist_ok=True)
    # Checkpoint format: append-only JSONL (one completed subtree per line).  The
    # previous code rewrote the entire state.json after EVERY subtree, which is
    # O(n^2) in I/O -- ~40 GB of writes for the 20,466-subtree d=6 run, and several
    # hundred GB once the per-subtree diagnostics are included.  state.json is
    # still written once at the end as a consolidated snapshot, and is still read
    # on startup if no JSONL exists, so old checkpoints resume unchanged.
    state = {}
    if STATE_JSONL.exists():
        for line in STATE_JSONL.read_text().splitlines():
            if not line.strip():
                continue
            try:
                rec = json.loads(line)
            except Exception:
                continue
            state[rec.pop("key")] = rec
    elif STATE.exists():
        state = json.loads(STATE.read_text())
    if state:
        print(f"!! RESUMING {len(state)} cached subtrees from {OUTDIR}. "
              f"Cached subtrees are NOT recomputed. For a fresh run pass "
              f"out=<new-dir>.", flush=True)

    # d=6: put the known realizer (type 379) first, so the anchor is checked
    # before the long tail.  Other dimensions have no such anchor.
    anchor = 379 if D == 6 else None
    tids = sorted(set(tids) - set(ESSELMANN_3FREE_KILLED) - set(LEMMA_KILLED),
                  key=lambda t: (t != anchor, t))
    queue = [(tid, ()) for tid in tids]

    # Precompute Aut(type) for every type up front, before the pool forks, so
    # workers inherit the cache via copy-on-write instead of recomputing VF2
    # per task.
    for tid in tids:
        _aut_for(tid)

    wave = 0
    while queue:
        todo = [(tid, pfx) for tid, pfx in queue
                if _key(tid, pfx) not in state]
        skipped = len(queue) - len(todo)
        print(f"--- wave {wave}: {len(todo)} subtree tasks ({skipped} cached), "
              f"nproc={nproc} ---", flush=True)
        next_queue = []
        if todo:
            with mp.get_context("fork").Pool(processes=nproc) as pool:
                for tid, pfx, keys, res, exh, ecount, el, err, diag in \
                        pool.imap_unordered(work, todo):
                    k = _key(tid, pfx)
                    if err:
                        print(f"tid {tid} pfx={pfx}: ERROR {err}", flush=True)
                        state[k] = {"error": err, "exhausted": False, "sec": el,
                                    "diag": diag}
                    else:
                        tag = ("EXHAUSTED" if exh else
                               f"TIMEOUT->refine d{len(pfx)+1}")
                        if keys:
                            tag = "**REALIZES** " + tag
                            fn = REALIZER_DIR / f"tid{tid}_{'_'.join(map(str,pfx)) or 'root'}.json"
                            fn.write_text(json.dumps(res, indent=1, default=str))
                        print(f"tid {tid} pfx={pfx}: distinct={len(keys)} "
                              f"enum={ecount} [{tag}] in {el}s", flush=True)
                        state[k] = {"keys": keys, "exhausted": exh,
                                    "enum_count": ecount, "sec": el,
                                    "diag": diag}
                    with STATE_JSONL.open("a") as _fh:
                        _fh.write(json.dumps({"key": k, **state[k]},
                                             default=str) + "\n")
        # refinement pass over the whole queue (incl. cached entries)
        for tid, pfx in queue:
            st = state.get(_key(tid, pfx), {})
            if st.get("exhausted") or st.get("error"):
                continue
            if len(pfx) >= MAX_DEPTH:
                print(f"tid {tid} pfx={pfx}: UNRESOLVED at max depth", flush=True)
                continue
            next_queue.extend((tid, tuple(pfx) + (li,)) for li in range(N_LAB))
        queue = next_queue
        wave += 1

    # final per-type verdicts.  A subtree is covered if it was exhausted itself
    # or (having timed out) all 10 of its children are covered — the refinement
    # tree partitions the search exactly, so root coverage = rigorous verdict.
    def covered(tid, pfx):
        st = state.get(_key(tid, pfx))
        if st is None or st.get("error"):
            return False
        if st.get("exhausted"):
            return True
        if len(pfx) >= MAX_DEPTH:
            return False
        return all(covered(tid, tuple(pfx) + (li,)) for li in range(N_LAB))

    verdict = {}
    for tid in LEMMA_KILLED:
        _f = _fl.get(str(tid), {})
        why = ("Lemma 4.2: a facet is disjoint from %d others, so it would be a "
               "compact Coxeter 5-polytope with %d facets, i.e. a simplex or "
               "smaller, and compact hyperbolic Coxeter simplices exist only for "
               "d <= 4 (Lanner)" % (_f.get("max_dashed_degree", 0),
                                    NF - 1 - _f.get("max_dashed_degree", 0))
               if _f.get("max_dashed_degree", 0) >= 3 else
               "Lemma 4.3: the type has a minimal non-face of size %d, which would "
               "require a Lanner subdiagram of that order; none exists above 5" % D)
        verdict[tid] = {"distinct": 0, "rigorous": True, "subtrees": 0, "by": why}
        print(f"VERDICT tid {tid}: distinct=0 RIGOROUS ({why[:46]}...)", flush=True)
    for tid in ESSELMANN_3FREE_KILLED:
        verdict[tid] = {"distinct": 0, "rigorous": True, "subtrees": 0,
                        "by": "Esselmann-3free-bound (Esselmann 1994, Lemma 6.7; "
                              "= Burcroff arXiv:2201.03437 Lemma 10.6): all "
                              "missing faces have size 2, so the polytope would "
                              "be 3-free and would need >= 2d = 12 facets > 10"}
        print(f"VERDICT tid {tid}: distinct=0 RIGOROUS (Esselmann 3-free "
              f"2d-facet bound)", flush=True)
    for tid in tids:
        ks = [s for key, s in state.items() if key.startswith(f"{tid}|")]
        keys = sorted(set(k for s in ks for k in s.get("keys", [])))
        rigorous = bool(ks) and covered(tid, ())
        verdict[tid] = {"distinct": len(keys), "rigorous": rigorous,
                        "subtrees": len(ks)}
        print(f"VERDICT tid {tid}: distinct={len(keys)} "
              f"{'RIGOROUS' if rigorous else 'NOT-RIGOROUS'}", flush=True)
    STATE.write_text(json.dumps(state, indent=1, default=str))
    (OUTDIR / "verdicts.json").write_text(json.dumps(verdict, indent=1))
    realizing = sorted(t for t, v in verdict.items() if v["distinct"] > 0)
    nonrig = sorted(t for t, v in verdict.items() if not v["rigorous"])
    print(f"\nRIGOROUS-DONE. realizing={realizing} not-rigorous={nonrig}", flush=True)


if __name__ == "__main__":
    main()
