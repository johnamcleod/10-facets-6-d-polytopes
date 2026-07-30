#!/usr/bin/env python3
"""Is the bounded-box floating-point fallback ever reachable in the d=6 run?

Deciding a labelling normally pins the ultraparallel weights by the cascade of
Section 5.3: repeatedly find an 8x8 minor of G whose index set contains exactly
one not-yet-pinned dashed edge (the minor is then a quadratic in that one
weight), or, failing that, a pair of minors over the same two unknowns (Sylvester
resultant).  Only if NEITHER is available does the code fall back to a
multistart L-BFGS-B search over the bounded box x_e in [1.001, 1000] -- the one
place where an emptiness verdict depends on a bounded weight range.

Whether such a minor exists depends only on the dashed-edge structure of the
combinatorial type, not on the labels or the weight values: the minor index set
is fixed, and "contains exactly one unpinned dashed edge" is a statement about
sets.  So reachability of the fallback is decidable per type, statically.

This script runs that decision for all 52 searched d=6 types.

Expected output:  NONE -- the fallback is structurally unreachable in d=6.

Caveat, stated precisely: this rules out the fallback being reached because the
cascade runs out of pinnable edges.  The cascade can also return "inconclusive"
for value-dependent reasons -- a minor that vanishes identically in its unknown,
a resultant that vanishes identically, or more than 256 leaves -- which this
static analysis cannot exclude.  What it does establish is that the *structural*
stall, the reason the fallback exists at all and the one that fires for d=4
type 15 (whose dashed edges form a perfect matching), never occurs here.

Run:  python3 paper/checks/cascade_reachability.py
"""
import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pipeline.stage4_gram import _build_minor_index   # noqa: E402

N, D = 10, 6


def cascade_stalls(dashed):
    """True iff the cascade runs out of pinnable edges for this dashed-edge set."""
    minors = _build_minor_index(dashed, N, D)
    unknown = set(dashed)
    while unknown:
        single = next((u[0] for S, de in minors
                       for u in [[e for e in de if e in unknown]]
                       if len(u) == 1), None)
        if single is not None:
            unknown.discard(single)
            continue
        pairs = defaultdict(list)
        for S, de in minors:
            u = tuple(sorted(e for e in de if e in unknown))
            if len(u) == 2:
                pairs[u].append(S)
        pair = next((p for p, ss in pairs.items() if len(ss) >= 2), None)
        if pair is None:
            return True
        unknown -= set(pair)
    return False


def structural_facts(types):
    """Two stronger statements than the per-survivor result: they hold over EVERY
    generated type, so they do not depend on the facet filter.

      (i) no d=6 type's dashed edges form a perfect matching;
     (ii) no d=6 type fails to have an 8-minor isolating a single dashed edge at
          the start of the cascade -- the condition that forces the pair path.

    Both matter because the one demonstrated false-negative region in dimension 4
    is the 4-cube type, whose dashed edges DO form a perfect matching and which is
    therefore forced onto the pair-resultant branch."""
    pm, forced = [], []
    for tid, t in sorted(types.items()):
        dashed = [tuple(sorted(m)) for m in t["missing_faces"] if len(m) == 2]
        if not dashed:
            continue
        covered = {x for e in dashed for x in e}
        if len(covered) == 2 * len(dashed) == N:
            pm.append(tid)
        ds = set(dashed)
        if not any(len([e for e in de if e in ds]) == 1
                   for _, de in _build_minor_index(dashed, N, D)):
            forced.append(tid)
    print(f"over all {len(types)} generated d=6 types:")
    print(f"   dashed edges form a perfect matching:            "
          f"{pm if pm else 'NONE'}")
    print(f"   no 8-minor isolates a single dashed edge:        "
          f"{forced if forced else 'NONE'}")
    print("   (the d=4 4-cube type G4 has both properties; no d=6 type has either)")
    return not pm and not forced


def main():
    types = {t["type_id"]: t
             for t in json.load(open(ROOT / "runs/d6_n10/stage2/types.json"))}
    surv = json.load(open(ROOT / "runs/d6_n10/facet_profile_survivors.json"))["survivors"]
    three_free = [t for t in surv
                  if {len(m) for m in types[t]["missing_faces"]} == {2}]
    searched = [t for t in surv if t not in three_free]

    stalling = []
    for tid in searched:
        dashed = [tuple(sorted(m)) for m in types[tid]["missing_faces"]
                  if len(m) == 2]
        if cascade_stalls(dashed):
            stalling.append(tid)

    strong = structural_facts(types)
    print()
    print(f"searched d=6 types: {len(searched)} "
          f"({len(three_free)} 3-free types excluded by theorem)")
    print("types where the cascade structurally stalls "
          "(bounded-box fallback reachable):",
          stalling if stalling else "NONE")
    print()
    print("RESULT:", "the [1.001, 1000] fallback is structurally unreachable "
          "for every searched d=6 type"
          if not stalling else f"reachable for {len(stalling)} types")
    return 0 if (not stalling and strong) else 1


if __name__ == "__main__":
    sys.exit(main())
