#!/usr/bin/env python3
"""Exact refutation of the labellings WITHOUT a wildcard.

The companion of checks/wild_exact_certify.py, and the second half of the
screen caveat.  A labelling with no wildcard is decided by the pinning cascade of
Section 5.3, which is algebraically structured but evaluated in double precision:
its minors, quadratic discriminants and root tests all sit against tolerances of
1e-7 to 1e-12.  A rejection from it is therefore a computational claim, not a
proof, and this script replaces those claims by exact certificates.

WHY THIS CASE NEEDS MORE THAN THE WILDCARD CASE.  For a wildcard labelling it
suffices to show that no point of the domain has rank Gr <= d+1: the wild entries
are free, so if the rank condition is unsatisfiable the labelling is dead.  For a
labelling with no wildcard the dashed weights are DETERMINED by the rank
condition, so rank-deficient points of the domain generally do exist and no rank
minor can refute the labelling.  What fails at those points is the signature.
Accordingly the certifier prunes a box when EITHER

    some (d+2)-principal minor is nonzero throughout it   (rank >= d+2), OR
    some principal submatrix has two negative eigenvalues throughout it
    (superhyperbolic, so the signature cannot be (d,1)),

the second read off interval-evaluated leading principal minors by Jacobi's rule.
Both are exact: coefficients over Q(sqrt2,sqrt3,sqrt5), intervals rounded
outward, the ultraparallel directions compactified by x_e = 1/t_e so no bound on
the weights is assumed.  See pipeline/utils/exact_certify.py.

THE ANCHOR.  Exactly one of these labellings IS realizable -- P^B6 itself.  It
must come through unrefuted, and the script fails loudly if it does not.  That is
the one test this whole file could fail in the dangerous direction.

Run:  python3 checks/plain_exact_certify.py [--workers=N]
"""
import collections
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.utils.exact_certify import certify_instance_full   # noqa: E402

# src=/dst= select another dump (e.g. the lemmas=valid run); paths are relative
# to the repository root.  Defaults are the artifacts of the 11-type run.
_ARG = {a.split("=", 1)[0]: a.split("=", 1)[1] for a in sys.argv[1:]
        if "=" in a and not a.startswith("--")}
SRC = ROOT / _ARG.get("src", "runs/d6_n10/plain_instances.jsonl")
DST = ROOT / _ARG.get("dst", "runs/d6_n10/plain_certificates.json")

MAX_MINORS = 8
MAX_INERTIA = 40
MAX_BOXES = 200_000


def _complete(src):
    meta = src.with_suffix(src.suffix + ".meta.json")
    if not (src.exists() and meta.exists()):
        return False
    m = json.loads(meta.read_text())
    n = sum(1 for line in src.read_text().splitlines() if line.strip())
    return bool(m.get("complete")) and m.get("records") == n


def _one(rec):
    t0 = time.time()
    ok, stats = certify_instance_full(rec, max_minors=MAX_MINORS,
                                      max_inertia=MAX_INERTIA,
                                      max_boxes=MAX_BOXES)
    return {
        "tag": rec.get("tag"),
        "n": rec["n"],
        "d": rec["d"],
        "dashed": rec["dashed"],
        "wild": [],
        "ordinary": rec["ordinary"],
        "cascade_decision": rec.get("cascade_decision"),
        "certified_unrealizable": bool(ok),
        "route": stats.get("route"),
        "boxes": stats.get("boxes"),
        "witness": stats.get("witness"),
        "undecided_reason": None if ok else stats.get("reason"),
        "seconds": round(time.time() - t0, 2),
    }


def main():
    workers = next((int(a.split("=")[1]) for a in sys.argv[1:]
                    if a.startswith("--workers=")), max(1, mp.cpu_count() - 2))
    if not SRC.exists():
        print(f"MISSING {SRC.relative_to(ROOT)} -- run  python3 run_wild_dump.py",
              file=sys.stderr)
        return 2
    if not _complete(SRC):
        print(f"{SRC.relative_to(ROOT)} is INCOMPLETE (no matching completion "
              f"marker) -- let run_wild_dump.py finish", file=sys.stderr)
        return 2

    recs = [json.loads(l) for l in SRC.read_text().splitlines() if l.strip()]
    t0 = time.time()
    if workers > 1 and len(recs) > 20:
        with mp.Pool(workers) as pool:
            out = list(pool.imap(_one, recs, chunksize=4))
    else:
        out = [_one(r) for r in recs]
    elapsed = time.time() - t0

    # The cascade's own verdict partitions the instances.  "candidate" is a
    # labelling the cascade passed on to exact certification.  Everything else must
    # be refuted exactly: "reject" is a cascade rejection, and "stalled" is a
    # labelling on which the cascade stalled for a value-dependent reason (not the
    # structural stall that Proposition 5.2 excludes) and which the bounded-box
    # fallback then rejected -- a numerical rejection, so it needs a certificate
    # just as much.  (Until 2026-09-28 every non-"reject" decision was counted as
    # a candidate; the 11-type run had no stalls, so the two readings agreed.)
    by_dec = collections.Counter(r["cascade_decision"] for r in out)
    routes = collections.Counter(r["route"] for r in out if r["certified_unrealizable"])
    rejected = [r for r in out if r["cascade_decision"] != "candidate"]
    candidates = [r for r in out if r["cascade_decision"] == "candidate"]
    cert_rejected = sum(1 for r in rejected if r["certified_unrealizable"])
    wrongly = [r for r in candidates if r["certified_unrealizable"]]

    summary = {
        "case": "d=6, n=10, labellings without a wildcard",
        "source": str(SRC.relative_to(ROOT)),
        "instances": len(out),
        "cascade_decisions": dict(by_dec),
        "cascade_rejected": len(rejected),
        "cascade_stalled": by_dec.get("stalled", 0),
        "exactly_certified_unrealizable": cert_rejected,
        "undecided_among_rejected": len(rejected) - cert_rejected,
        "realizable_candidates": len(candidates),
        "wrongly_certified_among_candidates": len(wrongly),
        "routes": dict(routes),
        "seconds": round(elapsed, 1),
        "max_minors": MAX_MINORS, "max_inertia": MAX_INERTIA,
        "max_boxes": MAX_BOXES,
    }
    DST.write_text(json.dumps({"summary": summary, "instances": out}, indent=1)
                   + "\n")

    print(f"labellings without a wildcard: {len(out)}, {elapsed:.1f}s "
          f"({workers} workers)")
    print(f"    cascade decisions              : {dict(by_dec)}")
    print(f"    to refute (rejected or stalled): {len(rejected)}")
    print(f"    of those exactly certified     : {cert_rejected}")
    print(f"    undecided                      : {len(rejected) - cert_rejected}")
    print(f"    candidates (should be realizable): {len(candidates)}")
    print(f"    of those wrongly refuted       : {len(wrongly)}")
    print(f"    certificate routes             : {dict(routes)}")
    print(f"    artifact -> {DST.relative_to(ROOT)}")

    problems = []
    if wrongly:
        problems.append(f"{len(wrongly)} labellings that the cascade passed as "
                        f"candidates were REFUTED exactly -- check them against "
                        f"the search's own exact acceptance verdict")
    if len(rejected) - cert_rejected:
        problems.append(f"{len(rejected) - cert_rejected} cascade rejections lack "
                        f"an exact certificate")
    print()
    print("RESULT:", "every non-wildcard rejection is exactly certified, and the "
          "realizable labelling survives" if not problems
          else "PROBLEMS: " + "; ".join(problems))
    return 0 if not problems else 1


if __name__ == "__main__":
    sys.exit(main())
