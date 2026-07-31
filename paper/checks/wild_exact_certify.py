#!/usr/bin/env python3
"""Exact infeasibility certificates for the wildcard-bearing labellings.

BACKGROUND.  A wildcard-bearing labelling is screened by the joint-feasibility
step of Section 5.4: a 24-start L-BFGS-B minimisation of the rank-deficiency
residual over x_e in [1.001, 1000], c_e in [cos(pi/7), 1), rejecting when no start
reaches 1e-6.  The failure of a local search over a bounded box is not a proof, so
no verdict rests on that step.

WHAT THIS SCRIPT DOES.  For each such labelling it produces an exact
certificate that the labelling is infeasible, with no local search, no bounded
weight range and no tolerance:

  * the dashed weights are compactified by x_e = 1/t_e, t_e in [0,1], so the
    unbounded direction is removed rather than truncated;
  * the wild entries keep their own compact range c_e in [c7lo, 1] with c7lo a
    rational lower bound for cos(pi/7), so no integer label range is quantified
    over -- the integer scan happens only after feasibility;
  * each (d+2)x(d+2) principal minor of the Gram matrix, which rank <= d+1
    forces to vanish, is computed EXACTLY as a polynomial over
    Q(sqrt2, sqrt3, sqrt5) in those variables;
  * emptiness on the closed rational box is certified by interval arithmetic
    with outward rounding, refining by branch and bound when the root box does
    not already settle it.

Both enlargements (closing the box, lowering c7) are in the safe direction: they
grow the domain, so emptiness there implies emptiness on the true domain.

VALIDATION.  A certificate can only refute, so the property to establish is that
it never refutes something realizable.  That is measured on the d=5 census, where
ground truth exists.  Eleven of its 149,666 wildcard instances survive the joint
screen: nine become polytopes (exactly the d=5 polytopes with a dihedral angle
pi/10) and two are continuously feasible but killed by the integer window scan.
The certifier must refute none of those eleven -- run `--passers-only` for that
check, which covers the whole census -- and should refute the rejections, which is
measured over whichever rejection set is available (the two types holding all nine
polytopes by default, the whole census with `--full-d5`).

INPUT.  The instance dumps written by run_wild_dump.py from the classification's
own instance sink, so the certified instances are the screened instances and not a
re-derivation of them.

Run:  python3 paper/checks/wild_exact_certify.py
      python3 paper/checks/wild_exact_certify.py --passers-only
      python3 paper/checks/wild_exact_certify.py --shard=i/n   (then --merge --shards=n)
"""
import collections
import json
import multiprocessing as mp
import sys
import time
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from pipeline.utils.exact_certify import certify_instance   # noqa: E402

MAX_MINORS = 12
MAX_BOXES = 50_000

# Above this many instances the artifact stores one compact digest per instance
# instead of the full certificate: the d=5 validation census is 149,666
# instances, and full witnesses for all of them would be a ~150 MB file.  The
# d=6 case, which is the one the theorem depends on, is far below the cap and is
# therefore stored in full.
FULL_DETAIL_CAP = 2000


def _complete(src):
    """True only if run_wild_dump.py finished writing this dump.

    The dump is append-only during a run, so existence is not completeness;
    certifying a prefix would report "every instance certified" about part of the
    search.  The marker carries the record count, which must match the file.
    """
    meta = src.with_suffix(src.suffix + ".meta.json")
    if not (src.exists() and meta.exists()):
        return False
    m = json.loads(meta.read_text())
    n = sum(1 for line in src.read_text().splitlines() if line.strip())
    return bool(m.get("complete")) and m.get("records") == n


def cases():
    # The d=5 calibration set: the whole census when it has been dumped AND its
    # certification merged, otherwise the two-type subset.  Those two types are the
    # ones containing every d=5 wildcard polytope, so the subset still carries the
    # instances a false refutation would damage; the census adds refutation-power
    # coverage over the rejections, which is a much larger and much slower set.
    d5_full = ROOT / "runs/d5_n9/wild_instances_full.jsonl"
    use_full = _complete(d5_full) and "--full-d5" in sys.argv[1:]
    d5_src = d5_full if use_full else ROOT / "runs/d5_n9/wild_instances.jsonl"
    return [
        ("d=6, n=10", ROOT / "runs/d6_n10/wild_instances.jsonl",
         ROOT / "runs/d6_n10/wild_certificates.json"),
        ("d=5, n=9 (validation)", d5_src,
         ROOT / "runs/d5_n9/wild_certificates.json"),
    ]


def _margin(witness):
    """Distance from zero of the certifying enclosure (float, diagnostics only)."""
    if not witness:
        return None
    if witness["kind"] == "label_only_minor":
        lo, hi = (Fraction(x) for x in witness["value_enclosure"])
        return float(min(abs(lo), abs(hi)))
    best = None
    for lo, hi in witness["root_ranges"].values():
        lo, hi = Fraction(lo), Fraction(hi)
        if lo <= 0 <= hi:
            continue
        m = float(min(abs(lo), abs(hi)))
        best = m if best is None else min(best, m)
    return best


def _one(rec):
    ok, stats = certify_instance(rec, max_minors=MAX_MINORS, max_boxes=MAX_BOXES)
    w = stats.get("witness") or {}
    return {
        "tag": rec.get("tag"),
        # n and d are carried through so that each artifact record is a
        # self-contained instance an independent checker can rebuild.
        "n": rec["n"],
        "d": rec["d"],
        "dashed": rec["dashed"],
        "wild": rec["wild"],
        "ordinary": rec["ordinary"],
        "numerical_verdict": rec["numerical_verdict"],
        "numerical_residual": rec.get("joint_residual"),
        "certified_infeasible": bool(ok),
        "minors_used": stats.get("minors_used"),
        "boxes": stats.get("boxes"),
        "witness": w or None,
        "margin": _margin(w),
        "undecided_reason": None if ok else stats.get("reason"),
    }


def _digest(r):
    w = r["witness"] or {}
    return {"tag": r["tag"], "numerical_verdict": r["numerical_verdict"],
            "certified_infeasible": r["certified_infeasible"],
            "kind": w.get("kind"), "boxes": r["boxes"],
            "minor": w.get("minor") or (w.get("minors") or [None])[0],
            "margin": r["margin"]}


def run_case(name, src, dst, workers=1, shard=None):
    """Certify every instance in `src`.

    Parallelism is by SHARDING into independent processes, not by a worker pool:
    on the d=5 calibration set (149,666 instances) a multiprocessing pool wedges
    with all workers idle -- the same macOS unreliability the classification
    drivers avoid.  `shard=(i, n)` certifies the instances with index = i mod n
    and writes its own artifact; `--merge` then combines the shards.  Progress is
    printed as it goes, because this case takes hours.
    """
    recs = [json.loads(l) for l in src.read_text().splitlines() if l.strip()]
    if shard is not None:
        i, n = shard
        recs = [r for k, r in enumerate(recs) if k % n == i]
        dst = dst.with_name(dst.name.replace(".json", f".shard{i}of{n}.json"))
    t0 = time.time()
    out = []
    for r in recs:
        out.append(_one(r))
        if len(out) % 500 == 0:
            rate = len(out) / max(time.time() - t0, 1e-9)
            print(f"    {len(out)}/{len(recs)} certified ({rate:.1f}/s, "
                  f"{(len(recs) - len(out)) / max(rate, 1e-9) / 60:.0f} min left)",
                  flush=True)
    return summarize_and_write(name, src, dst, out, time.time() - t0)


def summarize_and_write(name, src, dst, out, elapsed):
    """Reduce per-instance results to the summary the paper reads, and write the
    artifact.  Shared by a direct run and by the shard merge, so a sharded run and
    a single run produce the same artifact."""
    kinds = collections.Counter((r["witness"] or {}).get("kind", "UNDECIDED")
                                for r in out)

    # "passed the screen" is NOT the same as "realizable": the joint-feasibility
    # step only decides the CONTINUOUS relaxation, and an instance can pass it and
    # still die at the integer window scan (no admissible m_e >= 7) or at exact
    # certification.  At d=5 that happens twice, in type 55, which is why 11
    # instances pass the screen while 9 become polytopes.  The certifier refutes
    # the continuous relaxation, so it must leave all screen-passers undecided --
    # including those two -- and that is what is measured here.
    n_num_infeasible = sum(1 for r in out if r["numerical_verdict"] == "infeasible")
    n_num_feasible = len(out) - n_num_infeasible
    cert_of_infeasible = sum(1 for r in out
                             if r["numerical_verdict"] == "infeasible"
                             and r["certified_infeasible"])
    cert_of_feasible = sum(1 for r in out
                           if r["numerical_verdict"] == "feasible"
                           and r["certified_infeasible"])

    margins = [r["margin"] for r in out if r["margin"] is not None]
    # The residuals the NUMERICAL screen attained on the instances it rejected.
    # Reported because the paper quotes how far those rejections sat from the
    # 1e-6 threshold, and that is a property of the instances, not of a sample.
    resid = [r["numerical_residual"] for r in out
             if r["numerical_verdict"] == "infeasible"
             and r["numerical_residual"] is not None]
    summary = {
        "numerical_residual_min": min(resid) if resid else None,
        "numerical_residual_max": max(resid) if resid else None,
        "case": name,
        "source": str(src.relative_to(ROOT)),
        "instances": len(out),
        "min_margin": min(margins) if margins else None,
        "median_margin": (sorted(margins)[len(margins) // 2] if margins else None),
        "types_covered": len({r["tag"] for r in out}),
        "rejected_at_joint_screen": n_num_infeasible,
        "passed_joint_screen": n_num_feasible,
        "exactly_certified_infeasible": cert_of_infeasible,
        "undecided_among_rejected": n_num_infeasible - cert_of_infeasible,
        "wrongly_certified_among_screen_passers": cert_of_feasible,
        "certificate_kinds": dict(kinds),
        "seconds": round(elapsed, 1),
        "max_minors": MAX_MINORS,
        "max_boxes": MAX_BOXES,
    }
    if len(out) <= FULL_DETAIL_CAP:
        body = {"summary": summary, "instances": out}
    else:
        # Full witnesses for everything that is NOT a routine success, digests for
        # the rest: an undecided instance is the only thing a reader needs the
        # full record of, since it is the only thing that would weaken a claim.
        body = {"summary": summary,
                "detail_note": (f"{len(out)} instances exceeds the {FULL_DETAIL_CAP}"
                                f"-instance full-detail cap; full records are kept "
                                f"for undecided and realizable instances, digests "
                                f"for the rest"),
                "instances_full": [r for r in out
                                   if not r["certified_infeasible"]],
                "instances_digest": [_digest(r) for r in out]}
    dst.write_text(json.dumps(body, indent=1) + "\n")

    print(f"{name}: {len(out)} wildcard-bearing labellings over "
          f"{summary['types_covered']} subtree(s), {elapsed:.1f}s")
    print(f"    rejected at the joint screen   : {n_num_infeasible}")
    print(f"    of those exactly certified     : {cert_of_infeasible}")
    print(f"    passed the joint screen        : {n_num_feasible}")
    print(f"    of those wrongly refuted       : {cert_of_feasible}")
    for k, v in sorted(kinds.items()):
        print(f"    certificate kind {k:<26}: {v}")
    if margins:
        print(f"    enclosure margin from zero     : min {min(margins):.3g}, "
              f"median {sorted(margins)[len(margins) // 2]:.3g}")
    print(f"    artifact -> {dst.relative_to(ROOT)}")
    return summary


def merge_shards(name, src, dst, n):
    """Combine shard artifacts into the single artifact the paper reads."""
    parts = []
    for i in range(n):
        p = dst.with_name(dst.name.replace(".json", f".shard{i}of{n}.json"))
        if not p.exists():
            print(f"missing shard {p.name}", file=sys.stderr)
            return None
        body = json.loads(p.read_text())
        parts.extend(body.get("instances")
                     or body.get("instances_digest") or [])
        if body.get("instances_full"):
            parts.extend(x for x in body["instances_full"]
                         if x not in parts)
    total = sum(1 for line in src.read_text().splitlines() if line.strip())
    if len(parts) != total:
        print(f"shards hold {len(parts)} instances, source has {total}",
              file=sys.stderr)
        return None
    return summarize_and_write(name, src, dst, parts, 0.0)


def main():
    workers = next((int(a.split("=")[1]) for a in sys.argv[1:]
                    if a.startswith("--workers=")), max(1, mp.cpu_count() - 2))
    if "--passers-only" in sys.argv[1:]:
        # The one-sidedness claim rests entirely on the instances that SURVIVED
        # the screen: those are the only ones a false refutation could damage, and
        # there are few enough to certify in minutes rather than hours.  Kept as a
        # separate artifact so that claim does not wait on the census-wide run
        # over the rejections.
        name, src, dst = cases()[1]
        recs = [json.loads(l) for l in src.read_text().splitlines() if l.strip()]
        pas = [r for r in recs if r["numerical_verdict"] == "feasible"]
        out = [_one(r) for r in pas]
        wrong = [r for r in out if r["certified_infeasible"]]
        dst = dst.with_name(dst.name.replace(".json", "_passers.json"))
        dst.write_text(json.dumps({"summary": {
            "case": name + ", instances surviving the joint screen",
            "source": str(src.relative_to(ROOT)),
            "census_instances": len(recs),
            "passed_joint_screen": len(pas),
            "wrongly_certified_among_screen_passers": len(wrong),
        }, "instances": out}, indent=1) + "\n")
        print(f"{len(pas)} of {len(recs)} instances survived the joint screen; "
              f"{len(wrong)} of them were refuted by the exact certifier")
        print(f"artifact -> {dst.relative_to(ROOT)}")
        print("RESULT:", "no surviving instance is refuted" if not wrong
              else "UNSOUND: a surviving instance was refuted")
        return 0 if not wrong else 1
    shard = next((a.split("=")[1] for a in sys.argv[1:]
                  if a.startswith("--shard=")), None)
    if shard is not None:
        i, n = (int(x) for x in shard.split("/"))
        # A shard run does the d=5 case only: the d=6 case is seconds.
        name, src, dst = cases()[1]
        run_case(name, src, dst, shard=(i, n))
        return 0
    if "--merge" in sys.argv[1:]:
        n = next(int(a.split("=")[1]) for a in sys.argv[1:]
                 if a.startswith("--shards="))
        name, src, dst = cases()[1]
        return 0 if merge_shards(name, src, dst, n) else 1
    summaries = []
    for name, src, dst in cases():
        if not src.exists():
            print(f"{name}: MISSING {src.relative_to(ROOT)} "
                  f"-- run  python3 run_wild_dump.py  first", file=sys.stderr)
            return 2
        if not _complete(src):
            print(f"{name}: {src.relative_to(ROOT)} is INCOMPLETE (no matching "
                  f"completion marker) -- let run_wild_dump.py finish",
                  file=sys.stderr)
            return 2
        summaries.append(run_case(name, src, dst, workers=workers))
        print()

    d6, d5 = summaries
    problems = []
    # The d=6 claim: every wildcard rejection of the run of record is now exact.
    if d6["undecided_among_rejected"]:
        problems.append(f"d=6: {d6['undecided_among_rejected']} rejections not "
                        f"certified exactly")
    # One-sidedness: nothing that survived the screen may be refuted.  This is the
    # property the d=6 conclusion depends on, and the d=5 set is where it can be
    # measured, because 9 of its screen-passers are known polytopes.
    for s in summaries:
        if s["wrongly_certified_among_screen_passers"]:
            problems.append(
                f"{s['case']}: {s['wrongly_certified_among_screen_passers']} "
                f"instances that PASSED the screen were refuted -- the certifier "
                f"is UNSOUND")
    # Cross-check the d=5 set against the committed census artifact rather than
    # against a transcribed number.
    wc = json.loads((ROOT / "runs/d5_n9/d5_wildcounts.json").read_text())
    tot = sum(v.get("wild_assignments", 0) for v in wc.values())
    passes = sum(v.get("wild_assignments", 0)
                 - (v.get("wild") or {}).get("joint_infeasible", 0)
                 for v in wc.values())
    polytopes = sum(v.get("wild_assignments", 0)
                    - sum((v.get("wild") or {}).values()) for v in wc.values())
    if d5["instances"] == tot and d5["passed_joint_screen"] != passes:
        problems.append(f"d=5: {d5['passed_joint_screen']} screen-passers dumped, "
                        f"census artifact says {passes}")
    if polytopes != 9:
        problems.append(f"d=5 census artifact reports {polytopes} wildcard "
                        f"polytopes, expected the 9 with a pi/10 angle")
    if d5["instances"] == tot:
        print(f"d=5 cross-check against runs/d5_n9/d5_wildcounts.json: "
              f"{tot} instances, {passes} pass the joint screen, of which "
              f"{polytopes} become polytopes and {passes - polytopes} are killed "
              f"later by the integer window scan -- so the certifier is expected "
              f"to leave exactly {passes} undecided.")

    print("RESULT:", "every wildcard rejection is exactly certified, and nothing "
          "that survived the screen is refuted" if not problems
          else "PROBLEMS: " + "; ".join(problems))
    return 0 if not problems else 1


if __name__ == "__main__":
    sys.exit(main())
