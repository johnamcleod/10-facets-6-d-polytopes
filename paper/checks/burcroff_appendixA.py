#!/usr/bin/env python3
"""Reproduce Remark 3.5 of the paper: three corrections to [Burcroff 2024,
Appendix A], the list of 111 candidate combinatorial types H_1..H_111 of simple
5-polytopes with 9 facets having at least two pairs of disjoint facets.

Findings the paper asserts, each checkable here:

  (i)   H_31 and H_32 carry IDENTICAL missing-face lists, so the list has 110
        distinct entries, not 111.
  (ii)  H_7 and H_35 are not the missing-face systems of simple polytopes: in a
        simple d-polytope every (d-1)-subset of facets that is a face lies in
        exactly two vertices, and these violate that.
  (iii) One combinatorial type present both in the Ma-Zheng ground truth
        (data/ground_truth/5d9m.txt) and in our generator's output -- our tid89
        -- is absent from the 111.

None of this affects Burcroff's classification: H_7 and H_35 are recorded there
as realizing no polytope, and tid89 realizes none in our run either.  It matters
only because it shows our generator's 109 types are a strict SUPERSET of the
valid published candidates, which is the direction completeness requires.

INPUT.  This script needs the missing-face lists of H_1..H_111.  They are
transcribed below from arXiv:2201.03437, Appendix A.  The transcription was
produced by machine text extraction and then CONFIRMED with a second, independent
extractor (PyMuPDF and poppler's pdftotext agree entry for entry), so it is not
a rendering artefact of one tool.  To re-extract from scratch:

    pdftotext -layout 2201.03437.pdf out.txt
    # then regex  H(\\d+)\\s+((?:\\d+[,\\s]+)*\\d+)\\s+(\\d+)  over the Appendix A block

Run:  python3 paper/checks/burcroff_appendixA.py
"""
import itertools
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

N, D = 9, 5

# The three entries the paper singles out, transcribed from Appendix A.
H = {
    7:  "01,02,03,12456,1457,267,268,38,45678",
    8:  "01,02,03,12456,1457,268,37,38,45678",     # control: a valid neighbour
    31: "01,02,03,124,1456,278,356,378,45678",
    32: "01,02,03,124,1456,278,356,378,45678",
    35: "01,02,034,15,256,267,3478,5678",
    1:  "01,02,03,12,14,25678,34,35678,45678",     # control: H1 realizes 22
}


def parse(s):
    return [frozenset(int(c) for c in tok.strip()) for tok in s.split(",")]


def vertices(mf):
    """Vertices = the d-subsets of facets containing no missing face."""
    def is_face(S):
        return not any(m <= S for m in mf)
    return [frozenset(c) for c in itertools.combinations(range(N), D)
            if is_face(frozenset(c))]


def ridge_defects(mf):
    """(#vertices, multiplicity histogram of the (d-1)-subsets that are faces).

    In a simple d-polytope each such subset is an edge and lies in exactly two
    vertices, so any multiplicity other than 2 rules out simple polytopality."""
    V = vertices(mf)
    ridge = Counter()
    for v in V:
        for r in itertools.combinations(sorted(v), D - 1):
            ridge[frozenset(r)] += 1
    return len(V), Counter(ridge.values())


def main():
    ok = True

    print("(i) duplicate entry")
    same = parse(H[31]) == parse(H[32])
    print(f"    H_31 and H_32 have identical missing-face lists: {same}")
    print(f"    -> Appendix A has {'110' if same else '111'} distinct entries")
    ok = ok and same

    print()
    print("(ii) non-polytopal entries "
          "(a simple 5-polytope needs every ridge in exactly 2 vertices)")
    for lab in (1, 8, 7, 35):
        nv, hist = ridge_defects(parse(H[lab]))
        bad = sum(c for m, c in hist.items() if m != 2)
        verdict = "simple polytope" if bad == 0 else "NOT a simple polytope"
        print(f"    H_{lab:<3d}: {nv:3d} vertices, ridge multiplicities "
              f"{dict(sorted(hist.items()))}  -> {verdict}")
    ok = ok and all(sum(c for m, c in ridge_defects(parse(H[l]))[1].items()
                        if m != 2) > 0 for l in (7, 35))
    ok = ok and all(sum(c for m, c in ridge_defects(parse(H[l]))[1].items()
                        if m != 2) == 0 for l in (1, 8))

    print()
    print("(iii) a valid type absent from Appendix A")
    p = ROOT / "runs/d5_n9/stage2/types.json"
    if not p.exists():
        print("    SKIPPED (runs/d5_n9/stage2/types.json not present)")
    else:
        t = next(x for x in json.load(open(p)) if x["type_id"] == 89)
        mf = [frozenset(m) for m in t["missing_faces"]]
        nv, hist = ridge_defects(mf)
        good = set(hist) == {2}
        print(f"    our tid89: missing faces "
              f"{','.join(''.join(map(str, sorted(m))) for m in mf)}")
        print(f"    {nv} vertices, ridge multiplicities "
              f"{dict(sorted(hist.items()))}  -> "
              f"{'simple polytope' if good else 'NOT a simple polytope'}")
        print("    (present in data/ground_truth/5d9m.txt and in our 109; "
              "absent from H_1..H_111)")
        ok = ok and good

    print()
    print("RESULT:", "all three corrections reproduce" if ok else "MISMATCH")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
