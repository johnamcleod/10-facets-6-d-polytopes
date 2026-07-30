#!/usr/bin/env python3
"""Is the polytope we classify the one drawn in [Burcroff 2024, Figure 5]?

Burcroff's Figure 5 is captioned "The Coxeter diagram of the compact Coxeter
6-polytope with 10 facets constructed by Bugaenko".  Its edge multiplicities are
read off the PDF's vector drawing rather than by eye: each ordinary edge is drawn
as a stack of strokes of decreasing width, and by her own conventions (Section 4
of that paper) a single, double and triple line means m = 3, 4 and 5
respectively, while a dashed line means the facets diverge.  Transcribing the
line segments from the figure on page 37 of arXiv:2201.03437 gives, naming her
nodes C,D,E,F,G,A,B,I,J,H by position (left-to-right, top-to-bottom):

    triple (m=5) : C-D, H-F, H-J
    double (m=4) : D-E, D-I
    single (m=3) : D-A, A-B, B-G, F-G, J-G
    dashed       : E-I, E-F, I-J

This script builds that labelled graph and the labelled graph of the diagram in
our paper, and searches all 10! relabellings for an isomorphism preserving every
edge label ("dashed" being its own label, since the dashed weights are determined
by the rank condition and are not drawn).

Expected output: ISOMORPHIC: True, with the explicit node correspondence.

The SAME drawing appears earlier as Figure 5 of Burcroff's MSc thesis
(etheses.durham.ac.uk/14202, page index 57 of the PDF).  Its 21 line segments are
the paper figure's segments up to a rescaling and translation -- 3 triple edges at
widths (2.54, 1.81, 0.36), 2 double at (1.45, 0.72), 5 single and 3 dashed at
0.29 -- so the thesis and the paper carry one drawing between them, not two
independent ones.  That matters for provenance: see REVIEW_NOTES.md section 1.

To re-derive the segment list instead of trusting the transcription above:

    import fitz
    page = fitz.open("2201.03437.pdf")[36]
    for g in page.get_drawings():
        for it in g["items"]:
            if it[0] == "l":
                print(it[1], it[2], g.get("width"), g.get("dashes"))

Segments sharing endpoints and differing only in width are one multiple edge:
three strokes = triple, two = double, one = single, dash pattern = dashed.

Run:  python3 paper/checks/burcroff_fig5_isomorphism.py
"""
import itertools
import sys

import sympy as sp
from sympy import Rational, S, sqrt

N = 10
DASH = "DASHED"


def entry(m):
    return {2: S.Zero, 3: Rational(-1, 2), 4: -sqrt(2) / 2,
            5: -(1 + sqrt(5)) / 4}[m]


# ---- the diagram of our paper (nodes 1..10) -------------------------------
OURS_ORD = [(1, 2, 5), (2, 3, 3), (3, 4, 3), (4, 5, 3), (5, 6, 3),
            (9, 10, 5), (2, 7, 4), (2, 8, 4), (5, 9, 3), (6, 10, 5)]
OURS_DASH = [(6, 7), (7, 8), (8, 9)]

# ---- Burcroff, Figure 5, transcribed from the PDF vector drawing ----------
BUR = {"C": 0, "D": 1, "E": 2, "F": 3, "G": 4,
       "A": 5, "B": 6, "I": 7, "J": 8, "H": 9}
BUR_5 = [("C", "D"), ("H", "F"), ("H", "J")]
BUR_4 = [("D", "E"), ("D", "I")]
BUR_3 = [("D", "A"), ("A", "B"), ("B", "G"), ("F", "G"), ("J", "G")]
BUR_DASH = [("E", "I"), ("E", "F"), ("I", "J")]


def labelled(ordinary, dashed, index=None):
    M = [[S.Zero] * N for _ in range(N)]
    for a, b, m in ordinary:
        i, j = (index[a], index[b]) if index else (a - 1, b - 1)
        M[i][j] = M[j][i] = entry(m)
    for a, b in dashed:
        i, j = (index[a], index[b]) if index else (a - 1, b - 1)
        M[i][j] = M[j][i] = DASH
    return M


def edge_labels(M):
    return sorted(str(M[i][j]) for i in range(N) for j in range(i + 1, N)
                  if M[i][j] != S.Zero)


def main():
    A = labelled(OURS_ORD, OURS_DASH)
    B = labelled([(a, b, 5) for a, b in BUR_5]
                 + [(a, b, 4) for a, b in BUR_4]
                 + [(a, b, 3) for a, b in BUR_3],
                 BUR_DASH, index=BUR)

    same = edge_labels(A) == edge_labels(B)
    print("our diagram      : %d edges, labels %s"
          % (len(edge_labels(A)), edge_labels(A)))
    print("Burcroff Fig. 5  : %d edges, labels %s"
          % (len(edge_labels(B)), edge_labels(B)))
    print("edge-label multisets agree:", same)

    found = None
    for perm in itertools.permutations(range(N)):
        if all(A[i][j] == B[perm[i]][perm[j]]
               for i in range(N) for j in range(i + 1, N)):
            found = perm
            break

    print()
    print("ISOMORPHIC:", found is not None)
    if found is not None:
        inv = {v: k for k, v in BUR.items()}
        print("our node -> Burcroff Fig. 5 node:",
              {i + 1: inv[found[i]] for i in range(N)})
        w67, w78 = 2 * sqrt(2) + sqrt(10), 17 + 8 * sqrt(5)
        G = sp.eye(N)
        for a, b, m in OURS_ORD:
            G[a - 1, b - 1] = G[b - 1, a - 1] = entry(m)
        for (a, b), w in {(6, 7): w67, (7, 8): w78, (8, 9): w67}.items():
            G[a - 1, b - 1] = G[b - 1, a - 1] = -w
        print("our Gram matrix exact rank:", G.rank(), "(expected 7)")

    print()
    print("RESULT:", "our polytope IS the one drawn in Burcroff Figure 5"
          if (found is not None and same) else "NOT isomorphic")
    return 0 if (found is not None and same) else 1


if __name__ == "__main__":
    sys.exit(main())
