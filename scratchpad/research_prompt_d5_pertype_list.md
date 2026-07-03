# Research task (for a fresh Claude instance with web access / deep-research)

## Goal
Recover the **per-combinatorial-type breakdown** of the 51 compact hyperbolic Coxeter
**5-polytopes with 9 facets**, so we can localize a discrepancy in our own pipeline
(we reproduce only 12 of the 51 and need to know *which* combinatorial types account for
the missing ~39).

## Exactly what we need (in priority order)
1. **Which of the 322 combinatorial types (equivalently, which of the 109 types with ≥2
   disjoint facet pairs) realize at least one compact polytope, and how many each.**
   The ideal answer is a table: `combinatorial-type label → number of compact polytopes`.
2. If (1) isn't tabulated, then **the list of all 51 polytopes with, for each, an
   identifier of its combinatorial type** (any label that can be mapped to a type: the
   Fukuda–Miyata–Moriyama order-type index 1..322, the paper's P_k label, the vertex/facet
   incidence "vertex flag", or the number of disjoint facet pairs k).
3. Failing both, **the distribution of the 51 by number of disjoint facet pairs k**
   (how many of the 51 have k = 2, 3, 4, 5, 6). And the **maximum number of polytopes a
   single combinatorial type produces**.
4. Any **explicit Coxeter diagrams / Gram data** for individual d=5 polytopes (we can feed
   these straight through our pipeline to see where they are dropped) — especially for the
   low-k (k=2,3) types.

## Where to look (specific sources)
- **Primary:** Ma, Zheng, *Five-dimensional compact hyperbolic Coxeter polytopes with nine
  facets*, Transformation Groups (2023), **arXiv:2203.16049**. Get the **full PDF**, not
  the abstract. The results are in its **Section 7 ("Validation and Results")** and any
  appendix — that is where the complete list of Coxeter diagrams / hyperbolic lengths of
  Theorem 1.1 lives. Table 2 groups the 109 feasible types by disjoint-pair count.
- **Cross-check:** A. Burcroff, *Near classification of compact hyperbolic Coxeter
  d-polytopes with d+4 facets…*, Eur. J. Combin. 120 (2024) 103957, **arXiv:2201.03437**;
  and her Durham MSc thesis (2021), **etheses.dur.ac.uk/14202** — independently classified
  d=5 (she and Ma–Zheng agree on 51). Her tables may give the per-type list too.
- **Data repo:** `github.com/GeoTopChristy/HCPdm` — check for any final polytope list (not
  just the `output/P9_*` round-3 *candidate* files, which are pre-final-solve and, in our
  local clone, append-polluted). Look in the repo, its releases, and the paper's ancillary
  arXiv files (arxiv.org/e-print/2203.16049 or the "Other formats" / ancillary-files link).

## Context that makes the answer usable
- The 322 order types come from Fukuda–Miyata–Moriyama (oriented-matroid enumeration of
  5-polytopes with 9 facets). "k" = number of disjoint (non-adjacent) facet pairs =
  number of size-2 missing faces. Felikson–Tumarkin: a compact d+4 polytope needs k ≥ 2,
  so only the **109** types with k ≥ 2 can realize; the paper confirms this.
- We have verified our generator produces exactly these 109 types and that our k = 4,5,6
  counts (15, 7, 3 types) match the paper's Table 2 exactly. Our k=2/k=3 split is 50/34;
  the paper's Table 2 appears to be ~54/30 — **please confirm the exact k=2 and k=3 type
  counts and, if possible, the type labels in each group**, since a mismatch there would
  itself be a lead.

## Output format wanted
- A markdown table `type_label | k | #polytopes` covering every realizing type (or as much
  as the sources give), plus the k-distribution of the 51, plus the single largest per-type
  count. Cite the exact source (paper section/table number, page, or repo file) for each
  fact. Flag anything that is inferred vs directly stated. Note explicitly if Section 7's
  full list could not be obtained and what the blocker was.
