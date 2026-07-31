# Compact hyperbolic Coxeter 6-polytopes with 10 facets

Code, data and machine-checked certificates for the paper

> **Compact hyperbolic Coxeter 6-polytopes with 10 facets, and the completion of
> the d+4 classification** (`paper/paper.tex`)

**Result.** Up to isometry there is exactly one compact hyperbolic Coxeter
6-polytope with 10 facets — the polytope drawn in [Burcroff 2024, Fig. 5], which
our pipeline recovers independently from raw combinatorics. Combined with Felikson–Tumarkin (none for `d ≥ 8`,
unique for `d = 7`) and Burcroff / Ma–Zheng (`d = 4`: 348, `d = 5`: 51), this
completes the classification of compact hyperbolic Coxeter `d`-polytopes with
`d+4` facets in every dimension.

Please read [`REVIEW_NOTES.md`](REVIEW_NOTES.md) alongside the paper: it lists
every place where the text's claims outrun what the data strictly establish.

---

## 0. Script menu — which entrypoint does what

The repository accumulated one driver per phase of the work, and several are
superseded. This is the map; everything below is elaborated in §1 and §3.

**Verify the published results (no re-run, no pipeline import).**

| script | what it does | runtime |
|---|---|---|
| `verify_polytope.py` | the whole result in one command: exact signature of P^B6 over Q(√2,√5), no parabolic subdiagram, weights > 1, CoxIter cross-check | ≈ 3 min |
| `paper/make_tables.py` | re-derives **every** number the paper quotes from the run certificates, recomputes the exhaustion coverage recursion, and fails if any verdict rests on a truncation | seconds |
| `verify_diagram.py` | same checks for an *arbitrary* Coxeter diagram you supply | seconds |

**Re-run a census.** All three dimensions use the same driver; `d=` selects the
dimension and `out=` names a fresh state directory (without it, cached subtrees
are resumed and nothing is recomputed).

| command | what it produces | cost |
|---|---|---|
| `run_survivors_rigorous.py all <nproc> wildcard` | the **d=6** classification of record → `runs/d6_n10/<out>/` | 12 CPU-h |
| `run_survivors_rigorous.py all <nproc> wildcard d=4 out=…` | the **d=4** census, 348 polytopes over 30 types | hours |
| `validate_d5_wildcard.py <nproc> out=…` | the **d=5** census, 51 polytopes over 109 types | 1.4 CPU-h |
| `run_full_pipeline_d6.py` | Stages 2–4 end to end, starting from the order-type database (needs `otypes10.b16`) | hours |
| `apply_facet_profile_filter.py` | the facet filter, 387 types → 54 survivors | minutes |
| `python -m pipeline.validate_coverage --d 4 / --d 5` | the type generator reproduces the published 30 / 109 types exactly | minutes |
| `run_d4.py`, `run_d4_parallel.py`, `run_blockpaste_parallel.py` | earlier d=4 drivers (single-process, pooled, and Ma–Zheng block-pasting) | hours |
| `run_regression.py` | d=4 → 348 and d=5 → 51 as a regression gate | hours |

**Exact refutations (the screen replaced by certificates).**

| command | what it does | cost |
|---|---|---|
| `run_wild_dump.py` | re-runs the d=6 subtree with the instance sinks on, writing every labelling that reached a screen: 406 wildcard-bearing and 546 without | ≈ 18 min |
| `run_wild_dump.py d=5 out=…` | the same for the d=5 census, as the calibration set | ≈ 80 min |
| `paper/checks/wild_exact_certify.py` | exact infeasibility certificate for each **wildcard** labelling, plus the d=5 refutation-power measurement | seconds |
| `paper/checks/wild_exact_certify.py --passers-only` | census-wide one-sidedness: the 11 d=5 instances that survive the screen must **not** be refuted | ≈ 6 min |
| `paper/checks/plain_exact_certify.py` | exact certificate for each **non-wildcard** labelling | seconds |
| `pytest tests/test_exact_certify.py` | the certifiers' own tests, including that they do **not** refute P^B6 | seconds |

**Claim-by-claim checks.** `paper/checks/*.py`, one verdict line each — see §1.3
for the table of which claim each substantiates.

**Superseded, retained for provenance.** These produced intermediate results and
are not the route to any number in the paper: `run_survivors_fulllabel.py` (its
"complete" tag was a wall-clock heuristic), `run_stage4_d6.py`,
`run_stage4_d6_full.py`, `run_stage4_d6_noext.py`, `run_d6_ext0.py`,
`run_uncertain_final.py` with `analyze_final_results.py`, `run_targeted_stage4.py`,
`run_incomplete_backtrack.py`, `run_lowk_parallel.py`, `check_exhaustion.py`,
`apply_f5_filter.py`, `run_d5_stage4.py`, `run_d5_campaign.py`,
`validate_d5_fulllabel.py`, `run_truth_stage4.py` (fed published types in while the
generator was broken), `harvest_rejects.py`, `bench_stage1.py`, `probe_cost.py`,
`mz_solve_p8_17.py` (now driven by `paper/checks/mz_p8_17_crosscheck.py`).

---

## 1. How to verify the results

Everything needed to check the paper is in this repository: the verification
scripts, the claim-by-claim checks, and the run certificates the numbers are
derived from. None of it requires re-running the classification, and the two
commands in §1.1 and §1.2 do not import the pipeline at all.

### 1.1 The polytope (≈ 3 minutes)

```bash
python3 -m venv .venv && . .venv/bin/activate
pip install numpy scipy sympy mpmath      # pandas additionally for §3.4-3.5

python3 verify_polytope.py
```

Expected output ends with:

```
[1] EXACT inertia over Q(sqrt2, sqrt5) (congruence elimination)
    positive=6  negative=1  zero=3   -> signature (6,1), rank 7, nullity 3   OK
[2] EXACT rank and algebraic relation
    rank = 7 (expected 7);  w(7,8) - (w(6,7)^2 - 1) = 0   OK
[3] Ultraparallel weights > 1                                                OK
[4] Parabolic subdiagram scan (all 2^10-1 subsets)          none found       OK
[5] CoxIter    Cocompact: yes   Dimension: 6
               f-vector: (31, 93, 125, 95, 42, 10, 1)
               Euler characteristic: -67/288000
               Covolume: pi^3 * 67/540000                                    OK
[6] Agreement with the classification run's realizer records                 OK

  RESULT: ALL CHECKS PASSED
```

Step 5 is skipped unless CoxIter is built:

```bash
cd scratchpad/CoxIter && mkdir -p build && cd build && cmake .. && make
```

The polytope itself: facets `1..10`, ordinary edges `(i, j, m)` meaning dihedral
angle `π/m` (all unlisted intersecting pairs orthogonal)

```
(1,2,5) (2,3,3) (3,4,3) (4,5,3) (5,6,3) (9,10,5) (2,7,4) (2,8,4) (5,9,3) (6,10,5)
```

and ultraparallel ("dotted") pairs `(6,7)`, `(7,8)`, `(8,9)` with weights

```
w(6,7) = w(8,9) = 2√2 + √10 = √2(2+√5) ≈ 5.990705     minpoly x⁴ − 36x² + 4
w(7,8)          = 17 + 8√5           ≈ 34.888544      minpoly x² − 34x − 31
w(7,8) = w(6,7)² − 1        field of definition Q(√2, √5)
```

### 1.2 The classification: totals, coverage and the absence of truncations

```bash
python3 paper/make_tables.py
```

This reads the run certificates and re-derives every number the paper quotes. It
does not trust the stored verdicts: it recomputes the coverage recursion from the
subtree certificates and **fails** if any of the searched types is not closed, if
any verdict rests on a truncation at the refinement depth cap, on a wildcard
window that reached its scan edge, or on a per-assignment deadline, or if any
labelling that reached a screen lacks an exact refutation certificate. Expected
tail:

```
  NumTypes           387        NumSurvivors       54
  NumRequired        11         NumThreeFree       2      (types 159, 329)
  NumSubtrees        71         CPUHours           12
  NumAssignments     952        MaxDepth           2
  WildCertified      406        PlainCertified     545
  DFiveTotal         51         DFourFound         348
  RealizingType      379

all paper invariants re-checked OK
```

(`NumRequired` is 11 because Lemmas 4.2–4.4 removed types the earlier runs
searched; the 54 facet-filter survivors and the 387 generated types are unchanged.
`WildCertified` + `PlainCertified` + 1 accepted = the 952 labellings that reached a
screen.)

It also regenerates `paper/tables/pertype.tex` (Table 1) and
`paper/tables/summary_nums.tex`; those two files are the only route by which run
data reaches the paper, so no number in it is transcribed by hand.

### 1.3 Claim-by-claim checks

Each script substantiates one claim of the paper and prints a verdict line.

```bash
for f in paper/checks/*.py; do echo "== $f"; python3 "$f" | tail -2; done
```

| script | claim |
|---|---|
| `d4_census_reconciliation.py` | the d=4 census agrees with Burcroff Appendix A type by type (348) |
| `burcroff_appendixA.py` | three corrections to Burcroff's 111-entry d=5 candidate list |
| `burcroff_fig5_isomorphism.py` | our polytope is isomorphic to Burcroff Fig. 5, edge for edge and label for label |
| `lanner_label_bound.py` | re-derives Lannér's order-4 and order-5 diagrams (9 and 5) from the definition |
| `wildcard_soundness.py` | the m ≥ 7 wildcard predicate is insensitive to the value of m |
| `forward_check_margins.py` | forward-checking eigenvalue tolerances carry 10³–10⁶ margin |
| `cascade_reachability.py` | the bounded-box fallback is structurally unreachable at d=6; the pair-resultant path is never entered |
| `screen_margins.py` | margins at the cascade's decision thresholds |
| `mz_p8_17_crosscheck.py` | our d=4 solutions all appear in Ma–Zheng's own published candidate list |

Two notes on running them. `mz_p8_17_crosscheck.py` compares against Ma–Zheng's
published intermediate data, which is third-party material and is not redistributed
here; it prints `SKIPPED` with cloning instructions unless
`scratchpad/HCPdm` is present. `screen_margins.py` is a measurement rather than a
test and takes several minutes.

### 1.4 The certificates themselves

| path | contents |
|---|---|
| `runs/d6_n10/stage2/types.json` | the 387 combinatorial types, with missing faces and vertex sets |
| `runs/d6_n10/facet_profile_survivors.json` | the 54 survivors of the facet filter |
| `runs/d6_n10/d6_final/` | the d=6 classification of record: 71 subtree certificates over the 11 required types |
| `runs/d6_n10/d6_final/realizers/` | the realizer records for the unique polytope |
| `runs/d6_n10/d6_tangency/`, `d6_rest/` | the superseded two-part run over 52 types (1,156 subtrees, 173 CPU-h), retained for provenance; `make_tables.py` falls back to it only if `d6_final/` is absent |
| `runs/d6_n10/wild_instances.jsonl`, `plain_instances.jsonl` | every labelling that reached a screen, 406 + 546, with the verdict each was rejected on |
| `runs/d6_n10/wild_certificates.json`, `plain_certificates.json` | the exact refutation certificate for each of them |
| `runs/d5_n9/wild_certificates.json` | the d=5 calibration: the certifier refutes the rejections and none of the 9 realizable instances |
| `runs/d5_n9/d5_discfix.json` | the d=5 census, 51 polytopes over 109 types |
| `runs/d4_n8/survivors_discfix/` | the d=4 census, 348 polytopes over 30 types |

A subtree certificate records its prefix, whether the subtree was exhausted, how
many labellings it enumerated, the polytopes found, and per-gate diagnostics. A
type is decided when its root is exhausted, or when all six children are
recursively decided — the recursion `make_tables.py` recomputes.

### 1.5 The companion note on the polytope

`paper/polytope.tex` compiles standalone (`pdflatex polytope.tex`) and gives the
Coxeter diagram, the exact Gram matrix over Q(√2, √5), the automorphism computation
behind "exactly one up to isometry", the identification with Burcroff Fig. 5, and
the polytope's position relative to the Felikson–Tumarkin class.

---

## 2. Environment

Everything in the paper was produced on a single **Apple M1, 8 cores, 16 GB RAM,
macOS 26.5.2**, with **Python 3.10.0** and `numpy`, `scipy`, `sympy`, `mpmath`,
`pandas`, plus **CoxIter** built from source
(`scratchpad/CoxIter`, github.com/rgugliel/CoxIter).

---

## 3. Re-running the classification

### 3.1 Combinatorial types (Stage 1–2) — 387 types

Already committed as `runs/d6_n10/stage2/types.json` (387 entries). Regenerating
them requires the AAK order-type database file `otypes10.b16`
(14,309,547 records, ~570 MB) from

<http://www.ist.tugraz.at/aichholzer/research/rp/triangulations/ordertypes/>

converted once to the chirotope cache with the C parser in `pipeline/c/`. This
step takes several hours; the committed output makes it optional.

### 3.2 Generator validation — d=4 (30/30) and d=5 (109/109)

```bash
python3 -m pipeline.validate_coverage --d 4
python3 -m pipeline.validate_coverage --d 5
```

(These default to `runs/d{d}_n{n}/stage2`, the authoritative type lists. Older
`.../stage3` directories in the tree are stale pre-2026-06 artifacts; passing
`--stage3 runs/d4_n8/stage3` will report a spurious 27/30.)

Expected:

```
d=4:  Ground-truth 30 distinct;  ours 30 distinct
      COVERED 30/30   MISSING 0   SPURIOUS 0

d=5:  ours all: [(2,50),(3,34),(4,15),(5,7),(6,3)]   (= 109, all p ≥ 2)
      SPURIOUS 0;  the only MISSING truth types are the p = 0 and p = 1 ones,
      which Felikson–Tumarkin exclude for n = d+4.
```

### 3.3 Facet filter — 387 → 54 (≈ 20 minutes)

```bash
python3 apply_facet_profile_filter.py
```

Expected tail:

```
survivors (54): [8, 12, 17, 34, ..., 378, 379, 382]
killed: 333
matches hardcoded 54-list: True
```

Rewrites `runs/d6_n10/facet_profile_survivors.json`.

### 3.4 The exhaustive d=6 search — 54 → 1 (≈ 2,300 CPU-hours)

```bash
caffeinate -i python3 run_survivors_rigorous.py all 7 wildcard
```

Resumable: completed subtrees are checkpointed to
`runs/d6_n10/survivors_wildcard/state.json` and skipped on restart. Running it
against the committed state file reproduces the verdicts immediately without
recomputation. Expected final line:

```
RIGOROUS-DONE. realizing=[379] not-rigorous=[]
```

Types 159 and 329 are excluded by theorem (Esselmann's bound for 3-free
polytopes; see `REVIEW_NOTES.md` §2) and are not searched.

### 3.4a Other dimensions through the same driver

The driver takes an optional `d=D` argument (default 6, so every command above is
unaffected). For `D != 6` it reads `runs/dD_n{D+4}/stage2/types.json`, keeps its
own state directory, and takes the candidate list from that dimension's
`facet_profile_survivors.json` if present, else all generated types:

```bash
python3 run_survivors_rigorous.py all 7 wildcard d=4
```

Two dimension-aware details:

* **The Esselmann 3-free elimination applies only when `n < 2d`.** A 3-free
  compact Coxeter `d`-polytope needs at least `2d` facets, with equality only for
  the `d`-cube. For d=6, `10 < 12`, so 3-free types are killed. For d=4,
  `8 = 2·4`, so a 3-free type could be the 4-cube — which does admit compact
  Coxeter structures — and those types are **searched, not killed**.
* The known-realizer-first ordering (type 379) is a d=6 anchor and is skipped
  elsewhere.

**d=4 has not been run this way yet** — the published 338 came from the fixed
`{2,…,10,12}` alphabet, not the wildcard path. See `REVIEW_NOTES.md` §4.

### 3.5 The d=5 census through the identical code path — 51/51 (≈ 7 CPU-hours)

This is the paper's main soundness anchor.

```bash
python3 validate_d5_wildcard.py 7
```

Expected final line:

```
D5-WILDCARD-DONE. total=51 (expect 51)  problem-types=NONE
```

Output: `runs/d5_n9/wildcard_validation.json` — 109 types, every one
`exhausted: true` and `unbounded: false`, realizing counts
`{0: 22, 6: 18, 19: 6, 29: 3, 5: 1, 63: 1}`.

Repeat it with Burcroff's Lemma 5.5(b) low-weight caps enabled — the one solver
flag the d=6 run sets that the baseline above does not — so that the anchor
exercises the exact d=6 flag combination (≈ 2.9 CPU-hours):

```bash
python3 validate_d5_wildcard.py 7 55b
```

Same expected final line; output goes to
`runs/d5_n9/wildcard_validation_55b.json`.

### 3.6 Supporting checks

```bash
python3 paper/checks/wildcard_soundness.py        # m=7 proxy lossless; 0 discrepancies
python3 paper/checks/lanner_label_bound.py        # re-derives Lannér's 9 and 5 diagrams
python3 paper/checks/cascade_reachability.py      # bounded-box fallback unreachable
python3 paper/checks/forward_check_margins.py     # forward-check tolerances, exhaustive
python3 paper/checks/screen_margins.py            # cascade tolerances, per-type sample
python3 paper/checks/burcroff_appendixA.py        # the three Appendix A corrections
python3 paper/checks/burcroff_fig5_isomorphism.py # our diagram ≅ Burcroff Fig. 5
python3 paper/checks/wild_exact_certify.py        # exact refutation of every wildcard labelling
python3 paper/checks/plain_exact_certify.py       # exact refutation of every non-wildcard labelling
python3 paper/make_tables.py                      # regenerates the paper's tables
```

The two `*_exact_certify.py` scripts read instance dumps; regenerate those first
with

```bash
python3 run_wild_dump.py            # d=6: re-runs subtree 379|0,0 (~18 min), 406 + 546 instances
python3 run_wild_dump.py d=5 out=runs/d5_n9/wild_instances_full.jsonl   # the d=5 calibration census
```

One run writes both sinks: `runs/d6_n10/wild_instances.jsonl` (the 406
wildcard-bearing labellings) and `runs/d6_n10/plain_instances.jsonl` (the 546
without), each with a `.meta.json` completion marker that the certifiers require —
so a dump still being written cannot be certified and reported as complete. The
re-run must reproduce the run of record exactly (952 labellings, 406
wildcard-bearing, 1 realizer, exhausted); the dump is instrumentation only.
The unit tests of the certifiers, including the check that they do **not** refute
P^B6, are `tests/test_exact_certify.py`.

All except `screen_margins.py` run in seconds to a few minutes.
`screen_margins.py` re-runs whole types under instrumentation, so pass it a
short list of cheap type ids (it defaults to the sub-CPU-hour survivors).

`make_tables.py` also re-asserts every headline number in the paper and prints
`all paper invariants re-checked OK`.

---

## 4. Building the paper

```bash
cd paper
pdflatex paper && bibtex paper && pdflatex paper && pdflatex paper
```

Compiles clean (no warnings, no undefined references) with TeX Live 2026.

---

## 5. Where the data live

| path | contents |
|---|---|
| `runs/d6_n10/stage2/types.json` | the 387 combinatorial types (missing faces, `p`, a realizing affine Gale diagram) |
| `runs/d6_n10/facet_profile_survivors.json` | the 54 survivors of the facet filter, plus the witness that killed each of the 333 others |
| `runs/d6_n10/survivors_wildcard/state.json` | 20,466 subtree certificates: prefix, `exhausted`, labellings enumerated, seconds |
| `runs/d6_n10/survivors_wildcard/verdicts.json` | per-type verdict (`distinct`, `rigorous`, `subtrees`) |
| `runs/d6_n10/survivors_wildcard/realizers/` | the 12 distinct Gram configurations of type 379 with exact minimal polynomials |
| `runs/d5_n9/wildcard_validation.json` | the d=5 revalidation, 109 types |
| `runs/d5_n9/stage2/types.json` | the 109 d=5 combinatorial types with `p ≥ 2` |
| `runs/d4_n8/blockpaste/` | the d=4 per-type counts (338 total — **not** the published 348; see below) |
| `data/ground_truth/{4d8m,5d9m}.txt` | published combinatorial censuses used by `validate_coverage` |
| `paper/tables/` | LaTeX tables regenerated from all of the above |

## 6. Key source files

| path | role |
|---|---|
| `pipeline/utils/gale_exact.py` | exact-rational affine Gale face criterion |
| `pipeline/stage2_gale.py` | order types → Gale diagrams → combinatorial types, exact dedup |
| `pipeline/stage4_gram.py` | the search: enumeration, forward checking, wildcards, cascade screen, exact certification |
| `pipeline/utils/automorphisms.py` | VF2 automorphism group, orbit symmetry breaking |
| `pipeline/utils/coxiter.py` | CoxIter adapter |
| `apply_facet_profile_filter.py` | the 387 → 54 facet filter |
| `run_survivors_rigorous.py` | driver: prefix partitioning, exhaustion bookkeeping, verdicts |
| `validate_d5_wildcard.py` | d=5 revalidation through the identical code path |
| `verify_polytope.py` | standalone one-command verification of the final Gram matrix |


### 6.1 Where to find each part of the accompanying paper

The paper deliberately names no source files, so that its argument does not depend
on this repository. This table is the map from its claims to the code that produced
them. Section numbers refer to the paper in the companion repository
(`10-facets-6-d-polytopes-paper`).

| paper | claim | implementation |
|---|---|---|
| §3.1–3.3, Prop. 3.4 | 387 candidate combinatorial types from the order-type database via affine Gale duality | `pipeline/stage1_order_types.py` (reads `data/aak/otypes10.b16`), `pipeline/stage2_gale.py`; the exact face criterion is `pipeline/utils/gale_exact.py` |
| §3.4 | the generator reproduces 30 types at d=4 and 109 at d=5 exactly | `python -m pipeline.validate_coverage --d 4` and `--d 5` |
| §4, Lem. 4.3 | facet-admissibility filter, 387 → 54 survivors | `apply_facet_profile_filter.py` |
| §4, Lem. 4.6 / Cor. 4.7 | Esselmann's 2d bound kills the two 3-free types | `ESSELMANN_3FREE_KILLED` in `run_survivors_rigorous.py` |
| §4, Lem. 4.8 | Burcroff 5.5(b) low-weight edge caps | `burcroff_55b_low_weight_edges` in `pipeline/stage4_gram.py` |
| §4.1, Lem. 4.10 | wildcard predicate is insensitive to m ≥ 7 | `paper/checks/wildcard_soundness.py` |
| §4.1 | Lannér's order-4 and order-5 diagrams re-derived (9 and 5) | `paper/checks/lanner_label_bound.py` |
| §5.1 | enumeration with forward checking | `enumerate_labels_backtrack` in `pipeline/stage4_gram.py` |
| §5.2 | orbit symmetry breaking | `pipeline/utils/automorphisms.py` |
| §5.3 | the cascade screen; the tangency tolerance | `_structured_screen`, `_quad_roots_gt1`, `_DISC_RTOL` in `pipeline/stage4_gram.py` |
| §5.3, Lem. 5.4 / Prop. 5.5 | the bounded-box fallback is structurally unreachable at d=6; the pair path is never entered | `paper/checks/cascade_reachability.py` |
| §5.4 | wildcard range analysis and integer window scan | `_solve_wild_assignment` in `pipeline/stage4_gram.py` |
| §5.3, Prop. 5.3 | exact certificates for all 545 rejected non-wildcard labellings | driver `paper/checks/plain_exact_certify.py` |
| §5.4, Prop. 5.4 | exact infeasibility certificates for all 406 wildcard labellings, and the d=5 calibration | `pipeline/utils/exact_field.py`, `pipeline/utils/exact_certify.py`, driver `paper/checks/wild_exact_certify.py`, instances from `run_wild_dump.py`, tests `tests/test_exact_certify.py` |
| §5.5 | exact certification of an accepted labelling | `_refine_mpmath`, `_recognize_minpoly_and_verify`, `_inertia_mpmath` in `pipeline/stage4_gram.py` |
| §5.6 | exhaustion certificates and the coverage recursion | `process_type_stage4` (`exhausted`, `wild_unbounded`) and the driver `run_survivors_rigorous.py`; the recursion is recomputed by `paper/make_tables.py` in the paper repository |
| §6.1, eq. (1) | the Gram matrix, its exact signature and the absence of parabolic subdiagrams | `verify_polytope.py` (paper repository) |
| §7.1 | d=5 census reproduced at 51 | `validate_d5_wildcard.py 7 out=<dir>` |
| §7.2 | d=4 census reproduced at 348 | `run_survivors_rigorous.py all 7 wildcard d=4 out=<dir>`; the type-by-type comparison against Burcroff Appendix A is `paper/checks/d4_census_reconciliation.py` |
| §8(iv) | forward-checking tolerance margins | `paper/checks/forward_check_margins.py` |
| App. B | hardware, timings and the BLAS-pinning measurement | `probe_cost.py`; pin with `OMP_NUM_THREADS=1` and equivalents |
| App. C | branch counters and stall taxonomy | `probe_cost.py`, `paper/checks/screen_margins.py` |

Three checks import this pipeline and therefore live here rather than in the paper
repository: `paper/checks/cascade_reachability.py`,
`paper/checks/screen_margins.py`, `paper/checks/mz_p8_17_crosscheck.py`. The other
six are standalone and ship with the paper.

Material deliberately kept out of the paper — provenance arguments, the history of
the three defects corrected during development, and answers to questions a reader
may reasonably ask — is in `AUDIT_RESPONSES.md` in the paper repository.

---

## 7. Known limitations

Stated fully in §8 of the paper and in `REVIEW_NOTES.md`. The two that matter
most:

1. **The d=4 census is reproduced only to 338 of 348**, and via a different
   code path (fixed alphabet `{2,…,10,12}`) from the d=6 run. The d=5 census
   *is* reproduced exactly (51/51) through the identical d=6 code path.
2. **Emptiness verdicts are not exact certificates.** Acceptances are exact
   (100-digit certification + CoxIter), but a labelling is *rejected* by a
   screen evaluated in double precision, whose fallback branch searches
   ultraparallel weights only in the bounded box `[1.001, 1000]`.

## 8. Provenance of this code

Most of the code in this repository — the pipeline, the screening and
certification routines, the drivers, and the verification scripts — was written
by **Claude, a large language model developed by Anthropic**, working
interactively under the repository owner's direction. The owner set the
objectives, chose the mathematical approach, reviewed the code and output, and is
responsible for the results.

This is disclosed because it should affect how you check the work, not whether
you believe it. The design goal throughout has been to make the code's
correctness *checkable without trusting the code*:

* every reduction of the search space is a theorem, proved in the paper or cited
  by number — see `paper/paper.tex` §4;
* every accepted polytope is certified in exact arithmetic and re-checked by
  **CoxIter**, an independent third-party program;
* the pipeline reproduces the published d=5 census (51/51) through the identical
  code path used for d=6 — see §3.5 above;
* the final answer can be verified in one command that does not use the search at
  all — `python3 verify_polytope.py`.

A real defect *was* found in this code (see `REVIEW_NOTES.md` §3b-quater): the
screen's pair-resultant branch could refute realizable candidates depending on
facet numbering, which cost 11 of the 12 published 4-cubes. It was caught by
comparison against the independently published d=4 census, and it is confined to
a branch that dimensions 5 and 6 provably never reach. Treat that as a reason to
lean on the independent checks, and on `REVIEW_NOTES.md`, rather than on the
code's authorship.

## 9. Licence and citation

Please cite the paper and this repository. Third-party components retain their
own licences: CoxIter (`scratchpad/CoxIter`) is by Rafael Guglielmetti; the
Ma–Zheng data under `scratchpad/HCPdm` is by Jiming Ma and Fangting Zheng; the
order-type database is by Oswin Aichholzer et al.
