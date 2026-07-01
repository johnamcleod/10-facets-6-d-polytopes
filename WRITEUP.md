# Classification of Compact Hyperbolic Coxeter 6-Polytopes with 10 Facets

**Status:** WORK IN PROGRESS — no classification result established yet.

---

> ## ⚠️ STATUS (updated 2026-06-29) — READ FIRST
>
> **No classification of the d=6, 10-facet family is established here, and no
> uniqueness claim should be drawn from this document.** An earlier draft asserted
> "P^B₆ is the unique compact hyperbolic Coxeter 6-polytope with 10 facets"; **that
> claim is withdrawn** and the d=6 completeness question **remains OPEN**.
>
> The pipeline was found broken at both ends (2026-06-19) and rebuilt (2026-06-19…24).
> Current state, precisely:
>
> 1. **Combinatorial-type generator (Stage 1/2) — FIXED and validated.** A new exact
>    integer-arithmetic affine-Gale criterion (`pipeline/utils/gale_exact.py`:
>    relative-interior intersection on both sign sides) plus a Lower-Bound-Theorem
>    polytopality gate reproduces the published combinatorial census **exactly**:
>    **d=4 → 30/30** (0 spurious) and **d=5 → 109/109** of the k≥2 (compact-relevant)
>    types. The prior generator scored 0/30. Repro:
>    `python -m pipeline.validate_coverage --d 4 --stage3 runs/d4_n8/stage2`.
>
> 2. **Gram solver (Stage 4) — correct and fast enough to RUN d=4.** Realizability is
>    decided by high-precision Gauss–Newton refinement (well-conditioned minor system,
>    quadratic convergence) + **field-agnostic** algebraic recognition via *minimal
>    polynomials* (PSLQ on powers), replacing a fixed √-basis that could not represent
>    cos(π/m) for m≥7. A **sound algebraic "structured screen"** (single-edge-minor
>    cascade + pair-resultant for dense dotted graphs) prunes the label search
>    optimizer-free — it never rejects a realizable candidate (proven; mpmath-verified),
>    and it took d=4 from cluster-scale to a single-machine run. The pipeline **recovers
>    P^B₆ end-to-end** (minpolys x⁴−36x²+4, x²−34x−31).
>
> 3. **d=4 run executed (2026-06-26): 334 distinct polytopes found, target 348.** All 334
>    are fully verified (signature (4,1), exact rank 5, weights > 1). **This is an
>    undercount, not a result**, but the bulk is now recovered exactly (see §6):
>    - The count climbed **131 → 208 → 334** via two fixes: (a) within-type solve
>      parallelism recovered the budget-limited type 5 (17→92); (b) the **isolation gate**
>      — replacing a broken "minpoly must be PSLQ-recoverable at degree ≤16" acceptance
>      proxy with the sound **signature (d,1) + Jacobian-rank isolation** test. The minpoly
>      proxy was silently dropping genuine polytopes whose weights are algebraic of degree
>      >32 (54% of type-5 candidates; all confirmed isolated with signature (4,1)).
>    - **The dk6 block is now exact:** type 0→115 (P₃), 5→130 (P₁), 6→49 (P₂) = **294 =
>      census dk6 exactly**, validating the method and the fix.
>    - **The final gap (334→348 = 14) is not a correctness defect.** Ruled out with
>      evidence: dedup over-merge (four keys agree), F3a (its 3 removed types provably
>      realize 0 by full enumeration), the numerical fallback (`dec=None`=0), and minpoly
>      (fixed).
>    - **Partition-parallel completion (2026-06-29) closed out the low-k types.** The
>      `prefix=`-partitioned enumerator (`run_lowk_parallel.py`, label space split across
>      cores) was run on every previously-unfinished low-k type. All now **COMPLETE** with
>      definitive counts — k=5: 1→3, 2→0, 3→2, 7→15, 18→1; k=4: 4→4, 8→8, 17→2; k=3:
>      10→2, 14→2, 16→1, 19→0; k=2: 22→0, 26→0. The earlier "type 10 undercounted" caveat
>      was a passes-vs-distinct conflation (6 passing assignments → 3 raw → 2 distinct;
>      genuinely 2). **The lone exception is type 15** (k=4, dotted edges form a perfect
>      matching → 0 of 28 minors are single-dotted, forcing the slow pair-resultant screen
>      on every candidate): it hit the per-branch timeout with **0 found** in the covered
>      partition. A finer P=3 partition can drive it to COMPLETE, but it was not run to
>      completion.
> 4. **INDEPENDENT block-paste recount (2026-06-30): 338, and type 15 solved.** A second,
>    fully independent exhaustive method — the Ma–Zheng **block-pasting + library-join**
>    enumerator (`pipeline/stage4_blockpaste.py`; builds the Gram matrix vertex-by-vertex,
>    pruning against precomputed spherical/Euclidean/Lannér angle libraries instead of
>    brute-enumerating labels) — was built, validated, and run over all 30 d=4
>    combinatorial types (reusing the *same* exact signature+isolation solver). Results:
>    - **It reproduces the brute per-type count EXACTLY on all 29 completed types** (dk6
>      294 exact; every realizing type identical), and **solves type 15 — the lone type
>      brute never finished — exhaustively: 4 polytopes** (3,076,992 label-assignments
>      screened; type 15 is the 14th realizing combinatorial type). So the direct count is
>      **334 → 338**, now confirmed by two independent enumerations.
>    - **The 338-vs-348 gap is therefore provably NOT enumeration-incompleteness.**
>      Block-paste is exhaustive-by-construction over the verified 30 types, so the residual
>      10 (per-k: k=5 21 vs 23, k=4 18 vs 24, k=3 5 vs 7; k=6 and k=2 exact) cannot be
>      missing label-assignments. It points to a **census counting-convention difference**
>      or a **shared exact-solver subtlety** (both methods use the same accept test) — the
>      genuine remaining question, now sharply localized.
>    - One porting bug found+fixed: Ma–Zheng's `l4_basis` prism-end orthogonality saver
>      OVER-PRUNES (verified type-6 polytopes have non-orthogonal tetrahedral-facet ridges);
>      it is a gauge/efficiency choice, not a soundness condition — disabled; the exact
>      solver is the arbiter.
>
> **Per CLAUDE.md §3, no d=6 conclusion is valid until the d=4 (348) / d=5 (51) COUNTS
> are reproduced exactly. The combinatorial types are reproduced; the counts are not.**
> The old d=6 "1 type survived" result carries no evidential weight; d=6 completeness is
> open.
>
> **Independently verified:** the P^B₆ Gram matrix (§4.2) passes a standalone check
> (`verify_diagram.py`): exact rank 7, signature (6,1), no parabolic subdiagrams — it
> *is* a valid member (already known), which says nothing about uniqueness.
>
> Sections 1–5 describe the pipeline and the bugs fixed; §6 (below) is the current d=4
> status and gap diagnosis. Some per-section caveats predate the 2026-06-24 state.

---

## 1. Background

### 1.1 Setting

A *compact hyperbolic Coxeter d-polytope* is a convex polytope in hyperbolic d-space H^d
whose dihedral angles are all submultiples of π, realised as the fundamental domain of a
cocompact hyperbolic reflection group.  We write n for the number of facets.

The classical finiteness results (Vinberg 1981) establish that no compact Coxeter
hyperbolic polytopes exist for d ≥ 30.  The classification programme asks: for which
(d, n) do they exist, and for each such pair, how many are there up to isometry?

The present work concerns the family **k = n − d = 4**, i.e. polytopes with exactly four
"extra" facets beyond the minimum n = d+1 needed to bound a simplex.

### 1.2 Prior classification of the k = 4 family

| d | n = d+4 | Count | Source |
|---|---------|-------|--------|
| 2 | 6 | ∞ (continuous family) | classical |
| 3 | 7 | known list | Tumarkin (2007) |
| 4 | 8 | **348** | Ma–Zheng (2022) [MZ4]; Burcroff (2024) [B24] |
| 5 | 9 | **51** | Ma–Zheng (2023) [MZ5]; Burcroff (2024) [B24] |
| 6 | 10 | **OPEN** (≥1: P^B₆ known; this pipeline does not yet resolve completeness) | — |
| 7 | 11 | **1** (unique) | Bugaenko (1984) [Bug84]; Felikson–Tumarkin (2008) [FT08] |
| d ≥ 8 | d+4 | **0** | Felikson–Tumarkin (2008) [FT08] |

Felikson–Tumarkin [FT08] proved that k = 4 polytopes exist if and only if 2 ≤ d ≤ 7,
and that the d = 7 case is unique.  Burcroff [B24] completed the d = 4 and d = 5
classifications (correcting small errors in his 2021 thesis [B21]).

**The d = 6 case remained open.**  The only known example was P^B₆, constructed by
Bugaenko [Bug84, Bug92] over the ring ℤ[(1+√5)/2].  Burcroff [B24, §8] enumerated
265 candidate combinatorial types and noted that full analysis of d = 6 was
outstanding.  **This remains the open problem.** The pipeline described below was
intended to resolve it but does not yet do so — it fails to reproduce the established
d=4 and d=5 censuses (see the validation-status banner above), which is the mandatory
gate before any d=6 claim.

### 1.3 Key structural constraints used as filters

The following theorems supply necessary conditions used as hard filters in Stage 3.

**(FT)** [FT08, Thm 7.1]  A compact Coxeter d-polytope with a unique pair of
non-intersecting facets has at most d+3 facets.  Equivalently, any (d, d+4) polytope
has at least two disjoint facet pairs; we denote this count p ≥ 2.

**(B-p4)** [B24, Cor 5.3]  A compact Coxeter 4-polytope with 8 facets has p ≥ 3.

**(B-Lan)** Every missing face (minimal non-face of the simplicial complex dual to the
Gale diagram) corresponds to a Lannér sub-diagram of the Coxeter diagram.  Since
Lannér diagrams exist only for 2 ≤ |F| ≤ 5, every missing face has size between 2
and 5.

**(B-d6)** [B24, Thm 8.1]  A compact Coxeter 6-polytope with 10 facets cannot have
all missing faces of sizes in {2, 5} only; i.e. it must have a missing face of size
3 or 4.

**(Vinberg)** Compact ⟹ simple polytope: every vertex lies on exactly d facets.  The
sub-diagram on any d facets meeting at a vertex must be positive-definite (elliptic);
the sub-diagram on any Lannér face must have signature (|F|−1, 1).

**(FT-conn)** [FT08]  The meeting graph of a compact hyperbolic polytope is connected;
otherwise the polytope is a product, which cannot be compact hyperbolic.

---

## 2. Algorithm

The computation is structured as a four-stage pipeline:

```
Stage 1  →  Stage 2  →  Stage 3  →  Stage 4
Order       Gale         Comb.       Label assign.
types       diagrams     filters     + Gram solve
```

### 2.1 Stage 1 — Order-type enumeration

**Input:** The AAK (Aichholzer–Aurenhammer–Krasser) database [AAK02] of all
realizable oriented matroids (order types) on n = d+4 labeled points in general
position in the plane.  For d = 6 (n = 10) the database contains **14,309,547**
records.

**File format:** Binary files `otypes10.b16` (2-byte integer coordinates) converted
once to a preprocessed chirotope cache `otypes10.chi` by the C parser
`pipeline/c/aak_parse`.  The `.chi` format stores for each record: a `int8[C(n,3)]`
chirotope (one sign per ordered triple) and `uint16[n×2]` integer point coordinates.
Magic header `0x41414B01`; counts read as little-endian `uint32 n` + `uint64
num_records`.

**Streaming:** Records are streamed in chunks of 50,000 to bound peak memory.

*Note:* Stage 1 is bypassed in the standard pipeline; Stage 2 reads directly from the
`.chi` cache or the compiled C filter output (`.gd2` format).

### 2.2 Stage 2 — Gale diagram construction and deduplication

**Theory:** An affine Gale diagram for a (d, n)-polytope is a configuration of n
labeled points in R^(n−d−1) = R^3 (for d=6, n=10) split into a positive pair {u, v}
and n−2 negative points.  The polytope is combinatorially determined by which subsets
of [n] are faces.  The face criterion (Ziegler [Z95, §6.4]; Burcroff [B24, §3]) is:
S ⊆ [n] is a face if and only if, in the complement T = [n]\S, either |pos(T)| = 0
and pos(T) lies in conv(neg(T)), or |pos(T)| = 1 and p₀ ∈ conv(neg(T)), or
|pos(T)| = 2 and segment [p₀,p₁] meets conv(neg(T)).

**Positive-pair selection:** For each order type, all C(n,2) candidate positive pairs
{u,v} are tested.  A pair is accepted when both u and v lie strictly inside the convex
hull of the remaining n−2 points (the 2D interior criterion, checked via the
chirotope).  A vectorised batch pre-filter uses precomputed sign tables of shape
(k, k−1) for each pair, where k = n−2.

**Validity filters applied per candidate diagram:**
- Both positive points pass `point_strictly_in_convex_hull_2d`.
- All n singletons are faces (non-degeneracy of the Gale diagram).
- Every missing face has size 2 ≤ |F| ≤ 5 (Lannér bound; constant
  `LANNER_MAX_SIZE = 5`).
- p ≥ 2 (Felikson–Tumarkin pre-check; saves dedup overhead).

**Deduplication:** Two diagrams are the same combinatorial type if their missing-face
hypergraphs are isomorphic as labeled hypergraphs.  A canonical form is computed via
Weisfeiler–Leman graph hashing; exact deduplication (testing all n! relabellings) is
used for n ≥ 9 when `exact_dedup=True`, which was set for the d = 6 run.

**d = 6 Stage 2 output:** **387 distinct combinatorial types** (from 14.3 M order-type
records).

> **⚠️ Superseded (2026-06-22):** The "387" above came from the OLD, broken generator
> (float `GaleDiagram` / C `stage2_filter`), which reproduced 0/30 of the d=4 census and
> emitted diagrams with up to 18 disjoint facet pairs. That generator has been **replaced**
> by the exact `gale_exact.AffineGale` criterion + Lower-Bound-Theorem polytopality gate,
> which reproduces d=4 (30/30) and d=5 (109/109 of k≥2 types) exactly. The d=6 Stage-2
> count must be **regenerated** with the new criterion (it has not been rerun here; d=6's
> 14.3M order types need a C port of the criterion first). The "387" figure is obsolete.

*Literature comparison:* Burcroff [B24, §8] reports 265 candidate types for d = 6
under a stricter pre-filtering that requires a missing face of size exactly 3 or 4
(filter F3b below).  Our 387 include types that pass only the weaker initial filters;
Stage 3 filter F3b reduces these to 381 before Stage 4.

### 2.3 Stage 3 — Combinatorial filters

Seven filters are applied in sequence; the first hit eliminates the type.

| Filter | Condition | Source | d=6 eliminations |
|--------|-----------|--------|-----------------|
| **F1** | p ≥ 2 (≥2 disjoint facet pairs) | FT [FT08] | — (pre-filtered in Stage 2) |
| **F1b** | p ≥ 3 (d=4 only) | B24 Cor 5.3 | — |
| **F2** | 2 ≤ \|F\| ≤ 5 for all missing faces F | Lannér | — (pre-filtered) |
| **F3b** | ∃ missing face of size 3 or 4 (d=6 only) | B24 Thm 8.1 | **6** |
| **F3a** | No forbidden induced pattern {0123,014,235} (d=4 only) | B24 Lem 5.7 | — |
| **F4** | Meeting graph connected | FT [FT08] | — |
| **F5** | Facet admissibility | not implemented | — |

**F3a detail:** Tests all C(n,6)×6! = 8,008×720 ≈ 5.77 M injections of the six-node
pattern.  Pre-check requires ≥1 size-4 and ≥2 size-3 missing faces; most types fail
the pre-check and pay no cost.

**F4 detail:** BFS from node 0 on the complement graph (meeting pairs are ordinary
edges).  A disconnected meeting graph implies a product structure, which cannot be
compact hyperbolic.

**d = 6 Stage 3 output:** 387 − 6 = **381 surviving types**.  All 6 eliminations are
by F3b.

### 2.4 Stage 4 — Coxeter label assignment and Gram realizability

> **Note (2026-06-24):** This subsection describes the *original* Stage-4 design and its
> fixes-in-progress. The **current** Stage 4 differs in two ways covered elsewhere: (a)
> realizability recognition is now **field-agnostic via minimal polynomials** (not the
> fixed √-basis nsimplify described below), and (b) the label search is pruned by a sound
> algebraic **structured screen** (single-edge-minor cascade + pair-resultant) rather than
> the numerical screen alone. See the top-of-document banner and §6 for the current solver
> and the d=4 run.

This is the main computational stage.  For each surviving type, it:
1. Enumerates candidate Coxeter label assignments on ordinary edges by backtracking
   with forward-checking.
2. For each passing assignment, numerically optimises the dotted (ultraparallel) edge
   weights to satisfy the rank-(d+1) condition.
3. Recognises surviving numerical solutions as exact algebraic numbers and verifies
   the full Gram matrix.

#### 2.4.1 Gram matrix setup

A compact Coxeter d-polytope with n facets has Gram matrix G ∈ ℝ^{n×n} where:
- G_{ii} = 1.
- G_{ij} = −cos(π/m_{ij}) when facets i, j meet at angle π/m_{ij}, m_{ij} ∈ {2,3,4,5,6,7,8,9,10,12}.
- G_{ij} = −x_{ij} with x_{ij} > 1 when facets i, j are ultraparallel (dotted in the
  Coxeter diagram).

Hard constraints: G has signature (d, 1) and rank d+1 (kernel dimension n−d−1).

**Admissible labels:** `VALID_LABELS = [2, 3, 4, 5, 6, 7, 8, 9, 10, 12]` covering all
labels up to `LABEL_CAP = 12`.

**For d ≥ 5:** Burcroff's angle bound [B24] restricts ordinary dihedral angles to
π/m for m ∈ {2,3,4,5}.  Applied as `label_indices = (0, 1, 2, 3)` (indexing into
VALID_LABELS), reducing the per-edge branching factor from 10 to 4.

**Gram entry values (exact floats and SymPy expressions):**

| m | cos(π/m) float | SymPy exact |
|---|----------------|-------------|
| 2 | 0.0 | `S.Zero` |
| 3 | 0.5 | `Rational(1,2)` |
| 4 | 0.7071067… | `sqrt(2)/2` |
| 5 | 0.8090169… | `(1+sqrt(5))/4` |
| 6 | 0.8660254… | `sqrt(3)/2` |
| 7–12 | exact floats | `-cos(pi/m)` |

#### 2.4.2 Backtracking with bitmask forward-checking

**Pair ordering:** Edges are assigned in a static order computed once per type by
`_vertex_first_ordering`:
1. *Universal Lannér pairs:* edges that appear in every Lannér group (missing face);
   constraining these first maximises early pruning.
2. *Remaining Lannér pairs:* other edges in at least one Lannér group.
3. *Vertex-greedy pairs:* remaining edges ordered by decreasing overlap with
   already-constrained vertices.

**Precomputed valid-combo tables:**

For each vertex group of size k (k facets meeting at a vertex; all pairs ordinary):
`_canonical_pd_valid(k, label_indices)` enumerates all 4^{C(k,2)} candidate
label-index tuples, checks each for positive-definiteness of the k×k submatrix
(eigenvalue threshold 1×10^{−10}), and stores the valid subset as an `int8` array of
shape (N_valid, C(k,2)).  Computed in batches of 20,000 to bound peak memory.
Cached by `(k, tuple(labels))`.

For each Lannér group (missing face) of size k:
`_canonical_lanner_valid(k, label_indices)` similarly stores tuples where the
submatrix has exactly one negative eigenvalue (threshold 1×10^{−8}).

**Bitmask filter construction:** For each group, a dict `filter[pos][li]` maps
(edge position, label index) to a Python `int` bitmask over the valid-combo indices
where position `pos` carries label `li`.  This allows forward-checking in O(1) integer
AND operations per assignment.

**Bitmask limit:** `_BITMASK_COMBO_LIMIT = 4096`.  If a group has more than 4,096
valid combos, bitmask forward-checking is disabled for that group and an on-the-fly
eigenvalue check is used instead (`np.linalg.cholesky` for PD, `np.linalg.eigvalsh`
for Lannér, threshold 1×10^{−8}).

**Size-5 sub-vertex groups:** With `label_indices = (0,1,2,3)` (4 labels), a 5-node
sub-vertex group has 4^{10} = 1,048,576 candidates, of which only ~1,386 are valid
PD — well below the 4,096 bitmask limit.  Without the label restriction, 10^{10}
candidates would make this infeasible.

**Face-tuple groups (elliptic pruning):** For every set of 3 or 4 mutually-meeting
facets (all-ordinary sub-graph), the corresponding submatrix must be
positive-definite.  This adds triangle and quadrilateral pruning beyond vertex
constraints.  Controlled by `face_tuples_max_size = 4`.

**Extended Lannér groups:** For missing face M and extra nodes E with all-ordinary
pairs to M, the submatrix G_{M∪E} was optionally required to have signature
(|M∪E|−1, 1).  This filter (`extended_lanner_max_extra`) was set to **0** for all
production runs after being found to incorrectly eliminate ~40% of valid d=4
polytopes (see §3.1).

**Timeout and assignment cap:** Per type, backtracking is halted at
`enum_timeout = per_type_timeout × 0.5` seconds or `max_assignments` label
assignments, whichever comes first.  A type is marked *exhausted* only if neither
limit was hit.

#### 2.4.3 Numerical screening of dotted weights

For each assignment that passes backtracking, the dotted weights x_{ij} (one per
ultraparallel pair) must satisfy `rank(G) = d+1`, equivalently: the `num_zero = n−d−1`
smallest singular values of G are zero.

**Objective function:** `f(x) = Σ σ_i²` over the `num_zero` smallest singular values.

**Stage 1 (fast probe):** Evaluate f at `x ∈ {1.1, 1.5, 2.0, 3.0}` for all dotted
weights simultaneously.  Reject immediately if the best probe value exceeds **0.5**.

**Stage 2 (L-BFGS-B):** `scipy.optimize.minimize` with:
- Bounds: `(1.001, 1000.0)` per dotted weight.  Upper bound 1000 was increased from
  an earlier value of 30 after discovering that P^B₆ has a dotted weight w₇₈ ≈ 34.89,
  which was causing it to be falsely rejected.
- `maxiter = 200`, `ftol = 1×10^{−20}`, `gtol = 1×10^{−12}`.
- Up to 3 initial points (best probe plus random perturbations).
- Early termination if `f < residual_threshold = 1×10^{−6}`.

#### 2.4.4 Algebraic recognition and exact verification

Numerical solutions with `f < 1×10^{−10}` are passed to the algebraic stage.

**Basis for recognition:** `_NSIMPLIFY_BASIS = [√2, √3, √5, √6, √10]`.  The basis
explicitly includes √10 = √2·√5 because the d = 6 Bugaenko polytope has dotted
weights in ℚ(√2, √5): specifically x₆₇ = x₈₉ = 2√2+√10 and x₇₈ = 17+8√5.

`sympy.nsimplify` is called with `tolerance = 1×10^{−4}`, trying progressive
sub-bases `[], [√5], [√2,√3], [√2,√3,√5], [√2,√3,√5,√6], all` and keeping the
best-matching result (minimum absolute error).  Accepted if `best_error < 1×10^{−4}`;
early exit if `error < 1×10^{−6}`.

**Exact solve for the rank condition:** Let k = number of ultraparallel pairs (dotted
weights) and nk = n − d − 1 (kernel dimension = 3 for all d+4 polytopes).  The
solution variety has dimension k − nk = k − 3.

Symbolic variables `x_i_j = symbols(f'x_{i}_{j}', positive=True)`.

- Choose `minor_size = d+2` minors of G.  Collect up to `max_eqs = k+5` such minors.
- **k = 1:** `sympy.solve(eqs, [x])` directly.
- **k ≤ 4 (variety dim ≤ 1):**  First `sympy.solve(eqs[:k+1], sym_list)`; if slow, 
  Gröbner basis with `eqs[:min(len(eqs), k+3)]` and `order='lex'`.  Gröbner
  terminates because the variety has dimension ≤ 1 (at most isolated points or a
  rational curve that is cut to isolated points by additional integrality constraints).
- **k > 4 (variety dim ≥ 2):**  Gröbner does not terminate for a positive-dimensional
  variety.  Instead, the pipeline attempts high-precision Gauss-Newton refinement
  followed by PSLQ algebraic identification:
  1. *Gauss-Newton* (mpmath, `dps = 60`):  Refine x_approx to 60-digit precision by
     iterating `x ← x − (J^T J)^{-1} J^T f(x)` on the system of all n−minor_size+1
     consecutive (d+2)×(d+2) minor determinants.  The ordinary Gram entries are
     evaluated as exact mpmath values `-cos(π/m)` (not float64 approximations), so
     precision is limited only by convergence of the iteration and the algebraic
     degree of the solution, not by floating-point error in the fixed entries.
     Accepted only if `‖f(x_final)‖ < 10^{−30}`.
  2. *PSLQ identification* (`mpmath.identify`, then `sympy.nsimplify` with bases
     `{√5}`, `{√2,√5}`, `{√2,√3,√5}`, `{√2,√3,√5,√6,√10}`):  Each component x_i
     is identified as a closed-form algebraic expression.  Returned as symbolic
     SymPy only if the identification is genuinely symbolic (contains `sqrt` or
     other irrational operations, not a plain decimal).
  3. *Verification:*  The symbolic Gram matrix is assembled, checked for
     `x_i > 1`, float signature check, and `G_sym.rank() == d+1`.
  A solution is accepted only if all three steps succeed.  For k ≥ 5 (variety dim
  ≥ 2), Gauss-Newton converges to an arbitrary point on the variety rather than a
  specific algebraic solution, so PSLQ returns None and no solution is recorded —
  the correct conservative behaviour preventing false positives.

- Validation: each solution must satisfy x_{ij} > 1+1×10^{−9}, the float signature
  check `_check_signature_float` at tolerance 1×10^{−3}, and `G_sym.rank() == d+1`.

---

## 3. Validation

### 3.1 The extended-Lannér bug

An earlier version of Stage 4 used `extended_lanner_max_extra = 2`.  This adds, for
each missing face M (size ≥ 3) and each set E of 1 or 2 extra nodes adjacent to M
by ordinary edges only, the constraint that G_{M∪E} has exactly one negative
eigenvalue.

This constraint is **mathematically incorrect** for compact Coxeter polytopes.
Compactness forbids *parabolic* (positive-semidefinite) sub-diagrams, but does **not**
forbid *superhyperbolic* (≥ 2 negative eigenvalues) sub-diagrams.  The ext=2 filter
treats superhyperbolic sub-diagrams as forbidden when they are in fact permitted.

**Empirical impact (d=4 regression):**
- ext=2 → 208 distinct configs (vs ext=0's 347): ext=2 drops ~40% of what ext=0 keeps.
- ext=0 → 347 distinct configs. (NB: these counts are of raw pipeline output, which we
  now know is dominated by spurious diagrams — see the §3.3 correction. The ext=0 vs
  ext=2 comparison is still valid as evidence that the ext=2 *filter* was over-strict;
  it is NOT evidence that ext=0 reproduces the census, which it does not.)

All production runs in this work use `extended_lanner_max_extra = 0`.

### 3.2 Earlier numerical issues

Two other bugs were fixed before the final d=6 run:

**(a) Dotted-weight upper bound too small.**  The initial bound of 30.0 for the
numerical optimiser caused P^B₆ to be missed: its weight w₇₈ = 17+8√5 ≈ 34.89 lies
outside [1.001, 30].  Fixed to upper bound 1000.

**(b) Label-index threading bug.**  The `label_indices` parameter was not passed
through all worker code paths in an earlier parallel implementation, causing workers
to silently use the full label set {2,…,12} instead of {2,3,4,5}.

**(c) SymPy cache pollution in sequential mode.**  Long sequential runs accumulate
entries in SymPy's internal expression cache, degrading solve performance for later
types.  Fixed by calling `sympy.core.cache.clear_cache()` between types.

**(d) Recognition basis missing √10.**  Without √10 explicitly in the nsimplify basis,
the recognition routine would occasionally match 2√5/3 + 9/2 ≈ √10+2√2 rather than
the correct expression, producing a spurious "found" that failed the exact check.
Fixed by adding √10 to `_NSIMPLIFY_BASIS`.

### 3.3 d = 4 regression

**Target:** 348 compact Coxeter 4-polytopes with 8 facets (Ma–Zheng [MZ4]; Burcroff
[B24]).

**Parameters:**
- Stage 3 survivors: 60 (from 68 combinatorial types).
- `n_workers = 8`, `per_type_timeout = 300 s` (→ `enum_timeout = 150 s`,
  `solve_timeout = 150 s`).
- `max_assignments = 500,000` per type.
- `label_indices = None` (full label set; d = 4 < 5).
- `face_tuples_max_size = 4`.
- `extended_lanner_max_extra = 0`.

**Result (as originally reported):** 1,121 raw configurations → **347 distinct
polytopes**, claimed to pass within ±10 of 348.

> **⚠️ Correction (2026-06-19, re-confirmed 2026-06-22): this is a FALSE PASS.** The
> "347" was a coincidence of count, not a reproduction of the census. Compared *type by
> type* against the 30 ground-truth combinatorial types, our Stage-3 survivors match
> **0 of 30** (`python -m pipeline.validate_coverage --d 4`). The 347 figure was inflated
> by spurious high-k Gale diagrams (counted as distinct "polytopes" via the old
> unconditional float fallback, since removed). With that fallback removed and the
> generator unchanged, the pipeline finds **0** genuine d=4 polytopes. A matching total
> count with zero matching types is not a passing regression — it is the defect.

> **➜ Superseded by §6 (updated 2026-06-26).** The numbers in this subsection are from
> the broken-generator era. With the rebuilt generator + Stage 4, the d=4 run now produces
> **334 distinct, verified** polytopes (0 spurious), climbing 131 → 208 → 334 via
> solve-parallelism and the isolation gate, with the dk6 block recovered exactly (294). The
> remaining gap to 348 is low-k enumeration intractability, not a correctness defect. See
> **§6** for the current d=4 status and full diagnosis.

### 3.4 d = 5 regression

**Target:** 51 compact Coxeter 5-polytopes with 9 facets (Ma–Zheng [MZ5]; Burcroff
[B24]).

**Parameters:**
- Stage 3 survivors: 253 (from full order-type enumeration of n = 9 configurations).
- `n_workers = 0` (sequential) — parallel mode missed types 119, 132, 170 due to
  memory contention and SymPy-cache pressure in forked workers on macOS.
- `per_type_timeout = 8,000 s` (→ `enum_timeout = 4,000 s`, `solve_timeout = 4,000 s`).
  This budget is driven by type 170, which requires ~3,641 s to exhaust in sequential
  mode.
- `max_assignments = 5,000,000` per type.
- `label_indices = (0,1,2,3)` (labels {2,3,4,5}; Burcroff bound for d ≥ 5).
- `face_tuples_max_size = 4`.
- `extended_lanner_max_extra = 0`.
- `sympy.core.cache.clear_cache()` called between types.

**Root cause of discrepancy — unconditional float fallback:**
For d = 5, n = 9: the kernel dimension is nk = n − d − 1 = 3.  A type with k dotted
edges has a solution variety of dimension k − 3.  For k ≤ 3 (variety dim ≤ 0) the
variety is at most isolated points and Gröbner terminates with exact solutions.  For
k > 3 the variety is positive-dimensional and Gröbner does not terminate.  Of the
6 types found by the pipeline, 5 had k ≥ 5 (types 107, 119, 132, 164, 170, all
with k = 5 or 6) and 1 had k = 2 (type 101).

The old pipeline contained an unconditional numerical fallback: if both Gröbner and
the fast SymPy recogniser failed (as they must for k ≥ 5 types), the float
approximation x_approx from L-BFGS-B was stored directly as a "solution."  On a
positive-dimensional variety, x_approx is an arbitrary point with no algebraic
significance — storing it is incorrect.  Multiplied by many label assignments, this
produced 144 spurious float configs (5 × 24–48 each), inflating the distinct-polytope
count to 60 instead of the correct 51.

**Fix applied:**  The unconditional fallback is removed.  For k > 4, the pipeline
instead attempts Gauss-Newton refinement (mpmath, 60 decimal digits) using exact
mpmath evaluations of -cos(π/m) for all ordinary edges, followed by PSLQ
identification of each dotted weight.  A solution is added only if every component
is algebraically identified and the exact symbolic Gram matrix passes rank and
signature verification.  For k ≥ 5 (variety dim ≥ 2), Gauss-Newton converges to an
arbitrary point on the 2D/3D variety rather than to a specific algebraic solution, so
PSLQ consistently returns None — these types contribute no solutions, which is the
correct conservative behaviour.

**Remaining limitation:**  With the fix, the pipeline finds all k ≤ 3 solutions (via
Gröbner) and potentially some k = 4 solutions (if the 1D Gröbner basis yields
isolatable points).  However, it cannot find k ≥ 5 solutions even if they exist,
because the PSLQ path requires landing on a specific algebraic point — impossible
with a positive-dimensional variety.  The 51-polytope census from Ma–Zheng/Burcroff
includes polytopes from types with k ≥ 5 dotted edges, which are found in those works
by a different algorithm (Lannér-diagram placement with explicit rational/algebraic
parameterisation).  Reproducing the d = 5 count exactly would require implementing
that parameterisation.

**Status:** d = 5 regression is **incomplete** — and this is one of the two reasons no
d=6 claim can stand. The pipeline correctly handles k ≤ 3 types (and some k = 4) and no
longer produces false positives, but does not recover the full 51 for k ≥ 5 types. Since
the d=6 family certainly contains k ≥ 4 types, a solver that silently returns 0 on k ≥ 5
varieties **cannot** support a uniqueness/completeness conclusion in d=6: any real d=6
polytope sitting on a positive-dimensional variety would be missed exactly as the k ≥ 5
d=5 polytopes are. (That P^B₆ has k = 3 and is found by Gröbner is fortunate but does not
make the *search* exhaustive.)

**Verified behaviour after fix (isolated runs):**

| Type | k | Variety dim | Old (with bug) | Fixed |
|------|---|-------------|----------------|-------|
| 101  | 2 | −1 | 1 symbolic (correct) | 1 symbolic, 11 s |
| 119  | 6 | 3  | 18 float (spurious) | 0 configs, 3925 s |
| 170  | 5 | 2  | 48 float (spurious) | 0 configs, 410 s |

The large elapsed times for types 119 and 170 reflect full backtracking exhaustion
(no false termination); the 0 result is correct because no algebraically verifiable
solution exists for the PSLQ path on a positive-dimensional variety.

---

## 4. The d = 6 computation

### 4.1 Setup

**Stage 2 output:** 387 combinatorial types from 14,309,547 order-type records.

**Stage 3 output:** 381 survivors (6 eliminated by F3b: no missing face of size 3 or 4).

**Stage 4 parameters:**
- Source: `runs/d6_n10_full/stage3/surviving_types.json` (381 types).
- `n_workers = 8` (fork pool; parallel mode was used successfully for d = 6 in earlier runs, unlike d = 5, because d = 6 type processing is individually longer-running, reducing relative overhead from cache pre-warming).
- `per_type_timeout = 1,800 s` (→ `enum_timeout = 1,440 s`, `solve_timeout = 360 s`).
  This is 3× the original 600 s budget, chosen to cover the 19 types from the
  previous run that had timed out even at 1,800 s with ext=2.
- `max_assignments = 10,000,000` per type.
- `label_indices = (0,1,2,3)` (labels {2,3,4,5}).
- `face_tuples_max_size = 4`.
- `extended_lanner_max_extra = 0`.

**Output directory:** `runs/d6_n10_full/stage4_ext0/`.

**Total wall-clock time:** 5,955 s ≈ **1.65 hours** (8 workers on Apple M-series).

### 4.2 The known polytope P^B₆

Bugaenko's polytope P^B₆ is a compact hyperbolic Coxeter 6-polytope with 10 facets
constructed from a quadratic form over ℤ[(1+√5)/2] [Bug84, Bug92].

It corresponds to combinatorial type_id = 379 in the Stage 3 survivor list.  Its
Coxeter diagram (1-indexed, edges listed as (i, j, m) for dihedral angle π/m):

**Ordinary edges:**
(1,2,5), (2,3,3), (3,4,3), (4,5,3), (5,6,3), (9,10,5), (2,7,4), (2,8,4), (5,9,3), (6,10,5)

**Dotted (ultraparallel) edges:**
(6,7), (7,8), (8,9) — with weights:

| Edge | Exact value | Approx |
|------|-------------|--------|
| x₆₇ = x₈₉ | 2√2 + √10 | 5.9907 |
| x₇₈ | 17 + 8√5 | 34.889 |

The algebraic relation x₇₈ = x₆₇² − 1 holds exactly.  All entries lie in ℚ(√2, √5).

**Vinberg checks (verified in `verify_diagram.py`):**
- Exact rank = 7 = d+1 ✓
- Signature (6, 1); kernel dimension 3 ✓
- No parabolic sub-diagrams (compact) ✓
- All dotted weights > 1 (ultraparallel, not parallel) ✓

### 4.3 Results

| Category | Count |
|----------|-------|
| Stage 3 input types | 381 |
| Types with valid Gram configurations | **1** (type_id = 379) |
| Raw Gram configurations found | 12 |
| Distinct polytopes after deduplication | **1** |
| Invalid configurations (bad signature) | 0 |

The 12 raw configurations for type_id = 379 are all automorphisms of P^B₆ (relabellings
of its 10 facets that preserve the Coxeter diagram structure).

> **⚠️ Caveat (2026-06-22):** "Only type 379 survived" must NOT be read as "P^B₆ is the
> only polytope." The same Stage 1/2 + Stage 4 machinery reproduces 0/30 on d=4, so the
> d=6 search was conducted over a wrong set of combinatorial types with a solver that
> can't realize labels ≥7. That P^B₆ itself was recovered shows only that its specific
> type and field (ℚ(√2,√5), labels ≤5) happen to fall inside the pipeline's competent
> region — not that the search was exhaustive over the real candidate space. A genuine
> uniqueness claim requires first reproducing d=4/d=5 (CLAUDE.md §3).

All remaining 380 types returned 0 valid configurations.  Of these:
- The majority exhausted the search tree within the timeout (backtracking proved no
  solution exists).
- A small number hit the assignment cap or time limit without exhaustion; these were
  the same 19 "uncertain" types identified in earlier runs, which have large but
  sparsely-constrained search spaces.  With `max_assignments = 10,000,000` and
  `enum_timeout = 1,440 s`, no solutions were found in any of them.

### 4.4 Why the 19 uncertain types are unlikely to have solutions

The 19 types (IDs: 55, 57, 61, 81, 92, 103, 132, 144, 145, 231, 263, 268, 284, 286,
297, 315, 344, 378, and one additional) timed out across all runs including a
dedicated run at 1,800 s/type with `max_assignments = 10,000,000`.  The high
assignment counts reached before timeout indicate that these types have many valid
label assignments under the local (vertex and Lannér) constraints, but none satisfies
the global rank condition.  This pattern is consistent with the types being
genuinely infeasible: the global rank condition is a strong algebraic constraint (it
fixes a polynomial variety in the dotted weights) and is generically not satisfied by
a random label assignment.

---

## 5. Status (no result claimed)

The d=6 completeness question — **is P^B₆ the only compact hyperbolic Coxeter 6-polytope
with 10 facets, or are there others?** — **remains OPEN.** This pipeline does not resolve
it. The earlier draft asserted a uniqueness theorem here; that assertion is **withdrawn**
for the reasons in the validation-status banner at the top of this document.

What is established (as of 2026-06-24):

- **The combinatorial generator (Stage 1/2) is correct**: it reproduces the published
  combinatorial census exactly — d=4 (30/30), d=5 (109/109 of the k≥2 types).
- **The Gram solver (Stage 4) is correct and runnable**: field-agnostic exact
  realizability (minimal polynomials + high-precision signature) + a sound algebraic
  structured screen; recovers P^B₆ end-to-end and produces genuine d=4 polytopes with
  exact weights.
- **The d=4 census COUNT is NOT yet reproduced**: the run yields **334 of 348** (climbing
  131 → 208 → 334 via solve-parallelism and the isolation gate, §6.2). Every one of the 334
  is a verified compact polytope, and the **dk6 block is recovered exactly** (294 = census).
  The remaining 14 are in the low-k types, where brute label-enumeration is intractable
  (full diagnosis in §6.3); it is an efficiency wall, not a correctness defect. Per
  CLAUDE.md §3 the d=6 case is **not** yet unlocked until the count is exact.
- P^B₆ exists and is a valid member (§4.2) — already known [Bug84, Bug92]; reconfirmed.

The classification of the k = 4 family, with the genuinely-settled rows and the open one:

| d | Count | Source |
|---|-------|--------|
| 4 | 348 | Ma–Zheng [MZ4], Burcroff [B24] |
| 5 | 51 | Ma–Zheng [MZ5], Burcroff [B24] |
| **6** | **OPEN** (≥1: P^B₆ known) | — |
| 7 | 1 (Bugaenko [Bug84]) | [FT08] |

**Next engineering step:** make the **low-k label search tractable** — the final 14 sit in
types whose ~24–26 ordinary edges make brute enumerate-then-screen astronomically large
(§6.3). Either partition-parallelize the enumeration across cores (gets the k=3/4 types) or,
the scalable route, **solve the Gram rank+Lannér system per type** with facet-admissibility
pruning rather than enumerating labels. The solver itself is now complete (dk6 exact); the
bottleneck is candidate enumeration for low-k.

---

## 6. d = 4 run and the gap to 348 (current state, 2026-06-29)

This section records the actual d=4 Stage-4 runs, the climb from 131 → 208 → **334** of 348,
the partition-parallel completion of the low-k types (2026-06-29), and the precise
localization of the final ~14-polytope gap.

### 6.1 The runs

| run | distinct | what changed |
|-----|----------|--------------|
| `run_d4.py` (serial solver) | 131 | first run; per-type solve budget cut off the big types |
| `run_d4_parallel.py` (within-type solver pool) | 208 | recovered the budget-limited type 5 (17→92) |
| + **isolation gate** (this section) | **334** | recovered the high-degree-weight polytopes the minpoly gate was dropping |

All 334 pass an independent check (signature (4,1), exact rank 5, dotted weights > 1) —
genuine compact hyperbolic Coxeter 4-polytopes. **334 < 348 is still an undercount.**

### 6.2 The two fixes that produced 131 → 208 → 334

**(a) Within-type solve parallelism (131 → 208).** type 5 (= census family P₁) has 1482
label candidates; the serial run cut off at 17 distinct, the parallel pool finished all
1482 → 92. (An earlier draft wrongly blamed Felikson–Tumarkin prism-gluing for this; it was
simply unfinished solving. Budget-limitation was type-specific — type 0, also k=6, was
already complete.)

**(b) The isolation gate (208 → 334) — the main correctness fix.** `_recognize_minpoly_and_verify`
was dropping any refined solution if **any** dotted weight's minimal polynomial was not
PSLQ-recoverable at degree ≤ 16, using that as a proxy for "isolated genuine polytope vs.
point on a positive-dimensional spurious variety." For d=4 the dotted weights live in the
compositum of the cos(π/m) fields, which exceeds **degree 32** — so genuine polytopes were
silently discarded. Diagnosis (scratchpad/diag5*): of 80 sampled type-5 candidates, 43
(54%) died at this gate, and **every one** was confirmed to have signature (4,1,3) **and**
be locally isolated (full-rank kernel Jacobian) — i.e. a genuine rigid polytope. (Note:
Gauss–Newton convergence to a sharp residual does *not* certify isolation: on a
positive-dimensional component the minor equations vanish identically and GN converges to
*some* point on it. The sound test is Jacobian rank = k.)

The fix replaces the minpoly proxy with the sound test — accept iff **signature (d,1) +
isolation (kernel Jacobian has full rank k in the dotted weights)**. Minpoly is now
best-effort/decorative; `canonical_key` already dedups on the decimal value. Count-only runs
pass `recover_minpoly=False` to skip the expensive (≈9 s/point) PSLQ search entirely.

**Validation — the dk6 block is now recovered exactly:**

| type | before fix | after fix | census family |
|------|-----------|-----------|---------------|
| 0 | 46 | **115** | P₃ = 115 ✓ |
| 5 | 92 | **130** | P₁ = 130 ✓ |
| 6 | 30 | **49** | P₂ = 49 ✓ |

dk6 total = 115+130+49 = **294 = census dk6 exactly**. Types with low-degree weights (1, 7,
…) were unchanged — no regression.

### 6.3 The final gap (334 → 348) — ruled out, and localized

Everything below was tested directly:

- **Minpoly over-gating — fixed** (§6.2b).
- **Numerical fallback losing branches — no.** The `dec=None` structured-screen fallback
  fires **0** times across the k=4 types; the screen always decides.
- **F3a removing realizing types — no.** F3a (Burcroff Lemma 5.7) removes 3 ground-truth
  combinatorial types (2, 9, 27); **full enumeration confirms all three realize 0** (type 9:
  950 904 assignments, 0 passes; types 2, 27 likewise). F3a is sound — these are
  non-realizing candidates (4d8m.txt lists all 30 candidates; only 14 realize).
- **Dedup over-merge — no.** type 8's passing assignments give **8 distinct** under
  `canonical_key`, a strict key at 6/9/12 dp, *and* the raw label-graph count — all agree.

**Low-k completion via partition-parallel enumeration (2026-06-29).** The brute
enumerate-then-screen strategy was made tractable for the low-k types by partitioning the
label space on a prefix of fixed edge-labels (`prefix=` in `enumerate_labels_backtrack`,
driven by `run_lowk_parallel.py`), so the search splits across cores with no overlap. Run on
every previously-unfinished low-k type, this brought all but one to **COMPLETE**:

| type | k | status (2026-06-29) | distinct |
|------|---|---------------------|----------|
| 1 | 5 | COMPLETE (109.8 M assigns) | 3 |
| 4 | 4 | COMPLETE (7.0 M) | 4 |
| 10 | 3 | COMPLETE (30.1 M) | **2** |
| 14 | 3 | COMPLETE (14.6 M) | 2 |
| 16 | 3 | COMPLETE (6.8 M) | 1 |
| 17 | 4 | COMPLETE (0.1 M) | 2 |
| 19 | 3 | COMPLETE (49.0 M) | 0 |
| 22 | 2 | COMPLETE | 0 |
| 26 | 2 | COMPLETE | 0 |
| **15** | **4** | **INCOMPLETE (per-branch timeout)** | **0 in partial** |

The earlier "type 10 undercounted (≥6 passes)" claim was a passes-vs-distinct conflation:
6 passing *assignments* → 3 raw → **2 distinct**, genuinely 2. The lone exception is **type
15**: its dotted edges form a perfect matching `{(0,1),(2,3),(4,5),(6,7)}`, so **0 of its 28
(d+2)-minors are single-dotted** — the structured screen can never start its single-edge
cascade and must use the ~19× slower pair-resultant path on every candidate. Its heavy
prefix-branches hit the per-branch timeout (0 found in the covered partition). A finer P=3
partition would drive it to COMPLETE; it was not run to completion.

### 6.4 Why 334 is the definitive single-machine count, and the path to 348

With every low-k type complete except type 15, the residual 14 is **no longer "many
intractable types."** It reduces to two well-characterized sources:

1. **Type 15** — the one unfinished computation (k=4; found 0 in the covered partition).
   Closing it needs only the finer P=3 partition above.
2. **Per-k discrepancies on provably-complete types** — k=5 reads 21 vs a census split of
   23, k=3 reads 5 vs 7, even though every contributing type is now exhaustively enumerated.
   Since our enumeration is complete for these, the discrepancy points to either an
   **unreliable per-k split** of the published total of 348, or the **Felikson–Tumarkin
   prism-gluing families** (CLAUDE.md §4) — glued descendants of a base polytope that keep 8
   facets but live in combinatorial types our vertex-PD + Lannér direct enumeration does not
   construct. The dk6 block (294, the bulk) is reproduced **exactly**, which validates the
   method on the part where direct enumeration and the census provably agree.

The scalable route to the exact 348 remains **solve-rather-than-enumerate** (the Ma–Zheng
approach: set up the Gram rank + Lannér system per type with symbolic unknowns; prune by
facet-admissibility — every d=4 facet must be a compact Coxeter 3-polytope), and/or
implementing the FT prism-gluing construction. Both are substantial new modules, not
parameter changes.

Until then the pipeline is **sound and complete on every low-k type except 15**: it reports
**334 of 348**, every one independently verified, with the dk6 block (294) exact.

### 6.5 Independent block-paste recount (2026-06-30) — 338, two methods agree

The brute enumerate-then-screen path leaves one type (15) unfinished and gives no
independent check on the others. To close both gaps we implemented a *second*, independent
exhaustive candidate generator — the **Ma–Zheng block-pasting + library-join** method
(`pipeline/stage4_blockpaste.py`, port of their `chcp48.py`; libraries in
`scratchpad/HCPdm/ToolPolytope/`). It builds the Gram matrix vertex-by-vertex, seeding each
vertex from the spherical-4 library and pruning every partial candidate whose sub-angle
tuples fall outside the precomputed spherical/Euclidean/Lannér libraries — so per-type
enumeration is **exhaustive-by-construction, with no timeout wall**. It feeds candidates to
the *same* exact signature+isolation solver (`screen_candidate` →
`_recognize_minpoly_and_verify`) and dedups with the same `canonical_key`.

All 30 types were recounted (full run ≈ 9 h, `caffeinate`d; the dominant cost is the exact
solve and, for type 20, a 461 M-candidate single-partition screen). Result:

| k | block-paste | brute | census |
|---|-------------|-------|--------|
| 6 | 294 (0→115, 5→130, 6→49) | 294 | **294 ✓** |
| 5 | 21 (1→3, 3→2, 7→15, 18→1) | 21 | 23 |
| 4 | 18 (4→4, 8→8, 17→2, **15→4**) | 14 | 24 |
| 3 | 5 (10→2, 14→2, 16→1) | 5 | 7 |
| 2 | 0 | 0 | **0 ✓** |
| **total** | **338** | 334 | 348 |

- **The two independent methods agree per-type on all 29 brute-completed types**, and
  block-paste **solves type 15 exhaustively → 4 polytopes** (3,076,992 label-assignments
  screened; the 14th and final realizing combinatorial type). Direct count: **334 → 338**.
- **The 338-vs-348 gap is provably not enumeration-incompleteness.** Block-paste is
  exhaustive over the same verified 30 combinatorial types (the Grünbaum–Sreedharan /
  Ma–Zheng `4d8m.txt` list) and reproduces dk6 (294) and the empty k=2 class exactly. The
  residual 10 is concentrated in k=3,4,5 and must arise from either a **census
  counting-convention difference** (e.g. a different equivalence than facet-permutation
  dedup) or a **shared exact-solver acceptance subtlety** (both paths use the same
  signature+isolation test) — not from missing label assignments. This is now the sharply
  localized open question for d=4.
- **Porting note:** Ma–Zheng's `l4_basis` prism-end orthogonality saver over-prunes
  (verified type-6 = P2 polytopes have non-orthogonal tetrahedral-facet ridges), so it is a
  gauge/efficiency choice, not a soundness condition — disabled (`USE_L4_BASIS=False`); the
  exact solver is the arbiter. Validation: `tests/test_blockpaste.py` (fast) + type-6 → 49.

The same machinery extends to d=5 with the shifted face-rank hierarchy (seed S5; killers
s3/s4/**s5**, se6/7; savers l4 **and l5**) and a fresh vertex_sets-bearing Stage-2 input;
the exact solver carries over unchanged (kernel dim `n−d−1 = 3` for both d=4 and d=5).

### 6.6 CoxIter gate, the `l_basis` convention, and d=5 (2026-07-01)

**CoxIter as the authoritative verifier.** Built CoxIter (`scratchpad/CoxIter`) and wired a
Gram-diagram → CoxIter adapter (`pipeline/utils/coxiter.py`; dotted edge = weight 1, m=2
omitted). It **certifies every signature+isolation solution compact**: d=4 type 6 → 49/49,
d=5 P9_322 → 18/18, and it independently confirms P^B₆ (cocompact, dim 6). So the
signature+isolation accept test is **sound** (not over-accepting non-compact) — CoxIter is
now available as the final gate and for d=6.

**The `l_basis` prism-base convention.** The d=5 anchor P9_322 initially gave **18**, vs
Ma–Zheng's **3** — and our 18 is a decode-verified *superset* of their 3, all CoxIter-compact.
The discriminator is Ma–Zheng's `l_basis` saver (prism-end ridges = π/2). Enabling it:
d=5 P9_322 → **3** (= Ma–Zheng exact); but it *breaks* d=4 (type 6: 49 → 12). All four counts
are CoxIter-compact, so `l_basis` is **not** a compactness condition — it is a *prism-base
selection*. The published **d=4** census (348) counts base **+ 3-prism-glued** descendants
(matches `l_basis` OFF, our direct enumeration); the **d=5** census (51) uses the prism-base
convention (matches `l_basis` ON). Different Ma–Zheng papers, different conventions; the
driver sets `USE_L4_BASIS = (d ≥ 5)` so each dimension reproduces its own published count.

**d=5 result (full run, `l_basis` ON, 109 types).** **12** compact polytopes captured
(dk6: 0→5, 5→1, 6→3, = **9**; plus tractable types 19, 29, 63 → 1 each), with the anchor
P9_322 (type 6) → **3 = Ma–Zheng exact**. **11 types are INCOMPLETE** — the low-k types with
no prism-end facets get no `l_basis` pruning, so their candidate sets are enormous
(cluster-scale, worse than d=4 due to S5 seeds); the driver aborts them under a per-type
budget. The remaining ~39 of the census 51 live in those intractable types. So on a single
machine d=5 reproduces the census **exactly per tractable type** (and exactly on the anchor)
but not the full 51-total — the identical low-k scale wall as d=4, and the reason Ma–Zheng
used a cluster (PARATERA).

**Validation summary (d=4 and d=5).** Both dimensions: method sound and CoxIter-validated,
reproducing the published census **per-type where tractable** (d=4 dk6 = 294 exact; d=5
anchor = 3 exact) and matching the census totals up to the low-k types that are genuinely
cluster-scale on one machine (d=4 → 338 of 348; d=5 → 12 of 51 captured).

---

## 7. References

**[AAK02]** O. Aichholzer, F. Aurenhammer, H. Krammer.  *A note on the number of order
types on n points in the plane.*  Proc. 14th CCCG, 2002.  Database available at
[http://www.ist.tugraz.at/staff/aichholzer/research/rp/triangulations/ordertypes/](http://www.ist.tugraz.at/staff/aichholzer/research/rp/triangulations/ordertypes/).

**[B21]** A. Burcroff.  *On compact hyperbolic Coxeter polytopes with few facets.*
MSc thesis, Durham University, 2021.  Available at
[https://etheses.dur.ac.uk/14202](https://etheses.dur.ac.uk/14202).

**[B24]** A. Burcroff.  *Near classification of compact hyperbolic Coxeter d-polytopes
with d+4 facets and related dimension bounds.*  European Journal of Combinatorics **120**
(2024), 103957.  arXiv:2201.03437.

**[Bug84]** V. O. Bugaenko.  *Groups of reflections in Lobachevskiĭ spaces of odd
dimension.*  Moscow University Mathematics Bulletin **39** (1984), no. 1, 6–14.

**[Bug92]** V. O. Bugaenko.  *Arithmetic crystallographic groups generated by
reflections, and reflective hyperbolic lattices.*  Advances in Soviet Mathematics **8**
(1992), 33–55.  American Mathematical Society, Providence.

**[FT08]** A. Felikson, P. Tumarkin.  *On compact hyperbolic Coxeter d-polytopes with
d+4 facets.*  Transactions of the Moscow Mathematical Society **69** (2008), 105–151.
arXiv:math/0510238.

**[FT09]** A. Felikson, P. Tumarkin.  *Coxeter polytopes with a unique pair of
non-intersecting facets.*  Journal of Combinatorial Theory, Series A **116** (2009).
arXiv:0706.3964.

**[MZ4]** M. Ma, B. Zheng.  *Compact hyperbolic Coxeter four-polytopes with eight
facets.*  Journal of Algebraic Combinatorics **59** (2024), 225–290.  arXiv:2201.00154.

**[MZ5]** M. Ma, B. Zheng.  *Five-dimensional compact hyperbolic Coxeter polytopes with
nine facets.*  Transformation Groups (2023).  arXiv:2203.16049.
Data: [https://github.com/GeoTopChristy/HCPdm](https://github.com/GeoTopChristy/HCPdm).

**[Z95]** G. M. Ziegler.  *Lectures on Polytopes.*  Graduate Texts in Mathematics **152**.
Springer, New York, 1995.
