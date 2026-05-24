# Classification of Compact Hyperbolic Coxeter 6-Polytopes with 10 Facets — Implementation Specification

**Target for:** Claude Code (autonomous multi-stage research-software build)
**Mathematical goal:** Produce the complete, verified list (possibly empty) of compact hyperbolic Coxeter 6-polytopes with 10 facets — i.e. the `d=6, n=d+4` case — closing the last open dimension in the `n ≤ d+4` classification program.

---

## 0. Context and what "done" means

This is the final unresolved case in a program that is otherwise complete:

- Felikson–Tumarkin proved no compact hyperbolic Coxeter `d`-polytope with `d+4` facets exists for `d ≥ 8`, with the unique `d=7` example (Bugaenko). 
- `d=4` (348 polytopes) and `d=5` (51 polytopes) were settled independently by Burcroff (*Near Classification…*, arXiv:2201.03437) and Ma–Zheng (arXiv:2201.00154, arXiv:2203.16049), using **affine Gale diagrams enumerated via planar point-set order types**. The two independent methods agreeing is what gives confidence in these delicate lists.
- `d=6, n=10` is explicitly flagged by Burcroff as the only remaining dimension where new `d+4` polytopes may arise.

**Definition of success.** A reproducible pipeline that:
1. Enumerates every combinatorial type of simple 6-polytope with 10 facets that *could* carry a compact hyperbolic Coxeter structure (via order types), with each type produced exactly once (no isomorphic duplicates).
2. For each surviving type, determines all admissible Coxeter diagrams and decides Gram-matrix realizability in `H^6` (signature `(6,1)`, rank 7) using **exact arithmetic**.
3. Emits a final list of Coxeter diagrams + Gram matrices, each independently verified by **CoxIter** (cocompact + dimension 6), with a machine-checkable certificate per polytope.
4. Cross-checks the count against an independent re-derivation path before claiming completeness.

**Non-negotiable constraints.** Exact arithmetic throughout the realizability stage (no floating point in any step that can rule a polytope in or out). Every elimination must cite the lemma it invokes and log the witness (the missing-face set, the Lannér subdiagram, the failed minor, etc.). Completeness claims require the independent cross-check in Stage 6 to agree.

---

## 1. Mathematical objects and conventions

Fix `d = 6`, `n = 10`, so `k := n - d = 4` and the Gale dimension is `n - d - 2 = 2` (planar Gale diagrams). Index facets by `I = {0,1,...,9}`.

**Gram matrix.** For a compact Coxeter polytope with outward unit facet normals `e_i` in `R^{6,1}`, the Gram matrix `G = (⟨e_i, e_j⟩)` is symmetric `10×10` with:
- `G_ii = 1`.
- For facets meeting at dihedral angle `π/m_ij` (`m_ij ≥ 2` integer): `G_ij = -cos(π/m_ij)`.
- For facets that **do not** intersect (disjoint, "dotted edge"): `G_ij ≤ -1` (`= -1` parallel/ideal — excluded in compact case so `< -1`; `< -1` ultraparallel).
- Signature `(6,1)`, hence **rank exactly 7**. This rank-7 condition on a 10×10 matrix is the master arithmetic constraint: it forces 3 independent linear/algebraic dependencies among rows.

**Coxeter diagram `Σ`.** Nodes = facets. Between nodes `i,j`:
- no edge if `m_ij = 2` (orthogonal);
- edge labelled `m_ij` (or with `m_ij - 2` for the standard weight convention) if `3 ≤ m_ij < ∞`; usual graphical conventions: unlabelled = 3, bold/labelled 4, etc.;
- **dotted edge** if facets disjoint (`G_ij < -1`), carrying a real weight `> 1` (i.e. `-G_ij`).

**Missing face.** A subset `M ⊆ I` is a *missing face* if the facets in `M` have empty common intersection in `P` but every proper subset has nonempty intersection — equivalently a minimal non-face. The set of missing faces is the **missing-face list** and determines the combinatorial type (together with the nerve). Crucially: a missing face of size `s` corresponds to a **Lannér subdiagram** on those `s` nodes in any compact realization (the corresponding sub-Coxeter-system is a compact hyperbolic simplex group / finite-covolume diverging configuration). Disjoint pairs of facets are exactly the **size-2 missing faces** ("dotted edges").

**Disjoint-facet count `p`.** Number of pairs of non-intersecting facets. Two key facts that scope the search:
- (Felikson–Tumarkin, p=1 theorem) A compact polytope with exactly one disjoint pair has `n ≤ d+3`. Since here `n = d+4`, **`p ≥ 2`**.
- (Felikson–Tumarkin essential-polytopes) For `p ≤ n-d-2 = 2` the family is finite and finitely listable. We sit exactly at this boundary, so `p ∈ {2, 3, ...}` but finiteness is guaranteed.

---

## 2. The master theorem the pipeline is built on

**(Burcroff, Thm 3.5; Ma–Zheng equivalent.)** Every compact Coxeter `d`-polytope with `d+4` facets admits an **affine Gale diagram** consisting of an arrangement `A ⊆ R^2` of `d+4` points in general position, with a choice of **exactly two** of the points (lying in the interior of `conv(A)`) designated *positive* and the rest *negative*. The **combinatorial type of `P` is completely determined by the order type (chirotope) of the configuration** together with the positive/negative sign vector.

This is the linchpin. It converts "enumerate combinatorial types of simple 6-polytopes with 10 facets" — hopeless directly — into "enumerate planar order types on 10 points, then choose 2 interior points as positive." Order types quotient out coordinate redundancy automatically: **each combinatorial type is generated exactly once**, which is precisely the non-redundancy requirement.

**Realizability warning (must be respected, not ignored).** Polytopes with `d+4` facets are at the threshold where **Mnëv universality** appears for the realization space. This means: (a) an abstract order type / oriented matroid may be *non-realizable* over the reals, and (b) even a realizable combinatorial type may have a realization space that is not "nice." Consequences for us:
- The enumeration must distinguish **realizable** planar order types from abstract (non-realizable) ones. Use realizable order-type data only, or certify realizability.
- The Gram-realizability stage (Stage 4) is a real-semialgebraic feasibility problem and must be treated with exact methods (Gröbner / cylindrical algebraic decomposition / exact eigenvalue-signature tests), not numerically.

---

## 3. Pipeline overview (stages, each a separate verifiable artifact)

```
Stage 1  Enumerate realizable planar order types on 10 points        -> order_types/
Stage 2  Build affine Gale diagrams (choose 2 interior +)            -> gale/
         => combinatorial types of simple 6-polytopes, 10 facets
Stage 3  Combinatorial filtering (missing faces, p>=2, Lannér,        -> types_surviving/
         elimination lemmas)  -- kills most types
Stage 4  Coxeter-diagram assignment + EXACT Gram realizability        -> diagrams/, gram/
Stage 5  Verification with CoxIter + independent cert                 -> verified/
Stage 6  Independent cross-check path + final reconciliation          -> RESULT/
```

Each stage reads the previous stage's on-disk artifacts and writes its own, with a manifest (JSON) recording counts, parameters, code version hash, and per-item provenance. No stage recomputes a previous stage silently.

---

## 4. Stage 1 — Realizable planar order types on 10 points

**Objective.** A canonical, deduplicated list of all *realizable* order types of 10 points in general position in the plane, each with an explicit integer-coordinate realization (so downstream stages have exact coordinates).

**Data source decision (in priority order):**
1. **Preferred:** Aichholzer–Aurenhammer–Krasser order-type database. The complete realizable order-type database extends to **n = 11** (the 10-point and 11-point sets are fully enumerated; the 10-point database stores explicit small-integer coordinate realizations, historically 16-bit coordinates). Acquire the 10-point file; it is the authoritative realizable list and side-steps the realizability problem entirely for Stage 1.
2. **Fallback / cross-check generator:** **NumPSLA** (Rote, arXiv:2503.02336) generates abstract order types for up to 12 points from scratch and complements the AAK database. Use it to *regenerate* the 10-point realizable list independently as a consistency check on the database ingestion, and as the engine if database acquisition fails.
3. **Oriented-matroid route (degenerate-aware backup):** Finschi–Fukuda / Miyata et al. databases of rank-3 oriented matroids on ≤12 elements, filtered to realizable uniform ones. Heavier; use only if 1 and 2 are unavailable.

**Why 10 is tractable:** the realizable order-type count for 10 points is large but fully enumerated and storable; 11 points is the current practical edge (≈2.3 billion, ~100 GB), so 10 is comfortably inside the feasible regime. **Do not attempt to brute-force coordinates.** Ingest/generate order types.

**Implementation notes.**
- Represent each order type by its **chirotope**: the map `χ: (i,j,k) ↦ sign(det[p_j - p_i, p_k - p_i])` over all triples. Store canonically (lexicographically minimal chirotope over relabelings, or the database's RevLex index).
- Keep an explicit exact-integer coordinate realization per order type (from the database, or recovered via the generator). These coordinates feed Stage 2.
- **Acceptance test for Stage 1:** the count of realizable 10-point order types matches the published database cardinality; spot-check 20 random entries by recomputing all `C(10,3)=120` triple orientations from stored coordinates and confirming they match the stored chirotope.

**General position note.** The master theorem requires `A` in general position (no 3 collinear). Realizable *uniform* rank-3 oriented matroids = general-position planar point sets, which is exactly the AAK database content. Degenerate configurations are excluded by the theorem's hypothesis — confirm none leak in.

---

## 5. Stage 2 — Affine Gale diagrams → combinatorial types

**Objective.** From each order type, produce every valid affine Gale diagram for a *simple* 6-polytope with 10 facets, and from it the combinatorial data (vertex–facet incidence / missing-face list).

**Construction (per order type `A` with coordinates):**
1. Choose an unordered pair `{u,v} ⊂ A` to be the **positive** points; the other 8 are negative. By the theorem, both positive points must lie in the **interior of `conv(A_-)`** (Burcroff's lemma: a positive point not interior to the negative hull cannot give a polytope). Enforce this interiority as a fast pre-filter using exact orientation tests.
2. The resulting signed configuration is an affine Gale diagram. Recover the face lattice / cofaces via the standard Gale duality dictionary:
   - A subset `S ⊆ I` of facets forms a **face** of `P` iff the complementary points (the Gale-dual points indexed by `I \ S`) **positively span** (their relative interior of the positive hull contains the origin in the linear Gale picture; in the affine signed picture, the standard "the positive points among the complement are captured by the convex hull of the negative ones" coface criterion). Implement the precise affine-Gale coface test exactly as in Burcroff §3 / Ziegler *Lectures on Polytopes* Ch. 6. **This must be coded against the textbook criterion, with unit tests on `d=4` and `d=5` reproducing the known 348 and 51 counts.**
   - **Missing faces** are the minimal non-cofaces; compute them directly.
3. Enforce **simplicity** of `P` (every vertex on exactly `d=6` facets). Equivalently a condition on the Gale diagram (general position of the signed configuration). Non-simple types are *not* automatically excluded for compact Coxeter polytopes in general, BUT: by Felikson–Tumarkin the relevant compact `d+4` polytopes are simple — **verify this reduction explicitly from the literature before discarding non-simple types**; if simplicity is not fully guaranteed in `d=6`, retain non-simple types in a separate bucket rather than dropping them. (See Risk R3.)

**Deduplication.** Two (order type, positive-pair) choices can give isomorphic polytopes. Canonicalize each combinatorial type by computing a canonical form of its **missing-face hypergraph** (e.g. via `nauty`/`bliss` on the bipartite facet–missingface incidence, or canonical relabeling of the face lattice). Store by canonical key; collisions are the same type. This is where "constructed once" is enforced at the combinatorial level.

**Acceptance test for Stage 2:** Re-run the whole Stage 1→2 path for `d=4,n=8` and `d=5,n=9` (just change parameters) and reproduce **exactly 348 and 51** final polytopes after Stages 3–5, and the intermediate combinatorial-type counts reported by Burcroff/Ma–Zheng. This regression test is the single most important correctness gate in the project — do not proceed to trust `d=6` output until `d=4,5` reproduce.

---

## 6. Stage 3 — Combinatorial filtering (the heavy pruning)

Apply, in increasing order of cost, eliminations that remove combinatorial types before any Gram work. Every elimination logs `{type_id, rule, witness}`.

**F1 — Disjoint-pair count.** Compute `p` = number of size-2 missing faces. Discard any type with `p < 2` (Felikson–Tumarkin p=1 theorem forbids `p=1` at `n=d+4`; `p=0` gives simplices/Esselmann only). Keep `p ≥ 2`.

**F2 — Lannér-compatibility of missing faces.** Each missing face must be realizable as a Lannér diagram (compact hyperbolic simplex) on its node set. Lannér diagrams exist only in ranks `2..5` (dimensions of compact simplices ≤4, i.e. simplex on ≤5 facets). Therefore **every missing face has size between 2 and 5**; a missing face of size ≥6 is fatal — discard the type. (A size-`s` missing face needs a rank-`s` Lannér subdiagram; none exist for `s ≥ 6`.)

**F3 — Pairwise Lannér disjointness obstruction.** Burcroff's Lemma 5.7-type result: certain configurations of induced missing-face lists are non-realizable. Concretely, **two Lannér subdiagrams (missing faces) cannot be vertex-disjoint** in a compact Coxeter polytope of this class without forcing a contradiction in specific patterns; and specific small induced missing-face lists (e.g. the pattern isomorphic to `{0123, 014, 235}` used to kill `d=4` types G22–G24) are forbidden. Encode the known forbidden induced-sublist patterns from Burcroff §6 as a **pattern library**; discard any type containing a forbidden induced pattern. (Port the exact list from arXiv:2201.03437 §5–6; treat it as data, extensible.)

**F4 — "No two Lannér subdiagrams disjoint" / connectivity of the dotted structure.** Implement the global condition (FT) that the disjoint-facet structure cannot decompose the diagram in forbidden ways. Log each kill with the specific separating set.

**F5 — Rank feasibility (combinatorial precheck).** The eventual Gram matrix must have rank 7. A purely combinatorial necessary condition: the diagram cannot contain an induced **elliptic (finite-type) subdiagram of rank > 7** that would over-determine the form, nor a parabolic/affine subdiagram pattern incompatible with signature `(6,1)`. Implement the Vinberg sign/spectral necessary conditions at the combinatorial level where they don't need weights yet.

**Output.** A (hopefully small — `d=4` had ~30 combinatorial types reaching this point, `d=5` fewer) set of surviving combinatorial types, each with its missing-face list, disjoint-pair set, and the *abstract* Coxeter diagram skeleton (which edges are forced absent (orthogonal), which are dotted (disjoint), and which are "ordinary" edges with as-yet-unknown labels).

**Acceptance test for Stage 3:** On `d=4`, reproduce Burcroff's surviving combinatorial-type list (the `G1..G21` survivors after killing `G22..G30`, matching her Corollaries 6.2/6.3). Mismatch ⇒ stop and fix the lemma encodings.

---

## 7. Stage 4 — Coxeter labels + exact Gram realizability

For each surviving combinatorial type with its abstract diagram, decide which integer dihedral labels on the ordinary edges yield a Gram matrix of signature `(6,1)`. **This is the mathematically hardest stage; treat it adversarially.**

### 7.1 Bounding the labels (make the search finite and small)

- **Angle set.** Dihedral angles are `π/m`, `m ∈ {2,3,4,5,...}`. Compactness + Lannér constraints sharply bound large labels. Apply Burcroff's **low-weight lemma** (an ordinary edge inside the relevant subdiagrams must have low weight) and Esselmann/FT lemmas to cap `m`. Establish a proven per-edge cap (typically `m ≤ 5`, occasionally a specific higher value); **prove or cite the cap for each edge class rather than assuming `m ≤ 5` globally.** Log the justification.
- **Dotted-edge weights** are real (`> 1`), not free integers; they are solved for, not enumerated (see 7.3).
- Result: each surviving type yields a **finite** set of candidate integer-label assignments on ordinary edges (a small product set after caps and local elliptic constraints on each vertex figure — every facet's neighborhood must be a finite/elliptic Coxeter diagram because vertices are genuine, i.e. each rank-6 vertex subdiagram is positive definite).

### 7.2 Local elliptic constraints (cheap pruning before global solve)

For every vertex of `P` (set of 6 facets meeting at a point), the corresponding `6×6` principal submatrix of `G` must be **positive definite** (a spherical/finite Coxeter group — the link of the vertex is a spherical simplex). Enumerate vertices from Stage 2 incidence data; for each candidate label assignment, test all vertex submatrices for positive definiteness (exact: all leading principal minors `> 0`, computed over the relevant number field). Most label assignments die here. This is the discrete analogue of "every facet of the polytope is itself a spherical Coxeter polytope."

### 7.3 The global realizability problem (exact)

After local pruning, for each surviving labelled diagram we must decide existence of the dotted-edge weights making `G` have signature `(6,1)` and rank `7`, with all dotted weights `> 1`.

Set up `G` symbolically: known entries from labels (`-cos(π/m)`, algebraic numbers in `Q(cos π/m)` — i.e. in a real cyclotomic field; for `m ∈{2,3,4,5}` the field is `Q`, `Q`, `Q(√2)`, `Q(√5)` respectively), and unknowns `x_ab = -G_ab > 1` on dotted edges.

**Constraints:**
1. `rank(G) = 7` ⇔ all `8×8` minors vanish (gives polynomial equations in the `x_ab`).
2. Signature `(6,1)`: exactly one negative eigenvalue. Combined with rank 7 and the positive-definite vertex blocks, this is checked via an exact symmetric-matrix inertia computation (e.g. exact `LDL^T` / leading-minor sign sequence in the number field) once the `x_ab` are pinned.
3. `x_ab > 1` for each dotted edge (real inequality).
4. Vinberg's **local determinant** criteria: use the recursive local-determinant relations to reduce the rank/signature conditions to a smaller solvable system and to prune impossible sign patterns early.

**Solving method (exact, in order of preference):**
- The rank conditions (vanishing `8×8` minors) typically pin the unknowns to **finitely many algebraic values** (often the dotted weights are forced to specific values like `2cos(π/m)`-type or quadratic irrationals). Compute the solution variety with an **exact Gröbner basis** over `Q` (work in the compositum number field for the `cos` entries, or rationalize via minimal polynomials). Use a CAS with exact arithmetic: **Sage/Singular/Macaulay2** for Gröbner; **Sage/PARI** number fields.
- For each real solution branch, test the inequalities `x_ab > 1` and the signature exactly. Real solutions: isolate with exact real-root methods (CAD via **QEPCAD-B**/**Maple** if needed, or `Sage`'s exact real fields). Mnëv universality means we cannot assume the solution set is a point — handle positive-dimensional components by extracting whether any real point satisfies all inequalities (a real-feasibility question = existential theory of the reals; CAD or critical-point methods).
- **No floating point may decide feasibility.** Floating point may be used only to *guide* (e.g. guess which branch), never to accept/reject.

**Output.** For each combinatorial type: the finite set of fully-determined Gram matrices (entries as exact algebraic numbers) of signature `(6,1)`, rank 7, satisfying all constraints — i.e. candidate compact hyperbolic Coxeter 6-polytopes.

**Acceptance test for Stage 4:** reproduce the exact Gram matrices / Coxeter diagrams of the known `d=5,n=9` polytopes (all 51) and a sample of `d=4,n=8` (the 348), entry-for-entry, via this same solver path.

---

## 8. Stage 5 — Verification with CoxIter and certificates

For every candidate Gram matrix/diagram from Stage 4:

1. **CoxIter** (Guglielmetti) — the community-standard checker. Feed the Coxeter diagram; confirm:
   - dimension `= 6`,
   - **cocompact = true** (compactness),
   - finite covolume,
   - compute the f-vector (must be the f-vector of a simple 6-polytope with 10 facets), Euler characteristic, signature, arithmeticity, growth series. Install from source (`https://github.com/rgugliel/CoxIter`); it links PARI/GP and is exact for these tests.
2. **Independent in-house certificate** (do not trust a single tool): recompute signature/inertia of `G` exactly; verify each vertex block positive definite; verify each dotted entry `< -1`; verify rank 7; verify the abstract diagram's f-vector matches the Stage 2 combinatorial type. A polytope is **accepted** only if CoxIter and the in-house certificate agree on every property.
3. **Isometry/duplicate check:** two diagrams may give isometric polytopes. Canonicalize accepted polytopes by graph-isomorphism of the labelled (including dotted-weight) Coxeter diagram; report distinct isometry classes.

**Output.** `verified/` — the accepted polytopes with full data: Coxeter diagram, exact Gram matrix, f-vector, CoxIter report, in-house certificate, arithmeticity.

---

## 9. Stage 6 — Independent cross-check and completeness argument

The `d=4,5` results were trusted only because **two independent methods agreed**. Replicate that discipline:

- **Path A (primary):** the order-type→Gale pipeline above.
- **Path B (independent combinatorial generation):** generate simple 6-polytopes with 10 facets by an unrelated route and intersect the combinatorial-type list with Path A's. Options: (i) the **gift-wrapping / beneath-beyond** inductive enumeration of simple `d`-polytopes with `n` facets (Koranne-style, or via `polymake`'s polytope enumeration filtered to simple, 10 facets, dim 6); (ii) oriented-matroid enumeration (Finschi–Fukuda rank-`(n-d)=4`? — note the *linear* Gale picture has rank `n-d-1=3`; reconcile the two dual encodings and confirm they agree). Reconcile Path B's combinatorial types with Stage 3 survivors; any discrepancy must be explained.
- **Parameter-sweep regression:** the entire pipeline runs on `d=4` and `d=5` as automated tests in CI and must reproduce `348` and `51`. Lock these as golden tests.
- **Completeness statement.** Document the chain of theorems that makes the order-type enumeration *exhaustive* (master theorem 3.5 ⇒ every such polytope appears as some order-type+positive-pair), so that "we examined all order types on 10 points" ⇒ "we examined all combinatorial types" ⇒ (after exact realizability) "we found all polytopes." State each external theorem invoked with citation.

**Final deliverable `RESULT/`:**
- `polytopes.json` — the verified list (possibly empty) with all certificates.
- `report.md` — human-readable summary: count, each polytope's diagram + Gram matrix + invariants, the elimination census (how many types died at each filter and why), and the completeness argument.
- `provenance.json` — code version, data-source versions/hashes, tool versions, per-stage counts.
- If the list is **empty**, that is itself the headline result (no compact hyperbolic Coxeter 6-polytope with 10 facets exists) and must carry the full completeness argument; if **non-empty**, each polytope is a new example near the top of the known range.

---

## 10. Software/environment

| Need | Tool | Notes |
|---|---|---|
| Order-type data (n=10) | AAK order-type database | authoritative realizable list w/ integer coords |
| Order-type generation / cross-check | NumPSLA (Rote) | n≤12; independent regenerate |
| OM backup | Finschi–Fukuda / Miyata DBs | degenerate-aware, heavier |
| Polytope combinatorics, Gale duality | `polymake` | face lattices, Gale, simplicity, f-vectors |
| Graph canonical form / iso | `nauty`/`bliss` | dedup combinatorial types & final diagrams |
| Exact algebra / Gröbner | Sage + Singular/Macaulay2 | number fields `Q(cos π/m)`, minor ideals |
| Exact real feasibility | QEPCAD-B / Maple CAD / Sage | only if positive-dim solution components arise |
| Coxeter verification | CoxIter (github rgugliel/CoxIter) | cocompactness, dim, invariants; links PARI |
| Orchestration | Python | stage drivers, manifests, CI golden tests |

All exact-arithmetic stages: **no float decides inclusion/exclusion.** Pin every tool version in `provenance.json`.

---

## 11. Risk register (read before coding)

- **R1 — Realizability (Mnëv).** At `d+4`, realization spaces can be arbitrarily complicated. *Mitigation:* use realizable order-type data in Stage 1 (side-steps OM non-realizability); in Stage 4 treat Gram realizability as exact existential-theory-of-the-reals, never numeric; handle positive-dimensional solution components explicitly.
- **R2 — Correctness of Gale coface code.** Easy to get the affine-Gale coface/missing-face criterion subtly wrong. *Mitigation:* unit-test against `d=4`(348) and `d=5`(51) before touching `d=6`. This is the gate.
- **R3 — Simplicity assumption.** Confirm from FT/Burcroff whether *all* compact `d+4` polytopes in `d=6` are simple. If not fully guaranteed, retain non-simple types in a separate bucket and run the Gale dictionary in its general (non-simple) form rather than silently dropping them.
- **R4 — Label-cap soundness.** A wrong (too-tight) cap on `m_ij` could miss a polytope. *Mitigation:* prove/cite each cap; where uncertain, set the cap generously and let the exact solver discard — never cap below a proven bound for speed.
- **R5 — Scale.** 10-point order types are numerous; Stage 2's pair-choice multiplies by `C(10,2)=45` then interiority-filters. *Mitigation:* push the cheap interiority + `p≥2` + missing-face-size filters as early as possible (before canonicalization), parallelize Stage 2/3 embarrassingly, store survivors only.
- **R6 — Trusting one tool.** *Mitigation:* dual verification (CoxIter + in-house) and dual generation (Path A/B). Disagreement halts the claim.
- **R7 — Number-field plumbing.** Mixing `cos(π/5)` (needs `√5`) with `cos(π/4)` (needs `√2`) in one matrix needs a common field. *Mitigation:* build the compositum up front per diagram; keep minimal polynomials; test inertia in that field exactly.

---

## 12. Suggested build order for Claude Code

1. Scaffold repo, manifests, CI. Implement the **`d=4`/`d=5` regression harness first** with the parameter `d` abstract.
2. Stage 1 ingest (AAK) + Stage 1 cross-gen (NumPSLA) for n = d+4 with `d=4` (n=8) — small, fast — and verify against the known order-type counts.
3. Stage 2 Gale dictionary in `polymake`/own code; **make `d=4` reproduce its combinatorial-type census.**
4. Stage 3 filters with the Burcroff pattern library; **reproduce `d=4` survivors (G22–G30 killed).**
5. Stage 4 exact solver; **reproduce all 348 (`d=4`) then all 51 (`d=5`) Gram matrices.** Only now is the machinery trusted.
6. Run `d=5` end-to-end as a second golden test.
7. **Run `d=6, n=10`.** Then Stage 6 Path B independent generation and reconcile.
8. Write `RESULT/report.md` with the elimination census and completeness argument.

**Do not run `d=6` for record until steps 2–6 reproduce the published `d=4` and `d=5` numbers exactly.**

---

## 13. Key references (for the implementer to pull and encode)

- Burcroff, *Near Classification of Compact Hyperbolic Coxeter d-Polytopes with d+4 Facets…*, arXiv:2201.03437 — **the method paper**: Thm 3.5 (Gale/order-type), §5–6 elimination lemmas (low-weight, forbidden missing-face patterns), `d=4,5` results.
- Burcroff, MSc thesis, Durham (2021), etheses.dur.ac.uk/14202 — fuller proofs; the `d=6`-must-have-a-… proposition.
- Ma–Zheng, arXiv:2201.00154 (`d=4`) and arXiv:2203.16049 (`d=5`) — independent method for cross-checking logic.
- Felikson–Tumarkin, arXiv:math/0510238 (`d≥8` impossibility, `d=7` unique), arXiv:0706.3964 (p=1 ⇒ n≤d+3), arXiv:0906.4111 (essential polytopes, finiteness for p≤n-d-2).
- Vinberg, *Absence…* (local determinants, signature criteria) — Stage 4 inertia tooling.
- Aichholzer–Aurenhammer–Krasser, *Enumerating order types…*, Order 19 (2002) + n=11 extension — Stage 1 data.
- Rote, *NumPSLA*, arXiv:2503.02336 — Stage 1 generation/cross-check.
- Guglielmetti, **CoxIter**, github.com/rgugliel/CoxIter — Stage 5 verification.
- Ziegler, *Lectures on Polytopes*, Ch. 6 — Gale duality reference for Stage 2.
