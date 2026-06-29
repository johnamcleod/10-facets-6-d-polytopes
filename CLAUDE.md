# Project Brief: Compact Hyperbolic Coxeter 6-Polytopes with 10 Facets

**Audience:** Claude Code, operating the enumeration/classification pipeline.
**Read this first, every session. It overrides any stale assumption in the codebase or commit history about the project goal.**

---

## 0. The one thing that changes everything

The project was started on the premise that **no compact Coxeter hyperbolic 6-polytope with 10 facets exists**, and that the pipeline would produce a *nonexistence proof*.

**This premise is false.** Such a polytope is already known: **Bugaenko's polytope, P^B₆** (10 = d+4 = 6+4). It was constructed in the 1980s–90s from reflective hyperbolic lattices over the golden-ratio ring ℤ[(1+√5)/2].

Consequences, which are not optional:

- **A nonexistence proof is impossible.** Do not write code, proofs, or documentation that aim to show the (d=6, n=10) family is empty. It is not empty.
- **If the pipeline returns zero polytopes in (d=6, n=10), that is a BUG, not a result.** The pipeline is failing to reproduce a known object. Treat any empty output here as a defect to be diagnosed, never as evidence.
- **The real goal is COMPLETENESS / UNIQUENESS:** determine whether P^B₆ is the *only* compact Coxeter 6-polytope with 10 facets, or whether others exist. This is the genuine open problem — the last unresolved dimension in the entire d+4 family.

If a task, prompt, or piece of code conflicts with the above, stop and flag it rather than proceeding.

---

## 1. Where (d=6, n=10) sits in the literature

Classification of compact hyperbolic Coxeter d-polytopes by k = n − d (facets minus dimension):

- **k=1 (simplices):** Lannér (1950). Exist only for d = 2, 3, 4.
- **k=2:** Kaplinskaja (1974) + Esselmann (1996). Exist for d = 3, 4, 5.
- **k=3:** Esselmann (1994) + Tumarkin (2007). Exist for d = 2…6 and d = 8 (none in d=7); unique in d=8 (Bugaenko).
- **k=4 (THIS PROJECT):** Felikson–Tumarkin (2008, arXiv:math/0510238) proved they exist **iff 2 ≤ d ≤ 7**, sharp, with d=7 **unique** (11 facets, Bugaenko). Burcroff (2024, Eur. J. Combin. 120:103957; arXiv:2201.03437) and independently Ma–Zheng classified **d=4 (348 polytopes)** and **d=5 (51 polytopes)**.

**Dimension 6 is the single unfinished case in the k=4 family.** Burcroff: *"the only remaining dimension where new polytopes may arise is d = 6,"* where *"the only known polytope was constructed by Bugaenko."* Uniqueness in d=6 is **not** proven (unlike d=7) — that is the open question.

Key structural facts already established for the d=6, 10-facet case (use these as hard constraints, see §4):

- **(Felikson–Tumarkin)** Any compact polytope with exactly one pair of non-intersecting facets has ≤ d+3 facets. So a (d=6, n=10) polytope has **at least two pairs of disjoint facets.**
- **(Burcroff)** A compact Coxeter 6-polytope with 10 facets **must contain a missing face of size 3 or 4.**
- Burcroff's order-type enumeration yields **265 candidate combinatorial types** in d=6 (vs. 34 candidates → 14 realized in d=4; 186 candidates → 6 realized in d=5). These 265 were **not** fully analyzed — that analysis is the open work.

---

## 2. Mathematical constraints the pipeline must respect (invariants)

A compact hyperbolic Coxeter d-polytope with n facets, normals e₁…eₙ, has Gram matrix G = (gᵢⱼ):

- gᵢᵢ = 1.
- gᵢⱼ = −cos(π/mᵢⱼ) when facets i, j meet at dihedral angle π/mᵢⱼ (mᵢⱼ ∈ ℤ, ≥ 2).
- gᵢⱼ = −1 if the facets are parallel (not allowed in the compact case — see below).
- gᵢⱼ ≤ −cosh(ρ) < −1 if the facets diverge (ultraparallel) at distance ρ.

Hard invariants — assert these in code and check every candidate against them:

1. **Signature (d, 1).** G has exactly **one** negative eigenvalue and rank d+1. For d=6, n=10: signature (6, 1) with a 3-dimensional kernel. Negative inertia index > 1 ⇒ **superhyperbolic** ⇒ reject immediately. This is the single most common silent-failure point: a tolerance or sign error here makes valid polytopes look superhyperbolic and disappear.
2. **Compact ⇒ simple polytope** (Vinberg). Every vertex lies on exactly d facets.
3. **No parabolic subdiagrams** (compact ⇒ no ideal vertices/cusps). Parabolic subdiagram ⇒ reject.
4. **Elliptic subdiagram of order j ↔ (d−j)-face.** Face lattice is read off the diagram.
5. **Every facet is itself a compact hyperbolic Coxeter (d−1)-polytope.** For d=6 each facet must be a compact Coxeter **5-polytope** — and those are now fully classified. Use the d=5 list as an admissibility filter on facets.
6. **Missing face of size m ↔ Lannér subdiagram of order m.** Lannér diagrams exist only up to order 5, so **every missing face has size ≤ 5**; the diagram must contain a Lannér subdiagram of order < 5.
7. **Golden-ratio entries are expected.** P^B₆ lives over ℤ[(1+√5)/2]; entries like −cos(π/5) = −(1+√5)/4 will appear. **Do arithmetic symbolically or in the number field, not in floating point**, when testing signature/realizability — float rounding on √5 terms is a known false-negative source.

---

## 3. FIRST PRIORITY — validate against known anchors before any new work

Do not attempt the d=6 classification until the pipeline provably reproduces every known object. In order:

1. **Recover P^B₆.** Obtain its Coxeter diagram / Gram matrix (see §5), feed it through the pipeline AND through CoxIter, and confirm: dimension 6, cocompact (compact), 10 facets. If the pipeline rejects P^B₆, fix that before anything else — the bug that drops P^B₆ is almost certainly the bug producing the false "empty" result.
2. **Recover the unique d=7, 11-facet polytope** (Bugaenko / Felikson–Tumarkin).
3. **Reproduce the published d=4 (348) and d=5 (51) lists** exactly. Cross-check counts against both Burcroff (348/51) and Ma–Zheng. (Note: Burcroff's 2021 thesis listed 341/50; these were corrected to **348/51** after reconciliation with Ma–Zheng. Use 348/51.)

**Likely failure modes to check, in priority order:**

- **Combinatorial-type generator drops types.** The 265 d=6 candidate types come from point-set order types of 10 points (heavy). If the generator is incomplete, P^B₆'s type is silently absent and output is empty. Verify the generator emits all 265 types and that P^B₆'s type is among them.
- **Signature test too strict / float-based.** See §2.1 and §2.7. Re-run with exact arithmetic over ℚ(√5).
- **Facet-admissibility filter too aggressive**, wrongly excluding a valid 5-polytope facet.
- **Missing-face/Lannér placement** not covering size-3 and size-4 missing faces (both are required possibilities in d=6).

A passing benchmark (P^B₆ recovered, d=4/d=5 lists reproduced) is the gate to §4. Until then, the pipeline's "no examples" output carries zero evidential weight.

---

## 4. The genuine open problem — d=6 completeness (only after §3 passes)

Specialize the Burcroff / Ma–Zheng method to d=6:

1. Generate all **265** combinatorial types via point-set order types of 10 points → affine Gale diagrams (positive/negative planar point configuration; general position; exactly 2 positive points inside the convex hull of the negative points).
2. For each type, compute the **missing-face list**; discard types violating the necessary conditions: **≥ 2 disjoint facet pairs**, **a missing face of size 3 or 4**, and **every facet realizable as a compact Coxeter 5-polytope** from the complete d=5 classification.
3. For surviving types, place Lannér subdiagrams on missing faces; iterate admissible small dihedral angles (≤ π/6) constrained by the Lannér lists.
4. Solve the remaining Gram entries under the **rank/signature (6,1)** condition, **using exact arithmetic over ℚ(√5)**.
5. Discard superhyperbolic / parabolic / non-realizable solutions. Use **Vinberg local determinants** (factorization of det along a shared vertex/edge) to test superhyperbolicity fast without full eigen-decomposition.
6. **Validate every surviving Gram matrix with CoxIter** (dimension, cocompactness, facet count).

Expected — but unproven — outcome: only P^B₆ survives ⇒ uniqueness. If a *second* valid, CoxIter-confirmed diagram appears, the headline result flips to "≥ 2 members" and the task becomes full enumeration. Either way it resolves the last open d+4 case; both are publishable.

Optional accelerators: Alexandrov's ×₀-product superhyperbolicity results (Trans. AMS 376 (2023), 6989–7012; arXiv:2203.07248) and Felikson–Tumarkin lifting/face-recursion lemmas, to prune large families of glued-Lannér types a priori.

---

## 5. Where to get the P^B₆ diagram (for §3 step 1)

In order of usefulness:

- **Ma–Zheng data/code repo (most practical):** `github.com/GeoTopChristy/HCPdm` — ships intermediate data and Gram matrices in machine-readable form; fastest path to a matrix you can feed straight to CoxIter.
- **Primary source:** V. O. Bugaenko, *Arithmetic crystallographic groups generated by reflections, and reflective hyperbolic lattices*, Adv. Soviet Math. **8**, AMS (1992), 33–55. (Also: Bugaenko, Moscow Univ. Math. Bull. **39** (1984), 6–14, golden-ratio ring.) Explicit diagrams of the golden-ratio forms and their fundamental polytopes.
- **Secondary cross-check:** Burcroff, *On Compact Hyperbolic Coxeter Polytopes with Few Facets*, Durham MSc thesis (2021), `etheses.dur.ac.uk/14202` — discusses d=6 and gives diagram data/notation; cites the Bugaenko polytope as [10].

**Do NOT use** the Felikson–Tumarkin d+4 paper (arXiv:math/0510238) as a source for the *picture* — it covers only the d≥8 nonexistence and d=7 uniqueness and does not redraw the d=6 example. **Do NOT confuse** P^B₆ (10 facets) with "Bugaenko's 6-polytope with 34 facets" (a different, larger object via Allcock).

If transcribing from Bugaenko by hand: regenerate the Gram matrix, then confirm in CoxIter before trusting it. A single mislabeled dihedral angle (especially on the π/5 / golden-ratio edges) makes a valid polytope read as superhyperbolic and vanish.

---

## 6. Tooling

- **CoxIter** (`github.com/rgugliel/CoxIter`, R. Guglielmetti): given a Coxeter diagram, checks cocompactness/cofiniteness/arithmeticity and computes f-vector, Euler characteristic, signature, dimension, growth. **The standard verifier — every candidate must pass CoxIter.** Also ships a library of known-polytope graph files.
- **AlVin / VinAl** (Guglielmetti; Bogachev–Perepechko): implement Vinberg's algorithm for quadratic forms — the way to *re-derive* P^B₆ from its ℤ[(1+√5)/2]-form if you prefer construction over transcription.
- **Exact arithmetic over ℚ(√5):** use a symbolic/number-field backend (e.g. SageMath, or sympy with the field extension) for all signature/realizability tests. Floating point on √5 terms is banned for accept/reject decisions.

---

## 7. Definition-of-done

- [ ] §3 anchors all recovered: P^B₆ (d=6, 10 facets, compact via CoxIter), the d=7 unique polytope, and the d=4 (348) / d=5 (51) lists reproduced exactly.
- [ ] All 265 d=6 candidate combinatorial types generated; P^B₆'s type confirmed present.
- [ ] Each type processed through the §4 pipeline with exact ℚ(√5) arithmetic; every surviving Gram matrix CoxIter-validated.
- [ ] Result stated as uniqueness (only P^B₆) **or** an explicit list of all compact Coxeter 6-polytopes with 10 facets.
- [ ] No artifact anywhere in the repo claims or assumes the (d=6, n=10) family is empty.

---

## 8. Key references

- Felikson, Tumarkin, *On compact hyperbolic Coxeter d-polytopes with d+4 facets*, Trans. Moscow Math. Soc. 69 (2008), 105–151; arXiv:math/0510238.
- Felikson, Tumarkin, *Coxeter polytopes with a unique pair of non-intersecting facets*, J. Combin. Theory A 116 (2009); arXiv:0706.3964.
- Burcroff, *Near classification of compact hyperbolic Coxeter d-polytopes with d+4 facets and related dimension bounds*, Eur. J. Combin. 120 (2024) 103957; arXiv:2201.03437. (Thesis: Durham 2021, etheses.dur.ac.uk/14202.)
- Ma, Zheng, *Compact hyperbolic Coxeter four-polytopes with eight facets*, J. Algebraic Combin. 59 (2024), 225–290; arXiv:2201.00154. *...five-dimensional polytopes with nine facets*, Transform. Groups (2023); arXiv:2203.16049. Data: github.com/GeoTopChristy/HCPdm.
- Alexandrov, *Lannér diagrams and combinatorial properties of compact hyperbolic Coxeter polytopes*, Trans. AMS 376 (2023), 6989–7012; arXiv:2203.07248.
- Bugaenko (1984, 1992) — primary construction of P^B₆ (see §5).
- Vinberg — absence of compact Coxeter polytopes in d ≥ 30; local determinants; the algorithm.
