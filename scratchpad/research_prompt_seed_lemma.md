# Research task (for a fresh Claude instance with web access / deep-research)

## Why this matters
We reproduce compact hyperbolic Coxeter censuses (d=5: 51, exact) with a **two-phase method**:
1. SEED pass — force every prism-end facet orthogonal ("l_basis"/basis condition) → the
   *basis (seed)* polytopes. Cheap; eliminates most combinatorial types (0 seeds).
2. CENSUS pass — relax orthogonality on the seed-bearing types → all polytopes, via gluing
   compact simplicial prisms onto the orthogonal ends.

For a **uniqueness / non-existence** result in d=6 (is Bugaenko's P^B₆ the ONLY compact
Coxeter 6-polytope with 10 facets?), phase 1 must be *provably complete*: skipping a type on
the grounds that it has 0 seeds is only valid if **every realizable polytope of that type has
at least one orthogonal-prism-end (basis) representative.** We call this the SEED LEMMA:

    (0 basis seeds for a combinatorial type)  ⇒  (that type realizes 0 compact polytopes).

It held *empirically* for d=5 (the 6 seed-bearing types were exactly the 6 realizing types),
but we have not proven it or found it stated. If it is false even once in d=6, the pipeline
could silently miss a second polytope and report a FALSE "unique."

## Exactly what we need
1. **Is the seed lemma a theorem in Ma–Zheng?** Find the precise statement, in their
   d=5 paper (arXiv:2203.16049, esp. Section 5 and the "basis" / l5_basis definition, Table 13)
   and/or the d=4 paper (arXiv:2201.00154, chcp48 / l4_basis). Quote the lemma/proposition and
   its proof sketch. Does every compact polytope arise by gluing simplicial prisms onto a
   basis (orthogonal-end) polytope? Is the basis set proven COMPLETE (a complete set of seeds)?
2. **The exact prism-gluing construction.** State precisely: given a basis polytope with an
   orthogonal prism-end facet, which compact simplicial d-prisms may be glued, and how the
   glued Gram matrix / Coxeter diagram is formed. Enough detail to implement and to argue
   completeness. Cite the relevant Felikson–Tumarkin lemmas (arXiv:math/0510238, arXiv:0706.3964)
   the construction relies on.
3. **Does the lemma/construction depend on dimension or on the facet being a simplex?** In
   d=4 the basis condition OVER-prunes (real P2 polytopes have non-orthogonal tetrahedral-facet
   ridges — census is basis-OFF there), whereas in d=5 basis-ON gives the seeds and gluing
   gives the rest. Clarify why d=4 and d=5 differ and what the correct statement is for d=6.
4. **Any explicit statement about d=6 / P^B₆** being of prism type, or how the d+4 family's
   prism structure is used in Felikson–Tumarkin's d=7 uniqueness proof (arXiv:math/0510238) —
   since d=7 uniqueness is proven, its method is the template for a d=6 uniqueness argument.

## Sources
- Ma–Zheng d=5: arXiv:2203.16049 (full PDF, Section 5 + Table 13). d=4: arXiv:2201.00154.
- Felikson–Tumarkin d+4: arXiv:math/0510238 (has the d=7 UNIQUENESS proof — study its
  structure); "unique pair of non-intersecting facets": arXiv:0706.3964.
- Burcroff: arXiv:2201.03437 (+ Durham thesis etheses.dur.ac.uk/14202).

## Output
State whether the seed lemma is (a) proven in the literature — with citation and the argument,
(b) plausible but unstated, or (c) false / conditional. Give the prism-gluing construction
precisely enough to implement. Flag inferred vs directly-stated. This determines whether our
pipeline can ever PROVE d=6 uniqueness or only DISCOVER additional polytopes.
