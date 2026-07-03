# Prism / "l_basis" condition — derivation from Ma-Zheng source (2026-07-03)

Session goal (branch feat/prism-condition): derive the correct prism-end condition,
which the last commit (525a08a) flagged as the real open blocker for d=5 (and d=6).

## What the condition IS (read from HCPdm source — CONFIRMED)

Ma-Zheng name: `l4_basis` (chcp48, d=4) / `l5_basis` (chcp59, d=5). Definition is
IDENTICAL in both:

- Trigger: each **simplex facet** = facet with exactly `d` vertices
  (chcp: `l5_bf = [i : len(value_1[i])==5]`, value_1[i] = vertices on facet i; d=9 is
  the FACET count, dimension is 5; a (dim-1)-simplex facet has `dim` vertices).
- Edges: for simplex facet k, take ALL facets q sharing >=1 vertex with k (matrix entry
  != '0' and != 1), form edges {k,q}.
- Force each such edge to label **m=2 (pi/2, orthogonal)**: library `tl5_basis=[[2]]`,
  applied as a round-1 SAVER (keep only rows whose basis-edge column == 2).

Our `USE_L4_BASIS` (pipeline/stage4_blockpaste.py:189-203) computes the **byte-identical
edge set** — verified on d=4 type 6: both give
{(0,3),(0,5),(1,3),(1,5),(2,5),(3,4),(3,7),(5,6)}.

## When it is applied (CONFIRMED)

- d=4 chcp48: ALWAYS applied (loop bound hardcoded `-2` -> includes l4 AND l4_basis).
  No flag. Ma-Zheng's 348 census uses basis-ON.
- d=5 chcp59: applied under `flag` (loop bound `-flag-1`). flag=2 => includes l5_basis;
  flag=1 => skips it. README/default = flag=2 (basis YES). So d=5 census also basis-ON.

=> The condition is GLOBAL basis-ON for BOTH dimensions. It is NOT a per-facet or
dimension-dependent toggle. This CORRECTS the framing in WRITEUP 6.7 / commit 525a08a
("l_basis not a clean global ON/OFF ... needs Ma-Zheng's exact construction").

## The count semantics (CONFIRMED for the anchor)

chcp59 "data round 3: (N, ...)" -> N is the per-type census count. For P9_322 flag2 -> 3
= published census. Count is taken AFTER automorphism dedup (`per`, generated on the fly
via getper over S9) + 3 rounds of library killing/saving. So round-3 count == final
per-type census (combinatorial LSIE conditions are tight; no separate signature step
changes P9_322's 3).
chcp48 (d=4) is NOT self-contained: run48 hardcodes reading P8_17_per.txt and writing
output/P8_17/ -> only reproduces P8_17 as-is. chcp59 IS self-contained (generates per).

## The real discrepancy (localized)

- P9_322 anchor: OUR basis-ON = 3 = Ma-Zheng flag2. OUR basis-OFF = 18 = Ma-Zheng flag1.
  => our basis impl is FAITHFUL on the anchor.
- d=4 type 6 (P2, correct 49): our basis-OFF = 49 (=correct, basis is a no-op here since
  real P2 polytopes already satisfy it); our basis-ON = 12 (OVER-prunes). Edge set is
  identical to MZ, so the 12-vs-49 divergence is DOWNSTREAM of the edge set (label
  forcing / screening / dedup interaction) — a code bug on OUR side, NOT the condition.
  [NOTE: d=4 canonical run uses basis OFF (USE_L4_BASIS=(d>=5)) -> 338; the 12 only
   appears in scratchpad/test_lbasis.py which forces it ON for d=4.]

## RESOLUTION (2026-07-03, decisive — via stored polytopes, not the OOM batch)

The oracle batch (run Ma-Zheng chcp59 per intractable type on this machine) FAILED: the
intractable types (num=2, num=4, ...) OOM-kill after ~350-430s in Ma-Zheng's OWN round-1
pasting (empty captured output = SIGKILL). => Even Ma-Zheng's cheap combinatorial paste
blows up memory on ONE machine for these types; they used a cluster. So that oracle path
is compute-bound and abandoned.

Stronger, cheaper resolution from OUR stored CoxIter-verified polytopes:

**The basis (prism) condition is DIMENSION-DEPENDENT, and our code already implements it
correctly (`USE_L4_BASIS = (d>=5)`):**
- d=4 census P2 = 49 (basis-OFF). Checked the 49 (166 raw) stored polytopes' basis-edge
  labels: edges (3,7),(5,6) are pi/3 in 68/166 configs; (0,3),(0,5),(1,3),(1,5) pi/3 in
  21/166. Only 34/166 raw survive forcing basis->pi/2 (=>12 distinct). So REAL d=4
  polytopes violate basis-ON; the d=4 census is basis-OFF. Our basis-OFF=49 MATCHES.
- d=5 census P9_322 = 3 (basis-ON). basis-OFF gives 18 genuinely-distinct CoxIter-compact
  polytopes; census keeps 3 via basis. Our basis-ON=3 MATCHES.
=> basis is a genuine per-dimension SELECTION convention (for d=5 it picks 3 canonical of
   18 real; for d=4 it must NOT be applied). NOT gauge-dedup. Matches Ma-Zheng: chcp48
   (d=4) shipped in repo but NOT runnable here (missing Slis/Elis lib files); chcp59
   (d=5) self-contained, flag=2 default = basis-ON.

## CORRECTED CONCLUSION (supersedes WRITEUP 6.7 / commit 525a08a framing)

1. The prism/basis condition is DERIVED and our implementation is FAITHFUL: edge-set
   byte-identical to Ma-Zheng (static + runtime), and matches BOTH anchors
   (d=4 off->49, d=5 on->3). It is NOT an unresolved "needs Ma-Zheng's exact
   construction" blocker.
2. It IS dimension-dependent (d=4 OFF, d=5 ON); our `USE_L4_BASIS=(d>=5)` gate is right.
3. Therefore the d=5 12-vs-51 gap is NOT a wrong-basis problem. It is a per-type
   SOLVE/COMPUTE undercount on the ~11 intractable types (confirming d5-gap-diagnosis).
   Root cause: our per-type Gram-solve is far heavier than Ma-Zheng's paste+kill, and even
   that paste is cluster-scale for these types.

## STILL GENUINELY OPEN (the real remaining risk)
Is basis-ON uniformly correct across ALL d=5 types, or does it over-prune SOME d=5 types
the way it over-prunes ALL d=4 types? We only verified the P9_322 anchor. If basis-ON
over-prunes some d=5 types, a faster solver alone won't reach 51.
TEST NEEDED: get Ma-Zheng's actual 51 d=5 diagrams (paper / HCPdm output files) and check
whether every basis edge is pi/2 in them. If yes -> basis-ON uniform -> gap is pure
compute. If some basis edge != pi/2 -> basis is type-dependent within d=5.

## UNIFORMITY TEST RESULTS (2026-07-03, via injected chcp59 on tractable types)

Built scratchpad/mz_inject_run.py: injects OUR vertex_sets into chcp59 (bypasses num->line
lookup), runs flag2/flag1, reports round-3 count. Compared to our per-type basis-ON count.

- **basis-ON is FAITHFUL for d=5** (2nd + 3rd oracle anchors beyond P9_322):
  - tid0 (k=6): chcp59 flag2 = 5; OUR basis-ON = 5. MATCH. (flag1 OOM'd — basis makes it
    tractable.) chcp59 printed l5_basis(2)=[[7,...],[8,...]] = 2 simplex facets.
  - P9_322 (earlier): on=3=flag2. => our single-edge-pi/2 encoding reproduces Ma-Zheng's
    Lannér-set l5_basis on every tested type. basis-ON is NOT over-pruning d=5.
  - Confirmed on 59 "zero-paste" types: basis ON vs OFF give IDENTICAL paste counts, so
    basis is not the killer there.

- **The d=5 12-vs-51 undercount is in PASTE + SCREEN, not basis.** Our 109 types partition:
  6 found (=12 polytopes), 59 "zero-paste" (paste gives 0 candidates), 44 "screen-killed"
  (candidates exist, structured screen dec=False kills all). k-dist of found = {3:1,4:1,5:1,
  6:3}; only 3 types are k=6 (tids 0,5,6 -> 9). [WARNING earlier bug: len(paste_candidates())
  ==2 is the (cands,ordinary) TUPLE length, not candidate count -- always use cands,ord=... .]

- **But the zero-paste / screen-killed types appear to be LEGITIMATELY 0:**
  - tid11 (k=4, zero-paste): chcp59 flag2 = 0. OUR 0 is CORRECT.
  - tid38 (k=5, screen-killed): structured dec=False AND numerical floor res~4e-3 (no
    feasible rank-6 point) -- both agree infeasible. [chcp59 oracle: PENDING]
  If tid38 oracle = 0 too, then the ENTIRE gap lives in the 11 intractable types
  (compute), and our screen is SOUND for d=5. If tid38 oracle > 0, screen drops reals.

## ANSWER TO "is basis-ON uniform for d=5": YES.
basis-ON reproduces Ma-Zheng flag2 on every tractable type tested (P9_322->3, tid0->5),
does not over-prune (zero-paste types identical on/off), and Ma-Zheng use it uniformly
(global flag2 default). The census 51 IS defined with basis-ON. NOTE: basis-ON does drop
REAL compact polytopes as a CANONICAL SELECTION (P9_322: 18 real compact -> 3 kept) -- so
"51" counts prism-canonical representatives, not raw compact polytopes. Important framing
for d=6.

## RESOLVED (2026-07-03) — the 39 missing = prism-glued polytopes (Ma-Zheng per-type data)

Deep-research recovered Ma-Zheng per-type breakdown (arXiv:2203.16049v3 Table 13, Figs 6-9).
Only 6 of 109 types realize; each count splits basis(seed) -> final(census) via GLUING compact
simplicial 5-prisms onto orthogonal prism-end facets:
  P322 5->18(+13), P319 3->22(+19), P302 1->6(+5), P313 1->3(+2), P312 1->1, P284 1->1.
  basis total = 12, census = 51.
OUR pipeline reproduces the BASIS column EXACTLY: same 6 types, count multiset {5,3,1,1,1,1}=12
(our tids 0/6/5/19/29/63 -> 5/3/1/1/1/1). So the gap IS exactly the 39 prism-glued polytopes.
=> Missing step = PRISM-GLUING enumeration. l_basis ON forces prism-ends orthogonal (weight 2)
   = the SEED condition (12 seeds; campaign exhaustive+screen-sound confirms no more seeds). The
   39 have NON-orthogonal ends; Ma-Zheng generate them by gluing 5-prisms at orthogonal ends
   (l_basis OFF would brute-enumerate them = the branch that explodes). Gluing is the cheap
   structured route. Math/construction task, SHARED BY d=6 (P^B6 is prism-type).
CORRECTS earlier: "P9_322 basis-ON->3 = census" was a num-mislabel; P322 basis=5 (=our tid0),
   census=18. l_basis ON = seeds, NOT census. Generator k=2/3 split 50/34 confirmed correct
   (paper "54/30" was an extraction error).
NEXT: implement 5-prism gluing per basis seed at each orthogonal prism-end; validate finals
   {18,22,6,3,1,1}=51 via CoxIter. Ground truth: Table 13/16, Figs 6-9, Burcroff map Table 15.

## (superseded by RESOLVED above) CORRECTION — gap is upstream of the screen

Superseded the "compute wall" reading below. Verified since:
- Campaign was EXHAUSTIVE, not a timeout: runs/d5_n9/d5_campaign done with empty todo/wip,
  all 109 types COMPLETE; it SCREENED 11,967,720 candidates for tid1 (197,511,868 for tid16)
  to 0. So our completed answer is genuinely 12.
- Screen SOUND: 1500-candidate sample from tid1 all superhyperbolic (sig (2,7)x1499,(3,6)x1;
  min rank resid 1.07), zero feasible. Plus tid11 oracle 0, tid38/tid18 all superhyperbolic.
- Generator MATCHES Ma-Zheng: paper arXiv:2203.16049 Table 2 examines exactly 109 k>=2 types;
  we produce 109; k=4/5/6 counts (15/7/3) match EXACTLY. (k=2/3 split ours 50/34 vs paper
  ~54/30 — confirm.) So all 51 live in our 109.
- 51 proven independently (Ma-Zheng + Burcroff).
=> By elimination the 12-vs-51 gap is UPSTREAM of the (sound) screen: block-paste candidate
   GENERATION incompleteness or the exact SOLVE dropping realizers. A CORRECTNESS bug, not a
   scaling wall. The "port in-layer killing to scale" plan is therefore the WRONG fix.

REVISED NEXT STEP:
1. Get Ma-Zheng per-type realizing list (Section 7 of arXiv:2203.16049 / Burcroff thesis /
   HCPdm final data) -> which types hold the missing ~39. Prompt: research_prompt_d5_pertype_list.md
2. Paste-completeness audit: trace one known-real d=5 diagram (from a type we return 0 on)
   through paste->expand->screen->solve to find the drop stage. (tid0 paste already matches
   MZ count 5, so audit a MISSING type.)

## (superseded) FINAL SYNTHESIS (2026-07-03) — screen is sound; gap is compute in exploding types

Verified our screen is SOUND on every completable d=5 type tested (candidates it kills are
genuinely non-realizable). Feasibility = thorough multistart (scipy trf) search for dotted
weights w>1 giving a rank-6, signature-(5,1) Gram; done independent of our screen:
- tid11 (k=4, zero-paste): oracle chcp59 flag2 = 0; ours = 0. CORRECT.
- tid38 (k=5, screen-killed, 6 cands): all candidates best rank-resid ~1.4, best signature
  (2 neg,7 pos) = SUPERHYPERBOLIC. Correctly killed. [oracle still running/likely OOM]
- tid18 (k=5, screen-killed, 60 cands): ALL 60 rank-resid ~1.55, zero feasible. Correct.
So the 44 "screen-killed" + 59 "zero-paste" tractable types are LEGITIMATELY 0 where we can
complete them.

The missing ~39 census polytopes live in the EXPLODING-candidate types (the 11 campaign-
intractable ones): e.g. tid1 (k=5) blows the paste table to 12M+ candidates (>>solveable on
one machine). These are cluster-scale in BOTH pipelines (Ma-Zheng chcp59 OOMs on their
num=2, num=4 here too). We can't complete them -> our 0 there is a TIMEOUT, not a result.

Bottom line for the session goal:
- PRISM/BASIS CONDITION: fully derived, faithfully implemented, uniform & correct for d=5
  (matches oracle tid0=5, P9_322=3; not over-pruning). It is NOT the blocker. The last
  commit's "prism condition is an unresolved math blocker" is REFUTED.
- d=5 12-vs-51 GAP = pure COMPUTE on ~11 exploding-candidate types (re-confirms
  d5-gap-diagnosis with direct feasibility evidence that the OTHER types are correctly 0).
- To actually reach 51 (and unlock d=6): need a memory-bounded enumerator for the exploding
  types (Ma-Zheng's block-paste with their refined killing inserted EARLIER/in-layer, or a
  cluster). Not a math problem; an enumeration-scale problem.

## (superseded) earlier hypothesis via oracle batch

Run chcp59 flag2 (basis ON, canonical) + flag1 for the 11 "intractable" d=5 types
(nums 2,4,9,10,242,248,249,254,255,268,280) + P9_322 anchor. These 11 are exactly the
types our campaign got 0 on (base 12 came from the other ~98). 
map my_idx->num = {1:254,4:255,8:249,15:10,16:2,20:280,30:248,31:268,53:242,55:4,60:9}.

HYPOTHESIS: sum of flag2 over the 11 intractable ~= 39, so 12(ours,tractable) + 39 = 51.
That would prove: basis-ON is the correct global condition (no ambiguity); the d=5
12-vs-51 gap is purely OUR SOLVER/COMPUTE timing out on the 11 intractable types, NOT a
prism-condition problem. Ma-Zheng get finite counts on these because their combinatorial
pasting+killing is far cheaper than our per-type Gram solve.

If instead flag2 sum << 39, the condition story is more complex — revisit.

## NEXT AFTER BATCH
- Compare per-type flag2 vs our basis-ON counts.
- If hypothesis holds: reframe WRITEUP 6.7 (condition is settled global basis-ON; gap is
  compute on 11 types). Consider porting Ma-Zheng's cheap combinatorial paste+kill (no
  per-type Gram solve) as the enumerator — that's how they made d=5/d=6 tractable.
- Also fix our d=4 basis-ON over-prune bug (12 vs 49) since same code path feeds d=6.
