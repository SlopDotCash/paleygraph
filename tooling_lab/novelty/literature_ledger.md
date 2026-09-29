# Primary-source and search ledger

Search/retrieval date: **2026-09-05**. This is a bounded audit, not a systematic review or proof of absence. Paper dates below come from primary metadata where inspected. Most papers were reviewed at abstract/metadata level; CEGAR and Kemper's separating-invariants source also had full-text extraction available. No theorem from these papers is used to prove a Paley or prize result.

| ID | Primary source | What was verified and why it matters |
|---|---|---|
| L1 | Clarke et al., [Counterexample-guided Abstraction Refinement](https://web.stanford.edu/class/cs357/cegar.pdf) | The original methodology refines an abstraction using spurious counterexamples. This is direct prior art for the proposed research loop. |
| L2 | Kemper, [Separating Invariants](https://www.math.cit.tum.de/fileadmin/w00ccg/algebra/people/kemper/kemper.mega_hyper.pdf), source dated 2007-01-24 | Separating subsets are studied even for general function families; separating and generating are distinguished. A target-fiber auditor is not the invention of separating invariants. |
| L3 | Kemper–Lopatin–Reimers, [Separating invariants over finite fields](https://arxiv.org/abs/2011.07408), 2020-11-14 | Finite-field matrix-group separating invariants and degree bounds are explicit prior art. Our subset features need not be polynomial generators of such a ring. |
| L4 | Curto–Fialkow–Moeller, [The extremal truncated moment problem](https://arxiv.org/abs/math/0610882), 2006-10-28 | Moment positivity, variety constraints and consistency are distinct; low moments need not specify a unique representing measure. Actual-character realizability is an additional restriction. |
| L5 | Meka–Potechin–Wigderson, [Sum-of-squares lower bounds for planted clique](https://arxiv.org/abs/1503.06447), 2015-03-22 | Pseudoexpectations and moment-matching obstructions to low-degree algorithms are known. The random-graph result does not certify a deterministic Paley relaxation. |
| L6 | Kunisky–Yu, [A Degree 4 Sum-Of-Squares Lower Bound for the Clique Number of the Paley Graph](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.CCC.2023.30), published 2023-07-10 | The degree-four relaxation has value at least order p^(1/3). The paper also reports why random-graph pseudocalibration does not immediately transfer. This is a direct Paley tooling barrier. |
| L7 | Male, [Traffic distributions and independence](https://arxiv.org/abs/1111.4662), 2011-11-20; version 2018-03-08 | Traffic distributions enrich ordinary noncommutative moments; their independence rules have hypotheses on random families. Graph-indexed moment language is not new. |
| L8 | Ahn–Medarametla–Potechin, [Graph Matrices: Norm Bounds and Applications](https://arxiv.org/abs/1604.03423), initial 2016-04-12 | Shapes encode dependent random matrix entries; trace powers yield norm bounds. This is close to a graph-contraction compiler, but not a deterministic character-sum theorem. |
| L9 | Kunisky–Moore–Wein, [Tensor cumulants for statistical inference on invariant distributions](https://arxiv.org/abs/2404.18735), 2024-04-29 | Tensor-network invariants and cumulants form an explicit near-orthogonal polynomial basis for inference. This directly prevents claiming a new general tensor-cumulant calculus. |
| L10 | Collins–Gurau–Lionni, [Free cumulants and freeness for unitarily invariant random tensors](https://arxiv.org/abs/2410.00908), 2024-10-01 | Finite moment-cumulant transforms indexed by trace invariants exist. Also found: Buc-d'Alché–Lionni, [Properties of tensorial free cumulants](https://arxiv.org/abs/2605.01887), 2026-05-03, extending higher-order and product formulas. These hypotheses must not be silently imported into finite-field deterministic problems. |
| L11 | Hong et al., [Algorithm for computing μ-bases of univariate polynomials](https://arxiv.org/abs/1603.04813), initial 2016 | Algorithms compute polynomial syzygy module bases over arbitrary fields. Degree-filtered kernel computations are known tools. |
| L12 | Brakensiek–Gopi–Makam, [Generic Reed-Solomon Codes Achieve List-decoding Capacity](https://arxiv.org/abs/2206.05256), 2022-06-10 | Higher-order MDS controls intersections of column spans; generic RS achieves the relevant property. Generic is a material qualifier. |
| L13 | Brakensiek–Dhar–Gopi, [Generalized GM-MDS: Polynomial Codes are Higher Order MDS](https://arxiv.org/abs/2310.12888), 2023-10-19 | Polynomial and algebraic code generalizations expand the generic-rank toolkit. This does not discharge a prescribed special-field incidence pattern. |
| L14 | Gatermann–Parrilo, [Symmetry groups, semidefinite programs, and sums of squares](https://arxiv.org/abs/math/0211450), initial 2002 | Group representations reduce symmetric semidefinite problems. Symmetry-sector decomposition is established tooling. |
| L15 | Davies et al., [Advancing mathematics by guiding human intuition with AI](https://www.nature.com/articles/s41586-021-04086-x), 2021-12-01; Romera-Paredes et al., [Mathematical discoveries from program search with large language models](https://www.nature.com/articles/s41586-023-06924-6), 2023-12-14 online | Learning feature relationships and iterating executable proposals with evaluators are already published mathematical-discovery methods. Their presence supports evaluator-first tooling, not priority for the general loop. |
| L16 | Kunisky, [Spectral pseudorandomness and the road to improved clique number bounds for Paley graphs](https://arxiv.org/abs/2303.16475), 2023-03-29 | Character estimates are connected to weak spectral limits; extremal-edge convergence is separately conjectural. A finite bulk fit cannot certify the edge. |
| L17 | Xu, [Switching Graph Matrix Norm Bounds: From i.i.d. to Random Regular Graphs](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.CCC.2025.11), 2025 | A recent primary result transfers graph-matrix analysis beyond independent edges under specific random regular hypotheses. Constrained randomness has existing tooling; deterministic arithmetic transfer remains a separate task. |
| L18 | Diaconis–Sturmfels, [Algebraic Algorithms for Sampling from Conditional Distributions](https://www.math.ucdavis.edu/~deloera/MISC/LA-BIBLIO/trunk/Diaconis/Diaconis8.pdf), Annals of Statistics 26(1), 1998 | Integer moves preserving a linear sufficient statistic and Markov bases are established machinery. The prototype's primitive moment-preserving histogram move is a simple integer-kernel instance, not a new general mathematical invention. PDF retrieved directly on 2026-09-05; publication metadata confirmed on [Sturmfels's author page](https://math.berkeley.edu/~bernd/articles.html). |

## Literal search queries

Every query below was issued on 2026-09-05. Search hits from secondary indexes, Wikipedia, social media, and search-generated summaries were excluded as evidence for technical claims. A primary page or primary abstract is linked above for retained claims.

1. `Paley graph clique pseudorandomness sum of squares lower bounds planted clique Meka Potechin Wigderson`
2. `Reed Solomon codes higher order MDS GM MDS theorem matroid Brakensiek Gopi Makam`
3. `traffic probability cumulants graph polynomials random matrices Male`
4. `counterexample guided abstraction refinement automated mathematical conjecture discovery`
5. `site.arxiv.org traffic probability cumulants test graphs Male`
6. `site.arxiv.org Paley graph sum of squares pseudorandomness`
7. `site.nature.com Advancing mathematics guiding human intuition AI 2021 Davies knot representation theory`
8. `site.nature.com Mathematical discoveries from program search large language models FunSearch 2023`
9. `site.arxiv.org Gowers norms true complexity systems linear forms finite fields quadratic character`
10. `site.arxiv.org separating invariants finite groups Derksen Kemper separating invariants`
11. `site.arxiv.org truncated moment problem moments finitely atomic measures Curto Fialkow`
12. `site.arxiv.org "Paley" "cumulants"`
13. `site.arxiv.org "Reed-Solomon" "matroid"`
14. `"Paley" "cumulants"`
15. `"Reed-Solomon" "rank defect" matroid`
16. `"Paley" "separating invariants"`
17. `site.arxiv.org Kunisky Moore tensor cumulants finite free cumulants tensor`
18. `site.arxiv.org "Graph matrices: Norm bounds and applications"`
19. `site.arxiv.org Gatermann Parrilo symmetry groups semidefinite programs sums of squares`
20. `site.arxiv.org "mu-basis" "univariate" syzygy`
21. `"Diaconis" "Sturmfels" "Algebraic algorithms for sampling" 1998`

Follow-up direct retrieval: the Diaconis–Sturmfels paper L18 was opened from its primary-paper PDF URL supplied by the root implementation lane. This was a direct source read, not an additional search query.

The exact Paley/cumulant searches produced substantial irrelevant Paley–Wiener material. Their lack of a matching result is weak negative evidence. An author-page lead then revealed tensor cumulants by Kunisky–Moore–Wein, a substantially closer construction. This is an example of why a search with no literal phrase match cannot establish novelty.

The Gowers/true-complexity query found [Gowers–Wolf](https://arxiv.org/abs/0711.0185) and [general linear-form complexity](https://arxiv.org/abs/1403.7703). They were considered adjacent vocabulary, not used as a theorem or proposed sixth tool. Their fixed-system control should not be confused with uniform control over growing, arbitrary input families.

## Novelty dispositions

| Candidate | Known base | Proposed added combination | Disposition |
|---|---|---|---|
| Arithmetic target-fiber auditor | Separating invariants; moment problems; CEGAR | Actual-input realization, target leakage labels, finite target-range certificates, structural holdouts | Prototype-specific combination worth testing; global novelty unestablished |
| Constraint-preserving adversary | SOS/pseudoexpectations; CEGAR | Arithmetic-realization layer plus learned constraint selection | Known framework with a proposed domain interface |
| Signed contraction compiler | Traffic/graph matrices/tensor cumulants | Exact character zeros, anchor labels, signed collision corrections and conditioning ledger | Potential specialized algorithm; generic calculus already known and locally attempted |
| Degree/specialization dependency atlas | μ-bases; matroids; higher-order MDS | Syzygy-budget, characteristic and stack-realizability certificates in one record | Proposed tooling integration; bases and rank methods already known and locally attempted |
| Spectral observability audit | Representation reduction; Krylov spaces; Paley spectral route | Automated missing-sector and leakage certificates tied to arithmetic countermodels | Useful specialized diagnostic; no new spectral theorem |

## Limits and next audit steps

This pass did not exhaust MathSciNet/zbMATH, dissertations, code repositories, languages, citation chains, or private work. It did not prove theorems from abstracts or independently certify the historical local failure claims. A publication-level novelty claim would need exact algorithm pseudocode, comparison to its nearest mathematical equivalent, broader citation tracing, and independent expert review. A finite repository- or bounded-literature-relative absence claim can be checked; an unrestricted historical “provably never tried” claim cannot be established by this process.
