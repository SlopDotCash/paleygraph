# Exact pair-type twins with different sixth-kernel laws

This finite counterfactual isolates information discarded by the grouped all-input L2 instrument. Paley49 and Peisert49 have the same **aggregate inventory** of row-pair sign types, hence the same complete degree-six Johnson L2 spectrum, but different sixth-kernel distributions. Both graphs and their nonisomorphism are established constructions. The contribution here is the executable information-loss test and its exact certificates, not a new graph or a proof of the Paley conjecture.

The distinction is aggregation: the labeled table assigning a type histogram to each pair `(i,j)` is different. Its entries reveal the adjacency sign. The common object is the multiset of those histograms, with the same templates conditional on adjacency.

## Construction and validation

Work in `F7[X]/(X²+1)`, encoding `a+bX` by `a+7b`. The polynomial is irreducible because `−1` is a nonsquare in F7. The first primitive element found is `g=9=2+X`, of order48. Let `C_j={g^(4t+j):0≤t<12}`.

- Paley edges have nonzero differences in `C0∪C2`.
- Peisert edges have nonzero differences in `C0∪C1`.
- The exported trade removes294 undirected edges with differences in C2 and adds294 with differences in C1. This is an explicit cyclotomic trade, not a claimed Godsil–McKay switching sequence.

For each matrix S, put `S[x,y]=1` on edges, `−1` on nonedges, and0 on the diagonal. The tool checks every row and all2,401 ordered row pairs: zero row sums,24 positive and24 negative entries, symmetry, `S²=49I−J`, and strongly regular parameters `(49,24,11,12)`. The complete pair-type records and edge trade are exported. Across both graphs there are49 diagonal pairs,1,176 positive off-diagonal pairs and1,176 negative off-diagonal pairs, with matching type templates.

## What the equal second-order data fails to determine

For each six-element vertex set C, evaluate the exact integer

`E(C)=Σ_y ∏_{x∈C} S[y,x]`.

The compiled census visits all **13,983,816** six-sets in both graphs. It uses bit parity, with a separate direct-product validation path. For Paley this is the corresponding quadratic-character sixth kernel over F49; for Peisert the sign function comes from its two cyclotomic classes.

| Statistic over all six-sets | Paley49 | Peisert49 |
|---|---:|---:|
| Mean | −1/141 | −1/141 |
| Second moment | 2007/47 | 2007/47 |
| Third moment | 44555/35673 | 6145/3243 |
| Minimum | −29 | −29 |
| Number attaining minimum | 784 | 1176 |
| Maximum | 27 | 27 |
| Number attaining maximum | 3528 | 2352 |

The extrema are equal; the complete distributions are different. Their total variation distance is `793/11891`. Their third moments differ by `7680/11891`, so the first distributional moment distinguishing these twins is already the third. The exact degree0 through6 L2 energies agree individually, not just in their sum. They are saved for n=6,7,8.

The same arithmetic-labeled support

`C={0,1,X,6+5X,5+6X,6+6X}`

has values **19** and **−21**. This is the largest absolute pointwise difference,40, in the complete paired census. A pointwise difference alone would not distinguish relabelings; the different full histograms do. Row permutations preserve the sum defining E, while column permutations merely reindex its six-set arguments. Thus histogram inequality also excludes equivalence by independent row and column permutations.

For `T6(U)=Σ_{C⊆U,|C|=6}E(C)`, products `E(C)E(D)` depend, after summation over all C,D with fixed intersection size, only on the aggregate row-pair type inventory. This proves equality of all distance correlations and hence the whole Johnson L2 spectrum for every `6≤|U|≤43`; the displayed n values are direct exact evaluations. This does **not** show that pair-type information cannot yield any useful shared worst-case bound. It shows that this data does not determine the sixth-kernel law, tails or its arithmetic-labeled values. An instrument seeking those distinctions must preserve additional information.

## Reproduce and inspect

From the repository root:

```sh
python3 tooling_lab/round3/pair_type_twins/build_and_census.py
python3 tooling_lab/round3/pair_type_twins/controls.py
```

Requirements: Python standard library and a C++17 compiler (`clang++`). The Python driver imports the previously tested finite-field adapter and exact grouped coefficient/Johnson routines from rounds1 and2. The saved full census took2.66 seconds on this machine; timings are observations, not complexity guarantees. Enumerating all six-sets scales as `O(q^6)` and the bit implementation accepts at most63 vertices.

- [results.json](results.json): both matrices, cyclotomic classes, complete edge trade, exact spectra, marginal/joint histograms, moments and witnesses for every attained value.
- [pair_type_certificates.json](pair_type_certificates.json): every ordered pair's exact sign-type histogram.
- [control_results.json](control_results.json): exhaustive84-set GF9 Paley/Peisert and permutation controls, plus2,168 direct-product/bit-parity checks over GF49. The permutation control exhibits pointwise differences while retaining the same histogram.
- [Independent review](../novelty/twins_review.md): independently constructed GF49 arithmetic, all4,802 pair records, trade and witnesses. An ordinary-product census on1,712,304 six-sets containing0 reproduced **every marginal and joint bin**, using translation invariance to scale each count by49/6. It does not share the bit-parity census implementation.

GF49 is a square prime power; it is not a critical-size prime-field Paley instance. At q49 the nominal fourth-root size is2, below the degree-six diagnostic. This is a finite test of the information retained by an instrument. No asymptotic claim, prize certificate, graph-switching invention or universal novelty claim is made.

## Prior art and next tool boundary

The construction follows established Peisert graphs. Bhowmik and Barman give the cyclotomic definition and study clique counts in [their2022 paper](https://arxiv.org/abs/2205.03928). Sin studies the known Paley/Peisert spectral relationship in [The critical groups of the Peisert graphs](https://arxiv.org/abs/1606.00870). Brouwer and Van Maldeghem describe the two order49 conference graphs and known distinguishing clique data in [Strongly regular graphs, §10.18](https://homepages.cwi.nl/~aeb/math/srg/rk3/srgw.pdf). Thus even higher-order combinatorial differences between these graphs are established territory.

A bounded local search of existing Python and Markdown found a literature note about Paley/Peisert clique work, but no matching implemented pair-type-twin sixth-kernel census. That search is not proof of historical absence. The useful next interface is an explicit **information budget**: retain row triples or arithmetic labels, then test whether the added data distinguishes this certified pair. The third-moment discrepancy supplies a concrete acceptance test for such a refinement. It does not establish that a particular refinement would suffice for the full Paley problem.
