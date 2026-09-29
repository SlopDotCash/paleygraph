# Independent review of piece discovery

2026-09-05. No correctness defect was found in the bounded discovery algorithm or its successful list certificates. The reviewed core is `piece_discovery.py`, SHA256 `2c546e99c602906d9957b5d0c541b0d7ba15fd83359ac158f9d7bd703a1917a0`. The saved experiment output is SHA256 `e33f37635088e83a8897154e15fba4327e23195f0d4916d14d448b90cb54ed9a`.

[review_piece_discovery.py](review_piece_discovery.py) uses the frozen independent prime/GF9 arithmetic oracle, complete codeword enumeration, a full reduced Gaussian elimination distinct from production's forward solve, and direct polynomial multiplication. Its [evidence](review_piece_discovery.json) records:

- 928 inputs: all625 F5 words of length4 at k2,s3;135 F3 word/degree/threshold cases; all128 assignments of seven coordinates to two fixed affine GF9 pieces, plus40 arbitrary GF9 words.
- 1962 independent rank/consistency checks,930 exact dual inconsistency witnesses,1596 full polynomial branch identities, and805 complete factor-product checks.
- 793 certified static lists matched to all degree<k codewords;164 complete scalar-fiber comparisons,842 nodes, including140 whole-field cases.
- A complete codeword-support set-cover dynamic program verifies all928 exported cover-size lower bounds and805 claimed exact minimum cover sizes.
- 123 incomplete discoveries and12 verified covers without a strict-cap complete-list certificate are retained as failures of the requested next stage.
- Ten fresh export readbacks and five rejected corruptions: false minimum cover size, false symbolic node count, nonmonic relation, false maximal support, and false dual residual.

The oracle receives no planted labels from the production call. The GF9 exhaustive assignments are constructed by the review harness and then passed as words only. An accepted certificate is always compared against the entire tiny codebook, regardless of its hidden construction.

[review_piece_scale.py](review_piece_scale.py) independently checks both saved n1024,k64,F65537 outputs with direct field/polynomial arithmetic. It does not call the discovery solver, its factor-recovery routine, or the frozen list verifier. [Scale evidence](review_piece_scale.json) confirms:

| Discovered cover | Maximal supports | Exact minimum cover size | Static list at s410 | Complete scalar nodes |
|---|---|---:|---:|---:|
| Three polynomials |342,341,341|3|0|1|
| Two polynomials |512,512|2|2|131073|

For the large rank checks, the independent argument is stronger than another numerical solve. If a monic relation of Y-degree j≤r and weighted degree j(k−1) vanishes on the word, then substituting any one of the r discovered polynomials gives a univariate polynomial with more zeros than its degree. Each distinct `Y-h(X)` must divide the relation over F(X)[Y]. For j<r this is impossible; for j=r, monicity forces the exact product. The homogeneous ansatz has Y-degree<j, so it has zero kernel for every tested j≤r. This independently certifies the printed ranks64,191,381 as applicable, the lower-r inconsistencies, and the minimum-cover conclusions. The script also reconstructs all monic coefficients, checks2048 complete input equations, validates partitions and maximal supports, applies the strict root-count cap, and checks every symbolic scalar-fiber coefficient identity.

The mathematical separation is sound. A dual witness proves only that the specified monic affine system is inconsistent. A consistent selected solution need not be factorable, and the algorithm does not explore its nullspace. A Hensel jet is accepted only after the full polynomial identity holds. A covering set of recovered branches supplies a valid piece partition; completeness of the nearby codeword list follows separately from the strict root-count cap. Whole-field scalar behavior follows separately from the exact code anchor and linear transport. The output flags retain those distinctions.

Two independent negative controls expose actual limitations. Over F3, `(Y-X)(Y+X)(Y-1)` has three distinct polynomial branches but no simple split center anywhere in the base field. The routine correctly remains incomplete even with its full center budget. Over F5, the simple roots of `Y²-(1+X)` lift to local jets, but the k2 jets are not polynomial branches; the full-identity check correctly rejects them. Increasing the random splitting budget cannot repair the first obstruction.

The integer implementation is appropriate for its declared domains. NumPy forward-elimination products are bounded by `(p-1)²` under the explicit prime-size cutoff, and subtraction remains within int64. Certificate dot products use NumPy only when the entire accumulated sum is below2^63. GF9 tests use the generic exact-field route. This does not establish a scalable extension-field implementation at official prize parameters.

The current README's claims match the evidence: it states the label-free local interface improvement, preserves all16 saved hard-coset failures, describes the missing anchor and underdetermined relation spaces, and explicitly excludes the MCA nonjointness event. Its identifiability condition is sufficient rather than necessary. The review did not find a mathematical overclaim requiring a source edit.

This prototype fills the previously explicit local gap between supplied piece witnesses and verified scalar fibers. Its algebra is established polynomial reconstruction, finite-field factorization, Hensel lifting, root counting, and linear-code symmetry. The [prototype's primary/local source ledger](../piece_discovery/README.md) correctly identifies Guruswami–Sudan, Roth–Ruckenstein, Cantor–Zassenhaus, and already-local function-field/Hensel work. The useful contribution here is the implemented discovery-to-certificate interface and its measured failure boundaries; neither global historical novelty nor a prize proof is established.

Replay `python3 tooling_lab/round5/review/review_piece_discovery.py` and `python3 tooling_lab/round5/review/review_piece_scale.py`. The independent tiny suite took29.94s and the large algebraic replay3.84s in the recorded runs. All review artifacts stay in this directory; the earlier frozen rounds are unchanged.
