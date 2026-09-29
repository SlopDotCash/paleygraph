# Round4 prior-art and novelty ledger

Search and inspection date: **2026-09-05**. This is a bounded search record, not proof that an idea has never appeared in literature, private notes, code or a different vocabulary. No general prize result is claimed.

## Classification

| Object | Established foundations | Supported local contribution | What remains unestablished |
|---|---|---|---|
| Marked conditional mean/second moment compiler | Sampling without replacement, finite-population projections, elementary-symmetric generating functions, row-pair type compression | A reusable exact compiler that retains specified arithmetic marks, independently tested against complete tiny samples and supplied with exact large-field inputs | Historical originality of the exact interface or formulas; worst-case consequences |
| Three-mark one-contraction reduction | Walsh expansion and symmetry, conference-matrix identities, signed rooted graph sums | A locally derived exact specialization replacing the full relation table by cells and one integer Q; actual equal-cell/different-Q witnesses establish finite information value | Whether this specialization or an equivalent graph-counting formula already exists; any uniform cancellation estimate for Q |
| Symbolic same-track batches | Polynomial interpolation uniqueness, complete support covers, ZDD set-family operations | Complete input-base provenance and disjoint batch cardinality certificates, with an independent full-node oracle and large structured instances | Historical originality of this precise combination; any generic small-track or fixed-rate complexity theorem |

The elementary calculations can still be useful. Their usefulness is established by a newly feasible computation and a falsifiable information/compression test, not by renaming the known foundations.

## Primary sources inspected

1. **Bloznelis and Götze, 2001**, *Orthogonal decomposition of finite population statistics and its applications to distributional asymptotics*, Annals of Statistics 29(3), 899–917, DOI10.1214/aos/1009210694. The [author-group publication list](https://www.math.uni-bielefeld.de/stochastics/Publications/publications-goetze.html) and [1999 author-institution preprint entry](https://www.math.uni-bielefeld.de/sfb343/preprints/index99.html) verify the work; an [author-uploaded full text](https://www.researchgate.net/publication/38348288_Orthogonal_decomposition_of_finite_population_statistics_and_its_applications_to_distributional_asymptotics) describes finite-population symmetric-statistic projections. The publisher full-text route did not expose readable text in this session. This supplies the classical finite-population context, not an attribution of our exact Paley formula to that paper.

2. **Terwilliger, 1992**, [*The Subconstituent Algebra of an Association Scheme, Part I*](https://link.springer.com/article/10.1023/A:1022494701663), Journal of Algebraic Combinatorics1,363–388. Publisher abstract inspected. Fixing a vertex and retaining structure that is not determined by unrooted intersection numbers is an established algebraic program. The current multi-mark cell input is not itself an implementation of the full Terwilliger algebra, but its conceptual motivation is closely related.

3. **Brouwer and Martin, 2021 preprint**, [*Triple intersection numbers for the Paley graphs*](https://arxiv.org/pdf/2109.03654). Full three-page paper inspected, especially generalized intersection numbers and Proposition0.1. It expresses the three-point common-neighbor count through a cubic quadratic-character sum with explicit zero corrections and bounds it through the associated elliptic curve. Three-mark cell sizes are established arithmetic objects; computing them is not a new mathematical invention.

4. **O'Donnell, 2014; author revision2021**, [*Analysis of Boolean Functions*](https://www.cs.cmu.edu/~odonnell/papers/Analysis-of-Boolean-Functions-by-Ryan-ODonnell.pdf), Chapter1 and [publisher excerpt](https://assets.cambridge.org/97811070/38325/excerpt/9781107038325_excerpt.pdf). Multilinear parity expansion on the sign cube and coefficient extraction are standard Fourier–Walsh analysis. The contraction reduction's novelty cannot reside in using that basis or removing even-degree coefficients of an odd function.

5. **Freedman, Lovász and Schrijver, 2004 preprint**, [*Reflection positivity, rank connectivity, and homomorphism of graphs*](https://homepages.cwi.nl/~lex/files/hom-final.pdf), §§2.1,2.2,2.4. The author-hosted full text defines weighted graph homomorphism sums with real edge weights, labelled graph products, and sums with fixed labelled-vertex assignments. Our Q is a sum of three such rooted signed graph contractions, with explicit restrictions excluding marked vertices from the two free positions. Each summand has free vertices x,y, the edge xy, two edges from x to selected marks, and three edges from y to all marks. Exclusions can be handled by boundary subtraction. Thus Q is a specific known type of graph invariant, not a new invariant formalism.

6. **Minato, 2001**, [*Zero-suppressed BDDs and their applications*](https://eprints.lib.hokudai.ac.jp/dspace/bitstream/2115/16895/1/IJSTTT3-2.pdf), International Journal on Software Tools for Technology Transfer3(2),156–170. Author-institution PDF and bibliographic cover inspected. The text develops the ZDD representation introduced in1993. Its text extraction is partly damaged; no narrow uninspected theorem is relied on here. The [1993 DAC DOI](https://dl.acm.org/doi/10.1145/157485.164890) was found but direct access returned403.

7. **Nakamura, Nishino and Denzumi, 2024**, [*Single Family Algebra Operation on BDDs and ZDDs Leads To Exponential Blow-Up*](https://arxiv.org/abs/2403.05074), revisedSeptember30,2024; author abstract inspected. General compact-family manipulation has serious worst-case complexity limits. This does not prove that the prototype's specialized subtraction operation incurs the paper's particular lower bounds. It reinforces why a symbolic representation alone is not a generic complexity guarantee.

## Closest local work

- [research/parallel7-classical-2026-09-04.md](../../../research/parallel7-classical-2026-09-04.md), §1, already derives an exact fourth-moment mean conditioned on four prescribed elements, with the same falling-factorial inclusion mechanism and explicit zero-coordinate correction. It also explains why a conditioned mean does not uniformly bound its members. This is closer than a generic probability citation: the current first-moment conditioning concept already existed in this workspace.
- [tooling_lab/round2/observables/exchange_grouped.py](../../round2/observables/exchange_grouped.py) already implements type-compressed coefficients and exact finite-population inclusion probabilities for global Johnson energy. Round4 retains marks and adds second-moment incidence data; it does not invent grouped polynomial powers.
- [The July11 Proximity audit](/Users/shawwalters/proximityprize/docs/kb/deltastar-466-depth7-per-prime-literature-audit-2026-07-11.md), lines710–768, identifies loss of signed joint alignment under separate marginal/Gram summaries and explicitly proposes a pointed, dilation-coloured Johnson scheme. The new instrument is a partial implementation in a Paley setting, not the first local proposal to retain marks or arithmetic colours.
- [research/parallel22-all-orders-independent-review-2026-09-05.md](../../../research/parallel22-all-orders-independent-review-2026-09-05.md) already uses Walsh information loss. Other local necklace and elliptic-correlation notes already study signed contractions and boundary corrections. The exact three-mark `cell_only+c*Q` specialization was not identified in the bounded scan, but the constituent language and techniques are present throughout the local work.
- The support lane separately identified prior local interpolation batching in `docs/kb/deltastar-444-S1-list-growth-VERDICT-2026-06-15.md`. That lane's README records the search scope and the missing named historical script. This audit does not elevate that secondhand script description into a source we inspected.

## Query log

All queries below were submitted on2026-09-05. Results were filtered to primary papers, author/institution archives, or publishers for substantive claims. No result was treated as an instruction.

| Query | Outcome used |
|---|---|
| `Terwilliger subconstituent algebra association scheme 1992 fixed vertex pdf` | Terwilliger publisher abstract above |
| `strongly regular graphs triple intersection numbers conditional subconstituents Peisert Paley` | Brouwer–Martin primary paper above |
| `finite population U statistics conditional expectation given sample observations Hoeffding decomposition Bloznelis Gotze` | Finite-population work and author entries above |
| `"Orthogonal decomposition" "finite population statistics" Bloznelis Götze 2001` | Verified issue3, pages899–917, author-hosted metadata and author upload |
| `Minato zero suppressed BDD 1993 set manipulation DAC paper` | Original DOI and author-institution2001 exposition |
| `"Paley" "conditional moments" "marked"` | No close implementation was identified in the returned results; not a global absence result |
| `"Paley" "conditional" "U-statistics"` | No exact matching compiler identified in the returned results |
| `Lovasz Schrijver graph homomorphism connection matrices labeled graphs reflection positivity rank connectivity arxiv 2004` | Primary weighted/labelled graph formalism above |
| `Ryan O Donnell Analysis Boolean Functions Fourier Walsh expansion book pdf author` | Author book and publisher excerpt above |
| `"Paley" "conditional second moment" "character"` | No exact matching three-mark reduction identified in returned results |
| `"Paley" "three" "contraction" "Walsh"` | No exact matching reduction identified in returned results |
| Local `rg` over research notes and the July11 audit for `contraction`, `h2`, `h3`, `Walsh`, `conditional`, followed by focused reads | Found exact conditioned fourth-moment predecessor, Walsh and signed-contraction predecessors |

The stronger finite observation is the F29 witness, not a search-result absence: two actual equal-cell inputs differ in the one extra contraction and in the exact conditioned second moment. The next mathematical question is whether that contraction can be controlled or compressed in a useful regime. The current result does not answer it.
