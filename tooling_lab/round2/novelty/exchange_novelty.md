# Grouped exchange calculation: bounded correctness and novelty audit

Date: 2026-09-05. Verdict: the degree-six calculation is supported; its
foundations are established, and originality of the exact evaluated
energy shares remains unestablished. Call this a specialized exact
evaluation and a new local implementation, not an invented decomposition.

## Correctness and the repaired API issue

`review_exchange_grouped.py` independently verifies 22 character-row type
histograms, 550 coefficient extractions by multiplying one coordinate at
a time, and all 448 distinct row-pair histograms over F5, F13 and F17.
It also verifies 324 falling-factorial containment probabilities by
enumerating 35,033 ordered subset pairs. Evidence and source hashes are
in `exchange_grouped_review.json`.

For normalized distinct rows, remove x=0,1. The remaining character sums
of a, b and ab are each -1, and the number of coordinates is p-2. The
four linear equations give the advertised counts. For arbitrary distinct
row arguments y,z, affine normalization multiplies both signs by
chi(z-y). Raw histograms can consequently be sign-flipped. The requested
coefficient is unchanged: its total character degree is 2d, always even.
This justifies multiplying one normalized off-diagonal coefficient by
p(p-1), including when d is odd.

For fixed marked d-sets S,T with overlap r, let u=2d-r. Their r shared
points must enter the common region of a distance-j pair (C,D). Choose
a of the d-r first-only marked points and b of the d-r second-only ones
for the corresponding exclusive regions. The probability of that labeled
allocation is `(n-j)_(u-a-b)*(j)_a*(j)_b/(p)_u`. The binomial choices and
sum in the implementation therefore count every valid allocation once.
There is no missing factorial or unaccounted shared marked point.

One generic-degree assertion was incorrect: levels one and two need not
both vanish for odd d. It now applies only to even degrees. Three added
odd-degree comparisons match the earlier exact backend; at p=13,n=6,
level-two energy is 448/99 for d=3 and 64/77 for d=5. The grouped result
file was regenerated; all reported degree-six values are unchanged.

## Closest primary sources

- [Filmus, *Orthogonal basis for functions over a slice of the Boolean
  hypercube*, Section 4](https://arxiv.org/html/1406.0142): slice degree
  spaces, orthogonal norms, and the Johnson primitive idempotents are
  explicitly established. Initial preprint 2014; journal publication 2016.
- [Filmus–Mossel, *Harmonicity and invariance on slices of the Boolean
  cube* (2016)](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.CCC.2016.16):
  the harmonic representation and exchange-operator framework are known.
- [Bloznelis–Götze, *Orthogonal decomposition of finite population
  statistics* (2001)](https://doi.org/10.1214/aos/1009210694): the relevant
  sampling-without-replacement statistical framework predates this work.
  The authors' [1999 preprint record](https://www.math.uni-bielefeld.de/sfb343/preprints/index99.html)
  was verified; the publisher full-text and old preprint download did not
  render in this audit. No uninspected theorem from that paper is used.
- [Spielman, Paley graph lecture (2012), Lemma 5.7.1](https://www.cs.yale.edu/homes/spielman/561/2012/lect05-12.pdf):
  the exact two-row character correlation and balanced character rows
  underlying the six-type table are classical. The counts are order-two
  cyclotomic/intersection numbers.
- [Godsil, *Generalized Hamming schemes* (2010), Sections 3–4](https://arxiv.org/pdf/1011.1044):
  symmetrized tensor constructions and multivariate generating functions
  from finite type counts are established. This is close methodological
  prior art, not a claimed identification of the present target with a
  particular theorem there.

No inspected primary source gave the exact degree-six energy fractions
reported by this prototype. That bounded negative search does not make
them historically new. The calculation uses classical second-order
intersection data plus exact sampling probabilities and a known spectral
change of basis. It is plausible that an equivalent evaluation exists
under finite-population statistics, constant-weight codes, or association
schemes terminology; that possibility remains unresolved.

The algebra itself also limits what the result detects: any signed row
array with the same aggregated pair-type generating coefficients produces
the same overlap sums and degree-energy spectrum. Thus these all-input
L2 shares do not test arithmetic information absent from pair-type data,
or control rare exceptional sets. This is a direct consequence of the
implemented formula, not a claim from the cited literature.

## Prior local work and literal search log

The bounded local search covered Paley `research/` and
`/Users/shawwalters/proximityprize/docs/kb/`, besides current lab files.
Paley `parallel23-classical-independent-review-2026-09-05.md` already
derives exact Johnson-edge/Poincare identities; pass22 already discusses
Walsh information loss. More significantly, the proximity note
`deltastar-466-depth7-per-prime-literature-audit-2026-07-11.md`, lines
710–768, already contrasts Johnson-slice estimates with signed arithmetic
alignment, proves that forgetting a mark loses a phase-dependent
component, and proposes a pointed, dilation-coloured Johnson scheme.
This history rules out advertising the broad representation-theoretic
direction as new to this workspace. Its cited
[Filmus–O'Donnell–Wu paper](https://arxiv.org/abs/1809.03546) concerns
multislice log-Sobolev inequalities, not the current exact spectrum.

Web queries executed on the review date:

1. `quadratic character U statistics Hoeffding decomposition Johnson scheme Paley graph`
2. `quadratic character sums subsets Hoeffding decomposition variance Krawtchouk Eberlein`
3. `cyclotomic numbers order two Paley association scheme intersection numbers character sums`
4. `Filmus Mossel Harmonicity invariance slices Boolean cube arxiv Hoeffding`
5. `Bloznelis Gotze orthogonal decomposition finite population U statistics 2001`
6. `Paley graph character sums Hoeffding Johnson harmonic variance subsets`
7. `"quadratic character" "Hoeffding"`
8. `Godsil generalized Hamming schemes generating functions symmetric powers arxiv`
9. `"quadratic character" "U-statistics"`
10. `"Paley" "Hoeffding decomposition"`
11. `"Paley" "Eberlein"`

The narrow exact-phrase queries produced no matching primary calculation;
some results were irrelevant and were not used. Searches are not proof of
absence. Neither this review nor the resulting finite L2 identities prove
a Paley conjecture, a proximity-prize assertion, or global tool novelty.
