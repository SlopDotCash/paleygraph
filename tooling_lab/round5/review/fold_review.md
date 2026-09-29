# Independent review of trace folding and the seven-statistic adapter

2026-09-05. The reviewed symmetry reductions and coefficient formulas are sound. Both adapters match the independent complete-subset moment oracles. A validation gap in the new seven-statistic converter was identified and fixed before the final replay: attached records and attached power sums must agree, and family keys must be distinct Booleans.

Final reviewed hashes:

- `folded_moment.py`: `aa1cc9cf142e7b8f1e5630010d5c11c2d0051ff139dc5fc8c4c7858814dd32af`
- `residue_moment.py`: `8bd40e565abea4560c63b90c99de6449943172c7772c692cc57a469ac4d744a1`
- `verify_symmetry.py`: `b680c6efd25efeb6ca430488d677cae6436d607af569864c22a0154385505d1d`

[fold_review.py](fold_review.py) and [its exact evidence](fold_review.json) contain146 admissible conference-type inventories spanning all eight edge-sign patterns,1752 literal row-permutation/global-sign histogram checks, and24 literal column-polynomial symmetry checks that do not assume conference structure. The final adapters pass30 complete-subset comparisons for the19-evaluation interface and24 for the15-evaluation interface, including both GF49 twins at n6,7,8. Fifteen fold/unfold roundtrips recover the original normalized records exactly. The review also checks38 raw highest-degree identities,10,232 parameter-orbit instances, and eight malformed-input rejections. It replays full prime orbit censuses through1297 and samples the two largest primes plus every exceptional orbit; it does not claim another independent full convolution at those scales.

## Symmetry and information retained

For a distinct symmetric sign triangle, `theta=e01*e02*e12*tau` is unchanged by row permutation or global sign reversal: both the edge product and tau change sign under the latter. The monochromatic/mixed classification is also unchanged. For even d, the product of three row elementary coefficients is invariant under these operations, including before the union-inclusion functional.

Walsh inversion at canonical edges+++ gives theta≡15−q modulo8. At canonical edges++−, tau=−theta and integrality gives theta≡q−7 modulo8. These residues differ by4 for q≡1 modulo4. Therefore the **full theta histogram retains the family label**. It is incorrect to describe this complete scalar histogram as losing the triangle class.

Under the previously verified normalized-pair transitivity contract, the mixed histogram divides equally between the three normalized mixed edge patterns. Every mixed bin count is divisible by3; this supplies the exact inverse unfolding. The global factor remains q(q−1) per normalized parameter and is applied once. The S3 action is not implemented by dividing all orbit counts by6: harmonic parameters have orbit size3; equianharmonic parameters have size2 when p≡1 modulo3; all others have size6. The harmonic orbit is{−1,2,1/2}; the equianharmonic roots satisfy t²−t+1=0. Their collision in characteristic3 is outside the prime p≡1 modulo4 interface. No prime Möbius-action claim is transferred to Peisert; its literal row symmetries and separately checked semilinear normalization suffice.

These statements are conditional on the matrix and normalization contract. Powers-only inputs are trusted sufficient-statistic summaries, not independent certificates of graph realizability. The full-record paths now validate admissible support and attached moments. The current README states this distinction correctly.

## The highest coefficients and the smaller input budget

The analytic refinement from19 to15 evaluations is valid. Independently sign-gauge the three balanced rows so their triple product sum is theta. Even degree6 makes these row sign changes harmless for the target coefficient, although they need not preserve a surrounding conference matrix. For the formal column product, write `log F=G+theta L`.

Every theta-dependent monomial has positive odd exponent in each of the three target variables. Its lowest term is

```
t1*t2*t3*K(z),   K(z)=z−3z²+2z³.
```

At target(6,6,6), the theta6 coefficient can use only six lowest terms of L, so it is exactly `K(z)^6/720`. For theta5, the five factors must likewise use their lowest terms; any increase by two in a variable would exceed target6. The remaining(1,1,1) coefficient of exp(G) vanishes. Its singleton terms vanish by row balance, and G's own(1,1,1) term is the removed theta term. Every other decomposition of(1,1,1) includes a singleton. Consequently the theta5 coefficient is zero before inclusion, and each residual family polynomial is quartic.

Writing A for monochromatic powers and B for mixed powers, the known conference identities give A0=(q−5)/4, B0=3(q−1)/4, B1=3A1−6 and A2+B2=(q−3)(q+1). Hence the seven measured integers

```
A1,A2,A3,B3,A4,B4,A6+B6
```

suffice, together with q,n and the normalization contract. They are not the seven unlabelled powers theta0 through theta6. Five exact nodes determine each quartic remainder, a sixth checks the implementation, and three repeated-row templates give15 evaluations. The small-q fallback uses all admissible values and does not claim its interpolating polynomial outside that finite domain. The independent replay checks both exact raw leading coefficients at all19 union degrees and both families for q101, as well as the full tiny outputs.

## An exact identity explaining the n6/n7 transition

[fold_review_derivation.py](fold_review_derivation.py) derives the first class-difference polynomials, with [source-bound coefficients and checks](fold_review_derivation.json). Let Delta_k be the monochromatic-minus-mixed raw coefficient of `t1^6*t2^6*t3^6*z^k`, after canonical tau=theta and tau=−theta respectively. Then

```
Delta_6(q,theta)=0,

Delta_7(q,theta)=
 -q² +2q*theta² -32q*theta +44q
 -theta⁴/3 +(32/3)*theta³ -(134/3)*theta²
 +(1024/3)*theta -323.
```

The vanishing has a direct combinatorial explanation. Union size below6 cannot contain any selected6-set. At union size exactly6 all three selected6-sets must coincide, so their product is the sixth elementary coefficient of the pointwise triple product. That sign vector has q−3 nonzero entries and sum±theta; its degree-six Krawtchouk coefficient is even and independent of the triangle class. Thus the difference starts at union degree7.

After inclusion the entire family difference contains the factor

```
(n)_7/(q)_7,
```

with remaining expression `sum_(k>=7) Delta_k*(n−7)_(k−7)/(q−7)_(k−7)`, stopping at min(18,q). At n7 the difference is Delta7/binom(q,7); at n8 it is `8*Delta7/binom(q,7)+Delta8/binom(q,8)`. The stored Delta8 has theta4 coefficient16/3, while Delta7 has theta4 coefficient−1/3. The resulting theta4 coefficients are therefore

```
n7: -1/[3*binom(q,7)],
n8: (23−q)/[3*binom(q,8)].
```

The bivariate reconstruction has an explicit proof budget, rather than an empirical fit. Newton's recurrence expresses raw E_k through products of at most k sign power sums, each linear in q and tau. Hence total degree in(q,theta) is at most k; the separate logarithm argument gives theta degree at most6. The script evaluates enough independent q and theta nodes for these bounds, using160 nonnegative integral synthetic type inventories, and checks60 extra theta nodes and105 extra q evaluations. No graph realization is claimed for this interpolation grid. The underlying Newton backend already has a separate literal-column oracle review.

## What the rejected one-polynomial shortcut means

The distinct formal family polynomials at q101,n7 and n8 rule out one common degree-six polynomial on the full admissible synthetic-type domain. However, the actual Paley101 theta histogram has only seven support values, so one degree-six polynomial can interpolate that particular finite support. The full domain qualification matters.

The actual Paley1297 support has18 mixed and9 monochromatic theta values. A common polynomial of degree≤6 at n7 would have to equal the mixed polynomial on its18 points. It would then force the nonzero quartic difference above to vanish at all9 monochromatic points, which is impossible. The same argument applies at n8 because its theta4 coefficient is nonzero. This is a concrete obstruction to the one-polynomial interface on an actual parameter support. It is not a pair of actual graphs with matching ordinary moments and different third moments, nor a statement that theta itself loses the family.

The established ingredients remain row symmetry, Walsh inversion, finite-population moments, symmetric-polynomial identities and finite character weights. The useful new local specialization is the explicit retained-information budget and its executable adapters, including the exact union-degree transition and failure tests. No historical originality, uniform exceptional-set bound, or prize proof is established. The current source documentation keeps these limits visible.

Replay `python3 tooling_lab/round5/review/fold_review.py` and `/opt/miniconda3/bin/python tooling_lab/round5/review/fold_review_derivation.py` (the latter uses SymPy only to display exact polynomial expressions). Earlier third-moment and piece-discovery reviews remain unchanged.
