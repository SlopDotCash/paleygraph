# Independent coefficient-structure and critical-scaling review

The general degree-six contraction coefficient and its critical-size scaling pass independent review. The coefficient is derived from an exact generating-polynomial recurrence, not fitted from numerical values. The resulting absolute correction bound concerns one conditional second-moment term; it does not imply a relative error against the baseline, a bound on individual target values, or a prize result.

## Why seven binomial terms suffice

Fix the three nonzero marked signs of each of two rows and their mutual sign. Increasing q by4 adds one occurrence of each outside paired-column type `(a,b)∈{−1,1}²`. The marked factors and the two zero-site factors remain unchanged. Therefore the entire marked generating polynomial is multiplied by

`B(t,w,v)=∏_(a,b=±1) [1+v(at+bw+abtw)]`.

Independent expansion gives

`B = 1−2v²(t²+w²+t²w²)+8v³t²w²`

`    +v⁴[(t²w²−t²−w²)²−4t²w²]`.

Thus `B−1` has minimum v-degree2. The desired `[t^6 w^6]` coefficient uses at most12 outside columns. Terms of v-degree above12 cannot contribute, and multiplication never lowers an exponent. With `h=(q−17)/4`, the exact finite binomial identity reduces to

`B^h = Σ_(j=0)^6 binomial(h,j)(B−1)^j  mod v^13`.

This holds for every integer q≥17 with q≡1 mod4. Starting at17 ensures that every hypothetical three-mark paired histogram used in the Walsh expansion has nonnegative multiplicities: the smallest bulk multiplicity is3, and at most three marked columns are subtracted. The argument is about the coefficient compiler; an actual graph interpretation additionally requires a realized conference matrix.

The common factor B is independent of the marked patterns and of the mutual row sign. Hence it also multiplies their signed Walsh combination W. The reviewed coefficient uses the character `u0*u1*w0*w1*w2` and the gap `K_+−K_−`, averaged over64 pattern pairs. This matches the coefficient of the surviving `(2,3)` contraction.

The review independently constructs W17 by multiplying **literal columns**, including the marked columns, without the producer's grouped multinomial arithmetic, canonicalization or helper functions. All seven resulting binomial coefficient rows agree exactly with the producer's output. In this specialization the rows j=3 through6 vanish identically after the required coefficient extraction, so every union coefficient `G_k(q)` is in fact at most quadratic in q. This is an exact polynomial calculation, not an empirical degree guess.

## Degree and factorization

Write `r=n−3`. Replacing the coefficient of v^k by the inclusion probability `(r)_k/(q−3)_k` gives

`c(q,n,6)=Σ_(k=0)^12 G_k(q) (r)_k/(q−3)_k`.

The denominator is nonzero throughout q≥17. Clearing the common denominator `(q−3)_12` yields a polynomial of degree at most12 in r. Independent symbolic division verifies the exact factor

`(r)_4 (q−r−3)_4`.

These eight endpoint factors correspond to `n=3,4,5,6` and `n=q−3,q−2,q−1,q`. The remainder after division is zero, and the quotient has degree4 in r over the rational-function field Q(q). The complete quotient and source bindings are saved in [review_results.json](review_results.json).

The generic degree is12. One formal specialization in the stated arithmetic domain is exceptional: at q21, the degree drops to11 because `G12(q)=2(q−21)(q−19)`. No existence of an order21 conference matrix is claimed. The fixed-q degree12 factorizations at29,49,101 and1297 all agree exactly with the general formula.

An independent second literal-column reconstruction at q29 agrees with every `G_k(29)`. Fourteen exact coefficient evaluations at q17 and q29 agree, including boundary zeros. In particular,

`c(29,7,6)=−2/7475`,

matching the earlier actual-matrix moment witness. The numeric checks supplement the recurrence proof; they are not what establishes its general validity.

## Uniform critical-size correction

Let alpha>0 be fixed, let q tend to infinity through admissible orders, and let `n/q^(1/4)→alpha`. Then `r=n−3∼alpha q^(1/4)`, so `r=o(sqrt(q))`. In the producer's factorization

`c=−(r)_4(q−r−3)_4 P(q,r) / [2(q−3)_12]`,

the quartic satisfies `P(q,r)∼2q^5`: every other monomial is smaller under `r=o(sqrt(q))`. Consequently

`c(q,n,6)∼−alpha^4 q^(−2)`.

The symbolic check substitutes `q=t^4`, `r=alpha*t−3`; numerator and denominator have t-degrees40 and48, with leading-coefficient ratio `−alpha^4`. The same leading estimate applies when the integer n differs by `o(q^(1/4))` from the canonical path.

For any symmetric balanced conference sign matrix, put h2 and h3 equal to their three-mark sign products outside the marks and zero on the marks. Let a,b,c be the three signs joining the marked vertices, and let `L=ab+ac+bc∈{3,−1}`. The identities

`sum(h2)=−3−L`,

`||h2||²=3q−15−2L`,

`||h3||²=q−3`

follow from the off-diagonal Gram identity and the pointwise relation `h2²=3+2h2` on unmarked rows. Since `S²=qI−J` and S is symmetric, its operator norm is sqrt(q). Cauchy–Schwarz gives the exact uniform bound

`Q²=(h2^T S h3)² ≤ q(q−3)(3q−13) < 3q³`.

Therefore, uniformly over realized conference matrices and choices of the three marks,

`|cQ| ≤ [sqrt(3)*alpha^4+o(1)] q^(−1/2)`.

The o(1) depends on the prescribed n(q) scaling, not on the marked triple or matrix. This bounds the **absolute additive correction** in the conditional second moment. A relative statement would need a separate lower bound on the baseline term. No pointwise control of T6 follows from this conditional-average bound alone.

The norm identities and squared bound were also checked directly for all4,334 marked triples of Paley17 and Paley29. Exact finite squared correction bounds at q65537, q1000033 and q6700417 are recorded in [critical_scaling_review.json](critical_scaling_review.json). These computations support the algebraic proof rather than replacing it.

## Reproduction and scope

The review uses the existing SymPy1.14 installation:

```sh
/opt/miniconda3/bin/python tooling_lab/round4/coefficient_structure_review/review_coefficient_structure.py
/opt/miniconda3/bin/python tooling_lab/round4/coefficient_structure_review/review_critical_scaling.py
```

The first script imports no polynomial arithmetic from the producer. It independently expands256 signed marked pair polynomials across q17 and q29, reconstructs the general union coefficients, and checks the factorization. Its saved run took57.4 seconds. The second script binds to that reviewed formula and checks the critical scaling and uniform norm calculation.

This is an independent program and mathematical review, not human refereeing or Lean formalization. Earlier sources and results were left unchanged. The derivation and checks do not establish that the coefficient identity is historically new.
