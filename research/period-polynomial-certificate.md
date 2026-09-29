# A certified maximum for the resonant subgroup

For `p=6700417` and `H=⟨2⟩` of order 64, the exact maximum is attained
precisely at the frequency coset `H`:

\[
\boxed{M=\eta_1,\qquad
43.802482797626304198\le M\le43.802482797626304199.}
\]

Every other frequency coset has absolute period **strictly below 42**.
This sharpens the previous `43<M≤sqrt(1970)` result for this one field.
It does not prove the asymptotic Paley bound or resolve either prize
challenge. The proof below combines ordinary mathematics with exact
integer/rational certificates; the new computation is not Lean formalized.
No novelty claim is made for the polynomial method.

## All signed moments, independently reconstructed

Let `m=(p-1)/64=104694` and list the real periods as `x_1,…,x_m`, one
per multiplicative coset. Write `C_k(0)` for the number of ordered k-term
sums from `H` equal to zero. Orthogonality gives

\[
s_k:=\sum_{j=1}^m x_j^k=\frac{pC_k(0)-64^k}{64},\qquad 0\le k\le24.
\tag{1}
\]

In particular `s_0=m`, `s_1=-1`, `s_2=p-64`, and
`s_3=3p-64²`. The principal frequency has been removed in every moment.
Odd moments retain signs; they are not moments of absolute values.

The earlier carry computation already supplied all `C_k(0)` through
order 24. The new script independently checks every one of these values
using quotient convolution. If `f_r=1_H^{*r}`, symmetry under negation gives

\[
C_{2r}(0)=\sum_y f_r(y)^2,\qquad
C_{2r-1}(0)=\sum_y f_{r-1}(y)f_r(y).
\]

The functions are constant on nonzero multiplicative cosets, each of
size 64. Thus these inner products can be evaluated with 104695 entries
including zero. Only `f_0,…,f_12` are needed. Their convolution steps use
int64 arithmetic after checking an upper bound before every step; all
products and sums forming moments use arbitrary-size Python integers.
Total mass and agreement with the independently generated carry counts
are checked. The source and input hashes are saved with the results.

## The elementary certificate principle

For any real polynomial `q(X)=Σ_i c_i X^i` of degree at most 12, (1)
computes its exact squared norm over the period multiset:

\[
N(q):=\sum_jq(x_j)^2=\sum_{i,k}c_i c_k s_{i+k}.
\tag{2}
\]

Two elementary consequences suffice:

1. If `q(t)²>N(q)` throughout an interval, no period lies in that interval.
2. If `q(t)≥u>0` on an interval and `N(q)<2u²`, at most one coset period
   lies in that interval.

Both follow because each summand in (2) is nonnegative. In the second
case two periods would contribute at least `2u²`.

The coefficients below were constructed using orthogonal polynomials;
verification needs only (2), polynomial translation, and positivity of
rational coefficients. There is no reliance on approximate optimization.

## Constructing the three rational polynomials

Let `P_0,…,P_12` be the monic polynomials orthogonal for the inner product
`⟨f,g⟩=Σ_j f(x_j)g(x_j)`, with positive squared norms `h_r`.
The script computes them by the three-term recurrence and then directly
checks every entry of their Gram matrix using (1).
For a rational center `a`, define

\[
q_a(X)=\sum_{r=0}^{12}\frac{P_r(a)P_r(X)}{h_r}.
\]

Orthogonality gives `N(q_a)=q_a(a)`. All coefficients are rational.
The output contains every polynomial, norm, and positivity certificate.

**Upper endpoint.** Set `a=4381/100=43.81`. Every coefficient of
`q_a(a+Y)` is nonnegative, its constant coefficient is positive, and

\[
q_a(a)^2>N(q_a).
\]

Hence no period is at least 43.81. These are exact rational comparisons.

**Lower endpoint.** Use `q_{-26}(-X)` as a polynomial in the reflected
variable `X`. Every coefficient after translation `X=26+Y` is
nonnegative and its constant square exceeds its norm, evaluated using
the reflected signed moments `(-1)^k s_k`. Hence no period is at most -26.

**Uniqueness above 42.** For `q_{43}`, the certificate proves

\[
q_{43}(t)>3/5\quad(42\le t\le43.81),\qquad
N(q_{43})<18/25=2(3/5)^2.
\]

The interval inequality is verified in the Bernstein basis. Substituting
`t=42+(181/100)z`, all thirteen Bernstein coefficients of
`q_{43}(t)-3/5` are positive rational numbers. The degree-twelve Bernstein
basis functions `binom(12,j) z^j(1-z)^(12-j)` are nonnegative on `[0,1]`
and sum to 1. This proves positivity over the entire interval, not merely
at sampled points. The certificate principle now allows at most one
coset with period at least 42.

## The period at frequency one

A direct rational enclosure proves that `η_1` lies in the displayed
interval, and in particular exceeds 43. This also provides an independent
ordinary-mathematics lower proof, in addition to the earlier Lean result.

The computation bounds π using Machin's identity

\[
\pi=16\arctan(1/5)-4\arctan(1/239).
\]

One elementary verification sets `a=arctan(1/5)`, `b=arctan(1/239)`.
The double-angle formula gives `tan(4a)=120/119`, hence
`tan(4a-b)=1`. The alternating arctangent bounds put `4a-b` between
0 and `4/5<π/2`, so `4a-b=π/4`.
Each arctangent is bounded using its alternating power series, with
odd truncation index below and even truncation index above. The resulting
π interval is rounded outward to denominator `10^30`.

For each `h∈H`, put `b_h=min(h,p-h)` and bound
`cos(2π b_h/p)`. Both rational bounding angles lie below π, as checked
explicitly; cosine is decreasing there. At the upper angle, the cosine
series through index 21 gives a lower bound; at the lower angle, the
series through index 20 gives an upper bound. The remaining terms alternate
and decrease in magnitude, so these are rigorous bounds. Sum all 64 terms
and round the endpoints outward to denominator `10^18`.

No floating-point approximation to π, cosine, or a period is used in
these acceptance checks. The raw rational bounds are retained in the
output as well as the shorter decimal enclosure.

## Why this identifies the global maximum

The subgroup coset `H` has period `η_1>43`, while at most one coset has
period at least 42. Therefore every other positive period is below 42.
All negative periods are greater than -26 by the lower endpoint certificate.
Their absolute values are also below 42. This proves `M=η_1` and identifies
all maximizing nonzero frequencies as precisely the members of `H`.

## Relation to a uniform proof

This computation uses all signed moments through degree 24, with their
arithmetic correlations intact. It demonstrates a stronger finite use of
the moment data than taking a root of the last even moment. It does not
bound these data uniformly across fields.

The [recurrence-coefficient criterion](jacobi-coefficient-route.md) gives
one sufficient uniform estimate that would imply the target bound.
The required uniform estimate remains unproved. In particular the finite
certificate here must not be substituted for that estimate.

Reproduce the full certificate with:

```sh
python3 experiments/period_polynomial_bounds.py
```

It uses NumPy for the bounded integer convolutions and Python's standard
library for exact rational algebra. Output is
`results/period_polynomial_bounds.json`.
