# An integer polynomial for fourth-energy collisions

**Status: the Paley conjecture and uniform subgroup bound remain unproved.**
This note converts the fourth-energy excess into collisions among roots
of one integer polynomial. It gives a complete fixed-order classification
through n=32 and an exact identity over prime-power residue rings. The
general pointwise bound needed for growing n is not established. No
novelty or Lean-formalization claim is made.

## The concentration estimate checked first

Use the notation of [mixed-periods-and-shifted-energy.md](mixed-periods-and-shifted-energy.md):
`κ_s=C₀s`, `Σκ_s=n−1`, `D=Σκ_s²−(2n−3)≥0`, and
`T₀=Σ_s κ_s(κ_s−1)(κ_s−2)≤𝒳(H)`. Then

\[
 \sum_s\kappa_s(\kappa_s-2)=D-1,\qquad
 \sum_s\kappa_s(\kappa_s-2)^2=T_0-D+1.
\]

Weighted Cauchy–Schwarz therefore gives

\[
 (D-1)^2\le(n-1)(T_0-D+1).
 \tag{1}
\]

The matrix symmetries also give `𝒳≥3T₀−2(κ₀)₃`: sum the entries on
row zero, column zero, and the diagonal; only (0,0) is counted three
times. These exact constraints still incur a size-dependent loss.
The stronger literal assertion `D²≤𝒳` is false in the existing quartic
example p=6700417,n=64: D=12 and 𝒳=114. This does not refute an
unspecified constant C in `D²≤C𝒳`, and no such estimate is proved here.

## A polynomial over the integers

Fix a dyadic n≥4, put N=n/2−1, and let ζ=exp(2πi/n). Define

\[
 \lambda_j=(1+\zeta^j)^n\quad(1\le j\le N),\qquad
 P_n(Y)=\prod_{j=1}^{N}(Y-\lambda_j).
 \tag{2}
\]

The involution j↦−j preserves λ_j because ζ^(jn)=1. Consequently every
Galois automorphism permutes this set, which takes one representative
from each inverse pair except the roots ±1. The coefficients are rational
algebraic integers, hence integers. This also follows computationally from
the integral divisions in the Newton recurrence below.

The roots are real and distinct. Indeed,

\[
 \lambda_j=(-1)^j\bigl(2\cos(\pi j/n)\bigr)^n,
\]

whose absolute values are strictly decreasing for 1≤j≤N. They are
nonzero and all have absolute value below 2ⁿ. Thus the integer

\[
 A_n=\operatorname{disc}(P_n)P_n(2^n)
 \tag{3}
\]

is strictly positive. There is no zero determinant hidden in this
definition. For example

\[
 P_4(Y)=Y+4,\quad A_4=20,
\]

\[
 P_8(Y)=Y^3+120Y^2-2160Y-256,
 \quad A_8=2^{27}3^9\cdot5\cdot17^3\cdot41.
\]

No trigonometry is needed to construct P_n. Its power sums are

\[
 s_k=\sum_{j=1}^N\lambda_j^k
   =\frac12\left[n\sum_{a=0}^{k}\binom{nk}{na}-2^{nk}\right].
 \tag{4}
\]

Expand `(1+ζ^j)^(nk)` and sum over all n-th roots; only powers divisible
by n survive. Remove the root 1, whose contribution is 2^(nk), and the
root −1, whose contribution is zero, then divide by the inverse-pair
multiplicity two. Newton's identities
`e₀=1`, `k e_k=Σ_(i=1)^k (−1)^(i−1)e_(k−i)s_i` give the coefficients
of (2) with exact integer arithmetic.

## Reduction modulo a splitting prime

Let p≡1 modulo n be prime, and let h have order n in F_p*. Substituting
h for ζ in (2) gives the complete factorization of P_n modulo p. One
way to justify this substitution is to express each coefficient of the
product in Z[X]/(Φ_n(X)); its value at ζ is the corresponding integer,
so the polynomial identity is divisible by Φ_n. Here Φ_n=X^(n/2)+1.

The map u↦uⁿ on F_p* has kernel H=μ_n. Therefore `(1+h^a)^n` labels
the multiplicative coset of 1+h^a. Excluding h^(n/2)=−1, these labels
consist of the single value 2ⁿ from h⁰=1 and two copies of every root
of P_n. If e_b is the multiplicity of b in P_n modulo p, then

\[
 \kappa_b=2e_b+1_{b=2^n},\qquad
 \frac D4=\sum_b e_b(e_b-1)+e_{2^n}.
 \tag{5}
\]

This proves in particular that D is divisible by four. It also proves
the exact criterion

\[
 \boxed{E_2(\mu_n\subset\mathbb F_p)>3n^2-3n
        \quad\Longleftrightarrow\quad p\mid A_n.}
 \tag{6}
\]

Indeed the discriminant vanishes modulo p precisely when P_n has a
repeated root, and P_n(2ⁿ) vanishes precisely when the distinguished
value collides with a root. Both are exactly the positive contributions
in (5). Thus one complete factorization of A_n classifies all exceptional
primes p≡1 mod n, not only those below a chosen search cutoff.

## An exact valuation identity, including higher collisions

For r≥1, lift h uniquely to a root of Xⁿ−1 modulo pʳ, compatible as r
increases. Existence and uniqueness follow by Taylor expansion: the
derivative n h^(n−1) is invertible modulo p, so there is a unique
correction h↦h+t pʳ at the next step. Write H_r for the resulting order-n
subgroup of `(Z/pʳZ)*`. These are residue rings, not extension fields.

Let D_r be defined by the same label multiplicities modulo pʳ. The pair
sum argument is still valid over this ring: 1+u is a unit for every
u∈H_r except u=−1, whose sum is zero. Hence the additive energy over the
ring satisfies

\[
 E_2(H_r)=3n^2-3n+nD_r.
\]

The partitions refine as r increases, so D_r is nonincreasing and
nonnegative. Consistent lifted roots in Z_p give

\[
 \boxed{4v_p(A_n)=\sum_{r\ge1}D_r.}
 \tag{7}
\]

For the proof, the discriminant is the product of squared pairwise root
differences. Each unordered pair contributes twice its valuation, while
P_n(2ⁿ) contributes the valuations of the differences from 2ⁿ. A nonzero
difference has valuation equal to the number of precisions at which it
vanishes. Summing (5) over those precisions gives (7). All valuations
are finite because A_n≠0, so only finitely many D_r are nonzero.

Consequently `D₁≤4v_p(A_n)`, but equality need not hold. For n=16,p=17,
the exact sequence is `(D₁,D₂,D₃)=(196,12,0)` and `v₁₇(A₁₆)=52`.
Dropping the higher terms would incorrectly predict D₁=208.
The script checks the energy at each precision by direct pair sums,
independently of the polynomial multiplicities.

## Complete classification through order 32

The integer discriminants are computed by fraction-free elimination on
multiplication by P_n' in Z[Y]/P_n. Complete prime factorizations of A_n
are saved, with exact product checks. Trial division certifies the factors;
an unfactored remainder is never accepted as prime without verification.
The results for n=4,8,16 exactly match the independently archived
quadruple-norm enumeration.

| n | Number of exceptional primes p≡1 mod n | Largest such prime | Exceptions in n⁴/4≤p≤n⁴ |
|---:|---:|---:|:---:|
| 4 | 1 | 5 | none |
| 8 | 2 | 41 | none |
| 16 | 6 | 337 | none |
| 32 | 37 | 21523361 | none |

For n=32, all exceptional primes other than 21523361 are at most 194977;
the quartic window is [262144,1048576]. Thus, for every prime in that
window with p≡1 mod 32,

\[
 E_2(\mu_{32}\subset\mathbb F_p)=2976.
\]

This conclusion is certified by the complete factorization, rather than
the earlier exploratory prefix scan. Every exceptional prime in the
table also has its energy verified by direct pair counts and its full
valuation verified by the lifting identity. At the existing quartic
witnesses (p,n)=(6700417,64),(67403009,128), the values of v_p(A_n) are
3 and 6, with D₁=12 and 24 and D₂=0. These valuations can be computed
from the lifted roots without constructing or factoring the whole A_n.
The n=256 and n=512 circular examples have valuation zero.

## The unresolved quantitative issue

A bound `v_p(A_n)=O(n log n)` uniformly in the quartic window would
imply the necessary fourth-energy bound `E₂=O(n² log n)`. That valuation
bound is not proved, and is not asserted to be necessary for the target:
higher D_r can contribute to the valuation without increasing D₁.
Even the required fourth-energy bound would not settle logarithmic-depth
centered moments or the Paley conjecture.

The elementary height estimate does not suffice. Since every relevant
root difference and boundary difference has absolute value at most
2^(n+1),

\[
 A_n\le2^{(n+1)N^2},\qquad
 D_1\le\frac{4(n+1)N^2}{\log_2 p}.
\]

In the quartic window this upper-bound expression is at least the
trivial bound `(n−2)²` for every dyadic n≥16: use p≤n⁴ and
`n+1≥4log₂n`. Thus a sharper arithmetic argument is needed. The
polynomial criterion and small-order classification do not provide one.

Run `python3 experiments/kernel_discriminant.py`. Coefficients, complete
factorizations, every exceptional prime, direct energies, and lift data
are saved in `results/kernel_discriminant.json`. This is ordinary
mathematical reasoning with exact finite checks, not Lean formalization.

The subsequent [quadruple orbit argument](quadruple-orbits-and-cube.md)
proves a uniform power factorization of the odd part of A_n and isolates
the contribution from four distinct elements. It does not resolve the
quantitative issue above.
