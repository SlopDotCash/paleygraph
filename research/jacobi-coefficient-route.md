# A sufficient recurrence-coefficient bound

This is a conditional route to the subgroup square-root bound. The
implication below is proved; its uniform coefficient hypothesis is not.
It is not asserted to be easier than the original moment problem.

## Setup

For an even-order subgroup `H≤F_p*`, put `n=|H|`, `m=(p-1)/n≥3`,
and let `x_1,…,x_m` be its real nonprincipal periods, one per coset.
Use the unnormalized inner product

\[
\langle f,g\rangle=\sum_{i=1}^m f(x_i)g(x_i).
\]

The periods are distinct. Indeed, an equality between two coset sums
would give a rational polynomial vanishing at `ζ_p`, with constant
coefficient zero and degree at most `p-1`. The minimal polynomial
`1+X+⋯+X^(p-1)` forces that polynomial to vanish identically, so the
two cosets coincide. The minimal-polynomial fact follows, for example,
from Eisenstein applied after the substitution `X↦X+1`.

Consequently squared norms of nonzero polynomials of degree less than
`m` are positive. Define the monic orthogonal polynomials and their
three-term recurrence by

\[
P_0=1,\quad P_{-1}=0,\quad
P_{j+1}(X)=(X-\alpha_j)P_j(X)-\beta_jP_{j-1}(X),
\]

where `β_j=h_j/h_(j-1)>0` for `j≥1` and `h_j=⟨P_j,P_j⟩`.
At `j=0` the final term is zero. The recurrence follows by projecting
`XP_j` onto the orthogonal basis: coefficients below `P_(j-1)` vanish
because multiplication by X is symmetric for this inner product.

## The conditional theorem

Let `1≤d<m`. Suppose constants `A≥0`, `B>0` satisfy

\[
|\alpha_j|\le A\sqrt{n(j+1)}\quad(0\le j<d),\qquad
\beta_j\le Bnj\quad(1\le j\le d).
\tag{J}
\]

Then

\[
\boxed{M\le\left[A+\left(2+m^{1/(2d)}\right)\sqrt B\right]\sqrt{nd}.}
\tag{6}
\]

**Proof.** The zeros of `P_d` are the eigenvalues of the symmetric
tridiagonal matrix with diagonal `α_0,…,α_(d-1)` and off-diagonal
entries `sqrt(β_1),…,sqrt(β_(d-1))`. Expanding its characteristic
determinant gives exactly the displayed recurrence, proving the assertion.
For a unit vector, the absolute value of its quadratic form is at most
the largest absolute diagonal entry plus twice the largest off-diagonal
entry. Hence every zero lies in `[-L,L]`, where

\[
L=(A+2\sqrt B)\sqrt{nd}.
\]

Since `P_d` is monic, whenever `|t|>L`,
`|P_d(t)|≥(|t|-L)^d`. Also

\[
h_d=m\prod_{j=1}^d\beta_j
\le m(Bn)^d d!\le m(Bnd)^d.
\]

If any period had absolute value exceeding the right side of (6), its
single squared polynomial value would exceed `h_d`. This contradicts
`h_d=Σ_iP_d(x_i)²`. ∎

For `d=ceil(ln m)`, one has `m^(1/(2d))≤sqrt(e)` and `d≤2ln m`.
Then (6) becomes

\[
M\le\left[A+(2+\sqrt e)\sqrt B\right]\sqrt{2n\ln m}.
\]

If `n<d`, the trivial bound `M≤n` already gives `M≤sqrt(2n ln m)`.
Thus proving (J) with absolute constants throughout the specified
subgroup family would imply the requested square-root bound there.
The separate reduction to the sponsor's full MCA challenge is still needed.

## What the current data establish

For `p=6700417`, `n=64`, `d=12`, exact rational calculation gives

\[
\alpha_j^2\le64(j+1),\qquad \beta_j\le64j
\]

through all the required indices, so (J) holds there with `A=B=1`.
The additional finite inequality `|α_j|≤3(j+1)` also holds.
These are **finite checks**. They do
not establish (J) uniformly. The coarse bound (6) with these constants
is itself weaker than the trivial bound for this small example; the
sharp finite certificate uses the full individual polynomial coefficients.

The coefficients also explain why the earlier Gaussian fourth-moment
comparison fails without forcing an oversized second recurrence coefficient.
If `X` is a uniformly chosen coset period, then

\[
\mathbb EX=\alpha_0,\quad\operatorname{Var}(X)=\beta_1,\quad
\mathbb E(X-\alpha_0)^4
=\beta_1\left[\beta_1+\beta_2+(\alpha_1-\alpha_0)^2\right].
\tag{7}
\]

To verify the last formula, expand
`(X-α_0)²=P_2+(α_1-α_0)P_1+β_1P_0` and use orthogonality.
Here `β_2<2β_1` but the positive term involving `α_1-α_0` more than
compensates, yielding a central fourth moment greater than `3β_1²`.
This central comparison and the earlier raw coset-moment comparison use
slightly different normalizations because `α_0=-1/m`; both fail here.
Every assertion in this paragraph is checked exactly by
`experiments/period_polynomial_bounds.py`.

The coefficient condition (J) allows this skewness. Positivity of the
moment matrix alone supplies no uniform bound on its recurrence
coefficients, and the present argument proves no such bound.
The main conjecture remains open.

## Constant one fails in another quartic-window field

The finite observation `A=B=1` does not extend throughout the window.
Take `p=67403009`, `n=128`, and `H=⟨64701253⟩`. Exact trial division
certifies that p is prime, and modular exponentiation certifies the
generator's order 128. Here `n⁴/4≤p≤n⁴` and `m=526586`.

Label a nonzero multiplicative coset by the n-th power of any of its
members, and let `κ_c` count the elements of `(1+H)\{0}` in that coset.
The nonzero κ-values are one copy of 1, 57 copies of 2, and three
copies of 4. Thus `κ_H=0`, `Σκ_c²=277`, and

\[
E_2(H)=n^2+n\sum_c\kappa_c^2=51840.
\]

The signed moments of orders zero through four are
`m, -1, p-n, pκ_H-n², p(n+Σκ_c²)-n³`. They give

\[
\beta_2=\frac{76801470694776}{277291762225},\qquad
\beta_2-2n=\frac{5814779565176}{277291762225}>0.
\]

This disproves **B=1 in (J)** even at degree two inside the specified
window. It does not disprove (J) with another absolute B.

The extra energy has a concrete source. Every nondegenerate unordered
zero quadruple is an H-multiple of

\[
(1,1074964,1550267,64777777),
\]

whose entries sum to p and have generator exponents `(0,115,15,40)`.
There are exactly n such quadruples; each has four distinct entries and
24 orderings. Here nondegenerate means that no two entries are opposite.
The degenerate ordered count is `3n²-3n`, so the full count is
`3n²-3n+24n`. No zero triple exists, since `κ_H=0`.

The standard-library script
[`jacobi_low_order_obstruction.py`](../experiments/jacobi_low_order_obstruction.py)
checks energy by both pair-sum counts and complete zero-quadruple
enumeration. It checks β₂ independently by direct Gram–Schmidt, including
the Gram matrix. The exact output is
[`jacobi_low_order_obstruction.json`](../results/jacobi_low_order_obstruction.json).
This certificate has not been formalized in Lean.

## What the first two coefficient bounds would already require

Let X be a uniform coset period. Orthogonality and

\[
X^2=P_2+(\alpha_0+\alpha_1)P_1+(\alpha_0^2+\beta_1)P_0
\]

give the exact raw identity

\[
\mathbb E X^4=(\alpha_0^2+\beta_1)^2
+\beta_1(\alpha_0+\alpha_1)^2+\beta_1\beta_2.
\tag{8}
\]

If (J) holds through degree two, then `(1+√2)²<6` gives

\[
\mathbb E X^4\le(A^4+8A^2B+3B^2)n^2.
\]

Fourier orthogonality gives

\[
E_2(H)=\frac{n^4}{p}+\frac{p-1}{p}\mathbb E X^4.
\]

Thus uniform absolute A and B would imply `E₂(H)=O(n²+n⁴/p)`, and
in particular `O(n²)` in the quartic window. This is a stronger
low-moment requirement than the `O(n² log n)` bound obtained by combining
the requested spectral estimate with the second moment. No uniform proof
of this `O(n²)` energy estimate is supplied here.

One can weaken (J) to `|α_j|≤A√(nd)` and `β_j≤Bnd` over the same
indices. The proof of (6) still works, using `h_d≤m(Bnd)^d` directly.
Conversely, a bound `M≤C√(nd)` gives these weakened coefficient estimates
with `A=C`, `B=C²`: multiplication by X on the period measure has
operator norm M, so each diagonal or off-diagonal matrix entry has
absolute value at most M. Consequently this degree-scale variant is
equivalent to the desired estimate up to constants at logarithmic d.
It is a reformulation, not a new proved bound.
