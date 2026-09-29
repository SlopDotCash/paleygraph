# Positive products and logarithmic-depth moments

**Status: the uniform subgroup estimate, classical Paley conjecture, and
full prize reduction remain unproved.** This note extends the two-coset
recurrence to the moment depth needed by the spectral target. It proves
a sufficient product-moment condition, gives its exact arithmetic
formulation, and identifies a principal-frequency obstruction to a
stronger uncentered approach. The product-moment bound itself remains
unproved uniformly. No novelty or Lean-formalization claim is made.

## The positive part is enough

Fix a prime p and a dyadic subgroup H_N⊆F_p*. For every subgroup H_s
in its tower, put

\[
 \eta_s(b)=\sum_{h\in H_s}e_p(bh),\qquad b\in\Omega=\mathbb F_p^*.
\]

All L^r norms below use the **uniform probability measure on Ω**. For
s≥2 the periods are real, since −1∈H_s. At a step H_(2k)=K∪gK, with
|K|=k, write

\[
 X(b)=\eta_k(b),\quad Y(b)=\eta_k(bg),\quad P_k(b)=X(b)Y(b).
\]

Then `η_(2k)=X+Y`, and X,Y have the same distribution on Ω. For any
real x,y, writing t_+=max(t,0),

\[
 (x+y)^2\le\max(x^2,y^2)+3(xy)_+.
 \tag{1}
\]

If xy≤0, the square of the sum is at most the larger square. If xy>0,
assume |x|≥|y| and use `y²≤xy` in `x²+2xy+y²`. Thus negative
products require no estimate in (1).

For a fixed real q≥1, let

\[
 A_s(q)=\|\eta_s\|_{2q}^2,
 \qquad a_k(q)=\|(P_k)_+\|_q.
\]

Minkowski's inequality and equality of the child distributions give

\[
 \boxed{A_{2k}(q)\le2^{1/q}A_k(q)+3a_k(q).}
 \tag{2}
\]

Indeed the q-th power of the maximum of the child squares is at most
the sum of their q-th powers. This proof neither assumes independent
periods nor expands a high power using an uncontrolled binomial bound.

## An explicit sufficient condition at the required depth

Fix q≥4 and C≥1. Suppose that at every step of the target tower,

\[
 \boxed{\|(\eta_k(b)\eta_k(bg))_+\|_q\le Cqk.}
 \tag{PM+}
\]

Then

\[
 \boxed{A_s(q)\le4Cqs\quad\text{at every level }s.}
 \tag{3}
\]

The base group of order two has `|η₂|≤2`. For the induction step,
`2^(1/q)≤2^(1/4)<5/4`, since `(5/4)^4>2`. Hence (2) gives

\[
 A_{2k}(q)<\frac54(4Cqk)+3Cqk=4Cq(2k).
\]

No assumption is needed at the initial steps with k≤q, because the
trivial pointwise bound gives `|(P_k)_+|≤k²≤qk`. Thus the unproved
condition only starts at child orders larger than q.

Let `m=(p−1)/N≥e²` and choose

\[
 q=\max\bigl(4,\,2\lceil(\log m)/2\rceil\bigr).
\]

This is even, `q≥log m`, and `q≤2log m`. Every nonprincipal period is
repeated N times on Ω, so

\[
 M_N^2\le m^{1/q}A_N(q)
 \le4eCqN\le8eCN\log m.
 \tag{4}
\]

Equivalently, with the notation of [subgroup-target.md](subgroup-target.md),

\[
 Q_q(H_N)\le m(4CqN)^q.
\]

Thus a uniform (PM+) at this one chosen depth, along every needed tower,
would give the intended square-root bound with logarithmic loss. This
goes beyond the fourth-energy conclusion of the preceding note. The
even integer may be slightly larger than the originally chosen moment
order; its comparable logarithmic size gives the same sufficient
spectral conclusion directly through (4).

The whole-tower hypothesis still includes smaller subgroups outside the
fixed quartic window. An endpoint estimate alone need not control those
smaller groups: two large child periods can cancel. No equivalence with
the endpoint conjecture or the prize is asserted.

For nonuniform bounds, iterating (2) gives the more flexible sufficient
estimate

\[
 A_N(q)\le4(N/2)^{1/q}
   +3\sum_{k=2,4,\ldots,N/2}
      \left(\frac{N}{2k}\right)^{1/q}a_k(q).
 \tag{5}
\]

The weighted sum on the right can replace separate bounds at each
level. It too is unbounded here at the required scale.

## Exact balanced counts for a stronger condition

An integer-count criterion is available if the positive part in (PM+)
is replaced by the absolute product. For an even integer q, let

\[
 Z_{q,q}(K,gK)=
 \#\{(a_1,\ldots,a_q,c_1,\ldots,c_q)\in K^q\times(gK)^q:
                         \textstyle\sum_i a_i+\sum_i c_i=0\}.
\]

Character orthogonality, with the frequency b=0 removed, gives

\[
 \boxed{\|P_k\|_q^q
 =\frac{pZ_{q,q}-k^{2q}}{p-1}.}
 \tag{6}
\]

Here P_k is real and q is even, so `P_k^q=|P_k|^q`. For odd q the
same count gives a signed moment and cannot be used as the absolute
moment in (6). In particular `E_Ω P_k=−k²/(p−1)` is negative.

The sufficient arithmetic condition is therefore

\[
 \boxed{pZ_{q,q}-k^{2q}\le(p-1)(Cqk)^q.}
 \tag{CM}
\]

It implies (PM+) and hence (4). Only the one even logarithmic depth is
required; a statement for every moment order is not being assumed.
For q=2, Z_(2,2) is exactly the mixed energy B in
[dyadic-descent-and-mixed-energy.md](dyadic-descent-and-mixed-energy.md).

The absolute-product condition is potentially stronger. For example,
in F_17 with |K|=8 and |H|=16, direct exact moments give
`E P_k=−4`, `E P_k²=16`. Zero variance proves `P_k=−4` at every
nonzero frequency. Its positive part vanishes, whereas its absolute
L^q norm is four. The script verifies all eight signed moments in this
example, as well as the variance certificate.

## The arithmetic correlation still needing a bound

Let r_q(x) count ordered q-term sums from K, and set
`f_q(x)=r_q(x)−k^q/p`. Scaling by g exchanges the child cosets and
g²∈K, so

\[
 Z_{q,q}=\sum_x r_q(x)r_q(gx),\qquad
 Z_{q,q}-k^{2q}/p=\sum_x f_q(x)f_q(gx).
 \tag{7}
\]

Thus (CM) asks for control of an off-diagonal correlation of the
centered q-step additive walk. The principal term is explicit, rather
than hidden inside an uncentered count.

Cauchy–Schwarz alone only gives

\[
 Z_{q,q}-k^{2q}/p\le E_q(K)-k^{2q}/p,
 \quad\text{or}\quad \|P_k\|_q\le A_k(q).
\]

Putting that estimate into (2) gives the factor `3+2^(1/q)` on A_k,
which exceeds the factor two needed for a bound proportional to the
subgroup size. The new recurrence does not make this elementary
substitution strong enough. The exact difference identity is

\[
 \sum_x f_q(x)f_q(gx)
 =\sum_x f_q(x)^2
    -\frac12\sum_x\bigl(r_q(x)-r_q(gx)\bigr)^2.
\]

Some additional arithmetic control of this correlation, or directly of
the positive product, is still required.

## Why complex-root decoupling cannot bound the raw finite count

In the complex dyadic cyclotomic fields, g has degree two over the
field generated by K. A sum from K plus a sum from gK vanishes only
when both sums vanish. Write T_q(k) for the ordered zero-sum count of
q complex k-th roots of unity. The corresponding balanced count is
exactly

\[
 Z_{q,q}^{\mathbb C}=T_q(k)^2,
 \qquad
 T_{2r}(k)=(2r)![z^r]
       \left(\sum_{j\ge0}\frac{z^j}{(j!)^2}\right)^{k/2}.
 \tag{8}
\]

For odd q, T_q(k)=0. These intrinsic relations remain valid after
reduction modulo a splitting prime, so `Z_(q,q)≥T_q(k)²`.
But finite-field relations cannot be restricted to these intrinsic
ones. Already at q=6,

\[
 T_6(k)=15k^3-45k^2+40k\le15k^3.
\]

Since q=6 is even, all Fourier contributions in (6), before removing
b=0, are nonnegative. Therefore, whenever p≤(2k)^4,

\[
 Z_{6,6}\ge\frac{k^{12}}p\ge\frac{k^8}{16}.
\]

For every dyadic k≥64, `k²>16·15²`, and hence

\[
 \boxed{Z_{6,6}>T_6(k)^2.}
 \tag{9}
\]

This holds at **every eligible prime** in that range, in particular
throughout the quartic windows with parent order N=2k≥128. It does
not require any assumption about the number of such primes.

At k=64, the complex balanced count is 14065500160000, while the
principal contribution alone gives at least 17592186044416 using the
largest allowed p. Thus at least 3526685884416 extra ordered balanced
relations are forced. Their presence is compatible with (CM): its
left side subtracts exactly that principal contribution.

More generally the ratio of this lower bound to `225k^6` is `k²/3600`.
This parameter-dependent obstruction prevents treating the raw count
as the complex count plus a negligible error. It is not a refutation
of (CM), (PM+), or the spectral conjecture.

## Exact finite verification

The script computes r_q on multiplicative K-cosets. Multiplication by
g gives a fixed-point-free involution on the nonzero cosets, and the
dot product with the permuted counts gives Z_(q,q). Every convolution
step has an explicit int64 overflow bound; all moment products and the
subtraction `pZ−k^(2q)` use arbitrary-precision integers.

Two small towers, in F_17 and F_97, are checked independently by dense
additive convolution at every level. Including the target computation,
there are 80 exact signed product moments across nine computed levels.
The distinction between odd signed moments and even absolute moments
is retained in the saved data.

For p=6700417, N=64, the number of nonprincipal cosets is 104694 and
the smallest even integer at least its logarithm is q=12. Rational
Taylor bounds for exp(10) and exp(12) verify that choice. The three
steps with child orders 2,4,8 satisfy (CM) with C=1 by the pointwise
bound `k²≤12k`. Exact computations verify it at the remaining child
orders 16 and 32. Thus the whole tower satisfies the stated critical-
depth hypothesis in this one field. The already archived whole-group
q-th energy independently agrees with the resulting coarse moment
bound. This finite result does not establish uniformity.

The counts also show where complex decoupling fails in that example:

| Child order k | q | Finite balanced count | Complex balanced count |
|---:|---:|---:|---:|
| 16 | 4 | 518400 | 518400 |
| 16 | 6 | 2562534400 | 2556313600 |
| 32 | 2 | 1024 | 1024 |
| 32 | 3 | 26112 | 0 |
| 32 | 4 | 11677184 | 8856576 |

The positive centered third moment in the final child split and the
nonzero extra relations are retained; no independence approximation
is used. The separate forcing argument (9) is checked by exact integer
inequalities and at the known quartic prime p=67403009,N=128, without
claiming to enumerate its large twelve-term count.

Run `python3 experiments/mixed_high_moments.py`; results, exact counts,
overflow bounds, and source/input hashes are in
`results/mixed_high_moments.json`. The argument is self-contained and
uses ordinary inequalities and character orthogonality; no new external
theorem dependency is claimed.

The open obligation now reaches the necessary depth: bound the positive
product moments, or the stronger centered balanced correlation (CM),
uniformly through the relevant split prime-field towers. Neither the
finite checks nor complex-field separation provides that estimate.

## Operator representation and a limit of unsigned estimates

The [signed quotient note](signed-quotient-operators.md) represents the
product values by a sparse symmetric integer matrix R_H with the principal
frequency absent. Thus (PM+) is exactly
`(tr((R_H)_+^q)/m)^(1/q)≤Cqk`, where `m=(p−1)/(2k)`. In the quartic
window, its entrywise absolute matrix has spectral radius at least
`k²−k/4`. An unsigned operator-norm comparison is therefore almost trivial
in that regime. Choosing different coset representatives does not generate
independent signs: it is diagonal conjugation and preserves every closed
walk product. These facts clarify the remaining need for arithmetic signed
cancellation; they do not prove (PM+).
