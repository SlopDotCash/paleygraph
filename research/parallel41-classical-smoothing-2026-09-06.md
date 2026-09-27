# Pass 41: a positive smoothed sixth-moment bound and its recovery cost

Date: 2026-09-06. This concerns arbitrary-set Paley character discrepancy,
including the Sidon single-size test from passes 20 and 40. It proves a
uniform bound for a short average of translates of the actual character
convolution. It does not prove the corresponding bound for the original
convolution, its original quartic-star aggregate, the subgroup target, or
the separate prize reduction.

The positive result is: if 2≤n and n⁴≤p, averaging
t=⌈n^(1/5)⌉ independent uniform translates gives expected sixth moment
at most (325/32)pn³. On n=⌊p^(1/4)⌋ this uses only
t=p^(1/20+o(1)) translates. The exact error identities show why the
usual input-norm L² tolerance does not recover an original-moment estimate:
the random smoothing removes the original pairing in expectation, while
Cauchy–Schwarz on its error retains the |A|n>p threshold.

## 1. Actual character inputs and the bounds used

Let p be an odd prime, V⊆F_p have size n, and

\[
F(x)=\sum_{v\in V}\chi(x-v),\qquad
M_j=\sum_x F(x)^j.
\]

Odd M_j are signed moments. All norms use counting measure. Exact character
orthogonality gives

\[
M_1=0,\qquad M_2=pn-n^2,\qquad \|F\|_\infty\le n.   \tag{1}
\]

Expanding the fourth and sixth moments and applying the individual Weil
bound to each nonsquare polynomial gives

\[
M_4\le3pn^2+3\sqrt p\,n^4,\qquad
M_6\le15pn^3+5\sqrt p\,n^6.                         \tag{2}
\]

For clarity, a word with all multiplicities even has at most p as its
complete character sum. Pairing positions bounds the numbers of such words
by 3n² and 15n³. Every other monic word polynomial is a nonsquare
with at most four or six distinct roots, giving the two Weil terms. The
standard bound used here is stated in
[Volostnov, arXiv:1712.09355v1, Theorem 5, PDF p.4](https://arxiv.org/pdf/1712.09355v1),
checked in the primary text. There is no simultaneous-sign or independence
assumption about those individual polynomial sums.

If n⁴≤p, these imply

\[
M_4\le6pn^2,\quad M_6\le pn^3(15+5n),\quad
M_3^2\le M_2M_4\le6p^2n^3.                         \tag{3}
\]

The last inequality is Cauchy–Schwarz on F and F². None of (1)–(3)
requires V to be Sidon.

## 2. Exact sixth moment of a short random smoothing

Choose s₁,…,s_t independently and uniformly from F_p, with replacement,
and put

\[
H(x)=\frac1t\sum_{i=1}^tF(x-s_i)=F*\mu(x),\qquad
\mu=\frac1t\sum_{i=1}^t\delta_{s_i}.                 \tag{4}
\]

The independence here is of our deliberately sampled shifts. For a fixed x,
the variables F(x−s_i) are independent samples from the actual full-field
value distribution of F, whose moments are M_j/p. This does not assert
independence of the summands χ(x−v) defining F.

Since M₁=0, the only surviving index multiplicity patterns in a sixth
power are 6, 4+2, 3+3, and 2+2+2. Their coefficients yield the exact identity

\[
\boxed{\begin{aligned}
\mathbb E\|H\|_6^6
={}&\frac{M_6}{t^5}
+\frac{15(t-1)}{t^5}\frac{M_4M_2}{p}
+\frac{10(t-1)}{t^5}\frac{M_3^2}{p}\\
&+\frac{15(t-1)(t-2)}{t^5}\frac{M_2^3}{p^2}.
\end{aligned}}                                      \tag{5}
\]

The M₃² term must be retained. The identity is also valid for t=1 and
t=2, with the appropriate vanishing factors.

Combining (3), M₂≤pn, and (5) gives

\[
\boxed{\mathbb E\|H\|_6^6
\le pn^3\frac{15(t^2+7t-7)+5n}{t^5}.}               \tag{6}
\]

Take t=⌈n^(1/5)⌉. For n≥2, t≥2 and n≤t⁵. Also

\[
t^2+7t-7\le\frac{11}{4}t^2,
\]

because the difference is (7/4)(t−2)². Therefore

\[
\boxed{\mathbb E\|H\|_6^6
\le\left(5+\frac{165}{4t^3}\right)pn^3
\le\frac{325}{32}pn^3.}                             \tag{7}
\]

In particular there exists a shift multiset satisfying (7). The first
constant is 5+o(1) as n grows. This is a uniform positive estimate for the
smoothed character convolution, including every Sidon V on the critical
slice. Its proof explicitly retains the original M₆/t⁵ term; it does
not presuppose the desired original sixth-moment estimate.

## 3. The same smoothing meets the usual L² input scale

Let G=F/n, K=H/n. Any probability smoothing is an L² contraction, so
every choice of shifts in (4) satisfies

\[
\|G-K\|_2\le2\|G\|_2
=2\sqrt{\frac{p-n}{n}}.
\]

Since ||χ||₂²=p−1, the particular smoothing furnished by (7)
simultaneously has

\[
\|G-K\|_2\le\epsilon_2\|\chi\|_2,
\qquad
\epsilon_2=2\sqrt{\frac{p-n}{n(p-1)}}=O(n^{-1/2}).    \tag{8}
\]

Thus the positive sixth-moment bound and an input-norm L² approximation at
the natural scale are compatible. There is no claimed obstruction to
standard Croot–Sisask normalization for this same convolution. The issue is
what accuracy (8) supplies after pairing or passage to L⁶.

For the random construction there are sharper exact mean identities:

\[
\mathbb EH=0,\qquad
\mathbb E\|H\|_2^2=\frac{M_2}{t},\qquad
\boxed{\mathbb E\|F-H\|_2^2=(1+1/t)M_2.}            \tag{9}
\]

Independence and zero mean prove the second identity; expanding the square
and using EH=0 proves the third. The error is small at the input-norm
scale after division by n, but is approximately the entire L² norm of F.

## 4. Pairing with an arbitrary set: the exact exponent loss

Let A⊆F_p have size m and S=Σ_(a∈A)F(a). Suppose a smoothing satisfies
||H||₆⁶≤Cpn³ and ||(F−H)/n||₂≤ε₂√(p−1).
Hölder on H and Cauchy–Schwarz on the error give

\[
\boxed{\frac{|S|}{mn}
\le C^{1/6}\left(\frac{p}{mn^3}\right)^{1/6}
+\epsilon_2\sqrt{\frac{p-1}{m}}.}                   \tag{10}
\]

For the positive construction, ε₂=O(n^(−1/2)). The error certificate
in (10) is therefore O(√(p/(mn))). It gives a saving only when
mn exceeds p by a power; it does not improve the classical L² range.
On n=p^(1/4+o(1)), m=p^(α+o(1)), the two displayed powers are

\[
p^{(1/4-\alpha)/6+o(1)}
\quad\hbox{and}\quad
p^{(3/4-\alpha)/2+o(1)}.                            \tag{11}
\]

The smoothed term would allow α>1/4, but the available error bound
requires α>3/4. To certify a saving p^(−δ) through the L² error
term alone requires ε₂≤p^(−(1−α)/2−δ+o(1)). The larger
standard tolerance cannot be silently substituted for that accuracy.

There is an exact paired version of the recovery gap. Let

\[
\mathcal E=\sum_z r_{A-V}(z)^2,
\qquad \mathcal V=\mathcal E-\frac{m^2n^2}{p}\ge0.
\]

The complete character correlation identity gives

\[
\mathbb E_s\left(\sum_aF(a-s)\right)^2=\mathcal V.
\]

Indeed, expanding in a−v turns the complete sum over s into
p·1_(a−v=a′−v′)−1. For P_H=Σ_aH(a), independence now gives

\[
\boxed{\mathbb EP_H=0,\quad
\mathbb EP_H^2=\mathcal V/t,\quad
\mathbb E(S-P_H)=S,\quad
\mathbb E(S-P_H)^2=S^2+\mathcal V/t.}                \tag{12}
\]

If V is Sidon, its nonzero differences have multiplicity at most one, hence
\(\mathcal E\le mn+m(m-1)\). This is a positive estimate for the
smoothed pairing's variance. It does not estimate the original S: the error
has exactly that original pairing as its expectation.

## 5. A finite-support obstruction to this particular L² certificate

For Sidon V and an arbitrary probability measure μ with q=μ(0),
[pass 40, equation (15)](/Users/shawwalters/Desktop/paleygraph/research/parallel40-classical-structure-2026-09-06.md)
proved the exact bound

\[
\|F-F*\mu\|_2\ge\sqrt{\frac pn}(n-1)(1-q).          \tag{13}
\]

For the equal-weight t-shift smoothing (4), either every shift is zero,
or 1−q≥1/t. In the latter case the *Cauchy–Schwarz error certificate*
appearing after normalized pairing satisfies

\[
\boxed{\frac{\|F-H\|_2}{n\sqrt m}
\ge\frac{n-1}{n}\sqrt{\frac{p}{mnt^2}}.}            \tag{14}
\]

This is not a lower bound on the actual paired error, which can cancel.
It proves that controlling that error solely by its full-field L² norm
cannot certify a nontrivial bound when mnt²=o(p).

For the short smoothing in (7), n=p^(1/4+o(1)) and
t=p^(1/20+o(1)), so (14) diverges whenever m=p^(α+o(1))
with α<13/20. No nonidentity equal-weight smoothing of that length
can escape this certificate obstruction by choosing its shifts more cleverly.
If all shifts are zero then H=F, and establishing (7) for that choice would
require a direct bound on the original moment; the averaging argument does
not provide it. The diffuse random choice has the stronger typical error
size described by (9) and (11).

## 6. Moving the L² error into L⁶ does not recover the missing power

For any positive smoothing, D=(F−H)/n obeys ||D||∞≤2. If
||D||₂≤ε₂√(p−1), interpolation gives

\[
\|D\|_6^6\le\|D\|_\infty^4\|D\|_2^2
\le16\epsilon_2^2(p-1).                             \tag{15}
\]

To reach the target original scale ||F−H||₆⁶=O(pn³)
*through this inequality*, one would need ε₂=O(n^(−3/2)).
At the available ε₂=O(n^(−1/2)), (15) gives only
||F−H||₆⁶=O(pn⁵). On the critical slice its exponent is 9/4,
worse than the direct Weil bound O(p²) in (2). This is a limitation of
this interpolation estimate, not a lower bound on the true L⁶ error.

A genuine input-norm L⁶ estimate would have a different consequence. If

\[
\|(F-H)/n\|_6\le\epsilon_6(p-1)^{1/6},
\]

then (7) and Minkowski give

\[
\frac{\|F\|_6}{p^{1/6}\sqrt n}
\le C^{1/6}+\epsilon_6\sqrt n,                       \tag{16}
\]

with an inessential factor ((p−1)/p)^(1/6) omitted on the right.
Thus ε₆=O(n^(−1/2)) for a smoothing satisfying (7) would be
sufficient. This accurate L⁶ estimate is not supplied by (8).

For the specific random smoothing, Jensen makes the outstanding error
estimate especially explicit:

\[
\boxed{\mathbb E\|F-H\|_6^6\ge\|F\|_6^6,}          \tag{17}
\]

because EH(x)=0 for every x. Therefore proving the needed L⁶ error
bound in expectation for this construction already proves the original
uniform moment bound. No assertion is made that every other use of
almost-periodicity is impossible. The missing information is cancellation
in the error at the required dual or L⁶ scale, rather than the existence
of input-scale L² almost-periods.

## 7. Validation and precise scope

[The exact verifier](/Users/shawwalters/Desktop/paleygraph/experiments/parallel41_classical_smoothing.py)
enumerated 39,887 shift tuples in seven cases over p=17,29,97, with
t=1,2,3 as applicable. It checked (5), (6), (7), (9), (12), and
the Jensen inequality (17), keeping nonzero third moments and all zero
character entries. It uses exact integer arithmetic and rational fractions.
[The results](/Users/shawwalters/Desktop/paleygraph/results/parallel41_classical_smoothing_2026_09_06.json)
record all checks passing. The general results follow from the displayed
proofs, not the finite checks.

The positive theorem (7) is about H=χ*(1_V*μ), a weighted convolution.
It cannot be inserted as if it were M₆(V) into the original quartic-star
identity or the single-size sampling transfer. Equations (10)–(17) locate
the recovery gap for two concrete attempts: L² pairing and L²-to-L⁶
interpolation. This pass establishes no new uniform bound for the original
Sidon positive sixth aggregate or arbitrary-set Paley discrepancy.
