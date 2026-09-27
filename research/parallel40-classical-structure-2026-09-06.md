# Pass 40: a Sidon obstruction to the additive-structure bridge

Date: 2026-09-06. Scope: classical discrepancy for **arbitrary** subsets of a
prime field. This note supplies an exact applicability obstruction, not a new
character-sum upper bound. It does not use the dyadic subgroup target or supply
the separate Reed–Solomon prize reduction.

The useful conclusion is stronger than “a Sidon set has large doubling.” The
Sidon test sets permitted by the existing full-Paley reduction admit neither
many translations inside an envelope with small expansion nor nonzero small
L² almost-periods relative to the norm of their actual Paley character
convolutions. The latter statement does not obstruct the larger error tolerance
in the usual Croot–Sisask normalization; section 5 makes this distinction
explicit. Small subsets and modest enlargements do not repair the separate
doubling obstruction. Thus the existing small-doubling hypotheses cannot simply
be imposed on the reduced test sets.

## 1. The precise remaining hypothesis

Let p be an odd prime, χ the quadratic character with χ(0)=0, and

\[
F_V(x)=\sum_{v\in V}\chi(x-v),\qquad
M_{2r}(V)=\sum_{x\in\mathbb F_p}|F_V(x)|^{2r}.
\]

A Sidon set means that equality of two sums of two elements, with repetition
allowed, forces equality of the unordered pairs. This is the B₂ convention in
the local reductions.

One sufficient hypothesis already isolated in
[pass 20](/Users/shawwalters/Desktop/paleygraph/research/parallel20-inversion-moments-2026-09-05.md)
and [pass 5](/Users/shawwalters/Desktop/paleygraph/research/parallel5-classical-2026-09-04.md)
is the following restricted single-size statement. There exist unbounded fixed
integers r_j≥3, constants C_j, and β_j≥0 tending to zero such that, for all
sufficiently large primes and **every Sidon set** V of size
k=⌊p^(1/(r_j+1))⌋,

\[
M_{2r_j}(V)\le C_jp^{1+\beta_j}k^{r_j}.                 \tag{1}
\]

The local inversion-and-sign-removal transfer gives the same bound for every
k-set with factor 3^(2r_j). Sampling k-subsets of an arbitrary B, followed by
Hölder, gives, when |B|≥k,

\[
\frac{|\sum_{a\in A,b\in B}\chi(a-b)|}{|A||B|}
\le 3\left(\frac{C_jp^{1+\beta_j}}{|A|k^{r_j}}\right)^{1/(2r_j)}. \tag{2}
\]

For |A|,|B|>p^ε, choose j with 1/(r_j+1)+β_j<ε/2. The fixed
constant can then be absorbed into a threshold, giving saving
p^(−ε/(8r_j)). This summarizes the existing reduction; it is not a new proof
of (1). The inversion transfer's explicit field-size thresholds hold eventually
for these fixed r_j and h=2.

The cardinality restriction in (1) is essential. The forced-common-row
construction in pass 5 disproves unrestricted Gaussian moment bounds even for
sets with minimal lower additive energies. The almost-all sixth-moment bound
in [pass 23](/Users/shawwalters/Desktop/paleygraph/research/parallel23-classical-upper-2026-09-05.md)
also does not prove (1): it leaves exceptional Sidon sets, and the
[pass 24 review](/Users/shawwalters/Desktop/paleygraph/research/parallel24-classical-independent-review-2026-09-05.md)
does not remove them.

## 2. The structural theorem checked against this target

For a bridge through additive structure, the strongest directly relevant
variable-doubling result verified in this bounded lookup is
Schoen–Shkredov, *Character sums estimates and an application to a problem of
Balog*, [arXiv:2004.01885v1](https://arxiv.org/pdf/2004.01885v1), dated
4 April 2020. Publication metadata identifies Indiana Univ. Math. J. 71
(2022), 953–964, [DOI 10.1512/iumj.2022.71.8972](https://doi.org/10.1512/iumj.2022.71.8972).
The theorem statements below were checked in the arXiv primary text, not
assumed to be unchanged in an inaccessible published PDF. This is a comparison
of a specific structural method, not a claim to classify every current
character-sum theorem.

Writing |X+X|<K|X| and |Y+Y|<L|Y|, their Theorem 3 (PDF p.2) requires

\[
|X|>p^\delta,\quad |Y|>p^{1/3+\delta},\quad
|X||Y|^2>p^{1+\delta},\quad L\le p^{\delta/2},\quad
(\log K)^5\ll\delta^4\log p.
\]

It bounds the normalized binary character sum by
exp(−c(δ⁴ log p/(log K)²)^(1/3)). Theorem 14 (PDF p.9) has the
same bound and K condition, allows arbitrary Y, and instead requires
|X|>p^δ and |X|²|Y|³>p^(2+δ). The structural input, Lemma 13
(PDF pp.6–7), gives, for integers d≥2 and ℓ≥1,

\[
|Z|\ge\exp[-C\ell^3d^2(\log K)^2]|X|,
\qquad [d^\ell]\cdot Z\subseteq 2X-2X.              \tag{3}
\]

For fixed doubling constants, the complementary power-saving theorem of
[Volostnov, arXiv:1712.09355v1, Theorem 3, PDF p.3](https://arxiv.org/pdf/1712.09355v1)
applies when both sets exceed p^(1/3+δ) and both have bounded additive
doubling; it gives p^(−τ), with τ>0 depending on δ,K. The later
variable-K bound is valuable in a broader structural range, but its displayed
stretched-exponential saving is not a fixed power of p even for fixed K.

These are bounds for genuine binary character sums, not subgroup period
bounds. None asserts that arbitrary sets satisfy its doubling assumptions.

## 3. Exact obstruction for subsets and containers

Let V be Sidon, |V|=k≥2. Then

\[
|V+V|=\frac{k(k+1)}2,\qquad
r_{V-V}(t)\le1\quad(t\ne0),\qquad E^+(V)=2k^2-k.     \tag{4}
\]

For the difference assertion, two representations v₁−v₂=v₃−v₄≠0
give v₁+v₄=v₃+v₂. The Sidon property forces identical ordered
differences; the other matching would make the difference zero.

Every subset U⊆V is Sidon. If |U+U|≤K|U|, it follows that

\[
|U|\le2K-1.                                         \tag{5}
\]

Consequently no subset with polynomial size can have the subpolynomial
doubling required in the Schoen–Shkredov range. For fixed θ>0 and
k=p^(θ+o(1)), taking K≥(k+1)/2 makes (log K)^5 of order
(log p)^5, which violates (log K)^5≪δ⁴log p for any fixed
admissible δ. The lower bound for |Z| in (3) also drops below one
for every fixed d≥2, ℓ≥1 as k grows. It no longer guarantees a
nonzero amplifier.

Adding elements has an exact cost. If W contains V, |W|≤Dk, and
|W+W|≤K|W|, then

\[
KD\ge\frac{k+1}{2}.                                 \tag{6}
\]

Indeed V+V⊆W+W. More generally, for t=|V∩W|,

\[
\frac{t(t+1)}2\le K|W|.                             \tag{7}
\]

Thus a cover by m sets W_i, each of size at most Dk and doubling at most K,
requires

\[
k\le\sum_i|V\cap W_i|\le m\sqrt{2KDk},\qquad
m\ge\sqrt{\frac{k}{2KD}}.                           \tag{8}
\]

For k=p^(θ+o(1)), θ>0 fixed, the number of containers, their relative
sizes, and their doubling constants cannot all be p^o(1). This rules out a
subpolynomial-cost cover or enlargement as an automatic repair of the
structural hypotheses. It does not rule out a method able to pay a polynomial
loss and still prove an adequate character estimate.

There is also a cardinality obstruction before doubling is considered. On
the r≥3 sparse slice, k≤p^(1/4). When both variables have this size,
the Theorem 3 condition |Y|>p^(1/3+δ) fails, and the Theorem 14
condition |X|²|Y|³>p^(2+δ) fails because k⁵≤p^(5/4).
For an unbalanced pair a larger second variable might satisfy a size
condition; it still does not change (4)–(8) for the structured variable.

## 4. New exact bound on a translation envelope

For every nonempty T⊆F_p, put t=|T|. Then

\[
\boxed{|V+T|\ge\frac{k^2t}{k+t-1}.}                 \tag{9}
\]

Proof: the additive energy of V and T is

\[
E^+(V,T)=kt+\sum_{s\ne0}r_{V-V}(s)r_{T-T}(s)
\le kt+t(t-1).
\]

Cauchy–Schwarz applied to the representation function of V+T gives
(kt)²≤|V+T|E⁺(V,T), proving (9).

In particular, if |V+T|≤Lk for some L<k, then

\[
\boxed{t\le\frac{L(k-1)}{k-L}.}                     \tag{10}
\]

When k is polynomial in p and L=p^o(1), this is t≤(1+o(1))L.
An amplification that needs polynomially many distinct shifts while retaining
only subpolynomial enlargement of the summation set is therefore impossible
for these test sets. This is an obstruction to that geometric requirement,
not to all possible character-sum amplifications.

## 5. L² almost-periods: the normalization matters

All norms below use counting measure on F_p. Define Qf(x)=Σ_yχ(x−y)f(y).
The exact quadratic-character correlation identity is

\[
\sum_x\chi(x-u)\chi(x-v)=p\mathbf1_{u=v}-1,
\qquad
\|Qf\|_2^2=p\|f\|_2^2-\left|\sum_xf(x)\right|^2.    \tag{11}
\]

It holds for either congruence class of odd p modulo four. For u≠v,
the substitution t=(x-u)/(x-v), omitting x=v, transforms the sum to
Σ_(t≠1)χ(t)=−1; the diagonal case is p−1. Expanding the square
then proves the second identity.

Write ρ_V(s)=|V∩(V+s)|. Translation commutes with Q and the difference
of two translates of 1_V has mean zero. Hence

\[
\|F_V\|_2^2=pk-k^2,\qquad
\boxed{\|F_V(\cdot+s)-F_V\|_2^2=2p(k-\rho_V(s)).}   \tag{12}
\]

For every nonzero s, the second quantity is at least 2p(k−1).
Consequently

\[
\frac{\|F_V(\cdot+s)-F_V\|_2^2}{\|F_V\|_2^2}
\ge\frac{2p(k-1)}{k(p-k)}.                          \tag{13}
\]

For k→∞ with k=o(p), the lower bound tends to 2. Thus no nonzero
shift has an L² error tending to zero **relative to ||F_V||₂**. This
statement concerns that accuracy requirement, not every almost-periodicity
application to the same convolution.

In particular, the usual Croot–Sisask tolerance for the normalized convolution
G_V=χ*(1_V/k)=F_V/k is measured relative to ||χ||₂, rather than
||G_V||₂. Equation (12) gives, for every shift s, including zero,

\[
\|G_V(\cdot+s)-G_V\|_2^2
=\frac{2p(k-\rho_V(s))}{k^2}\le\frac{2p}{k}.         \tag{13a}
\]

Since ||χ||₂²=p−1, **all shifts** therefore satisfy
||G_V(·+s)−G_V||₂≤ε||χ||₂ whenever

\[
k\ge\frac{2p}{\epsilon^2(p-1)}.                     \tag{13b}
\]

With the slightly larger tolerance ε√p, the sufficient threshold is
exactly k≥2/ε². The factor p/(p−1) in (13b) retains the exact
χ(0)=0 convention. For any fixed ε>0, both thresholds hold eventually
on the sparse slice. Indeed
||G_V||₂²/||χ||₂²=(p−k)/(k(p−1)), so the input-norm tolerance
is larger by a factor asymptotic to √k than the same relative tolerance
measured against the convolution's norm.

Consequently (13) does not obstruct standard almost-periodicity at its usual
input-norm scale, even for this same convolution. Equations (6) and (10)
are separate obstructions to small-doubling containers and small translation
envelopes; their validity does not depend on this norm comparison.

There is a related exact obstruction to positive smoothing. Let μ be a
probability measure on F_p and q=μ(0). Since 0≤1_V*μ≤1 and both
functions have mass k,

\[
\|1_V-1_V*\mu\|_1
=2\left(k-\sum_s\mu(s)\rho_V(s)\right)
\ge2(k-1)(1-q).                                     \tag{14}
\]

On V the difference is nonnegative and its sum is at least (k−1)(1−q).
Cauchy–Schwarz on those k coordinates, followed by (11), yields

\[
\boxed{\|F_V-F_V*\mu\|_2^2
\ge \frac{p(k-1)^2}{k}(1-q)^2.}                     \tag{15}
\]

The application of (11) is exact because 1_V−1_V*μ has mean zero.
For smoothing error tending to zero relative to ||F_V||₂, (15) forces
q→1 on the sparse slice. It rules out claiming that diffuse positive
smoothing preserves the convolution at that relative accuracy. It does not
rule out averaging justified at the larger input-norm tolerance above.

## 6. Validation and the remaining mathematical step

The proofs of (4)–(15) are supplied above and do not import any moment
independence assumption. A separate bounded exact-integer calculation checked
all 538 Sidon sets of sizes 2–4 in F_7, F_11 and F_13. It checked
5,916 nonzero-translation identities and 42,628 set-pair/smoothing cases with
|T|=1 or 2, including (9), (14), (15), and the isometry in (11).
The finite checks are corroboration, not the basis of the general claims.
No builds, central edits, or persistent experiment files were produced.

The remaining task is to obtain cancellation for the Sidon slice (1), or
supply another transfer sufficient for arbitrary-set discrepancy. A proposed
bridge using small-doubling containers or small translation envelopes must
respect (6) and (10). A bridge additionally requiring L² error small relative
to the convolution itself must respect (13) and (15); standard Croot–Sisask
tolerances for the same convolution are not ruled out by those identities.
Merely supplying a stronger general theorem about sets with small doubling,
a low-additive-energy decomposition, or small subgroup periods does not
furnish the missing transfer. Neither this note nor the cited structural
theorems prove the remaining uniform estimate.
