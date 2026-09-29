# Exact identities for two-anchor necklace sums

**Status: neither Paley conjecture is proved.** This note evaluates an
infinite subfamily of the character sums arising in spectra of localized
Paley graphs. It concerns the classical quadratic-character problem, not
the thin dyadic additive-period bound or the prize benchmark.

For singleton labels drawn from two distinct anchors, every necklace in
which one anchor occurs once or twice reduces exactly to monochromatic
necklaces. In particular, for cycle length k and two occurrences,

\[
 |N_{k,\ell}|\le k^2p^{k/2}\qquad(p\equiv1\pmod4).
 \tag{1}
\]

This is a uniform theorem in k and p, using an existing Kloosterman-sum
estimate. It does not control all the mixed necklaces required for a
localization theorem. No novelty or Lean-formalization claim is made.

Subsequent progress in [planar-necklace-reductions.md](planar-necklace-reductions.md)
evaluates arbitrary two-block words and all binary singleton words of length
six. In particular, the two named three-occurrence examples below now have
uniform exact formulas. The general three-occurrence problem at larger
lengths and the complete localization problem remain open.

## Source boundary

[Kunisky, arXiv:2303.16475v1](https://arxiv.org/html/2303.16475v1),
Definition 1.13 and Conjecture 1.14, uses cyclic character sums with
nonempty labels. Its degree-one case and Appendix B.1 treat constant
labels; the full degree-two estimates and the spectral-edge conclusion
require more. Sections 1.4 and 7 separate weak spectral convergence from
control of extreme eigenvalues. We use that framework, without adopting
those open conclusions.

The analytic input here is
[Lu–Zheng–Zheng, arXiv:1305.3405v3](https://arxiv.org/html/1305.3405v3),
Lemma 2.1, equation (2.3), with tensor exponents (1,1). Remark 2.4 gives
the invariant dimension R=1. Its twisted Kloosterman bound, including the
floor in its constant, yields (4) below. The relevant statement and proof
were read in the versioned HTML. The deep sheaf-theoretic theorem is an
imported literature result, not reproved here. Both HTML sources are
archived and hashed in the source manifest.

## Definitions and the monochromatic input

Let p be prime with p≡1 mod 4. All matrices below include every element
of F_p unless explicitly restricted. Set

\[
 S_{xy}=\chi(x-y),\quad D_b=\operatorname{diag}(\chi(x-b)),
 \quad D=D_0,\quad T_j=\operatorname{tr}(DS)^j\quad(j\ge1),
 \quad t_j=T_j/(p-1),\quad t_0=1.
 \tag{2}
\]

The convention t₀=1 is intentional: it is the normalized trace of the
identity on F_p*, not the full p-dimensional identity. Scaling every
coordinate by a nonzero element shows that t_j is an integer. Equivalently,

\[
 t_j=\sum_{u_1\cdots u_j=1}\prod_{i=1}^j\chi(1-u_i),
 \qquad u_i\in F_p^*.
 \tag{3}
\]

We have t₁=0 and the useful bound

\[
 \boxed{|t_j|\le(j-1)p^{(j-1)/2}\quad(j\ge1).}
 \tag{4}
\]

For completeness, write Kl_j(a) for the unnormalized j-variable
Kloosterman sum on the product fiber a. The quadratic Gauss expansion
of each factor in (3), followed by grouping the auxiliary variables by
their product, gives the exact identity

\[
 t_j=p^{-j/2}\sum_{a\ne0}\chi(a)|\operatorname{Kl}_j(a)|^2.
 \tag{5}
\]

In the cited lemma, the twisted sum in (5) is at most
`floor(j−1/j) p^((2j−1)/2)`. This floor is j−1 for j≥2, giving (4).
For j=1 the sum is zero directly. This derives the stated constant from
the original source, rather than using a finite moment fit.

## One occurrence of the second anchor

Choose distinct anchors a,b. At each of k vertices put one factor
χ(x_i−a), except at a single vertex put χ(x_i−b); multiply by the
cyclic edge factors `∏χ(x_i−x_(i+1))` and sum over F_p^k. Call this
number N_k^(1). Then

\[
 \boxed{N_k^{(1)}=-t_k.}
 \tag{6}
\]

Translate a to zero. For b≠0, substituting x_i=b y_i contributes
χ(b)^(2k)=1, so the sum is independent of b. Summing it over all b
gives zero, since the second anchor occurs once. The b=0 term is T_k,
and all other p−1 terms are N_k^(1). This proves (6).

## Two occurrences and the exact formula

Now use b at two vertices separated by ℓ edges and a at the other
k−2 vertices, where `1≤ℓ≤k−1`. Put

\[
 r=k-\ell,\qquad s=\min(\ell,r),\qquad d=|\ell-r|.
\]

The resulting necklace N_(k,ℓ) is independent of the two distinct
anchor values. Its exact value is

\[
 \boxed{
 N_{k,\ell}=p\,t_\ell t_r-p^s t_d-t_k
       +2(-1)^k\sum_{j=1}^{s-1}p^j.
 }
 \tag{7}
\]

The sum is empty when s=1. The t₀ term matters when the two arc lengths
are equal. Formula (7) includes k=2: it gives N_(2,1)=1−p, as required
for the monochromatic two-vertex cycle at b.

### Kernels with the zero coordinate retained

Define the symmetric matrix

\[
 K_j=S(DS)^{j-1}\quad(j\ge1).
\]

Homogeneity and the monochromatic trace imply

\[
 (K_j)_{xx}=\chi(x)t_j\quad(x\ne0),\qquad (K_j)_{00}=0.
 \tag{8}
\]

Indeed every entry scales by χ(u) under `(x,y)↦(ux,uy)`, and
`T_j=Σ_x χ(x)(K_j)_(xx)`. Also

\[
 (K_j)_{0x}=(-1)^{j-1}\chi(x)\quad(x\ne0).
 \tag{9}
\]

To verify (9) and the trace product needed below, let U be the restriction
of S to F_p*, let D_* be the restriction of D, and put C=D_*U. If 1
and w=(χ(x))_(x≠0) denote the two orthogonal vectors on F_p*, the exact
correlation identity `S²=pI−J` gives

\[
 CC^T=C^TC=pI-11^T-ww^T,\quad
 C1=-1,\quad Cw=-w,\quad D_*CD_*=C^T.
 \tag{10}
\]

Here `C1=−1` means the negative all-ones vector. For example,
`U1=−w` follows from S1=0, and `Uw=−1` follows by restricting the
column `S²e_0=p e_0−1`. These identities also show that C is normal.

The restriction of K_j is D_*C^j. Its row at zero is w^T times
C^(j−1), proving (9). On the span of 1,w, C acts as −I. On the
orthogonal complement, CC^T=pI. Therefore, reducing the product
`(C^T)^ℓ C^r` on these two invariant spaces gives

\[
 \operatorname{tr}(K_\ell K_r)
   =p^s(p-1)t_d+2(-1)^k(p-p^s).
 \tag{11}
\]

More explicitly, the nonzero-coordinate trace is
`p^s((p−1)t_d−2(−1)^d)+2(−1)^k`; the pairs with exactly one zero
coordinate contribute `2(p−1)(−1)^k`. Their sum is (11). This
derivation retains the two exceptional directions and every zero-coordinate
term; neither may be replaced by a generic square-root eigenvalue.

### Averaging the anchor

With a=0, the necklace with second anchor b is
`tr(D_b K_ℓ D_b K_r)`. The elementary two-point identity is

\[
 \sum_b\chi(x-b)\chi(y-b)=p\,1_{x=y}-1.
\]

Consequently,

\[
 T_k+(p-1)N_{k,\ell}
 =p\sum_x(K_\ell)_{xx}(K_r)_{xx}
     -\operatorname{tr}(K_\ell K_r)
 =p(p-1)t_\ell t_r-\operatorname{tr}(K_\ell K_r).
\]

Insert (11), divide by p−1, and use
`(p^s−p)/(p−1)=Σ_(j=1)^(s−1)p^j` to obtain (7). ∎

## Quantitative consequence and its scope

For k≥2, the four terms in (7), divided in magnitude by p^(k/2),
are bounded respectively by

\[
 (\ell-1)(r-1),\quad k-1,\quad k-1,\quad \tfrac12.
\]

The second bound uses t₀=1 if d=0 and (4) otherwise. The last uses
p≥5 and `Σ_(j=1)^(s−1)p^j≤p^s/(p−1)`. Their sum is at most
`(k−2)²/4+2k−3/2≤k²`, proving (1).

Thus, for this family,

\[
 p^{-(k/2+1)}|N_{k,\ell}|\le k^2/p.
 \tag{12}
\]

This tends to zero for fixed k, and even uniformly for k=o(√p).
Swapping anchor names also covers patterns where all but one or two
vertices use b. Hence every pattern using only the singleton labels
{a},{b} is covered through length five. Constant singleton patterns are
covered by (4). The full degree-two problem also includes the label
{a,b}; those words have not been covered by this assertion.

## The next averaging step exposes a weighted cubic sum

If b occurs at three marked vertices separated by positive arc lengths
ℓ₁,ℓ₂,ℓ₃ with total k, the same argument gives the exact relation

\[
 T_k+(p-1)N^{(3)}=
 \sum_{x,y,z} C_3(x,y,z)
 (K_{\ell_1})_{xy}(K_{\ell_2})_{yz}(K_{\ell_3})_{zx},
 \tag{13}
\]

where

\[
 C_3(x,y,z)=\sum_b\chi((x-b)(y-b)(z-b)).
\]

For distinct x,y,z, this is a signed elliptic-curve trace. For a repeated
pair, `C₃(x,x,z)=−χ(z−x)`; for three equal entries it is zero.
Unlike the two-point correlation, it is not determined just by equality
of indices. The singleton patterns 000111 and 010101 are the first
length-six cases with three occurrences of both anchors.

The direct absolute-value approach still gives no power saving in (13).
Indeed (11) implies

\[
 \|K_j\|_F^2=(p-3)p^j+2p\le p^{j+1}.
\]

Using `|C₃|≤2√p` and
`Σ_(x,y,z)|A_xy B_yz C_zx|≤||A||_F||B||_F||C||_F` bounds the sum
on the right of (13) by `2p^(k/2+2)`. After division by p−1 this is
only O(p^(k/2+1)), the unsaved scale of the necklace conjecture.
This absolute-value argument supplies no saving. The later planar reductions
establish cancellation for all binary words of length six, while the general
weighted sum at larger lengths still needs a bound.

Even proving all fixed-degree necklace estimates would still leave the
spectral-edge and growing-localization issues. This note supplies no
arbitrary-small-set Paley bound and no thin-subgroup estimate. It is a
partial evaluation along an alternate route, not a replacement target.

## Exact verification

[The script](../experiments/localized_necklace_identities.py) uses Python
integers throughout, including its NumPy object matrices. It constructs
S and D from the Legendre symbol and verifies the kernels, exceptional
directions, and both exact formulas. A separate multiplicative convolution
computes the t_j; small literal sums over all coordinate tuples check
the original necklace definition. The cases are p=5,13,17,29,41,61,97,
through lengths 12 or 10. No spectral floating-point approximation is used.

[Results](../results/localized_necklace_identities.json) retain the trace
sequences, representative values, check counts, and source hashes. The
three-occurrence values remain exact finite data; the later uniform formulas
come from graph identities, not a fit to these values. Their anchor-averaging identity is independently
checked in three small fields, including repeated-index contributions.
