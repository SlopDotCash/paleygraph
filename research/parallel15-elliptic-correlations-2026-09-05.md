# Higher elliptic correlations and what their square-root scale misses

Status: a source-dependent bound for all translated products, with exact
degeneracies and boundary corrections; an abstract sign-kernel construction
shows why their orders alone do not force the spectral edge. No proof or
disproof of the Paley conjecture is claimed. No stronger clique or
asymptotic Paley operator bound follows from this pass. Independent
mathematical review remains outstanding.

Use the notation of the [elliptic model](parallel14-elliptic-model-2026-09-05.md):
p≡1 mod 4 prime, χ(r)=χ(r−1)=1, E:y²=x(x−1)(x−r), G=E(F_p),
H=G/G[2], L=|H|. Write f for the descended function on H. It is real,
even, zero exactly at 0, and ±1 elsewhere. Its lift to G is χ(y),
with value zero at O. In particular |H[2]|=4 and L=2|C|+4.

## 1. Every translated product, including repeated shifts

For a multiset of shifts t₁,…,t_k∈H, let O be its distinct shifts
of odd multiplicity and B its distinct shifts of positive even
multiplicity. Write s=|O|, e=|B|, and

\[
g_O(h)=\prod_{t\in O}f(h+t),\qquad
R(h)=\prod_{i=1}^k f(h+t_i).
\]

The empty product is 1. The following identity is exact:

\[
R(h)=g_O(h)-\sum_{b\in B}g_O(-b)\mathbf1_{h=-b}.       \tag{1}
\]

Indeed, an odd power of a number in {−1,0,1} equals that number;
a positive even power is 1 except at its zero. Distinct shifts give
distinct zeros, and O∩B is empty. Thus the coefficients g_O(−b)
in (1) are ±1, not uncomputed exceptional stalks.

For any character η of H put

\[
\widehat g_O(\eta)=\sum_{h\in H}g_O(h)\overline{\eta(h)}.
\]

If s>0, uniformly in the shifts, η, their orders, and k,

\[
|\widehat g_O(\eta)|\le s\sqrt p,
\quad
\widehat R(\eta)=\widehat g_O(\eta)
 -\sum_{b\in B}g_O(-b)\overline{\eta(-b)},
\quad |\widehat R(\eta)|\le s\sqrt p+e.                \tag{2}
\]

If s=0, there is instead the exact formula

\[
\widehat R(\eta)=L\mathbf1_{\eta=1}
                   -\sum_{b\in B}\overline{\eta(-b)}. \tag{3}
\]

Proof of the bound in (2). Lift each t∈O to T∈G. On E consider
the quadratic Kummer sheaf of the rational function
∏_{t∈O}y(P+T). Its odd divisor support is the union of the cosets
E[2]−T. These cosets are disjoint precisely because the shifts are
distinct in H. Each contributes four points with odd order, using
div(y)=(0,0)+(1,0)+(r,0)−3O. Thus there are exactly 4s nontrivial
tame punctures. No cancellation between distinct cosets is possible.

Tensor with the unramified Lang-character sheaf for the lift of η̄.
It has weight zero and does not change those local monodromies.
The resulting rank-one sheaf is geometrically nontrivial, with
H_c^0=H_c^2=0. The tame Euler characteristic on a genus-one curve
with 4s punctures is −4s, so dim H_c^1=4s. The trace and weight
bounds give a full sum of absolute value at most 4s√p. Periodicity
under G[2] makes that full sum exactly 4 ĝ_O(η). This proves (2).
Formula (1) then restores the even-multiplicity zero locations.

The cohomological inputs are the same primary Katz GKM statements
reviewed in pass fourteen: §§2.3.1–2.3.3, §3.6 and §4.3, printed
pp.32–33,40–41,59–60. The [source PDF](../sources/katz-gauss-kloosterman-monodromy.pdf)
and earlier audit are preserved. This pass applies those inputs to
4s disjoint punctures; it does not assert a new independent review
of the cohomological foundations or use a new unreviewed theorem.

## 2. The corresponding Fourier correlation and squared matrix

Let F_α=Σ_h f(h)ᾱ(h). Fourier inversion gives, for any k≥1,

\[
\widehat R(\eta)=L^{1-k}
 \sum_{\alpha_1\cdots\alpha_k=\eta}
       \prod_{j=1}^k F_{\alpha_j}\alpha_j(t_j).        \tag{4}
\]

Equations (1)–(3) therefore control these particular joint sums of
Fourier coefficients, with all diagonal multiplicities included.
They are not a bound for an arbitrary contraction of Fourier indices.

For M_{u,v}=f(u+v)f(u−v), its full square satisfies

\[
(M^2)_{u,v}=\sum_h f(h+u)f(h-u)f(h+v)f(h-v).         \tag{5}
\]

Apply (1)–(3) to the multiset {u,−u,v,−v}. This classifies every
degeneracy:

| Row pair | Exact value or bound for (M²)_uv |
|---|---|
| u=v∈H[2] | L−1 |
| u,v∈H[2], u≠v | L−2 |
| u∉H[2], v=±u | L−2 |
| Exactly one of u,v lies in H[2] | At most 2√p+1 in absolute value, with its single correction explicitly given by (1) |
| Both ordinary, v≠±u | At most 4√p in absolute value |

These are genuine uniform row-correlation estimates. Taking absolute
values of all off-diagonal entries still permits an accumulated error
of order L√p for the squared matrix. That route does not reach the
required operator scale. The following construction gives a stronger
reason to retain more than the orders of these bounds.

## 3. An abstract obstruction preserving the sign kernel and border

For every eligible p≥10000 and r, there exists an even function f′ on
the same H, again zero exactly at 0 and ±1 elsewhere, such that:

1. For every translated multiset and every character, (1) and (3)
   remain exact, and the bound in (2) holds with 50s√p+e in place
   of s√p+e. In particular |F′_η|≤50√p.
2. M′_{u,v}=f′(u+v)f′(u−v) has the same ±-fiber identifications,
   four translation symmetries, and full exceptional border as M.
3. On H₀=H\H[2], the centered normalized matrix

\[
 Z'=\frac{2}{\sqrt{7p}}\left(M'_{H_0}-J_{H_0}/\sqrt p\right)
\]

   has norm at least

\[
 \|Z'\|\ge\frac4{\sqrt7}\left(1-\frac2{\sqrt p}\right).
                                                               \tag{6}
\]

This lower bound stays strictly above 1 as p grows. The original
normalization is exact: passing from G to H divides both the kernel
and the all-ones matrix by a fiber factor of four. The ordinary
±-fibers then give exactly the old C matrix with additional zeros.

This is an obstruction for the abstract sign-kernel information and
the stated square-root *orders*. It does not preserve the sharp
constant 1 in the original Fourier bound or s in the product bound.
No original four-puncture Kummer realization is asserted for it. It
is not a Paley counterexample and does not preserve the ambient
Paley projection identity. It must
not be combined with the different projection-only construction in
pass twelve as though one example satisfied both sets of constraints.

### Construction and proof

The finite elliptic group has at most two cyclic factors, and hence

\[
H\simeq\mathbb Z/n\mathbb Z\times\mathbb Z/d\mathbb Z,
\qquad d\mid n,\qquad d,n\text{ even},\qquad nd=L.
\]

For completeness, its prime-primary subgroups have rank at most two
because the geometric ℓ-torsion on an elliptic curve has rank two
for ℓ≠p. Here p does not divide |G|: Hasse gives |G|<2p, whereas
16 divides |G| and p is odd, so |G| cannot equal p. The finite
abelian-group classification gives the asserted two cyclic factors.
Quotienting the rational full 2-torsion divides both even invariant
factors by two; full rational 4-torsion makes the resulting n,d even.

Put N=⌈√p⌉. Hasse gives

\[
L\ge(\sqrt p-1)^2/4\ge16N\qquad(p\ge10000).
\]

The last inequality follows from N≤√p+1 and z²−66z−63≥0 for
z≥100. The Hasse input was checked in the archived primary
[Bekker–Zarhin Theorem 4.1](https://arxiv.org/html/1702.02255v2#S4.Thmthm1).
It is an imported classical theorem, not a new point-counting result.

Choose

\[
k=\min(N,\lfloor n/4\rfloor),\quad \ell=\lceil N/k\rceil,
\quad A=\{1,\ldots,k\}\times\{0,\ldots,\ell-1\},\quad a=|A|.
\]

Here n≥√L≥8. If k=N then ℓ=1≤d. Otherwise k≥n/8 and
ℓ≤8N/n+1≤d/2+1≤d. Thus A is a well-defined rectangle,
N≤a<2N, A∩(−A)=∅, and A avoids H[2]. The first coordinate
lies strictly between 0 and n/2, which proves the last two assertions.

Let

\[
D=(A+A)\cup(A-A)\cup(-A-A).
\]

It is symmetric. Each of the three rectangle sum/difference sets
has at most (2k−1)min(d,2ℓ−1)≤4a elements, so |D|≤12a<24N.
Set f′(h)=1 on D\{0}, leave f′(0)=0, and keep f elsewhere.
Only some negative signs in D change; write Δ for their number.
Then Δ≤|D|<24N, and ‖f′−f‖₁=2Δ.

For any s distinct shifts, the two products can differ only in the
union of s translates of the modified set. Their absolute difference
is at most 2 pointwise, so its ℓ¹ norm is at most 2sΔ. Multiplying
by any character does not change that bound. Therefore (2) becomes

\[
s\sqrt p+2s\Delta+e
 \le s\sqrt p+48sN+e\le50s\sqrt p+e.                \tag{7}
\]

The last inequality uses N≤√p+1 and √p≥100, so 48N≤49√p.
Repeated-shift identities and the all-even case are unchanged exactly
because the zero set remains {0}. This proves the first assertion at
every order, without a restriction on the number of factors.

For distinct u,v∈A, neither u−v nor u+v is zero and both are in D.
Thus M′ restricted to A is J−I. Evenness implies identical rows
on u and −u, so the block on A∪(−A) consists of four copies of
J_a−I_a. The unit constant vector on that set has M′ Rayleigh
quotient 2(a−1), and the J/√p term contributes 2a/√p. Hence

\[
\langle v,Z'v\rangle
=\frac{4(a-1-a/\sqrt p)}{\sqrt{7p}}
\ge\frac4{\sqrt7}(1-2/\sqrt p),                     \tag{8}
\]

using a≥√p. This proves (6).

For t∈H[2], evenness gives M′_{t,u}=f′(t+u)², equal to 1
unless t=u, when it is 0. Thus every exceptional coupling and the
J₄−I₄ exceptional block remain exact. Also
M′_{u+t,v+t}=M′_{u,v} for every t∈H[2], and M′ has identical
rows on u,−u with zeros precisely when v=±u. These facts retain
all four symmetry blocks and the complete exceptional quotient from
pass fourteen, without assuming the action on ordinary fibers is free.

There is a concrete reason the ambient Paley identity cannot survive:
the quotient now contains a clique of size a among ordinary vertices,
and all three original anchors are adjacent to it. A zero-diagonal
±1 matrix containing that clique has norm at least a+2>√p. An
ambient matrix satisfying S²=pI−J would have norm √p. Thus the
new construction deliberately retains a different set of inputs from
the prior projection obstruction; the full arithmetic problem remains.

## 4. Verification and consequence for the next attempt

The [verifier](../experiments/parallel15_elliptic_correlations_2026_09_05.py)
checks parity reduction and every entry of (5) in five small elliptic
quotients. It certifies simultaneous character-twist bounds via exact
positive leading minors of s²pI−A_gA_gᵀ, where A_g is convolution
by g_O. This works even when g_O is not even and its Fourier values
are complex. It also checks (4) over exact cyclotomic rings.

The planted examples at p=10009 and p=65537 are reconstructed from
actual curve points. Product-group coordinates are verified through
both generating translations and a bijection of all points, avoiding
a quadratic-sized group table. The second case tests the balanced
two-dimensional rectangle, not just a long cyclic interval. Every
exceptional coupling, the planted quadratic form, a rational Rayleigh
lower bound above 1 and ℓ¹ perturbation certificates are checked.
See the [results](../results/parallel15_elliptic_correlations_2026_09_05.json).

The higher-correlation theorem is more information than the individual
Fourier bound. Its uniform square-root order, even together with the
sum-and-difference form and exact border, is insufficient for the
sharp edge. A further attempt needs the sharp constants and/or
additional arithmetic identities controlling how the correlations fit
together. The construction does not rule out such an argument. The
previous logarithmic energy range, subgroup and classical gaps, and
exact official prize bridge remain unchanged and unproved.
