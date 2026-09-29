# Multiple kernel seeds, and a symmetry obstruction to the old span

**Status: the old kernel span misses an entire symmetry sector, even
with arbitrarily many powers. Multiple seeds remove that particular
obstruction and admit uniform correlation bounds. The full restricted
Paley spectral bound remains unproved.** No novelty claim is made.

Throughout p≡1 mod 4 is prime, G=F_p*, χ(0)=0, and
K_r=S(DS)^(r−1), k_r(t)=K_r(1,t), as in the
[previous kernel proof](parallel5-anchor-aggregate-2026-09-04.md).

## 1. An actual missing spectral sector

The exact identity k_r(t^(-1))=χ(t)k_r(t) holds for every r≥1.
Indeed K_r is symmetric and K_r(x,y)=χ(x)k_r(y/x) for nonzero
x,y. Apply both facts to (x,y)=(t,1).
Define the unitary involution on functions on G by

\[
 (Jf)(t)=\chi(t)f(t^{-1}).
\]

Every k_r is in its +1 eigenspace. The fixed points of inversion are
1 and −1, both with χ=1; the other p−3 elements occur in pairs.
Consequently

\[
 \dim\ker(J-I)=(p+1)/2,\qquad \dim\ker(J+I)=(p-3)/2.
 \tag{1}
\]

Thus adding powers to the single-seed family can never fill the
ambient function space. On the common neighborhood
Q={t:χ(t)=χ(t−1)=1}, the involution is ordinary inversion and all
the restricted k_r are inversion-even. The restricted character
matrix S_Q commutes with inversion, since
χ(t^(-1)−u^(-1))=χ(tu)χ(t−u)=χ(t−u) on Q.

There is an exact example where the omitted sector contains the
largest eigenvalue. At p=13, Q={4,10}, and

\[
 S_Q=\begin{pmatrix}0&-1\\-1&0\end{pmatrix}.
\]

Every restricted k_r is a multiple of (1,1), with eigenvalue −1.
The orthogonal vector (1,−1) has eigenvalue +1. For the compressed
Paley projection P_Q=(I−J_2/13+S_Q/√13)/2, its eigenvalues are

\[
 \lambda_{\rm even}=\frac12-\frac1{2\sqrt{13}}-\frac1{13},
 \qquad
 \lambda_{\rm odd}=\frac12+\frac1{2\sqrt{13}}.
\]

The larger one is invisible to the old span. This refutes a universal
claim that its extremal eigenvector must lie in that span; it does
not establish an asymptotic obstruction for every prime.

## 2. Distinct seeds and their exact correlation dimensions

For a nonzero seed b define

\[
 k_{b,r}(t)=k_r(t/b),\qquad h_{b,r}(t)=p^{-(r-1)/2}k_{b,r}(t).
\]

At p=13 the rank-one seed b=3 restricts to (1,−1) on Q, so this
enlargement actually detects the omitted eigenvector in the example.
For a set T⊂G disjoint from the seeds b,c, put
w_T(t)=∏_(a∈T)χ(t−a), m=|T|. If T is nonempty, then for every
multiplicative character ρ and all ranks r,s≥1,

\[
 \boxed{\left|\sum_{t\in G}\rho(t)w_T(t)k_{b,r}(t)k_{c,s}(t)\right|
 \le\bigl(mrs+r+s-1_{b=c}\bigr)p^{(r+s-1)/2}.}       \tag{2}
\]

For T empty and b≠c the same formula holds with m=0. For T empty,
b=c and r≠s it is the preceding unequal-rank bound. The remaining
diagonal b=c,r=s has the exact correction

\[
 q_{b,r}(t)=k_{b,r}(t)^2+p^{r-1}1_{t=b}-p^{r-1},
 \quad\left|\sum_t\rho(t)q_{b,r}(t)\right|
 \le(2r-2)p^{r-1/2}.                                  \tag{3}
\]

Equation (3) is the earlier correction under t=bu; its Mellin
transform is multiplied by ρ(b). Its singular addition is at b,
not at one.

**Proof of (2).** The imported rank-r sheaf F_r has trace −k_r,
weight r−1, and a nontrivial tame pseudoreflection at one. Pull it
back by t↦t/b. Tensor the two actual pullbacks with the indicated
Kummer sheaves. The resulting rank-rs sheaf has the trace in (2),
including its actual stalks. All ramification is tame.

If T is nonempty, a point in T is regular for both pulled kernels
and has scalar quadratic inertia for the product. Thus it has no
local invariant or coinvariant there, and no global one. If T is
empty and b≠c, a possible invariant would give an isomorphism of
irreducible factors (so r=s). At b one side has a nontrivial
pseudoreflection and the other is unramified; L_ρ is also unramified
at b. This is impossible, including rank one. If only the ranks
differ, irreducibility and rank mismatch exclude it.

For b≠c, remove 0,∞,b,c and T. The open-curve H_c¹ dimension is
(m+2)rs. Restoring the actual stalks at b and c subtracts
(r−1)s+r(s−1), giving mrs+r+s. For b=c the corresponding count is
(m+1)rs−(r−1)(s−1)=mrs+r+s−1. Stalks at the anchors in T are
zero. The tensor injects into the middle extension of its generic
restriction, so it has no punctual H_c⁰. The trace formula and
upper H_c¹ weight r+s−1 prove the bound.

The geometric inputs are the already checked
[Katz Section 2](https://web.math.princeton.edu/~nmk/g2hyper62finalcorrected.pdf)
and [Katz-GKM curve machinery](https://web.math.princeton.edu/~nmk/Katz-GKM.pdf).
There is no assertion of independence of different sheaves beyond the
explicit singularity argument, and no rank<p restriction.

The disjointness from T is material. If r=s=1, b≠c and T={b,c},
the product k_(b,1)k_(c,1)w_T equals χ(bc) off b,c and zero there.
Its trivial-character sum is χ(bc)(p−3). A square-root estimate with
a fixed coefficient would be false. In an adjacency cell, rank-one
seeds at its own anchors similarly become constant multiples of one
another. They must not be included in the isometry below.

## 3. Uniform quadratic bounds for the enlarged family

Choose L distinct nonzero seeds B, ranks 1≤r≤R, and put

\[
 U_R=\sum_{r=1}^Rr^2,\quad
 C_{L,R}=\frac{LR(3R+1)}2-R-1,\quad
 \Delta_{L,R}=C_{L,R}/\sqrt p+2/p.
\]

Let z be any complex coefficient vector indexed by B×{1,...,R},
and write F_z(t)=Σ_(b,r)z_(b,r)h_(b,r)(t). From (2)–(3),

\[
 \left|\frac1p\sum_t\rho(t)|F_z(t)|^2
               -1_{\rho=1}\|z\|_2^2\right|
 \le\Delta_{L,R}\|z\|_2^2.                             \tag{4}
\]

Indeed the absolute error row sum at rank r is at most
[L(Rr+Σs)−R−1]/√p+2/p, maximized at r=R. The same bound holds
for column sums, so it bounds the operator norm even for complex ρ.
On the equal-seed diagonal the exact correction contributes
−(1_(ρ=1)+ρ(b))/p, with absolute value at most 2/p.

For nonempty T disjoint from all seeds, mrs+r+s−1_(b=c)≤(m+2)rs
and Cauchy–Schwarz on the coefficient indices gives

\[
 \left|\frac1p\sum_t\rho(t)w_T(t)|F_z(t)|^2\right|
 \le\frac{(m+2)LU_R}{\sqrt p}\|z\|_2^2.                \tag{5}
\]

Thus the gain is uniform over a growing number of seed locations as
well as ranks. It is not obtained by assuming that every seed shares
the old inversion symmetry.

## 4. Adjacency cells and the limits of the enlargement

Let A⊂G be a nonempty set of a anchors, disjoint from B, with prescribed
signs ε_x. Let C be the cell χ(t−x)=ε_x for all x∈A, with t∈G\A.
One may also prescribe χ(t)=ε_0; let e_0 be one if this extra condition
is imposed and zero otherwise. Put

\[
 \Gamma_a=a/2+2(1-2^{-a}),
 \quad
 E_{a,L,R}=
 \frac{2^{-a}C_{L,R}+\Gamma_a LU_R}{\sqrt p}
 +\frac{2^{1-a}+aLU_R/2}{p}.
\]

Then

\[
 \boxed{\left|\frac1p\sum_{t\in C}|F_z(t)|^2
              -2^{-(a+e_0)}\|z\|_2^2\right|
       \le E_{a,L,R}\|z\|_2^2.}                        \tag{6}
\]

Expand the a nonzero-anchor factors. The empty subset uses (4).
The average of |T|+2 over nonempty subsets is Γ_a, giving (5)'s
contribution. If the extra condition at zero is present, average
the two Mellin twists ρ=1 and χ; (4)–(5) apply to both with the
same errors and only the trivial twist contributes a main term.
At each nonzero anchor the soft indicator is at most 1/2. Removing
these points costs at most aLU_R||z||²/(2p), since
|h_(b,r)(t)|≤r and hence |F_z(t)|²≤LU_R||z||². This proves (6)
with every anchor correction retained.

For fixed σ≥0, ε>0 with σ+ε<1/2, the relative error in (6)
tends to zero uniformly when

\[
 L\le p^\sigma,\quad R=O(\log p),\quad
 a\le(1/2-\sigma-\epsilon)\log_2p.
\]

Nevertheless, taking all p−1 seeds does not give a full spectral
estimate from (6). Even at rank one the coefficient dimension is
then too large for an isometry on a proper cell. Explicitly, the
matrix F_(t,b)=h_(b,1)(t) equals S_(G,G)diag(χ(b)), and

\[
 FF^T=pI-\mathbf1\mathbf1^T-vv^T,\qquad v_t=\chi(t).
 \tag{7}
\]

Its eigenvalues are p on the orthogonal complement of {1,v}, and
one on those two directions. It is invertible. Thus a coefficient
vector can represent a delta function at any point outside the cell;
its cell mass is zero. No relative isometry error less than one can
hold on that full coefficient space for a nonempty proper cell.
This is distinct from a spectral estimate for S_C, which might still
hold by a different argument.

The old single-seed eigenvector claim is therefore refuted, and the
new seed family has a proved range. The remaining obligation is to
control the full restricted operator or its weighted word aggregate,
not to assume that this range covers every eigenvector.

## Verification

The [verifier](../experiments/parallel6_seeded_kernels_2026_09_04.py)
checks the involution, the exact p=13 extremal-vector example, the
rank-one full-seed Gram matrix, pairwise bounds including all Mellin
modes through convolution-operator certificates, quadratic combinations,
and adjacency cells with and without a prescribed zero-anchor sign.
The seed/anchor collision example is retained as a negative control.
[Results](../results/parallel6_seeded_kernels_2026_09_04.json) record
finite cases and hashes. No finite test certifies the imported
cohomology theorems or the unresolved full Paley estimate.
