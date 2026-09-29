# Pairwise hypergeometric kernels and a bounded aggregate across ranks

**Status: uniform pairwise bounds and a signed quadratic aggregate are
proved below. They do not bound the full necklace aggregate.** The
new equal-rank correction generalizes the previous symmetric-square
Legendre correction. The last section quantifies the remaining cost of
passing to another translated anchor. No novelty claim is made.

Throughout, p=1 mod 4 is prime, G=F_p^*, chi(0)=0, and

\[
 S_{xy}=\chi(x-y),\quad D=\operatorname{diag}(\chi(x)),\quad
 K_r=S(DS)^{r-1},\quad k_r(t)=(K_r)_{1t},\qquad r\ge1.
\]

The [preceding proof](parallel3-necklace-2026-09-04.md) identifies k_r
with the r-fold multiplicative convolution of chi(1-t). It constructs
a rank-r middle-extension sheaf F_r whose trace at F_p points is -k_r,
including at t=1. On U=P^1 minus {0,1,infinity}, it is pure of weight
r-1 and geometrically irreducible. Its monodromy at zero is one
unipotent Jordan block; at one it is a nontrivial pseudoreflection with
determinant chi^r. All three local monodromies are tame.

These inputs, including the duality formula used below, are from
[Katz, G2 and Hypergeometric Sheaves, Section 2, pp.3-5](https://web.math.princeton.edu/~nmk/g2hyper62finalcorrected.pdf).
Euler-Poincare, the curve trace formula, and the upper weight bound are
the [previously checked Katz-GKM input](https://web.math.princeton.edu/~nmk/Katz-GKM.pdf),
Sections 2.3 and 3.6. Those established theorems are imported.

## 1. Unequal ranks, with the actual tensor stalk

For any multiplicative character rho of G and distinct positive r,s,

\[
 \boxed{\left|\sum_{t\in G}\rho(t)k_r(t)k_s(t)\right|
 \le(r+s-1)p^{(r+s-1)/2}.} \tag{1}
\]

Use T=F_r tensor F_s tensor L_rho on G_m, tensoring the actual middle
extensions. Its generic rank is rs, while its stalk at one has dimension
(r-1)(s-1). A compactly supported H^2 on U would yield a homomorphism
between the irreducible rank-r factor and a twist of the dual rank-s
factor. Rank mismatch rules this out. The actual tensor embeds in the
middle extension of its restriction, so it has no punctual sections.
The open-closed exact sequence, as in the preceding proof, gives

\[
 \dim H_c^1(G_m,T)=rs-(r-1)(s-1)=r+s-1.
\]

Its weights are at most r+s-1. The trace formula proves (1) for every
rho, including the trivial character, with no restriction r,s<p.

## 2. The exact correction at equal rank

Define, for t in G,

\[
 q_r(t)=k_r(t)^2+p^{r-1}1_{t=1}-p^{r-1}.
\]

Then for every multiplicative rho, including the trivial character,

\[
 \boxed{\left|\sum_{t\in G}\rho(t)q_r(t)\right|
 \le(2r-2)p^{r-1/2}.} \tag{2}
\]

For r=1 the right side and q_1 itself are zero. The constant subtraction
and singular-point addition are both part of the formula.

Here are the geometric and linear-algebra details. Katz's duality and
additive-character-change formulas, applied to the lists of trivial
and quadratic characters, give H_r^vee=H_r(2r-1): any change of additive
character contributes no translation or constant twist, since the two
list lengths agree and chi(-1)=1. The constant Gauss twist in F_r has
eigenvalue G_p^(-r), with G_p^2=p. Therefore

\[
 F_r^\vee\cong F_r(r-1),\qquad
 F_r\otimes F_r\cong\operatorname{End}(F_r)(-(r-1)). \tag{3}
\]

Split End(F_r)=1 direct-sum End^0(F_r) using the trace projection divided
by r, over the characteristic-zero coefficient field. Let Q_r be the
middle extension to G_m of End^0(F_r)(-(r-1)) on U. It has generic
rank r^2-1, and its trace on U is k_r(t)^2-p^{r-1}.

We must calculate its trace at one, rather than substitute a generic
formula there. Write V for the local rank-r representation. Since its
inertia generator is a pseudoreflection, it is either diag(-1,1,...,1)
(r odd) or J_2 direct-sum 1^(r-2) (r even). In the latter case
J_2 tensor J_2 has two invariant vectors, each mixed J_2 tensor 1
block has one, and each 1 tensor 1 block has one. In both cases

\[
 \dim (V\otimes V)^{I_1}=(r-1)^2+1,
 \qquad \dim(V^{I_1}\otimes V^{I_1})=(r-1)^2. \tag{4}
\]

Under (3), the global invariant tensor is the identity endomorphism,
with Frobenius eigenvalue p^{r-1}. It is not in the smaller tensor of
invariants: as an endomorphism, every element of that smaller space
has image in the proper subspace V^(I_1), whereas the identity does not.
It therefore generates the one-dimensional quotient in (4). The stalk
trace of the middle extension of F_r tensor F_r is consequently
k_r(1)^2+p^{r-1}. Removing the global invariant line shows that the
trace of Q_r at one is k_r(1)^2. Thus its trace on all G is exactly
q_r(t), and its stalk dimension at one is (r-1)^2.

For any rho, Q_r tensor L_rho has no geometric invariant or coinvariant.
For rho=1 this is Schur's lemma and removal of the scalar endomorphisms.
For rho different from one, a rank-one constituent of End(F_r) would give
F_r isomorphic geometrically to its rho-twist. At zero their inertia
eigenvalues would be respectively all one and all rho. A nontrivial
Kummer character has nontrivial inertia at zero, so this is impossible.
Geometric semisimplicity follows from purity; alternatively the
coinvariant assertion follows by duality and the same argument.

This middle-extension sheaf has H_c^0=H_c^2=0. Its Euler characteristic
on G_m is -(r^2-1)+(r-1)^2=-(2r-2). Its H_c^1 therefore has dimension
2r-2 and weights at most 2r-1, proving (2).

## 3. Matrix forms and the large exceptional eigenvalue

On G define the real symmetric convolution matrices

\[
 (V_{r,s})_{xy}=\chi(xy)(K_r)_{xy}(K_s)_{xy}.
\]

Their convolution kernels are chi(t)k_r(t)k_s(t). Equation (1) gives

\[
 \|V_{r,s}\|\le(r+s-1)p^{(r+s-1)/2}\quad(r\ne s).
\]

With v=(chi(x))_(x in G), (2) gives the exact corrected operator bound

\[
 \boxed{H_r=V_{r,r}+p^{r-1}I-p^{r-1}vv^T,
 \qquad \|H_r\|\le(2r-2)p^{r-1/2}.} \tag{5}
\]

The exceptional mode can also be checked without sheaves. The Mellin
eigenvalues of C^*=(DS)_(G,G) are Jacobi sums J(psi,chi), of absolute
value sqrt(p) except for psi=1,chi, where both are -1. Parseval gives

\[
 \sum_{t\in G}k_r(t)^2
 =\frac{(p-3)p^r+2}{p-1}
 =p^r-2\sum_{a=0}^{r-1}p^a.
\]

It follows that

\[
 V_{r,r}v=\left(p^r-2\sum_{a=0}^{r-1}p^a\right)v,
 \qquad H_rv=-2\left(\sum_{a=0}^{r-2}p^a\right)v. \tag{6}
\]

The empty sum is zero. For r=2, (5) is precisely the previously proved
Legendre correction. Applying the unequal-rank bound blindly on the
diagonal would miss the term of size p^r in (6).

## 4. A signed quadratic aggregate across lengths

For R>=1, set h_r(t)=k_r(t)/p^((r-1)/2), 1<=r<=R, and define

\[
 G^\rho_{rs}=\frac1p\sum_{t\in G}\rho(t)h_r(t)h_s(t),\qquad
 c_R=\frac{3R^2-R-2}{2},\qquad
 \Delta_R=\frac{c_R}{\sqrt p}+\frac2p.
\]

Equations (1)-(2) give the uniform matrix bound

\[
 \boxed{\|G^\rho-1_{\rho=1}I_R\|\le\Delta_R.} \tag{7}
\]

For distinct indices the entry bound is (r+s-1)/sqrt(p). On the
diagonal, (2) gives

\[
 G^\rho_{rr}=1_{\rho=1}-\frac{1+1_{\rho=1}}p
       +\frac1{p^r}\sum_t\rho(t)q_r(t).
\]

The last term is at most (2r-2)/sqrt(p). The maximum absolute row sum
of the resulting error bound is c_R/sqrt(p)+2/p; its maximum column
sum is the same. This proves (7), even for a complex rho, when the
matrix need not be Hermitian.

In particular, simultaneously for every complex coefficient vector a,

\[
 \boxed{\left|\frac1p\sum_{t\in G}\rho(t)
       \left|\sum_{r=1}^R a_rh_r(t)\right|^2
       -1_{\rho=1}\sum_{r=1}^R|a_r|^2\right|
       \le\Delta_R\sum_{r=1}^R|a_r|^2.} \tag{8}
\]

This is a uniform aggregate, rather than one bound for each chosen
coefficient vector. It is asymptotically an isometry with vanishing
nontrivial Mellin biases when R=o(p^(1/4)), in particular when R grows
logarithmically. It involves one-parameter kernels indexed by their
lengths; it is not the sum over arbitrary necklace label patterns.

## 5. Why another anchor is still costly

For a function w on G, let w_hat(rho)=sum_t w(t) conjugate(rho(t)),
and let w_bar=(p-1)^(-1)sum_t w(t). Mellin inversion and (8) imply

\[
 \left|\frac1p\sum_t w(t)\left|\sum_r a_rh_r(t)\right|^2
           -\bar w\|a\|_2^2\right|
 \le\frac{\Delta_R}{p-1}\sum_\rho|\widehat w(\rho)|\,\|a\|_2^2.
 \tag{9}
\]

For w(t)=chi(1-t), these Mellin coefficients are Jacobi sums. Exactly
two have absolute value one and the other p-3 have absolute value
sqrt(p). The multiplier in (9) is therefore

\[
 \frac{(p-3)\sqrt p+2}{p-1}.
\]

For fixed R>=2 it cancels the sqrt(p) saving in (7), leaving only an
O(R^2) bound by this argument. This concerns the strength of this upper
estimate, not a lower bound on the actual weighted sum. Cancellation
among the Mellin coefficients remains possible and unproved. Thus the
uniform kernel aggregate is a useful input, but inserting another
translated anchor still needs a stronger argument.

## Verification

The [exact verifier](../experiments/parallel4_kernel_aggregate_2026_09_04.py)
checks the pair kernels, their exceptional eigenvalues, and the corrected
operator bounds by integer positive-definiteness certificates. It also
checks the local pseudoreflection tensor dimensions and that the identity
endomorphism supplies the missing one-dimensional stalk. Selected signed
aggregate identities and inequalities are checked with exact rational
quadratic-field arithmetic. [Results](../results/parallel4_kernel_aggregate_2026_09_04.json)
record the finite cases and hashes. The cohomological proof supplies the
uniform statement; finite tests alone do not do so.
