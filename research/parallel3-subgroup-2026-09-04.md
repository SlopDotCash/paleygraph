# Why unramified cyclotomic structure does not control the mixed-energy defect

**Status: no uniform mixed-energy bound or new infinite quartic family
is proved.** This pass establishes a uniform algebraic description of
the defect and a concrete obstruction to using unramifiedness, an
invertible derivative, or Galois symmetry to eliminate it. The remaining
pointwise arithmetic estimate is identified below. No prime scan or
integer factorization was performed in this pass.

## Two orders in the same number field

Let n=2k≥4 be dyadic, d=n/4, ζ a primitive n-th root, and
t=ζ+ζ⁻¹. Define C_d by

\[
 C_1(T)=T,\qquad C_{2d}(T)=C_d(T)^2-2.
\]

Its roots are \(2\cos((2j-1)\pi/(2d))\), 1≤j≤d. This is the
minimal polynomial of t: these are the d distinct Galois conjugates
of t. Set

\[
 A(T)=-(2+T)^k,\quad \alpha=A(t)=(1+\zeta)^n,
 \quad\mathcal O=\mathbb Z[t],\quad\mathcal A=\mathbb Z[\alpha].
\]

The equality for α follows from
\((1+\zeta)^2=\zeta(2+t)\) and ζ^k=−1. The minimal polynomial of
α is precisely R_n from
[the primitive-fiber note](parallel2-subgroup-2026-09-04.md).
Its conjugates are distinct, so α also generates this degree-d field.
Consequently \(\mathcal A\subseteq\mathcal O\) has finite index

\[
 I_n=[\mathcal O:\mathcal A]\in\mathbb Z_{>0}.
\]

This discussion requires no assertion that either displayed order is
globally maximal.

Let M_n be the integer matrix whose j-th column is the coefficient
vector of α^j in the basis 1,t,…,t^(d−1), for 0≤j<d. Then

\[
 \boxed{I_n=|\det M_n|,\qquad
 \operatorname{disc}(R_n)=I_n^2\,2^{d-1}d^d.} \tag{1}
\]

**Proof.** The determinant gives the index of the two integer lattices.
The trace Gram matrix for the α-power basis is M_nᵀ times the trace
Gram matrix for the t-power basis times M_n. Their determinants are
the corresponding polynomial discriminants.

It remains to calculate the discriminant of C_d. At
T_j=2cosθ_j, θ_j=(2j−1)π/(2d),
\(|C_d'(T_j)|=d/\sin\theta_j\). For d≥2,
C_d(2)=C_d(−2)=2, which gives
\(\prod_j(4-T_j^2)=4\), hence
\(\prod_j\sin\theta_j=2^{1-d}\). Therefore
\(\operatorname{disc}(C_d)=2^{d-1}d^d\). The d=1 case is immediate.
∎

In particular, every odd prime is prime to the ambient discriminant.
Repeated roots of R_n at an odd prime can nevertheless occur, through
the index I_n.

## The defect is a loss of rank of the order inclusion

Let p≡1 mod n be prime. C_d splits into distinct roots t₁,…,t_d in
F_p. Evaluation identifies

\[
 \mathcal O/p\mathcal O\cong\mathbb F_p^d.
\]

Let λ_i=A(t_i), and let c₁,…,c_r be the multiplicities of the distinct
λ-values. Evaluation of M_n yields the rows
\((1,\lambda_i,\ldots,\lambda_i^{d-1})\). There are r distinct rows,
which are independent by the Vandermonde determinant. Thus

\[
 \boxed{\operatorname{rank}_{\mathbb F_p}M_n=r,\qquad
 \dim\ker(M_n\bmod p)=d-r=\sum_i(c_i-1).} \tag{2}
\]

The map \(\mathcal A/p\mathcal A\to\mathcal O/p\mathcal O\)
need not be injective even though \(\mathcal A\to\mathcal O\)
is injective over the integers. Its image consists of coordinate
functions constant on each λ-fiber. Its kernel is the nilradical of
\(\mathbb F_p[Y]/(R_n)\), of dimension d−r.

This distinguishes two facts: the ambient cyclotomic algebra modulo p
is reduced, while the algebra defined by R_n modulo p can contain
nilpotents. Unramifiedness of the ambient order does not remove the
fiber collisions counted by B.

## Even the original polynomial derivative is always invertible

The **unreduced** polynomial A(T)=−(2+T)^k satisfies

\[
 A'(T)=-k(2+T)^{k-1},\qquad
 \boxed{|N_{\mathbb Q(t)/\mathbb Q}(A'(t))|
              =k^d2^{k-1}.} \tag{3}
\]

Indeed the absolute norm of 2+t is 2, by C_d(−2), including d=1.
The right side of (3) is a power of two. At every odd prime, A' is
therefore nonzero at every root of C_d. Nevertheless, two or more
distinct roots of C_d can have the same value under A. Invertibility
of the derivative does not prove injectivity on this finite set.

The word “unreduced” matters: differentiating the representative of
A modulo C_d gives a different polynomial derivative. No assertion
about that alternative derivative is made.

## A quartic counterexample to Galois-stable fibers

Take the already certified quartic field

\[
 p=67403009,\quad n=128,\quad d=32,\quad g=64701253.
\]

The element g has order 128. Write
\(\lambda_j=(1+g^j)^{128}\bmod p\). Direct modular arithmetic gives

\[
 \lambda_{13}=\lambda_{25}=12400204,
\]

but after applying the exponent permutation induced by the Galois
automorphism ζ↦ζ³,

\[
 \lambda_{39}=7533465\ne29566144=\lambda_{75}.
\]

Thus the fiber equivalence relation at this prime is not preserved by
the global Galois action. In particular, the fibers are not necessarily
cosets of a subgroup of that Galois group. Here they consist of 28
singletons and two doubletons; M_n has rank 30 modulo p. All roots
of C_d remain distinct and all values A'(t_i) remain nonzero.

The reason is that a Galois automorphism moves the prime ideal used
for reduction. It does not act on a chosen copy of F_p by fixing
that reduction map. Complete splitting provides no nontrivial field
automorphism of F_p that could implement this permutation.

This is a concrete obstruction to a prospective proof that assumes
the cyclic Galois group, unramified reduction, or invertible derivative
forces equal-sized or singleton fibers. None of those assumptions
controls the pointwise index defect.

## A precise remaining index estimate

Lift the distinct roots t_i to Z_p, which is possible uniquely by
Hensel's lemma, and put λ_i=A(t_i). From (1), since p is odd,

\[
 \boxed{v_p(I_n)=\sum_{i<j}v_p(\lambda_i-\lambda_j).} \tag{4}
\]

All differences are nonzero in characteristic zero. If c_(r,i) are
the multiplicities of the λ-values modulo p^r, then

\[
 v_p(I_n)=\sum_{r\ge1}\sum_i\binom{c_{(r,i)}}2.
\]

At the first precision the earlier fiber formula gives

\[
 \boxed{B-k^2=4k\sum_i\binom{c_{(1,i)}}2
                \le4k\,v_p(I_n).} \tag{5}
\]

Hence a uniform bound

\[
 v_p(I_{2k})=O(k\log k)
 \qquad\bigl((2k)^4/4\le p\le(2k)^4,\ p\equiv1\pmod{2k}\bigr)
 \tag{6}
\]

would imply the needed near-quadratic mixed energy at those endpoints.
Appropriate bounds at the smaller tower levels are additionally needed
for the earlier whole-tower implication. Estimate (6) is **not proved**;
it can be stronger than the energy obligation because it counts
collisions at every p-adic precision.

The elementary height bound does not close (6). All d real roots of
R_n lie in (−2^n,0), so

\[
 \log_2 I_n\le\frac{n d(d-1)}2.
\]

It gives only \(v_p(I_n)=O(k^3/\log p)\), and (5) then costs
O(k⁴/log p), worse than the elementary bound B≤k³/2 in the quartic
window. This is a limitation of that height estimate, not a proof
that more precise arithmetic information about I_n cannot help.

## Literature boundary

Do Duc, Leung, and Schmidt's
[Main Theorem 1](https://arxiv.org/html/1903.07314) bounds every
cyclotomic number for a subgroup of cardinality k by 3 when
\(p>(\sqrt{14})^{k/\operatorname{ord}_k(p)}\). Applied to the child
K and its cosets, this would give B≤3k². In our prime-field setting,
\(\operatorname{ord}_k(p)=1\). For dyadic k≥16,
\(14^{k/2}>16k^4\), so the theorem's hypothesis lies beyond the
entire working window p≤(2k)⁴. The inequality follows at k=16 by
integer comparison and persists on doubling, since the exponential
side gains a factor 14^(k/2), larger than the polynomial factor 16.
The paper's stronger theorem for prime k does not cover dyadic k≥4.
This is a verified failure of a quantitative bridge, not a claim that
the literature rules out the desired estimate.

## Exact verification

`experiments/parallel3_subgroup_2026_09_04.py` constructs C_d and the
matrix M_n with integer arithmetic. For n=4,8,16,32,64, it checks
(1) against an independent multiplication-matrix discriminant of R_n,
whose coefficients were previously obtained by Newton's identities.
It also checks (3) as a separate determinant identity. Bareiss
elimination verifies every division is exact.

Modular certificates cover the dense case (p,n)=(17,16) and the
quartic cases (67403009,128) and (17189277697,512). They check
ambient distinctness, invertible A', rank defect, the collision
fibers, and the explicit Galois counterexample above. The coincident
roots in all three selected cases separate modulo p², certifying
the exact corresponding index valuations without factoring I_n.

The result file records full hexadecimal indices and discriminants,
the modular certificates, and source hashes. These verify the
algebraic obstruction and its examples; they do not prove (6), the
balanced logarithmic-depth moment bound, or either Paley target.
