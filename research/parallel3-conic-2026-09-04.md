# Exact conic and symmetry models of the two-anchor matrix

**Status: exact changes of representation, with no new uniform spectral
bound.** The boundary block and off-diagonal Mellin terms below are part
of the problem. They cannot be discarded when bounding extreme eigenvalues.
This complements the [aggregate criterion](parallel2-spectral-transfer-2026-09-04.md).
All calculations below are elementary; no novelty claim is made.

Let p be prime, p=1 mod 4, and let chi(0)=0. Write

\[
 R=\{x\in\mathbb F_p:\chi(x)=\chi(x-1)=1\},\qquad
 S_R=(\chi(x-y))_{x,y\in R}.
\]

Thus R is the common neighborhood of the clique {0,1}. Set N=p-1.

## 1. The conic lift, including the four exceptional parameters

For t in G=F_p^*, define

\[
 x(t)=\left(\frac{t+t^{-1}}2\right)^2,
 \qquad f(t)=\chi(t^2-1).
\]

The identities y=(t+t^(-1))/2, z=(t-t^(-1))/2 give y^2-z^2=1.
For x in R, the four choices of square roots y^2=x, z^2=x-1
give four distinct parameters t=y+z. They form
{t,-t,t^(-1),-t^(-1)}. There are exactly two parameters above 1,
namely +/-1, and two above 0, namely the roots of t^2=-1.
Consequently x(G)=R union {0,1}, every ordinary fiber has size four,
and |R|=(p-5)/4.

For every t,u in G, including those exceptional parameters,

\[
 x(t)-x(u)=\frac{(t^2-u^2)(t^2u^2-1)}{4t^2u^2},
 \qquad
 M_{t,u}:=\chi(x(t)-x(u))=f(t/u)f(tu).
 \tag{1}
\]

Let T be G with the four exceptional parameters removed. The incidence
matrix F from T to R has F^T F=4I and M_(T,T)=F S_R F^T.
Therefore F/2 is an isometry and the nonzero part of the lifted
operator is exactly 4S_R; its orthogonal complement is in the kernel.
This wording permits additional zero eigenvalues of S_R.

The full matrix M has an equally explicit quotient. For each x in R
take the normalized fiber indicator 1_(x(t)=x)/2. Add

\[
 q_+=(1_{x(t)=0}+1_{x(t)=1})/2,\qquad
 q_-=(1_{x(t)=0}-1_{x(t)=1})/2.
\]

These |R|+2 vectors are orthonormal and span all fiber-constant vectors.
In this basis M is the following matrix, with zero on the complementary
space:

\[
 \boxed{
 \begin{pmatrix}
 4S_R&4\mathbf1&0\\
 4\mathbf1^T&2&0\\
 0&0&-2
 \end{pmatrix}.}
 \tag{2}
\]

Indeed, every ordinary vertex has character +1 against both anchors,
and the two anchors have character +1 against each other. Each ordinary
fiber has size four and each anchor fiber size two. These facts give
every entry of (2). In particular, deleting the four parameters removes
a coupled border of size one as well as the eigenvalue -2. The coupling
has norm 4sqrt(|R|), of the same order as the desired spectral scale;
its rank alone does not justify ignoring it.

## 2. Exact Mellin transform: two Jacobi products, not diagonalization

Choose a generator g of G, put zeta=exp(2pi i/N), and set

\[
 F_j=\sum_{r=0}^{N-1}f(g^r)\zeta^{-jr},\qquad
 e_a(g^r)=N^{-1/2}\zeta^{ar}.
\]

Indices are taken modulo N. Fourier expansion of each factor in (1)
and the sum over u give the exact matrix entry

\[
 \boxed{\langle e_a,M e_b\rangle
 =\frac1N\sum_{\substack{j\bmod N\\2j=a+b}} F_jF_{j-b}.}
 \tag{3}
\]

There are two terms if a+b is even and none otherwise. Since
f(-t)=f(t) and f(t^(-1))=f(t), F_j=0 for odd j and all F_j are real.
In particular (3) vanishes unless a,b are even and a=b mod 4.
These are blocks, not individual Mellin eigenvalues.

For rho(g)=zeta and J(alpha,chi)=sum_(v in G) alpha(v)chi(1-v),
substitution v=t^2 gives, including every exceptional character,

\[
 \boxed{F_{2j}=J(\rho^{-j},\chi)+J(\chi\rho^{-j},\chi).}
 \tag{4}
\]

Here both terms are needed: the number of preimages of a nonzero v is
1+chi(v). In particular F_0=-2, since each of J(1,chi) and
J(chi,chi) is -1. Formula (3) still couples different character indices.
The exact verifier records nonzero off-diagonal entries already at p=13.
Bounds on each Jacobi sum do not by themselves control the norm of this
growing matrix: for example, a d by d matrix of ones has entries one
and norm d.

The actual compression in the aggregate criterion also contains a rank-one
term. With U=F/2 on T,

\[
 P_{R,R}=\frac12I+
 U^T\left(\frac{M_{T,T}}{8\sqrt p}-\frac{J_T}{8p}\right)U.
 \tag{5}
\]

Thus neither the omitted parameters nor the constant character correction
can be removed from an attempted extreme-eigenvalue argument by appealing
only to a limiting spectral distribution.

## 3. Six exact automorphisms and the sizes of their blocks

The transformations x -> 1-x and x -> 1/x preserve R and S_R.
For the second assertion use
chi(1/x-1/y)=chi(y-x)chi(xy)=chi(x-y) on R.
They generate the usual six fractional-linear transformations permuting
{0,1,infinity}, an action of S_3 on R. The action can be nonfaithful
in small fields; the following character calculation includes those cases.

Let m=|R|, t=1 if chi(2)=1 and 0 otherwise, and s=2 if p=1 mod 3
and 0 otherwise. The permutation character has values m,t,s on the
identity, a transposition, and a 3-cycle respectively. A transposition
is conjugate to x -> 1-x, whose sole possible fixed point is 1/2.
A 3-cycle has fixed points x^2-x+1=0. For p=1 mod 4, these belong
to R precisely when p=1 mod 3: their order is six and p=1 mod 12.

It follows by averaging the three irreducible S_3 characters that the
trivial, sign, and standard multiplicities are

\[
 \boxed{d_+=\frac{m+3t+2s}{6},\qquad
 d_-=\frac{m-3t+2s}{6},\qquad
 d_{\rm st}=\frac{m-s}{3}.}
 \tag{6}
\]

The real symmetric operator S_R therefore decomposes into blocks of
sizes d_+, d_-, d_st, with every eigenvalue of the standard block repeated
twice. This follows from the isotypic decomposition and the fact that all
three S_3 irreducibles are absolutely irreducible over the reals. Empty
blocks are allowed. The centering matrix J_R acts only on the trivial
isotypic component, so the same decomposition applies to (5).

This reduces the size of an exact spectral computation. It does not
bound the norm of any of the three blocks as p grows. In particular,
neither (3) nor (6) supplies the subexponential aggregate estimate required
by the preceding pass.

## Verification

The [verifier](../experiments/parallel3_conic_2026_09_04.py) checks every
conic kernel entry in 14 prime fields, both quotient identities, all six
automorphisms, and the integer numerators of the isotypic projectors.
For the five fields p<=41 it independently evaluates every entry in (3)
in Z[zeta_(p-1)] by polynomial reduction modulo the cyclotomic polynomial;
it also checks (4) and all claimed zero blocks. No floating-point acceptance
condition is used. The [results](../results/parallel3_conic_2026_09_04.json)
record counts, block sizes, witnesses, and input hashes.
