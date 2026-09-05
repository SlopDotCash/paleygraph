# Cartesian Möbius expansion and the centered Weil trace

Status: an exact group-theoretic reformulation, uniform subgroup escape,
and a source-dependent operator expansion estimate. A concrete Weil
average shows why a small operator norm alone does not control the
required trace. No new bilinear cancellation, clique bound, or prize
result is claimed. The representation and expansion inputs are classical.

Throughout p≡1 mod4 is prime, A,B⊂F_p are nonempty, m=|A|, n=|B|,
χ(0)=0, and S(A,B)=Σ_{a∈A,b∈B}χ(a−b).

## 1. An exact Cartesian family in SL₂

Put u_s=[[1,s],[0,1]], w=[[0,−1],[1,0]], and

\[
g(s,t)=u_swu_t=\begin{pmatrix}s&st-1\\1&t\end{pmatrix}.
\]

For sets U,V of sizes m,n, let G(U,V)={g(s,t):s∈U,t∈V}.
This parametrization is injective. Its action on the projective line is
g(s,t)z=s−1/(z+t), with g(s,t)∞=s and g(s,t)(−t)=∞.

For fixed x,y∈P¹(F_p),

\[
\#\{g\in G(U,V):gx=y\}\le\max(m,n).                 \tag{1}
\]

If x,y are finite, the equation is (s−y)(t+x)=1, which determines
either parameter from the other, giving the sharper bound min(m,n).
If x=∞ or y=∞, it fixes s or t. The pair x=y=∞ has no solution.
For x,y∈F_{p²}\F_p, the same finite equation gives min(m,n), and
there is no zero denominator. Rational Möbius transformations preserve
the set of nonrational projective points.

Consequently every coset of a Borel subgroup meets G(U,V) in at
most max(m,n) points. Every coset of a split or nonsplit torus
normalizer meets it in at most twice that number: such a normalizer
fixes or interchanges two points over F_{p²}. The same conclusion
holds for a union of two cyclic cosets whose cyclic subgroup lies
in a Borel or a nonsplit torus.

Using the proper-subgroup classification in
[Lyamkin, Section 1.5, Theorem 10 and Lemma 5](https://www.mathnet.ru/php/archive.phtml?jrnid=sm&option_lang=eng&paperid=9707&wshow=paper),
the remaining exceptional subgroups have order at most 120. Thus

\[
\boxed{\ |G(U,V)\cap hH|\le\max\{2\max(m,n),120\}
\quad(H<\mathrm{SL}_2(F_p)\text{ proper}).\ }          \tag{2}
\]

Primality excludes subfield subgroups. For m,n≥p^ε, the uniform
probability measure μ on G(U,V) therefore satisfies
μ(hH)≤p^{−ε/2} for every such coset once p is sufficiently large.

The same primary source's Theorem 11, with the Fourier operator norm
defined in Section 1.1, then gives some κ(ε)>0 such that

\[
\sup_{\pi\ne1}\left\|\frac1{mn}\sum_{s\in U,t\in V}
\pi(g(s,t))\right\|\ll_\varepsilon p^{-\kappa(\varepsilon)}.
\tag{3}
\]

The supremum is over nontrivial irreducible unitary representations.
This is an application of the existing expansion theorem, not a new
expansion theorem. No numerical value of κ is extracted. The source
is available in live primary HTML; direct archival downloads returned
403, so an archived copy of that source is not claimed.

There is also a direct generation observation. If either parameter
set contains two distinct elements, quotients with the other parameter
fixed give a nonidentity u_s. Since p is prime, its powers give all
u_t. One element of G(U,V) then gives w, and these generate SL₂(F_p).
Generation alone is weaker than the quantitative coset estimate (2).

## 2. The exact energy of the family

Let r_U(d)=#{(s,s')∈U²:s'−s=d} and E_+(U)=Σ_d r_U(d)².
The multiplicative energy E(G)=Σ_h#{(g,g')∈G²:g^{-1}g'=h}² is

\[
\boxed{E(G(U,V))=n^2E_+(U)+m^2E_+(V)-m^2n^2.}       \tag{4}
\]

Indeed, with d=s'−s,

\[
g(s,t)^{-1}g(s',t')=
\begin{pmatrix}1+td&t'-t+tt'd\\-d&1-dt'\end{pmatrix}.
\]

When d≠0, the matrix recovers d,t,t' uniquely and has multiplicity
r_U(d). When d=0 it is u_{t'−t}, with multiplicity m r_V(t'−t).
These cases are disjoint and give (4). In particular
E(G)≤m²n²(m+n−1). The formula retains the actual Cartesian structure;
it does not assume the additive energies are small.

## 3. A normalized Weil representation and its exceptional trace

Let e_p(z)=exp(2πiz/p), and use the p-dimensional unitary model

\[
(U_tf)(x)=e_p(tx^2)f(x),\quad
(Ff)(x)=p^{-1/2}\sum_y e_p(-2xy)f(y).
\]

For g=[[a,b],[c,d]], its kernel is

\[
\rho(g)_{x,y}=
\begin{cases}
\chi(c)p^{-1/2}e_p((ax^2-2xy+dy^2)/c),&c\ne0,\\
\chi(a)e_p(abx^2)\,1_{y=ax},&c=0.
\end{cases}                                          \tag{5}
\]

In particular ρ(u_t)=U_t and ρ(w)=F. These normalizations can be
checked directly. The quadratic Gauss identity is
Σ_z e_p(qz²+ℓz)=√p χ(q)e_p(−ℓ²/(4q)) for q≠0, with the
usual p or 0 answer for q=0. Right multiplication by u_1 multiplies
the kernel by e_p(y²). Right multiplication by w is evaluated by
that Gauss identity. For c≠0 its inner quadratic has coefficient
d/c; if d=0 the resulting delta function gives the c=0 branch of
ρ(gw), and otherwise the scalar is χ(d)/√p. The c=0 branch gives
the c≠0 branch of ρ(gw) directly. Since u_1,w generate SL₂(F_p),
these identities prove the representation law. Every operator is
unitary by its diagonal, permutation and Fourier factorization.

The classical trace formula is also recorded in
[Thomas, Section 2](https://arxiv.org/html/math/0610644v3#S2).
Here the quadratic phases correspond to additive character e_p(2z),
so its one-dimensional Weil index is 1. Directly from (5),

\[
\operatorname{tr}\rho(g)=
\begin{cases}
\chi(a+d-2),&a+d\ne2,\\
\sqrt p\,\chi(c),&a+d=2,\ c\ne0,\\
\sqrt p\,\chi(b),&a=d=1,\ c=0,\ b\ne0,\\
p,&g=I.
\end{cases}                                          \tag{6}
\]

Thus it is incorrect to use χ(tr g−2) on the unipotent locus.

Apply (5) to U=A+2,V=−B and put
T_{A,B}=(mn)^{-1}Σ_{a,b}ρ(g(a+2,−b)). Then

\[
\boxed{\operatorname{tr}T_{A,B}
=\frac{S(A,B)}{mn}+\frac{\sqrt p\,|A\cap B|}{mn}.}    \tag{7}
\]

The second term is exact. For A=B with n=m≪√p it can be much
larger than one, even though |S(A,B)|/(mn)≤1. For disjoint sets
the correction vanishes. Neither case allows dropping the term
without checking it.

## 4. What the actual operator estimate controls

For α_A(z)=m^{-1}Σ_{a∈A}e_p(az), the full kernel is

\[
(T_{A,B})_{x,y}=p^{-1/2}e_p(2x^2-2xy)
                 \alpha_A(x^2)\alpha_B(-y^2).         \tag{8}
\]

In particular (T_{A,B})_{0,0}=p^{−1/2}, so
‖T_{A,B}‖≥p^{−1/2}. The generic inequality
|tr T|≤p‖T‖ therefore cannot give a bound below one for this
family, irrespective of the strength of an upper estimate on its
operator norm. An estimate on a centered or otherwise modified
operator would require a separate argument.

There is a sharper identity retaining the factorization. Set
s_A=S(A,A), s_B=S(B,B). Quadratic Gauss summation gives

\[
\|T_{A,B}\|_F^2
=\frac1p\left(\frac pm+\frac{\sqrt p\,s_A}{m^2}\right)
          \left(\frac pn+\frac{\sqrt p\,s_B}{n^2}\right).
\tag{9}
\]

Applying Cauchy–Schwarz directly to the diagonal in (8) yields

\[
|S(A,B)+\sqrt p\,|A\cap B||^2
\le(\sqrt p\,m+s_A)(\sqrt p\,n+s_B).                 \tag{10}
\]

The quantities on the right are nonnegative. At m,n=o(√p),
using only |s_A|≤m² and |s_B|≤n² gives the familiar square-root
scale √(p/(mn)) for the normalized expression. This does not improve
the arbitrary-set threshold. Equations (7)–(10) locate the required
cancellation inside a centered trace, rather than establishing it.

## 5. A concrete trace-one average with norm tending to zero

This distinction persists even within the actual Weil representation.
Let p≥13, g₀=[[3,−1],[1,0]], and average its conjugates:

\[
V_p=\frac1{|\mathrm{SL}_2(F_p)|}
       \sum_h\rho(hg_0h^{-1}).
\]

Let Rf(x)=f(−x), with parity projections P_±=(I±R)/2.
Their dimensions are (p±1)/2. The commutant of this representation
is exactly span{I,R}. Here is an elementary verification. Commutation
with every U_t forces a matrix entry K_{x,y} to vanish unless
x²=y². Write its two possible entries on a nonzero row as a_x,b_x
and its (0,0) entry as d. Commutation with F at row or column zero
gives a_x+b_x=d and a_x=a_{−x}. At nonzero x,y, it then gives
(a_x−a_y)(e_p(−2xy)−e_p(2xy))=0. The second factor is nonzero,
so a_x is constant and K is a linear combination of I and R.
This also proves the two parity representations are irreducible and
inequivalent; for p≥13 both are nontrivial.

The conjugate average commutes with the whole representation.
By (6), tr ρ(g₀)=1 and tr(Rρ(g₀))=tr ρ(−g₀)=χ(5).
The two traces determine the average exactly:

\[
\boxed{V_p=\frac{1+\chi(5)}{p+1}P_+
             +\frac{1-\chi(5)}{p-1}P_-,\qquad
\operatorname{tr}V_p=1,\quad
\|V_p\|=\frac2{p+\chi(5)}.}                          \tag{11}
\]

It is a positive semidefinite scalar multiple of one parity projection.
For every j≥1, tr(V_p^j)=((p+χ(5))/2)^{1−j}. Even very small
higher-power traces coexist with the first trace being exactly one.
The standard dimension loss in recovering the first trace is sharp.

The conjugacy-class measure in (11) is **not asserted to be a
Cartesian family G(A+2,−B)**. It is therefore not a Paley
counterexample. It proves a specific limitation of deducing the
needed trace estimate from operator norm or higher-power trace
control alone, even for this genuine representation and positive
probability averages.

## 6. The remaining obligation

Subgroup escape and operator expansion are available for every
polynomial-sized Cartesian family above. They do not yet give

\[
\left|\operatorname{tr}T_{A,B}
      -\frac{\sqrt p\,|A\cap B|}{mn}\right|
\ll_\varepsilon p^{-\delta(\varepsilon)}.             \tag{12}
\]

By (7), (12) is exactly the classical Paley cancellation target in
this normalization. Establishing it needs more than the generic
operator argument tested here; the Cartesian structure and the
exceptional trace correction must be retained. It is not a new
independent sufficient theorem that has already been proved.

The [exact verifier](../experiments/parallel18_mobius_trace_2026_09_05.py)
checks the energy, transporter bounds over F_p and F_{p²}, the
normalized representation and trace formulas using cyclotomic
arithmetic, and complete conjugacy-class average matrices. Finite
checks do not establish the imported asymptotic expansion theorem
or the missing centered trace bound. Independent mathematical review
is outstanding. The full Paley and official prize goals remain open.
