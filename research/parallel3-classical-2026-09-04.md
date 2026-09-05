# The directional elliptic pairing: Fourier duality and another failed norm route

This pass derives an exact additive Fourier transform of the elliptic
trace vector in terms of ordinary Kloosterman sums. It also proves that
a Fourier-L1 sufficient input fails, **even after any bounded number of
cross-ratio coordinates are deleted**. The counterexamples lie inside a
range where the usual Weil argument already proves the fourth-moment
bound. Thus this is a precise obstruction to a proof strategy, not a
counterexample to Paley or to (LM). The prime construction has no
polynomial upper bound in \(n\): this refutes the input quantified over
**all** \(p,B\), and does not refute a version restricted to
\(n\ge p^\varepsilon\) for a fixed positive \(\varepsilon\).
That restricted possibility remains relevant to Paley. No literature
novelty is claimed.

Use the notation and proved formulas of
[the preceding elliptic reduction](parallel2-classical-2026-09-04.md):
\(V_B(\lambda)\) is the signed count of distinct ordered quadruples
having cross-ratio \(\lambda\),

\[
c=\frac{(n)_4}{(p-2)(p-3)},\qquad W_B(\lambda)=V_B(\lambda)-cE(\lambda),
\quad\lambda\notin\{0,1\}.
\]

Extend \(W_B\) by zero at 0 and 1. In this note \(E\) denotes the
**full** vector \(E(\lambda)=\sum_x\chi(x(x-1)(x-\lambda))\), including
its two singular parameters. The earlier exact formula gives

\[
M_4(B)=R_{p,n,B}+\langle W_B,E\rangle,
\qquad R_{p,n,B}\le3pn^2.
\tag{0}
\]

The extra singular entries do not change this pairing, since \(W_B\)
vanishes there.

## 1. Exact Fourier transform into Kloosterman sums

Write \(e_p(z)=\exp(2\pi iz/p)\), and use the unnormalized transform
\(\widehat f(t)=\sum_xf(x)e_p(-tx)\). Define

\[
\tau=\sum_x\chi(x)e_p(x),\qquad
\mathrm{Kl}(a)=\sum_{u\ne0}e_p(u+a/u).
\]

Then

\[
\boxed{\widehat E(0)=0,\qquad
\widehat E(t)=\tau\chi(t)e_p(-t/2)\mathrm{Kl}(t^2/16)
\quad(t\ne0).}
\tag{1}
\]

**Proof.** Put \(f(x)=\chi(x(x-1))\). Interchanging sums gives
\(\widehat E(t)=\tau\chi(t)\widehat f(t)\) for nonzero \(t\), and
zero at \(t=0\). To compute \(\widehat f(t)\), use

\[
1+f(x)=|\{y:y^2=x(x-1)\}|.
\]

The affine curve is parametrized bijectively by \(u\in\mathbb F_p^*\):
\(u=2x-1+2y\), \(u^{-1}=2x-1-2y\). Consequently

\[
\widehat f(t)=
 \sum_{u\ne0}e_p\left(-\frac t4(u+u^{-1}+2)\right)
=e_p(-t/2)\mathrm{Kl}(t^2/16).
\]

Here the Fourier transform of the constant term vanishes because
\(t\ne0\). The substitution in the last equality is \(v=-tu/4\).
This proves (1), including all phases and signs.

The Kloosterman sum is real, by replacing \(u\) with \(-u\) in its
complex conjugate. Parseval therefore gives the exact directional formula

\[
\boxed{\langle W_B,E\rangle=
\frac{\overline\tau}{p}
\sum_{t\ne0}\chi(t)e_p(t/2)\mathrm{Kl}(t^2/16)\widehat W_B(t).}
\tag{2}
\]

This is an explicit oscillatory expression, not an independent estimate
of its magnitude.

The classical bound \(|\mathrm{Kl}(a)|\le2\sqrt p\), with this
unnormalized convention, is the \(n=2\) case of the bound stated in
Section 2 of [Lu–Zheng–Zheng, *On the distribution of Jacobi sums*](https://arxiv.org/html/1305.3405v3#S2).
That bound is imported, not proved here. Since \(|\tau|=\sqrt p\),
(0)–(2) imply

\[
M_4(B)\le3pn^2+2\sum_{t\ne0}|\widehat W_B(t)|.
\tag{3}
\]

Thus a uniform estimate

\[
\sum_{t\ne0}|\widehat W_B(t)|
 \le C\bigl(pn^2+n^4\log p\bigr)
\tag{FL1}
\]

would suffice for (LM) at fourth-moment depth. Unlike the previously
refuted coefficient-L2 input, this norm can sometimes be smaller than
its Cauchy–Schwarz upper bound \(p\|W_B\|_2\). Nevertheless, (FL1)
is also false; Sections 3–4 prove an unbounded family.

## 2. Two averaging operations do not supply cancellation

The order-six cross-ratio symmetry has matching signs on both vectors:

\[
\begin{aligned}
E(1-\lambda)&=\chi(-1)E(\lambda),&
V_B(1-\lambda)&=\chi(-1)V_B(\lambda),\\
E(1/\lambda)&=\chi(\lambda)E(\lambda),&
V_B(1/\lambda)&=\chi(\lambda)V_B(\lambda).
\end{aligned}
\tag{4}
\]

For \(E\), substitute \(x=1-y\) and \(x=y/\lambda\), respectively.
For \(V_B\), permute the ordered quadruple by swapping its first and
third entries, or its third and fourth entries, respectively. Those
permutations induce the displayed transformations and sign factors.
The same laws hold for \(W_B\). Therefore

\[
W_B(1-\lambda)E(1-\lambda)=W_B(\lambda)E(\lambda)
=W_B(1/\lambda)E(1/\lambda).
\]

The product is constant on each orbit. Averaging over those orbits
cannot create cancellation; exceptional shorter orbits simply have
fewer identical summands. Likewise, for any \(u\ne0,v\),
\(V_{uB+v}=V_B\), since cross-ratios are unchanged and the sign gains
the square factor \(u^2\). Averaging affine images of \(B\) leaves
the entire coefficient vector and pairing unchanged.

Nor can every nonzero Fourier coefficient of \(E\) be bounded at the
smaller \(\sqrt p\) scale. The full norm identity is

\[
\sum_\lambda E(\lambda)^2=p^2-2p-1.
\]

Together with (1), Parseval gives

\[
\sum_{t\ne0}\mathrm{Kl}(t^2/16)^2=p^2-2p-1,
\quad
\max_{t\ne0}|\widehat E(t)|^2
 \ge\frac{p(p^2-2p-1)}{p-1}.
\tag{5}
\]

Thus the \(p\)-sized Fourier scale in (3) is necessary. A uniform
improvement there cannot replace cancellation in the weighted sum.

## 3. Intervals have many quadruples at each of infinitely many fixed ratios

Fix a positive integer \(h\). There are
\(\Omega_h(n^2\log n)\) distinct ordered integer quadruples from
\(\{1,\ldots,n\}\) whose cross-ratio equals \(-h^2\).

Here is a quantified construction. Set
\(N=\lfloor n/4\rfloor\), \(C_h=h(h+1)\), and
\(Q=\lfloor N/C_h\rfloor\). Choose coprime integers \(1\le r<s\),
a positive integer \(\ell\) with \(\ell s^2\le Q\), and exclude
\(s=hr\). Take the following four offsets:

\[
a_0=0,\quad b_0=(h^2+1)\ell rs,\quad
c_0=\ell r(s-hr),\quad d_0=h\ell s(hr-s).
\tag{6}
\]

All offsets have absolute value at most \(C_h\ell s^2\le N\),
and they are pairwise distinct. Their cross-ratio is

\[
\frac{(c_0-b_0)d_0}{c_0(d_0-b_0)}=-h^2.
\]

Translate them by each \(a\in\{N+1,\ldots,n-N\}\). The resulting
quadruples lie in the interval and are all distinct. The first entry
recovers \(a\); the ratio

\[
\frac{c_0}{b_0}=\frac{1-h r/s}{h^2+1}
\]

recovers the reduced fraction \(r/s\), and then \(b_0\) recovers
\(\ell\).

Let \(\varphi\) be Euler's totient. The number of offset patterns is
exactly

\[
T_h(Q)=\sum_{s=2}^{\lfloor\sqrt Q\rfloor}
 \varphi(s)\left\lfloor\frac Q{s^2}\right\rfloor
-\begin{cases}\lfloor Q/h^2\rfloor,&h\ge2,\\0,&h=1.\end{cases}
\tag{7}
\]

The subtraction removes precisely the coprime pair \((r,s)=(1,h)\).
With \(\kappa=\prod_{\ell\ \mathrm{prime}}(1-\ell^{-2})>0\),

\[
\sum_{s\le R}\frac{\varphi(s)}{s^2}=\kappa\log R+O(1).
\]

For completeness, expand \(\varphi(s)/s=\sum_{d\mid s}\mu(d)/d\)
and interchange sums. The result is
\(\sum_{d\le R}\mu(d)d^{-2}H_{\lfloor R/d\rfloor}\).
Use \(H_m=\log m+O(1)\), absolute convergence of
\(\sum d^{-2}\log d\), and
\(\sum_{d>R}d^{-2}=O(R^{-1})\). The leading coefficient is
\(\sum\mu(d)/d^2=\kappa\); positivity follows from convergence of
the Euler product. Floor errors in (7) total \(O(Q)\), so

\[
T_h(Q)=\frac\kappa2Q\log Q+O_h(Q).
\]

Consequently the interval cross-ratio count obeys the explicit lower bound

\[
\boxed{N_B(-h^2)\ge(n-2N)T_h(Q)
=\left(\frac{\kappa}{16h(h+1)}+o_h(1)\right)n^2\log n.}
\tag{8}
\]

This is an integer configuration count, with no random-model assumption.

## 4. FL1 fails even after any fixed number of coordinate deletions

Fix \(K\ge0\), and take \(h=1,\ldots,K+1\) in (8). For each
sufficiently large \(n\), choose a prime \(p\ge n^8\), larger than
all the fixed ratios, for which all nonzero interval differences are
quadratic residues. Such arbitrarily large primes were constructed in
[the previous pass](parallel2-classical-2026-09-04.md#the-required-prime-family-exists)
using quadratic reciprocity and a cyclotomic Euclid argument. No new
prime-distribution assumption is imposed here.

For \(B=\{1,\ldots,n\}\subseteq\mathbb F_p\), all signs in
\(V_B\) are positive. The rational configurations above remain valid
modulo \(p\), their coordinates are distinct, and the \(K+1\) values
\(-h^2\) are distinct and avoid 0 and 1. By (8), each has

\[
V_B(-h^2)\ge c_Kn^2\log n
\]

for some positive \(c_K\) and all sufficiently large \(n\). Centering
changes each coordinate by at most \(cp\le2n^4/p\), using the trivial
bound \(|E(\lambda)|\le p\). Thus
\(W_B(-h^2)\ge(c_K+o(1))n^2\log n\).

Let \(R\) be **any** set of at most \(K\) cross-ratio coordinates,
and put \(Z=W_B\mathbf1_{\mathbb F_p\setminus R}\). At least one
of those \(K+1\) large coordinates, say \(\lambda_0\), survives.
Moreover,

\[
|\widehat Z(0)|\le\|Z\|_1
 \le\|V_B\|_1+c\|E\|_1\le3n^4
\]

for these primes: \(\|V_B\|_1=(n)_4\),
\(\|E\|_1\le p^2\), and \(c\le2n^4/p^2\).
Fourier inversion at the surviving coordinate gives

\[
\boxed{\sum_{t\ne0}|\widehat Z(t)|
 \ge p|Z(\lambda_0)|-|\widehat Z(0)|
 =\Omega_K(pn^2\log n).}
\tag{9}
\]

On the other hand,
\(pn^2+n^4\log p\le2pn^2\) when \(p\ge n^8\) and \(n\) is
sufficiently large. Thus no constant can make (FL1) hold uniformly,
even if one may first delete any \(K\) exceptional coordinates chosen
separately for each \(B\).

**Size-range limitation.** The prime construction supplies arbitrarily
large \(p\ge n^8\), but supplies no upper bound \(p\le n^C\) for any
fixed \(C\). Accordingly this argument does not put the examples in
any fixed range \(n\ge p^\varepsilon\). It disproves the all-prime,
all-set version of (FL1), including bounded coordinate deletions. A
version explicitly restricted to polynomial-sized sets is not refuted
here and could still be useful for Paley. Restricting the sufficient
input's range is a different repair from deleting coordinates.

This strengthened version matters because deleting finitely many
coordinates would otherwise be inexpensive for the moment estimate.
For \(p\ge11\) and \(n\le\sqrt p\), the contribution from any
\(K\) coordinates has absolute value at most \(3Kpn^2\): use
\(|V_B(\lambda)|\le n^3\), \(|E(\lambda)|\le2\sqrt p\), and
\(c\le2n^4/p^2\), to bound each summand by
\(2\sqrt p\,n^3+8n^4/p\le(2+8/p)pn^2\le3pn^2\).
In particular the weakened sufficient estimate would have been

\[
M_4(B)\le(3+3K)pn^2+2\sum_{t\ne0}|\widehat Z(t)|.
\]

Equation (9) shows that such a finite correction still cannot repair
the proposed Fourier norm argument.

There is no contradiction to (LM). Our counterexamples satisfy
\(n\le p^{1/8}\), and the ordinary termwise Weil estimate already
gives \(M_4(B)\le6pn^2\) throughout the larger range
\(n\le p^{1/4}\). The missing cancellation concerns the **direction**
of the pairing; both the coefficient-L2 and the present Fourier-L1
majorants demand more than the actual moment bound requires.

## 5. What remains usable

Equation (2) provides a concrete Kloosterman-weighted representation for
further analysis. Equations (4)–(5) rule out free cancellation from the
sixfold symmetry or a smaller uniform Fourier scale. Equations (8)–(9)
show that even finitely many exceptional configurations do not fix the
Fourier-L1 approach. A successful estimate would have to keep the
arithmetic relation between the Kloosterman phases and the coefficient
vector, treat a growing family of concentrated configurations, or prove
a norm bound only in a size range that excludes these examples.

No new uniform Paley bound for small arbitrary sets is claimed.

The accompanying script checks Fourier identities exactly in
\(\mathbb Z[z]/(1+z+\cdots+z^{p-1})\), not by floating-point evaluation.
It also checks 900 symmetry instances, 54 directional pairings and
affine invariances, 106 nonzero-frequency transforms, and 20 integer
configuration records (12 with an independent complete count). The
directional identities are checked after clearing denominators.
