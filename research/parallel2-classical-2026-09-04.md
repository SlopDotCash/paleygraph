# Fourth moments: an exact elliptic reduction and a failed L2 upper-bound route

This second bounded pass concerns upper-bound structure. It derives an
exact signed cross-ratio representation of the fourth moment and then
**disproves a natural sufficient L2 estimate**, including its centered
version. The failure occurs even for sets small enough that the ordinary
Weil argument already bounds their fourth moments. It therefore rules
out this particular proof route, not the logarithmic moment hypothesis
(LM) or the Paley conjecture. No claim of novelty in the literature is
made.

Throughout, \(p\ge5\) is prime, \(\chi(0)=0\), \(B\subseteq\mathbb F_p\),
and \(n=|B|\). Set

\[
F_B(x)=\sum_{b\in B}\chi(x-b),\quad
M_4(B)=\sum_xF_B(x)^4,\quad
A_B=\sum_{b\in B}F_B(b)^2.
\]

## 1. Exact cross-ratio coordinates

For \(\lambda\notin\{0,1\}\), define the signed Legendre elliptic trace

\[
E(\lambda)=\sum_{x\in\mathbb F_p}\chi(x(x-1)(x-\lambda)).
\]

For an ordered quadruple of distinct elements \((a,b,c,d)\), put

\[
\lambda(a,b,c,d)=
 \frac{(c-b)(d-a)}{(c-a)(d-b)},\qquad
\epsilon(a,b,c,d)=\chi((c-a)(d-b)).
\]

Then

\[
\boxed{
\sum_x\chi((x-a)(x-b)(x-c)(x-d))
 =\epsilon(a,b,c,d)E(\lambda(a,b,c,d))-1.}
\tag{1}
\]

**Proof, including the missing point.** Write
\(h=(c-b)/(c-a)\) and \(u=h(x-a)/(x-b)\). This fractional linear
map sends \(a,b,c,d\) to \(0,\infty,1,\lambda\). For \(x\ne b\),

\[
u(u-1)(u-\lambda)
=h(h-1)(h-\lambda)
 \frac{(x-a)(x-c)(x-d)}{(x-b)^3}.
\]

The character of the prefactor is
\(\chi(h(h-1)(h-\lambda))=\epsilon(a,b,c,d)\). The image of
\(\mathbb F_p\setminus\{b\}\) omits the finite point \(u=h\).
Its contribution, multiplied by \(\epsilon\), is exactly one. The
original term at \(x=b\) is zero, proving the subtraction in (1).

Let the signed coefficient vector be

\[
V_B(\lambda)=
\sum_{\substack{a,b,c,d\in B\ \text{distinct}\\
                \lambda(a,b,c,d)=\lambda}}
\epsilon(a,b,c,d).
\]

Combining (1) with the multiplicity partition from
[the earlier fourth-moment identity](moments-and-obstructions.md#3-an-exact-fourth-moment-decomposition)
gives

\[
\boxed{M_4(B)=p(3n^2-2n)-n^4+3n^2-3n-6A_B
                +\langle V_B,E\rangle.}
\tag{2}
\]

Both vectors in the inner product are indexed by
\(\mathbb F_p\setminus\{0,1\}\). This is a signed identity;
absolute values have not been applied to individual elliptic traces.

## 2. Exact centering and a sufficient upper-bound input

The elliptic vector has the exact squared norm

\[
\boxed{\|E\|_2^2=p^2-2p-3=(p-3)(p+1).}
\tag{3}
\]

Indeed, set \(f(x)=\chi(x(x-1))\). Then
\(E(\lambda)=\sum_xf(x)\chi(x-\lambda)\),
\(\sum f=-1\), and \(\sum f^2=p-2\). The two-point character
correlation gives \(\sum_{\lambda\in\mathbb F_p}E(\lambda)^2=p(p-2)-1\).
The omitted terms are \(E(0)=-\chi(-1)\) and \(E(1)=-1\), each
with square one.

There is also an exact full-set formula:

\[
\boxed{V_{\mathbb F_p}(\lambda)=p(p-1)E(\lambda).}
\tag{4}
\]

To prove it, fix an ordered pair \((a,b)\) and use affine invariance
to take \((a,b)=(0,1)\). The remaining quadruples of fixed cross-ratio
are parametrized by \(h\notin\{0,1,\lambda\}\), with
\(c=1/(1-h)\), \(d=\lambda/(\lambda-h)\). Their signs are
\(\chi(h(h-1)(h-\lambda))\), and summing these signs is precisely
\(E(\lambda)\). There are \(p(p-1)\) ordered choices of \((a,b)\).

Write \((n)_4=n(n-1)(n-2)(n-3)\), which is zero for
\(n\in\{0,1,2,3\}\), and define

\[
c_{p,n}=\frac{(n)_4}{(p-2)(p-3)},\qquad
W_B=V_B-c_{p,n}E.
\]

This is the exact random-set centering: if \(B\) is a uniformly chosen
\(n\)-element subset, every distinct ordered quadruple is included with
probability \((n)_4/(p)_4\), so (4) gives
\(\mathbb E[V_B]=c_{p,n}E\).

Equations (2)–(3) now give the exact identity

\[
M_4(B)=R_{p,n,B}+\langle W_B,E\rangle,
\tag{5}
\]

where

\[
R_{p,n,B}=
p(3n^2-2n)-6n^3+14n^2-9n-6A_B
 +\frac{3(n)_4}{p-2}.
\]

For every \(0\le n\le p\), \(R_{p,n,B}\le3pn^2\). For \(n\ge3\),
use \((n)_4/(p-2)\le n^3\), \(p\ge n\), and \(A_B\ge0\);
the correction to \(3pn^2\) is at most
\(-3n(n-1)(n-3)\le0\). The cases \(n=0,1,2\) follow directly.
Thus Cauchy–Schwarz yields

\[
\boxed{M_4(B)\le3pn^2+p\|W_B\|_2.}
\tag{6}
\]

In particular, the following input would suffice for the fourth-moment
instance of (LM):

\[
\|W_B\|_2\le C\left(n^2+\frac{n^4\log p}{p}\right)
\quad\text{uniformly in }p,B.
\tag{CS-input}
\]

This is an additional sufficient hypothesis suggested by applying
Cauchy–Schwarz to (5). It is not claimed equivalent to (LM). The next
section proves it false for every constant \(C\).

## 3. Affine copies disprove CS-input

Take odd \(n=2m-1\), with \(m\ge6\), and let
\(B=\{1,\ldots,n\}\subseteq\mathbb F_p\). Choose a prime
\(p\ge n^8\) such that every nonzero integer of absolute value less
than \(n\) is a quadratic residue modulo \(p\). Arbitrarily large
such primes exist, as justified below. Every sign in \(V_B\) is then
positive. In particular, \(V_B\) counts ordered cross-ratios without
cancellation.

There are \((m-1)(m-2)(m-3)\) ordered quadruples
\((1,b,c,d)\) with distinct \(b,c,d\in\{2,\ldots,m\}\). Translate
each of these quadruples by each \(s\in\{0,\ldots,m-1\}\). All the
resulting quadruples lie in \(B\), and the \(m\) translates of each
original quadruple have identical cross-ratio and sign. All generated
quadruples are distinct: their first coordinate identifies the
translation, and then the remaining coordinates identify the original.

If several original quadruples have the same cross-ratio, their
contributions add positively. Therefore

\[
\boxed{\|V_B\|_2^2\ge
 m^2(m-1)(m-2)(m-3)\ge n^5/256.}
\tag{7}
\]

For the last inequality, each \(m-i\ge m/2\), \(1\le i\le3\),
and \(m\ge n/2\). There is no probabilistic assumption about
cross-ratio distribution in this counting argument.

Meanwhile, for \(p\ge11\),

\[
c_{p,n}\|E\|_2\le \frac{2n^4}{p}.
\]

The reverse triangle inequality and (7) imply

\[
\|W_B\|_2\ge n^{5/2}/16-2n^4/p
             =\Omega(n^{5/2}).
\tag{8}
\]

But for \(p\ge n^8\), monotonicity of \(\log p/p\) gives
\(n^4\log p/p\le8\log n/n^4\le n^2\) for these \(n\).
Thus the right side of (CS-input) is at most \(2Cn^2\), whereas (8)
is larger for sufficiently large \(n\). This disproves (CS-input) for
every constant \(C\), including after the exact random-set centering.

### The required prime family exists

Put \(L=8\prod_{\ell\le n,\ \ell\ \mathrm{odd\ prime}}\ell\).
For \(p\equiv1\pmod L\), quadratic reciprocity gives
\(\chi(\ell)=1\) for each odd prime \(\ell\le n\), while
\(p\equiv1\pmod8\) gives \(\chi(2)=\chi(-1)=1\). Multiplicativity
then gives the required signs for all differences in the interval.

One does not need a quantitative prime-in-progressions theorem here.
For completeness, primes congruent to one modulo any fixed \(L>1\)
are arbitrarily large by an elementary cyclotomic argument. Given a
bound \(Y\), choose an integer \(x\ge3\) divisible by \(L\) and all
primes at most \(Y\), and take a prime divisor \(q\) of \(\Phi_L(x)\).
Because \(\Phi_L(0)=1\), this divisor does not divide \(x\), is greater
than \(Y\), and does not divide \(L\). Modulo \(q\), the polynomial
\(X^L-1\) is separable. Its cyclotomic factorization shows that a root
of \(\Phi_L\) has order exactly \(L\): a smaller order would put it
in another cyclotomic factor and produce a repeated root. Hence
\(L\mid q-1\). Also \(\Phi_L(x)>1\), since its complex roots have
absolute value one and \(x\ge3\), so a prime divisor exists.

For each \(n\), choose the resulting prime larger than \(n^8\).
This yields the unbounded family required to refute a uniform constant;
a finite numerical violation would not suffice.

## 4. What the failure means

The exact elliptic reduction retains useful information, but bounding its
signed coefficient vector by its complete L2 norm loses too much. The
large norm in (7) comes from many affine copies of the same configuration,
and does not imply those configurations align with the elliptic trace
vector. Even subtracting the exact random-set mean does not remove it.

This failure is separate from the prior forced-row obstruction. Indeed,
the standard termwise Weil estimate already gives

\[
M_4(B)\le3pn^2+3\sqrt p\,n^4\le6pn^2
\quad\text{whenever }n\le p^{1/4}.
\]

Our counterexample family satisfies the stronger restriction
\(n\le p^{1/8}\). Thus (LM) holds in this very range while
(CS-input) fails. The sufficient input is simply too strong. Any useful
continuation of (5) must preserve directional cancellation in
\(\langle W_B,E\rangle\), or handle affine configuration multiplicity
before taking absolute values; the uniform norm estimate above cannot
be its missing lemma.

No new worst-case small-set upper bound is proved here. This pass removes
a concrete candidate lemma from consideration and provides exact formulas
for a more selective investigation of the signed aggregate.

## 5. Verification

`experiments/parallel2_classical_2026_09_04.py` checks the fractional-linear
identity, elliptic norm, full-set vector, both fourth-moment formulas,
the centering, and the Cauchy–Schwarz inequality in exact integer/rational
arithmetic. It also counts interval cross-ratio collisions and checks
the zero/square conditions at explicit larger primes. Those latter
examples illustrate the collision mechanism; their field sizes are not
claimed to satisfy the asymptotic choice \(p\ge n^8\).

The result file records final source hashes. The asymptotic refutation
rests on the proof in Section 3, not on extrapolation from computations.
