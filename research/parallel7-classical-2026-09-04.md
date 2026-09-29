# Conditioning, robustness, and a countermodel to the variance route

**Status: the adversarial single-size estimate is still unproved.** This
pass derives an exact conditional fourth-moment average and proves a
specific obstruction to promoting the previous variance result to a
uniform bound. The countermodel below lives on the **actual single-size
slice** n=floor(p^(1/3)), for every sufficiently large prime. It is an
artificial statistic, not the fourth moment of another character and
not a counterexample to (SS), (CS), (LM), or Paley.

The countermodel has exactly the same mean and variance as M4, is
nonnegative and affine-invariant, is a polynomial of degree at most
four in the membership indicators, and has square-root bounds for its
squarefree quartic coefficients. Nevertheless its maximum divided by
pn^2 tends to infinity, and the large values persist throughout a
neighborhood of a planted set. Thus even this collection of additional
constraints does not repair the variance-to-maximum inference. The
missing input must use more of the character factorization or another
constraint not retained here. No novelty claim is made.

Throughout, p>=11 is prime, chi(0)=0,

\[
 F_B(x)=\sum_{b\in B}\chi(x-b),\quad
 M(B)=\sum_xF_B(x)^4,\quad A_B=\sum_{b\in B}F_B(b)^2.
\]

Expectations and variances without further qualification are over all
n-element subsets of F_p with equal probability. Write mu=E M and
sigma^2=Var M. The [preceding pass](parallel6-classical-2026-09-04.md)
proved their exact formulas. In the range n->infinity, n=o(p),

\[
 \mu=3pn^2(1+O(1/n+n/p)),\qquad
 \sigma^2=24pn^4(1+O(1/n+n/p)).                         \tag{1}
\]

## 1. An exact average conditioned on four prescribed elements

Fix a four-element set Q, and choose B uniformly among n-element
supersets of Q, where 4<=n<=p. Put

\[
 a_j=\frac{(n-4)_j}{(p-4)_j}\quad(1\le j\le4),
\]

and define

\[
\begin{aligned}
c_4&=1-4a_1+6a_2-4a_3+a_4,\\
c_A&=6a_1-18a_2+18a_3-6a_4,\\
c_2&=-4a_1+16a_2-20a_3+8a_4,\\
c_{tt}&=3a_2-6a_3+3a_4,\\
c_t&=a_1-7a_2+12a_3-6a_4.
\end{aligned}
\]

Then

\[
 \boxed{\mathbb E[M(B)\mid Q\subset B]
       =c_4M(Q)+c_A A_Q+C_{p,n},}                       \tag{2}
\]

where

\[
\begin{aligned}
C_{p,n}={}&[c_A(p-5)+c_2](4p-16)\\
 &+c_{tt}[p(p-5)^2+8(p-5)+4]
   +c_t[p(p-5)+4].
\end{aligned}
\]

In particular c_4=(p-n)_4/(p-4)_4. At n=4 this gives M(Q), and
at n=p it gives zero, as required.

**Proof.** Fix x and put s=F_Q(x), z=1_(x in Q), t=p-5+z.
The population chi(x-b), b outside Q, has first and third power
sums -s, and second and fourth power sums t. If T is its sum over
a uniformly sampled (n-4)-subset, expansion by multiplicity gives

\[
\begin{aligned}
\mathbb E T&=-a_1s,\\
\mathbb E T^2&=a_2s^2+(a_1-a_2)t,\\
\mathbb E T^3&=-a_3s^3-3(a_2-a_3)st
                 +(-a_1+3a_2-2a_3)s,\\
\mathbb E T^4&=a_4s^4+6(a_3-a_4)s^2t
 +(4a_2-12a_3+8a_4)s^2\\
 &\qquad +(3a_2-6a_3+3a_4)t^2
 +(a_1-7a_2+12a_3-6a_4)t.
\end{aligned}
\]

Consequently E(s+T)^4=c_4s^4+c_A s^2t+c_2s^2+c_tt t^2+c_t t.
The two-point character correlation gives sum_x s^2=4p-16;
sum_x s^2t=(p-5)(4p-16)+A_Q. Also sum_x t=p(p-5)+4 and
sum_x t^2=p(p-5)^2+8(p-5)+4. These identities prove (2), including
the corrections at the four zero coordinates.

Equation (2) computes a conditioned mean, not a uniform bound on
its members. It will let the verifier compute an exact covariance
without enumerating the entire large slice.

## 2. An affine-invariant planted quartic statistic

Now take n=floor(p^(1/3))>=4 and I={0,...,n-1} in F_p. Let H be
the family of all four-element sets which are affine images of some
four-element subset of I. Each edge occurs once, irrespective of how
many representations it has. Write h=|H| and

\[
 U(B)=|\{Q\in H:Q\subset B\}|,\quad
 b_j=\frac{(n)_j}{(p)_j},\quad
 Z(B)=12\sqrt p\,[U(B)-h b_4].                           \tag{3}
\]

This is a real polynomial of degree four on the slice; it is invariant
under all affine permutations of F_p and has mean zero. The coefficient
of each squarefree monomial of degree four is either 0 or 12 sqrt(p).
Equivalently, in the ordered distinct-quadruple convention for M,
the added coefficient is either 0 or sqrt(p)/2.

There are p(p-1) affine maps, so

\[
 h\le p(p-1)\binom n4\le p^2n^4/24.                    \tag{4}
\]

The following elementary variance bound is sufficient:

\[
 \boxed{\operatorname{Var}Z\le28n^6+24n^5
                  =o(pn^4).}                          \tag{5}
\]

**Proof of (5).** Let d(T) count edges of H containing T. Affine
two-transitivity gives

\[
 \sum_{|T|=1}d(T)^2=16h^2/p,\qquad
 \sum_{|T|=2}d(T)^2=36h^2/\binom p2.
\]

For triples, d(T)<=p-3 and sum d(T)=4h; for four-element sets the
sum of squares is h. Disjoint edges have nonpositive covariance
under sampling without replacement, because b_8<=b_4^2. Counting
overlapping ordered pairs, and allowing overcounting, therefore gives

\[
 \operatorname{Var}U\le
 b_7\frac{16h^2}{p}
 +b_6\frac{72h^2}{p(p-1)}
 +b_5\,4(p-3)h+b_4h.                                   \tag{6}
\]

Using b_j<=(n/p)^j, the sharper first inequality in (4) when needed,
and (n)_4<=n^4, the right side is at most

\[
 \frac{n^{15}}{36p^4}+\frac{n^{14}}{8p^4}
       +\frac{n^9}{6p^2}+\frac{n^8}{24p^2}.
\]

Multiply by 144p and use p>=n^3 to obtain (5). No independence of
overlapping edges is assumed.

On the other hand,

\[
 Z(I)=12\sqrt p[\binom n4-hb_4]
       =(1/2+o(1))\sqrt p\,n^4.                        \tag{7}
\]

Indeed U(I)=binom(n,4), and

\[
 \frac{hb_4}{\binom n4}
 \le\frac{(n)_4}{(p-2)(p-3)}=O(n^{-2}).
\]

Uniformly in B,

\[
 Z(B)\ge-12\sqrt p\,h b_4
      \ge-\frac{n^8}{2p^{3/2}}=-O(n^{7/2}).             \tag{8}
\]

Thus mu+Z is already nonnegative for all sufficiently large p, has
the correct mean and a stronger variance upper bound than M, but its
maximum is much too large for a constant multiple of pn^2.

## 3. Matching the exact variance too

The preceding construction can be calibrated to match the **exact**
variance of M without losing its other properties. This matters because
the earlier pass computed a variance identity, not only an upper bound.

Put F=M+Z and tau^2=Var F. The triangle inequality and reverse triangle
inequality for centered L2 norms give

\[
 |\tau-\sigma|\le\sqrt{\operatorname{Var}Z},\qquad
 \tau/\sigma=1+O(n^{-1/2}).
\]

In particular tau>0 for all sufficiently large primes. Define the
artificial statistic

\[
 \boxed{G(B)=\mu+\frac{\sigma}{\tau}\,[M(B)+Z(B)-\mu].}  \tag{9}
\]

It satisfies all of the following simultaneously on the precise SS
slice, for every sufficiently large prime:

1. E G=mu and Var G=sigma^2, exactly.
2. G is invariant under every affine permutation of F_p and is the
   restriction of a polynomial of degree at most four in 1_B(b).
3. G(B)>=0 for every n-element B.
4. Its squarefree quartic coefficient in ordered-tuple normalization
   is bounded in absolute value by 3 sqrt(p).
5. G(I)/(pn^2)>=(1/2-o(1))n^2/sqrt(p), which tends to infinity like
   p^(1/6).

For (3), the exact second moment and Cauchy-Schwarz give the uniform
lower bound

\[
 M(B)\ge\frac{[n(p-n)]^2}{p}=(1-o(1))pn^2.
\]

Equation (8) is o(pn^2), while sigma/tau=1+O(n^(-1/2)) and
mu=O(pn^2). Inserting these into (9) gives
G(B)>=(1-o(1))pn^2, uniformly. For (4), the ordered distinct-root
coefficient of M is

\[
 K(Q)=\sum_x\chi\Bigl(\prod_{q\in Q}(x-q)\Bigr),
\]

The [exact elliptic reduction](parallel2-classical-2026-09-04.md)
and the genus-one Weil bound give |K(Q)|<=2 sqrt(p)+1. Adding (3)
changes it by at most sqrt(p)/2. Since sigma/tau=1+O(n^(-1/2)),
the resulting bound is less than 3 sqrt(p) for sufficiently large p.
Thus the countermodel even satisfies the usual degree-four polynomial
Weil constant. No aggregate trace estimate is being imported.
For (5), use (7), nonnegativity of M, and mu=o(sqrt(p)n^4).

The construction is not claimed to satisfy all identities of character
sums. In particular, no character-factorization representation of G
is asserted; G is defined as a different statistic on subsets. Its point
is the precise logical separation: mean and variance identities, affine
symmetry, polynomial degree four, nonnegativity, and a constant-times-
sqrt(p) bound for each distinct-root coefficient do not imply (SS).
It does not rule out an argument using additional exact coefficient
relations, such as the common character-factorization identity.

## 4. Large values can persist on a whole neighborhood

This countermodel does not rely on an isolated set. Fix theta in (0,1).
For every n-set B with |B intersect I|>=theta n,

\[
 U(B)\ge\binom{\lceil\theta n\rceil}{4},\qquad
 G(B)\ge(\theta^4/2+o(1))\sqrt p\,n^4,                 \tag{10}
\]

uniformly in B. The same proof as (7)--(9) applies. The proportion of
such B is

\[
 \frac{\sum_{k=0}^{n-\lceil\theta n\rceil}
               \binom nk\binom{p-n}{k}}{\binom pn}
 =\exp[-(\theta+o(1))n\log(p/n)].                       \tag{11}
\]

To see (11), the summands increase in this range for large p, and
log binom(p,n)=n log(p/n)+O(n). The last numerator term has logarithm
(n-ceil(theta n)) log(p/n)+O(n). The remaining factor of at most n+1
terms does not change the displayed exponent. This proportion is
exponentially smaller than any inverse power of p. It is therefore
compatible with the exact variance, even though the large values
survive replacement of a fixed positive fraction of the planted set.

There is also an elementary robustness estimate for the **actual** M.
If B and C are n-sets with |B\C|=|C\B|=k, then

\[
 \boxed{|M(B)^{1/4}-M(C)^{1/4}|
             \le(8pk^3)^{1/4}.}                        \tag{12}
\]

Indeed F_B-F_C is a character convolution of a signed vector with
2k nonzero entries and total sum zero. Its squared L2 norm is exactly
2pk, while its absolute value is at most 2k. Hence its fourth moment
is at most 8pk^3, and the L4 triangle inequality proves (12). At a
putative value M(B)>=Lpn^2 this preserves at least one sixteenth of
that value for

\[
 k\le(Ln^2/128)^{1/3}.                                 \tag{13}
\]

At fixed L this certified radius is O(n^(2/3)), whose ball has
proportion exp[-(1+o(1))n log(p/n)]. Such a ball is far too small to
contradict the previous polynomial exceptional-probability bound.
Equation (10) shows that even replacing (13) by a fixed-fraction radius
would not resolve this counting problem alone.

## 5. Verification and remaining input

The companion script checks (2) by exact averaging over prescribed
small-field completions. It checks the exact L2 identity and the
fourth-power bound underlying (12) on bounded full families of pairs.
For fixed primes on the SS slice it constructs H exactly, verifies
its affine symmetry and uniform one- and two-point degrees, and
computes Var U by the exact intersection distribution of ordered
edge pairs. It obtains Cov(M,U) from (2), separately on the disjoint
affine orbits of H. Thus the calibration variance in (9) is represented
exactly as A+B sqrt(p), without enumerating all large n-subsets or
using floating-point evidence as proof.

These identities do not supply an upper estimate for adversarial
character moments. The actual missing task remains a uniform signed
quartic bound, or another structural input distinguishing genuine
character products from the countermodel. No hypothesis in (SS),
(CS), or (LM) has been established or refuted by this pass.

Final exact test counts and source hashes are recorded in
`results/parallel7_classical_2026_09_04.json`.
