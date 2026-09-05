# A growing quartic class with quadratic mixed energy

**Status: this pass proves an unconditional growing class, not a
worst-case bound or a Paley proof.** As the dyadic order N tends to
infinity, a proportion at least
\(1-\log(2)/36=0.9807459116\ldots\) of its eligible quartic primes
have exactly the intrinsic fourth energy at **every** dyadic level.
A density-one class has nearly intrinsic fourth energy at every level.
The prime-counting input below has been checked for powers of two; no
unproved assertion about primes in progressions is assumed. No novelty
or formalization claim is made.

## Statements and notation

All logarithms are natural. Let N=2^m≥4 and

\[
 \mathcal P_N=\{p\text{ prime}:N^4/4\le p\le N^4,\quad p\equiv1\pmod N\}.
\]

For a dyadic divisor s≥2 of N let H_s be the order-s subgroup of F_p*,
and set \(D_s(p)=E_2(H_s)-(3s^2-3s)\ge0\). For s=2k≥4 write
K=H_k, L_s=H_s\setminus K, and

\[
 B_s=\|1_K*1_{L_s}\|_2^2,\qquad
 T_s=\#\{(a,b,c,d)\in K^3\times L_s:a+b+c+d=0\}.
\]

These are nonnegative integer counts, with B_s≥k². The collision mass
of the primitive-half polynomial R_s is

\[
 C_s=\sum_i\binom{c_{s,i}}2
     =\frac{B_s-k^2}{4k}=\frac{B_s-k^2}{2s}\in\mathbb Z_{\ge0},
 \qquad \mathcal C_N=\sum_{4\le s\mid N}C_s.                 \tag{1}
\]

Here c_(s,i) are its distinct-root multiplicities after reduction, as
proved in [the primitive-fiber note](parallel2-subgroup-2026-09-04.md).
The displayed sums over divisors always mean dyadic divisors.

**Theorem 1 (density one, all levels).** For every real η>0, outside at
most

\[
 \boxed{\frac{(4N^2-4)\log2}{3\eta\log(N^4/4)}}              \tag{2}
\]

primes in \(\mathcal P_N\), one has simultaneously

\[
 \begin{split}
 &\sum_{2\le s\mid N}\frac{D_s(p)}{s^2}\le\eta,
 \qquad E_2(H_s)\le3s^2-3s+\eta s^2,\\
 &B_s\le(1+2\eta/3)k^2,\qquad
 \sum_{4\le s\mid N}\frac{C_s}{s}\le\eta/12,\qquad
 \mathcal C_N\le\eta N/12.                                \tag{3}
 \end{split}
\]

Moreover

\[
 |\mathcal P_N|\sim\frac{3N^3}{8\log N}.                    \tag{4}
\]

Consequently η=N^(-1/2) gives (3) on a proportion 1−O(N^(-1/2)) of
eligible primes, uniformly in the whole tower. More generally, for any
positive w_N→∞, η=w_N/N gives (3) outside a proportion O(1/w_N), with
total primitive collision mass at most w_N/12. Choosing w_N=log N gives
a density-one class with \(\mathcal C_N\le\log N/12\), and no
primitive collisions at any level s<12N/log N. The last assertion follows
from C_s≤s log N/(12N)<1 and integrality.

**Theorem 2 (an exact growing class).** Let \(\mathcal G_N\) consist of
primes in \(\mathcal P_N\) for which, at every dyadic level,

\[
 E_2(H_s)=3s^2-3s,\qquad B_s=k^2,\qquad T_s=0.              \tag{5}
\]

Then

\[
 \boxed{\liminf_{N\to\infty,\ N\text{ dyadic}}
      \frac{|\mathcal G_N|}{|\mathcal P_N|}
       \ge1-\frac{\log2}{36}=0.9807459116\ldots .}          \tag{6}
\]

In particular this class contains quartic primes for every sufficiently
large dyadic order. Every R_s is squarefree there. This supplies a
proved growing infinite class for a condition that the earlier
primitive-fiber pass had only checked finitely in this window. It does
not assert that the complementary set is nonempty, or bound its worst
member.

## The exact recurrence and its telescoping consequence

Splitting an ordered zero-sum quadruple by the number of its entries in
L_s gives

\[
 E_2(H_s)=2E_2(H_k)+6B_s+8T_s,\qquad
 D_s=2D_k+6(B_s-k^2)+8T_s.                                \tag{7}
\]

Multiplication by an element of L_s exchanges the cosets, explaining
the matching pure and 3+1 patterns. Since −1 belongs to K for k≥2,
the 2+2 count equals the displayed mixed additive energy. Also
B_s≥k² follows by summing r(x)²≥r(x) over the integer representation
counts r=1_K*1_(L_s). The intrinsic paired quadruples give D_k≥0.

Using (1), dividing (7) by s, and starting from D_2=0 yields

\[
 \boxed{D_N=12N\mathcal C_N+
                8\sum_{4\le s\mid N}\frac Ns T_s.}         \tag{8}
\]

Thus D_s≥12s C_s and D_s≥6(B_s−k²). Formula (8) also gives the
standalone average estimate

\[
 \frac1{|\mathcal P_N|}\sum_{p\in\mathcal P_N}\mathcal C_N(p)
       \le\frac{\log2}{18}+o(1),                          \tag{9}
\]

once (4) and the norm budget below are used. This is an average of a
nonnegative integer, not a pointwise upper bound. It already gives a
proportion at least 1−log(2)/18−o(1) with B_s=k² at all levels;
Theorem 2 proves the stronger exact-energy statement (6).

## The norm budget used here

For a dyadic order s≥2 set A_s=s^4−3s²+3s. The elementary norm
argument in [cyclotomic-prime-average.md](cyclotomic-prime-average.md)
gives the unconditional estimate

\[
 \sum_{p\equiv1\ (s)}D_s(p)\log p
 \le\frac{A_s}{2}\log\frac{4s^4}{A_s}
 \le s^4\log2.                                           \tag{10}
\]

For clarity, the argument does not assume a bound for any one prime.
Normalize a zero-sum quadruple to first entry 1 and put
f_a(X)=1+X^(a_2)+X^(a_3)+X^(a_4). In the degree-d field with
d=s/2 and Φ_s(X)=X^d+1, discard the intrinsic polynomials f_a(ζ)=0.
There are R=A_s/s remaining exponent triples. The determinant of
multiplication by f_a is a nonzero integer. At a splitting prime its
valuation is at least the number of primitive roots g for which f_a(g)=0,
by the nullity bound for an integer matrix modulo p.

At each complex embedding, orthogonality gives
\(\sum_a|f_a(\zeta)|^2=4s^3\), and intrinsic terms contribute zero.
AM–GM over the R nonzero values, followed by multiplication over the d
embeddings, bounds the product of absolute determinants by
\((4s^3/R)^{dR/2}\). Summing the primitive-root nullities over exponent
triples gives \(dD_s(p)/s\): every choice of primitive root enumerates
the same normalized quadruple count. Dividing the resulting logarithmic
inequality by d/s proves the first inequality in (10).

For the second, the function y↦y log(4/y) increases on 0<y≤1,
because its derivative is log(4/y)−1>0. Substitute y=A_s/s^4≤1.
Nonnegativity of D_s is essential. Only finitely many primes contribute
for each fixed s, since they divide a finite product of nonzero norms.

Sum (10), divided by s², over all dyadic s≤N and restrict primes to
\(\mathcal P_N\). Every prime in that set splits at every smaller
level and satisfies log p≥log(N^4/4). Since
\(\sum_{2\le s\mid N}s^2=(4N^2-4)/3\), Markov's inequality proves
(2). The rest of (3) follows from (7) and D_s≥12s C_s. This proof
works for every η>0, including η depending on N; its constants do not
depend on η. Likewise (8), (10), and log p≥4log N−log4 imply (9)
using (4).

## Why there are enough quartic primes for powers of two

The external input is Thorner–Zaman,
[Refinements to the prime number theorem for arithmetic progressions,
Corollary 3.1 and equation (3.2)](https://arxiv.org/html/2108.10878#S3.SS1).
It covers q≥2 with squarefree part d=rad(q). Its discussion removes
the possible exceptional real zero once q exceeds an effectively
specified constant depending only on d. Equation (3.2) then applies
when h/φ(q)≥x^(7/12+ε_0), with effective constants depending only on
d and a fixed ε_0>0. Thus it includes even prime powers, not just
powers of a fixed odd prime.

Specialize to

\[
 q=N=2^m,\quad d=2,\quad a=1,\quad
 x=N^4,\quad h=3N^4/4,\quad \varepsilon_0=1/12.
\]

For all sufficiently large N the exceptional zero is absent and

\[
 \frac h{\varphi(N)}=\frac32N^3
       \ge x^{7/12+1/12}=N^{8/3}.
\]

The source's error specializes to

\[
 \sum_{p\in\mathcal P_N}\log p
 =\frac32N^3\left[1+
  O\left(\exp\left\{-c\frac{(\log N)^{1/4}}
                               {(\log\log N)^{3/4}}\right\}\right)\right]
                                                                  \tag{11}
\]

for some absolute c>0; the implied constant is absolute along these
moduli. Indeed log(qx/h)=log(4N/3), the squarefree part is fixed, and
the other denominator term in (3.2) has size
(log N)^(3/7)(log log N)^(3/7), smaller than the displayed leading
denominator. Its error therefore tends to zero.

This is a single application to the full annulus (N^4/4,N^4]; it does
not subtract two estimates with an uncontrolled error. The lower
endpoint N^4/4 is a power of two greater than two, so the open/closed
endpoint makes no difference. Throughout the annulus,
log p=4log N+O(1), proving (4). The analytic ε_0 was fixed once;
the variable tolerance η in Theorem 1 is unrelated to it. The cited
effective constants imply an eventual threshold, but this note does
not claim a numerical value for that threshold.

## Repeated entries cost a negligible prime set

Define the integer boundary product

\[
 F_N=\prod_{j=1}^{N/2-1}
           \left(2^N-(1+\zeta_N^j)^N\right)=P_N(2^N).
                                                                  \tag{12}
\]

The labels represent the inverse pairs of roots other than ±1, and
the corresponding values are invariant under inversion. Galois
permutations therefore make F_N a rational algebraic integer. For
the displayed j,
\((1+\zeta_N^j)^N=(-1)^j(2\cos(\pi j/N))^N\) is real, with
absolute value strictly below 2^N. Consequently

\[
 0<F_N<2^{(N+1)(N/2-1)}.                                  \tag{13}
\]

If an extra zero-sum quadruple has a repeated entry, write it as
{a,a,b,c}. Set u=−b/a and v=−c/a. Then u,v∈H_N and u+v=2.
If u=v then u=v=1 and the quadruple is paired, contrary to its being
extra. Thus h=u/v lies in H_N\{±1}, and

\[
 v=2/(1+h),\qquad (1+h)^N=2^N.
\]

Reduction of (12) at this splitting prime now gives p|F_N. This also
covers a triple repeated entry. Conversely the same calculation shows
that a boundary root supplies an extra repeated-entry quadruple, but
only the forward implication is needed.

By (13), the number of eligible primes with p|F_N is at most

\[
 \frac{(N+1)(N/2-1)\log2}{\log(N^4/4)}
       =O(N^2/\log N)=o(|\mathcal P_N|).                   \tag{14}
\]

No factorization of F_N is required.

## Proof of the exact-class proportion

At a prime not dividing F_N, every extra zero-sum multiset has four
distinct entries. Its scalar stabilizer in H_N is trivial: a nontrivial
subgroup of a dyadic group contains −1, whereas a four-element multiset
invariant under negation is two opposite pairs and would be intrinsic.
Thus each orbit has N multisets and 24N ordered quadruples. It follows
that

\[
 D_N>0,\quad p\nmid F_N\quad\Longrightarrow\quad D_N\ge24N.
                                                                  \tag{15}
\]

Using (10), (14), and (15) gives the explicit bound

\[
 \#\{p\in\mathcal P_N:D_N>0\}
 \le\frac{(N+1)(N/2-1)\log2}{\log(N^4/4)}
       +\frac{N^3\log2}{24\log(N^4/4)}.                   \tag{16}
\]

Divide by (4). The first term tends to zero as a proportion and the
second tends to log(2)/36. Finally, if D_N=0, every nonnegative term
in (7) vanishes, recursively down to H_2. This is equivalent to (5),
and proves Theorem 2.

## What is and is not improved

These are uniform theorems for growing classes in precisely the quartic
window, with all smaller levels included even though those levels are
outside their own quartic windows. On the exact class, the cubic trivial
energy estimate improves to the exact quadratic value. This is not an
improvement to a worst-case Stepanov bound: the exceptional primes in
(2) or (16) can still have uncontrolled pointwise collision mass.

The density-one theorem, the positive-proportion exact class, and the
mean collision bound have different quantifiers. None proves that every
eligible prime is good. None gives the signed centered high moments at
logarithmic depth required by the subgroup spectral target. Indeed,
even in the exact fourth-energy class, higher zero-relation counts
cannot all remain intrinsic because the principal Fourier contribution
is already n^(2r)/p. No transfer from these fourth-moment results to
the classical Paley clique conjecture or the proximity prize is claimed.

## Exact bounded checks

The companion script passed on 16 fixed prime certificates and 68 dyadic
levels, including the three quartic resonances at parent orders 64, 128,
and 512. It compares literal pair sums, primitive-root
fibers, and the exact telescoping identity (8). It verifies the boundary
product test and the orbit divisibility in (15), including cases with
repeated-entry exceptions. It also checks small integer F_N values and
their height bounds, the normalized-tuple AM–GM and norm divisibility
certificates through order 16, and the rational exponent/margin algebra
in the prime-count specialization. These finite checks audit the
identities; the analytic asymptotic follows from the cited theorem,
not from a prime scan or a numerical fit.

The repeated-entry exception is necessary, not just a technical proviso:

| N | p | D_N | F_N mod p | Repeated-entry orbits u | Four-distinct orbits v |
|---:|---:|---:|---:|---:|---:|
| 64 | 6700417 | 768=12N | 0 | 1 | 0 |
| 512 | 17189277697 | 12288=24N | 11555655343 | 0 | 1 |

The first row would contradict a claim of universal 24N divisibility.
The second satisfies the hypothesis used in (15). Literal pair counts
also give B_64=1024, T_64=96 and B_512=67584, T_512=0.

Run the companion Python script
[parallel4_subgroup_2026_09_04.py](../experiments/parallel4_subgroup_2026_09_04.py).
The source, this note, and the input-certificate hashes are recorded in
[the result](../results/parallel4_subgroup_2026_09_04.json).
