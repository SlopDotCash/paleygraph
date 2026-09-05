# Unbalanced product-ratio fibers and a complete small-order classification

**Status: the general aggregate saving remains open. This pass gives an
exact intrinsic product-ratio ledger and a complete finite-class upper
bound for actual subgroups of orders 4, 8, and 16 in their quartic prime
windows. No asymptotic character-sum exponent improves.**

Let \(H\le\mathbb F_p^*\) have dyadic order \(n\ge4\), and suppose
\(n^4/4\le p\le n^4\). For every eligible prime in the indicated
three order classes, the result is

| n | Eligible primes | Opposite-free sixth count \(R_6(H)\) | \(E_3(H)\) | \(\sum_{\rho\in H}N_\rho\) |
| ---: | ---: | --- | --- | --- |
| 4 | 16 | 0 | 400 | 100 |
| 8 | 95 | 0 | 5120 | 640 |
| 16 | 579 | 480 at \(p=33713,37201,41521\); 0 otherwise | 51040 at those three primes; 50560 otherwise | 3190 at those three primes; 3160 otherwise |

Every opposite-free relation at the three exceptional primes is
product-unbalanced at **every** three-versus-three partition. It has
multiplicity pattern \((4,1,1)\). Thus this classification directly
controls the remainder left open by the balanced-subset estimate in
[pass 23](parallel23-subgroup-upper-2026-09-05.md), on a finite family.
The assertions for orders 4 and 8 follow from an analytic norm bound.
Order 16 additionally uses an exact, exhaustive cyclotomic-norm
classification, independently checked by full enumeration over all
579 eligible prime fields.

This is a small-order result with explicit arithmetic certificates,
not a theorem for an unbounded family of subgroup orders. No novelty,
human-refereeing, or formal-certification claim is made.

## 1. Orbit sums retain the unresolved sixth energy

Use the actual distribution and correlations from pass 23:

\[
 w(z)=\#\{(a,d)\in H^2:ad-a-d=z\},\qquad
 N_\rho=\sum_z w(z)w(z/\rho),\qquad B_\rho=nN_\rho.
\]

Here \(B_\rho\) counts ordered equal-sum triples \(x,y\in H^3\)
with \(\prod x_i=\rho\prod y_i\). For a nonzero coset \(H_j\),
put \(W_j=\sum_{z\in H_j}w(z)\). Grouping the dilation correlations
by their orbits gives exactly

\[
 \sum_{\rho\in H}N_\rho
 =n w(0)^2+\sum_j W_j^2.                             \tag{1}
\]

The zero orbit is fixed by all \(n\) dilations, explaining its factor
\(n\). More significantly, the orbit masses themselves are already
the actual nonzero triple-sum counts:

\[
 \boxed{W_j=r_3(x)=n\mathbf1_{j=0}+(C^2)_{j0}
       \quad(x\in H_j),\qquad w(0)=\kappa_0.}        \tag{2}
\]

For the first equality, divide \(ad-a-d\) by \(ad\in H\).
Counting its membership in \(H_j\), and using \(-1\in H\), is
equivalent to counting \(u,v\in H\) with \(1-u-v\in H_j\).
Convolution and the cyclotomic identity \(C_{ij}=C_{-i,j-i}\)
give \(n\mathbf1_{j=0}+(C\kappa)_j\), exactly the previously
derived formula for \(r_3(x)\) on \(H_j\). At zero, inversion of
\(ad=a+d\) gives \(a^{-1}+d^{-1}=1\), hence \(w(0)=r_2(1)=\kappa_0\).

Consequently (1) is precisely
\(E_3=n^2w(0)^2+n\sum_jW_j^2\), not a new bound. It makes clear
why treating orbit masses as a fresh incidence input would be circular
unless a separate estimate on them is supplied. Applying Cauchy on
each orbit still yields the previously recorded fourth-power/log
scale. This pass found no general improvement from that step.

## 2. An exact intrinsic ledger in every product-ratio fiber

Call a zero-sum six-word intrinsic when its multiset is a union of
three opposite pairs. After negating its last three coordinates,
write \(I_\rho\) for its contribution to \(B_\rho\), and let
\(H^{(2)}=\{h^2:h\in H\}\). Then

\[
 \boxed{
 I_\rho=
 \begin{cases}
 6n^3-9n^2+4n,&\rho=1,\\
 18n(n-2),&\rho\in H^{(2)}\setminus\{1\},\\
 0,&\rho\notin H^{(2)}.
 \end{cases}}                                        \tag{3}
\]

To verify this without positional overcounting, an intrinsic pair of
triples is either a permutation pair, or can be written as multisets

\[
 x=\{a,-a,z\},\qquad y=\{b,-b,z\}.
\]

The latter form has product ratio \(a^2/b^2\). If the ratio is one,
the two opposite value pairs coincide and the triples are
permutations. All permutation pairs therefore contribute exactly
\(36\binom n3+9n(n-1)+n=6n^3-9n^2+4n\) to \(I_1\).

For a fixed nonidentity square ratio, there are \(n/2\) choices of
the ordered opposite value pairs \(\{a,-a\},\{b,-b\}\); the second
is determined by the first and is different from it. For each such
choice, \(n-4\) values of \(z\) give \(6\cdot6=36\) ordered
triple pairs. The other four values give \(3\cdot6=18\) each. Thus

\[
 I_\rho=\frac n2[36(n-4)+4\cdot18]=18n(n-2).
\]

No nonsquare ratio can occur. Summing (3) gives exactly

\[
 \sum_\rho I_\rho=15n^3-45n^2+40n=T_6(n).             \tag{4}
\]

In general one must retain a third category: a nonintrinsic word
containing an opposite pair. If its fiber count is \(J_\rho\), and
\(R_\rho\) counts opposite-free words at the fixed product ratio,
then

\[
 nN_\rho=I_\rho+J_\rho+R_\rho,\qquad
 \sum_\rho J_\rho
 =(15n-60)(E_2-T_4)+60n(r_2(2)-1)-30n\mathbf1_{3\in H}.
                                                               \tag{5}
\]

The last identity is the existing
[opposite-pair decomposition](parallel21-subgroup-next-input-2026-09-05.md).
If \(E_2=T_4=3n^2-3n\), every \(J_\rho\) vanishes: deleting an
opposite pair leaves an intrinsic four-word. Then (3) subtracts all
intrinsic contributions fiber by fiber. In general, \(\sum_{\rho\ne1}
R_\rho\) can still include words balanced at another partition;
it must not be identified with the fully unbalanced remainder.

## 3. A dyadic norm bound that uses actual multiplicative closure

Let \(n=2^k\), \(d=n/2\), and let \(\zeta\) be a complex primitive
\(n\)-th root of unity. For an opposite-free multiset of \(s\) powers
of a generator of \(H\), reduce its exponent polynomial modulo
\(X^d+1\). It becomes

\[
 f(X)=\sum_{j=0}^{d-1}a_jX^j,
\]

where opposite-freeness ensures that the absolute values of the
nonzero coefficients are the original multiplicities. Divide by their
positive gcd to obtain primitive integer coefficients, and put
\(Q=\sum_j a_j^2\) for this divided polynomial. The division does
not change vanishing modulo any prime \(p>s\).

The polynomial \(X^d+1\) is irreducible over the rationals: shifting
by one gives an Eisenstein polynomial at 2, since \(d\) is a power
of two. In particular the nonzero polynomial \(f\), of degree less
than \(d\), has nonzero integer norm

\[
 \mathcal N(f)=\prod_{u\text{ odd}\bmod n} f(\zeta^u).
\]

For completeness, the Eisenstein step follows from the elementary
mod-2 test: any monic integer factorization after the shift would have
both factors reducing to powers of \(X\), forcing both constant
terms even and their product divisible by 4, although the constant
term is 2. The norm is also the integer determinant of multiplication
by \(f\) in \(\mathbb Z[X]/(X^d+1)\).

The primitive-root orthogonality in this dyadic case gives

\[
 \frac1d\sum_{u\text{ odd}\bmod n}|f(\zeta^u)|^2=Q.
\]

Indeed, the sum over odd \(u\) of \(\zeta^{u(i-j)}\) is zero
whenever \(0<|i-j|<d\). AM–GM therefore gives

\[
 0<|\mathcal N(f)|\le Q^{d/2}=Q^{n/4}.               \tag{6}
\]

If the original word vanishes in \(\mathbb F_p\), its generator is
a root of \(X^d+1\) modulo \(p\) at which \(f\) vanishes.
Reducing the determinant, or the norm product, modulo \(p\) shows
\(p\mid\mathcal N(f)\). Thus

\[
 \boxed{\text{an opposite-free relation forces }p\le Q^{n/4}.}
                                                               \tag{7}
\]

This is a necessary arithmetic condition, not an energy restatement.
For four entries, primitive content removal gives \(Q\le10\): the
largest case is multiplicities \((3,1)\). For six entries it gives
\(Q\le26\), with largest case \((5,1)\). A single supported value
has primitive coefficient one and causes no exception.

For \(n=4,8,16\), \(10^{n/4}<n^4/4\); hence every zero-sum
four-word in the quartic window is intrinsic and \(E_2=T_4\).
For \(n=4,8\), also \(26^{n/4}<n^4/4\), excluding every
opposite-free six-word. This proves the first two rows of the table
without a prime-field search.

## 4. Order 16: exact finite norm classification

At order 16 the analytic bound alone does not exclude all six-word
patterns. It immediately excludes all but
\((5,1),(4,1,1),(3,2,1),(3,1,1,1)\): after content removal, every
other pattern has \(Q\le10\). The remaining question is finite.

The certificate enumerates all opposite-free six-multisets of
exponents modulo 16 and takes their lexicographically least translation
with an occupied exponent moved to zero. There are 1,688 distinct
representatives. Their multiplicity counts are

| Pattern | Representatives | Possible quartic prime divisors of their norms |
| --- | ---: | --- |
| \((6)\) | 1 | none |
| \((5,1)\) | 14 | none |
| \((4,2)\) | 14 | none |
| \((4,1,1)\) | 84 | 33713, 37201, 41521 |
| \((3,3)\) | 7 | none |
| \((3,2,1)\) | 168 | none |
| \((3,1,1,1)\) | 280 | none |
| \((2,2,2)\) | 28 | none |
| \((2,2,1,1)\) | 420 | none |
| \((2,1,1,1,1)\) | 560 | none |
| \((1,1,1,1,1,1)\) | 112 | none |

“Possible quartic” means a prime divisor \(p\equiv1\pmod{16}\)
in \([16384,65536]\). Each of the three primes occurs in eight
translation representatives, with a single norm value:

| p | One exponent representative | Norm |
| ---: | --- | ---: |
| 33713 | \((0,0,0,0,1,10)\) | \(67426=2p\) |
| 37201 | \((0,0,0,0,1,14)\) | \(74402=2p\) |
| 41521 | \((0,0,0,0,1,4)\) | \(83042=2p\) |

The complete norm and factorization ledger is in the result artifact;
this table is the output of exhaustive exact arithmetic, not an
unproved symbolic classification. The norm calculation uses the
recursion

\[
 f(X)=e(X^2)+Xo(X^2),\qquad
 \operatorname{Norm}_{X^d+1}(f)
 =\operatorname{Norm}_{Y^{d/2}+1}\bigl(e(Y)^2-Yo(Y)^2\bigr),
\]

reducing the polynomial modulo \(Y^{d/2}+1\) at each step. Every
norm is separately checked against an exact rational determinant of
the multiplication matrix. Trial division checks the complete prime
factorization.

At each exceptional prime the unique normalized unordered pair with
sum \(-4\) is as follows:

| p | \(\{b,c\}\subset H\), with \(b+c=-4\pmod p\) |
| ---: | --- |
| 33713 | \(\{4508,29201\}\) |
| 37201 | \(\{17702,19495\}\) |
| 41521 | \(\{8512,33005\}\) |

These are distinct elements, neither is \(1\) or \(-1\), and they
are not opposites. Their membership is verified in the actual
order-16 subgroup. The four-term intrinsic result already proved
implies that there cannot be two different unordered pairs with the
same nonzero sum: their equality would be a nonintrinsic zero-sum
four-word. Thus every exceptional six-word is a permutation and
scaling of \((1,1,1,1,b,c)\), and its exact count is

\[
 R_6=16\cdot\frac{6!}{4!}=480.                       \tag{8}
\]

All 24 exceptional exponent representatives have odd exponent sum.
Therefore the product of the six entries belongs to
\(H\setminus H^{(2)}\). Translation of all six exponents does not
change that parity, nor does changing the primitive generator. Since
\(-1\in H^{(2)}\) at order 16, every three-versus-three ratio
\(-\prod_Ih_i/\prod_{I^c}h_i\) is also nonsquare in \(H\).
It cannot equal one. This proves that all 480 words are unbalanced
at every partition, not merely at the fixed split used by (3).

Equations (4), (5), and (8) now prove the exact order-16 energy and
aggregate values in the opening table. In particular \(E_3\le15n^3\)
on all three finite order classes, including their exceptional primes.

## 5. Verification and the remaining general input

The new [verifier](../experiments/parallel24_subgroup_unbalanced_2026_09_05.py)
completed with `PASS`; the
[results](../results/parallel24_subgroup_unbalanced_2026_09_05.json)
contain every norm representative and every eligible prime field.
There are 1,770 independently checked norm determinants across the
three orders and 690 actual quartic prime fields. Direct additive
convolution and separately weighted normalized six-word enumeration
agree in every field. The latter checks intrinsic and opposite-free
words separately, including 178,782 normalized multiset matches.

The verifier also checks 10,088 individual product-ratio correlations,
the intrinsic ledger (3) in every such fiber, and 26,674 orbit masses
against actual triple-sum counts. The three exceptional fields have
their full subgroups, concrete six-word witnesses, and residual fibers
recorded. All four nonzero residual fibers in each exceptional field
are nonsquare. Their masses are 144, 144, 96, and 96, summing to 480.
These exact checks use Python integers and rational arithmetic, with
no floating-point acceptance threshold.

For general growing \(n\), the norm bound (7) is exponential in
\(n\), whereas the target prime scale is polynomial in \(n\).
The small-order exclusion therefore does not extend automatically
to the requested asymptotic family. The unrestricted orbit estimate
in (1) remains the original sixth-energy problem, and the exact
intrinsic subtraction in (3) does not control the remaining positive
terms. No new general aggregate upper bound, density improvement,
worst-case cancellation exponent, or prize bridge is established.

## Input hashes

| Input | SHA-256 |
| --- | --- |
| `experiments/parallel24_subgroup_unbalanced_2026_09_05.py` | `3d0ffb8477da636b3466b97360e907ff554f54255beeb4951898b4464b96cca5` |
| `results/parallel24_subgroup_unbalanced_2026_09_05.json` | `c75846a83de374db55cdd41ef3f5a394a18af059314795f89cd2e605762e7336` |
| `research/parallel23-subgroup-upper-2026-09-05.md` | `1b97e491bc1c0b2cfe391388b0340d89d1f2a7168822db116513b63f4148bd41` |
| `research/parallel21-subgroup-next-input-2026-09-05.md` | `84005053f3954ee2d36eccf85b2786f1fe00d6bd956c76b23b1b40e82b991c22` |
| `research/mixed-periods-and-shifted-energy.md` | `de7d7202e55253955b62f73c973d2bdee42b48123eb344c875d5f0527a21f4aa` |
