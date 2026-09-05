# Dyadic descent and the mixed-energy obligation

**Status: no uniform fourth-energy bound or Paley proof is established.**
This note gives an exact recurrence, an explicit sufficient mixed-energy
estimate, and a descent theorem when the field of definition grows. The
descent theorem does not apply to the steps of the target prime-field
tower. It also gives restrictions on the integer factors from
[quadruple-orbits-and-cube.md](quadruple-orbits-and-cube.md).
No novelty or Lean-formalization claim is made.

## The exact energy recurrence

Let H be a dyadic multiplicative subgroup of order 2k, with k≥2. Let K
be its subgroup of order k and write L=gK, where g∈H\K. Both K and L
are closed under negation. Work in a finite field of odd characteristic,
or a compatible unramified prime-power residue ring where these roots
of unity are distinct modulo the characteristic. Define

\[
 \begin{aligned}
 B&=\#\{(a,b,c,d)\in K^2\times L^2:a+b+c+d=0\},\\
 T&=\#\{(a,b,c,d)\in K^3\times L:a+b+c+d=0\}.
 \end{aligned}
\]

Equivalently B is the mixed additive energy of K and L:
`B=Σ_x r_{K+L}(x)²`. Splitting the ordered zero-sum quadruples by their
number of entries in L gives

\[
 \boxed{E_2(H)=2E_2(K)+6B+8T.}
 \tag{1}
\]

There are two pure patterns, six patterns with two entries from each
coset, and eight patterns with three entries from one coset. Multiplying
by g exchanges the cosets because g²∈K, so the corresponding counts
agree. No independence assumption is used.

Writing r_A for the ordered pair-sum count in A, we have
`B=Σ_x r_K(x)r_L(x)`. Since `E₂(L)=E₂(K)`, Cauchy–Schwarz gives

\[
 k^2\le B\le E_2(K).
\]

Also `r_{K+L}(0)=0` because the cosets are disjoint, whereas r_K(0)=k.
Applying Cauchy–Schwarz only at nonzero sums therefore gives

\[
 \boxed{T^2\le (E_2(K)-k^2)B.}
 \tag{2}
\]

These inequalities alone do not close a quadratic-energy induction.
Substituting B≤E₂(K) in (1)–(2) leaves the upper-bound expression
`8E₂(K)+8sqrt(E₂(K)(E₂(K)−k²))`, whose coefficient exceeds the factor
four needed when the subgroup size doubles. This is a limitation of
that substitution, not an impossibility theorem for using the recurrence.

## A sufficient bound requiring only B

Fix a prime p, a dyadic subgroup H_N⊆F_p*, a number Λ≥1, and a
constant C≥1. Suppose that at **every step** K⊂H_(2k)⊆H_N one has

\[
 B(K,H_{2k}\setminus K)\le Ck^2\Lambda.
 \tag{MB}
\]

Then every subgroup H_s in this tower satisfies

\[
 \boxed{E_2(H_s)\le22C s^2\Lambda.}
 \tag{3}
\]

For the proof, the base group of order two has energy six. If
`E₂(K)≤22Ck²Λ`, then (2), with its negative term dropped, gives
`T≤sqrt(22) Ck²Λ`. Formula (1) yields

\[
 E_2(H_{2k})\le(44+6+8\sqrt{22})Ck^2\Lambda
 <88Ck^2\Lambda=22C(2k)^2\Lambda.
\]

The strict inequality follows from `sqrt(22)<19/4`. This proves the
induction. Thus a uniform (MB) with Λ=max(1,log p), at every level of
each target tower, would imply the necessary fourth-energy estimate
throughout the quartic window, where log p is comparable to log N.
No separate upper bound on T would be needed.

Conversely, uniform energy bounds at every level imply mixed-energy
bounds by B≤E₂(K). Hence these two *whole-tower* families of bounds are
equivalent up to constants. An energy bound just at the target endpoint
does not supply the appropriately normalized bound at all smaller
levels. The smaller levels also lie outside the fixed quartic window.
This distinction is part of the hypothesis, not an implicit reduction.
Neither (MB) nor the requisite endpoint energy bound has been proved.

## Which new quadruple orbits B counts

Use the orbit counts ε,u,v from the preceding note. A pure-coset
nontrivial orbit in H is exactly an old orbit from K. Every new orbit
has either two entries in each coset (balanced), or three in one coset
and one in the other (unbalanced). Write u_(22),v_(22),u_(31),v_(31)
for the numbers of new orbits of the two repeated-entry or four-distinct
types, respectively. The triple-entry type is always unbalanced; its
number is `ε(H)−ε(K)`. Counting orderings and the free scalar action gives

\[
 \boxed{\begin{aligned}
 B&=k^2+4k u_{22}+8k v_{22},\\
 T&=k\bigl(\varepsilon(H)-\varepsilon(K)+3u_{31}+6v_{31}\bigr).
 \end{aligned}}
 \tag{4}
\]

For example, a balanced four-distinct orbit has 2k multisets, and each
has four orderings in K²×L², giving 8k. An unbalanced four-distinct
orbit has k multisets whose majority coset is K, with six orderings of
the first three entries, giving 6k.

The repeated-entry balanced count is at most k/2−1, as also follows
from the polynomial description below. Consequently
`B≤3k²−4k+8k v_(22)`. A uniform bound on these **balanced** four-distinct
orbits at all levels would therefore suffice for (MB). The unbalanced
orbits could then be controlled through (2); they need not be bounded
separately. No balanced-orbit estimate of the needed strength is supplied.

## The integer factors for the new balanced orbits

For k≥4 let λ_i be the labeled roots of P_k from
[kernel-discriminant.md](kernel-discriminant.md). Set

\[
 J_k=\left|\prod_{i<j}(\lambda_i+\lambda_j)\right|,\qquad
 C_k=|P_k(-2^k)|.
 \tag{5}
\]

These are positive integers. The products are symmetric in the roots,
and no factor vanishes: the roots have distinct absolute values, all
strictly smaller than 2^k. If
`Q_k(Y)=∏_i(Y−λ_i²)`, then Q_k has integer coefficients and

\[
 \operatorname{disc}(Q_k)=\operatorname{disc}(P_k)J_k^2.
 \tag{6}
\]

This computes J_k without floating-point roots. The polynomial Q_k is
`(-1)^d P_k(sqrt(Y))P_k(-sqrt(Y))`, where d=deg P_k; the odd powers
cancel identically.

At any compatible prime-power precision, the k-th power labels cosets
of K and multiplying by g negates a label. The kernel multiplicities
are one copy of 2^k and two copies of each λ_i. Therefore

\[
 u_{22}=\#\{i:\lambda_i=-2^k\},\qquad
 v_{22}=\#\{i<j:\lambda_i=-\lambda_j\}.
 \tag{7}
\]

These are counts of the original labeled roots, so multiplicities are
retained. Formula (7) can also be recovered by substituting the kernel
counts in `B−k²=kΣ_b κ_bκ_(-b)` and comparing with (4).

Work in an unramified splitting extension at an arbitrary odd prime p.
Summing (7) over all precisions gives

\[
 v_p(C_k)=\sum_r u_{22,r},\qquad v_p(J_k)=\sum_r v_{22,r}.
\]

The old/new orbit partition then proves the integer divisibilities

\[
 \boxed{U_k\operatorname{odd}(C_k)\mid U_{2k},\qquad
 V_k\operatorname{odd}(J_k)\mid V_{2k}.}
 \tag{8}
\]

The p-adic valuations of the respective quotient factors are exactly
the sums of the new unbalanced orbit counts. This refines the previous
power factorization; it does not bound those valuations.

## Descent when the field of definition doubles

Write `d_n(p)=ord_n(p)`. Suppose `d_(2k)(p)=2d_k(p)`. Then the field
generated by H is a quadratic extension of the field generated by K,
and `{1,g}` is a basis. The same holds for their unramified rings of
integers and their reductions modulo p^r.

In a zero sum with entries in K∪gK, the sums of the two coefficients
of this basis must separately vanish. A three-and-one pattern is
impossible because its singleton coefficient is a unit. In a two-and-two
pattern each pair must be opposite. Hence

\[
 \boxed{B=k^2,\quad T=0,\quad
 E_2(H)=2E_2(K)+6k^2,\quad D(H)=D(K).}
 \tag{9}
\]

In fact all nontrivial multiset orbits descend: ε,u,v agree with those
for K at every precision. This is stronger than a first-level count.

A second useful case is p≡−1 modulo n. Frobenius sends every element
of H_n to its inverse. Applying it to a zero-sum quadruple gives a
zero inverse sum. Thus its first and third elementary symmetric
functions vanish, and its root polynomial is even. Unique factorization
over the residue field forces roots and their negatives to have equal
multiplicity. The quadruple consists of opposite pairs. The same
conclusion holds at every lifted precision since the roots have distinct
residues. In particular

\[
 p\equiv-1\pmod n\quad\Longrightarrow\quad p\nmid U_nV_n.
 \tag{10}
\]

The target has p≡1 modulo the *entire* subgroup order. Its tower
therefore stays inside F_p at every step: `d_s(p)=1`. The basis
separation used in (9) is unavailable there. Treating g∉K as though it
implied that g lies outside the field generated by K would be false.

## Prime support of the new factors

Let n=2^s≥8, and put `a=v₂(p−1)`, `b=v₂(p+1)`. The elementary order
formulas are

\[
 d_n(p)=
 \begin{cases}
 2^{\max(0,s-a)},&p\equiv1\pmod4,\\
 2^{\max(1,s-b)},&p\equiv3\pmod4.
 \end{cases}
\]

They follow by successively factoring p^(2^j)−1; after the first
squaring, each new factor has exactly one factor of two. Together with
(9)–(10), they show that every odd prime newly dividing U_n/U_(n/2)
or V_n/V_(n/2) must satisfy

\[
 \boxed{p\equiv1\ \text{or}\ n/2-1\pmod n.}
 \tag{11}
\]

More precisely, all p-adic valuations of U_n,V_n stabilize once n reaches

\[
 m(p)=\begin{cases}2^a,&p\equiv1\pmod4,\\2^{b+1},&p\equiv3\pmod4.
 \end{cases}
\]

Every later step has doubled field degree. Combining this with the
previous complete factorizations through order 32 gives a uniform
consequence: for every dyadic n≥32,

\[
 V_n=7\cdot17^5\cdot47\cdot79\cdot Z_n,
 \tag{12}
\]

where every prime factor of Z_n is congruent to 1 or −1 modulo 32.
To check completeness of this assertion, the other residue classes are:

| Class of the odd prime p | Stabilization order | Possible factors of V_n |
|---|---:|---|
| p≡3 mod 8 | 8 | none |
| p≡5 mod 8 | 4 | none |
| p≡9 mod 16 | 8 | none |
| p≡7 mod 16 | 16 | only 7, valuation 1 |
| p≡17 mod 32 | 16 | only 17, valuation 5 |
| p≡15 mod 32 | 32 | only 47 and 79, valuation 1 each |

The remaining classes modulo 32 are exactly 1 and −1. These facts use
complete lower-order factorizations and the proved stabilization, not
factorization of all larger V_n. In particular V_64 is still not fully
factored here. Restriction (11) does not exclude any prime already
satisfying the target condition p≡1 mod n.

## A norm cutoff and its quantitative limit

There is another necessary condition for p|V_n:

\[
 \boxed{p^{d_n(p)}\le2^{n/2}.}
 \tag{13}
\]

Indeed a four-distinct zero-sum orbit has no opposite pair. The sum f
of its four complex roots reduces in the basis
`1,ζ_n,...,ζ_n^(n/2−1)` to four distinct nonzero coefficients ±1.
Thus f is nonzero. Averaging |f|² over the n/2 primitive embeddings
gives exactly four: every cross term has mean zero because no exponent
difference is zero or n/2. Arithmetic–geometric mean yields
`|Norm(f)|≤2^(n/2)`. A prime ideal above p has residue degree d_n(p),
so a zero reduction forces p^d_n(p) to divide this nonzero integer norm.

In a splitting prime field, (13) rules out four-distinct nontrivial
orbits when p>2^(n/2). This covers the quartic windows through n=32,
but not any dyadic n≥64: there `2^(n/2)>n⁴`, already above the entire
window. It does not provide the missing growing-order estimate.

It does allow the induction (3) to start later: for every initial level
s with `2^(s/2)<p`, (6) of the preceding orbit note gives
`E₂(H_s)≤9s²`. One may therefore start at the largest such dyadic level
inside the target tower and require (MB) only above it. When this level
is below the endpoint, it has size comparable to log p. The remaining
steps between that size and the quartic-scale endpoint are still
uncontrolled.

For context, the publisher's abstract of
[Do Duc–Leung–Schmidt (2020)](https://alco.centre-mersenne.org/articles/10.5802/alco.86/)
reports a cyclotomic-number bound of three when
`p>(sqrt(14))^{n/d_n(p)}`. That exponential threshold likewise lies
above the working prime-field window for dyadic n≥8. Only the abstract
was inspected and archived; the full theorem is not a dependency of the
self-contained arguments above, and no current-best claim is made.

## Exact verification and remaining work

The script checks (1), (2), (4), and (7) by independent pair counts,
multiset-orbit enumeration, and polynomial-root labels. There are 44
field/precision cases, including ten cases where the field degree doubles.
The latter include quadratic extensions of prime fields and quadratic
extensions of quadratic fields, with exact root-order certificates.

The known quartic examples separate the two mixed terms:

| p | Order of H | Order of K | B | T | New four-distinct orbits |
|---:|---:|---:|---:|---:|---|
| 6700417 | 64 | 32 | 1024 | 96 | none |
| 67403009 | 128 | 64 | 4608 | 0 | one balanced orbit |

Thus neither T=0 nor B=k² is valid uniformly in the target window.
This refutes only those literal assertions; it leaves (MB) open.
Other finite checks outside the quartic window test descent and
stabilization and are not presented as target-regime spectral evidence.

Exact integer checks verify (8) through doubled order 64, full prime
support through order 32, and selected stabilized valuations at orders
32 and 64. Source hashes and all results are in
`results/dyadic_energy_descent.json`. Run
`python3 experiments/dyadic_energy_descent.py`.

The next quantitative obligation is (MB), or equivalently suitable
control of the balanced four-distinct orbits at every needed level.
Neither the field-degree descent nor the norm cutoff applies strongly
enough along the split prime-field tower. Even proving this
fourth-energy obligation would leave the high centered moments and
the full reduction to the prize unresolved.

The subsequent [positive-product analysis](positive-product-moments.md)
extends the two-coset split to logarithmic-depth moments. Its sufficient
positive-product hypothesis reaches the spectral target directly, while
its uniform arithmetic estimate remains open. The principal term in
the corresponding balanced count cannot be omitted.
