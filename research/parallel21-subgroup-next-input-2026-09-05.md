# Sixth-energy exceptions: an exact relation decomposition

**Status: the uniform square-root bound for the dyadic subgroup periods is
unproved. No cancellation exponent is improved here.** This lane isolates the
sixth-energy input into opposite-free four- and six-term relations. It proves
an exact decomposition, checks every order-64 sixth-energy exception listed
by the sigma audit, and constructs genuine symmetric sets in the same
quartic prime fields whose fourth energy is exactly intrinsic but whose
sixth energy grows on the fourth-power scale. These constructed sets are
not multiplicative subgroups, so they do not refute the subgroup target.

The useful next input is now explicit: control the six-term relations
containing no opposite pair, together with the already necessary fourth
energy. Multiplicative structure must enter that control. Fourth energy and
prime-field order alone do not suffice.

## 1. The exact decomposition

Let p>3 be prime, let S=-S be a nonempty subset of F_p* of size n, and put

\[
 r_2(x)=\#\{(u,v)\in S^2:u+v=x\},\qquad
 E_j(S)=\#\{(x,y)\in S^j\times S^j:\textstyle\sum_i x_i=\sum_i y_i\}.
\]

By negating the last j coordinates, E_j is also the ordered zero-sum
count of length 2j.
Write

\[
 T_4(n)=3n^2-3n,\qquad T_6(n)=15n^3-45n^2+40n.
\]

These count the zero-sum words whose multiplicity at u equals that at -u
for every opposite pair {u,-u}. They are intrinsic relations: they vanish
formally before any reduction modulo p. Let R_6(S) count the ordered
zero-sum six-tuples **containing no opposite pair**, and define

\[
 U(S)=\sum_{u\in S}\bigl(r_2(2u)-1\bigr),\qquad
 J(S)=\#\{u\in S:3u\in S\}.
\]

Then

\[
 \boxed{E_3(S)=T_6(n)+(15n-60)\bigl(E_2(S)-T_4(n)\bigr)
                 +60U(S)-30J(S)+R_6(S).}                 \tag{1}
\]

For a multiplicative subgroup H containing -1, put K=r_2(2) and
\(\epsilon_3=1_{3\in H}\). Multiplicative invariance gives the compressed
formula

\[
 \boxed{E_3(H)=T_6(n)+(15n-60)\bigl(E_2(H)-T_4(n)\bigr)
                 +60n(K-1)-30n\epsilon_3+R_6(H).}        \tag{2}
\]

No analytic estimate is used in either identity.

### Proof

Set m=n/2. Intrinsic words of length six have respectively three, two,
or one occupied opposite pairs. Their numbers are

\[
 720\binom m3+180m(m-1)+20m=T_6(n).
\]

Similarly, the intrinsic four-term count is
\(24\binom m2+6m=T_4(n)\).

A zero-sum four-tuple is either intrinsic or contains no opposite pair:
once an opposite pair is deleted, the two remaining entries must also be
opposite. A four-tuple with no opposite pair has multiplicity pattern
1111, 211, or 31. Pattern 22 forces an opposite pair; pattern 4 is
impossible in characteristic p>3. Denote the three ordered counts by
A_1111, A_211, A_31. Thus their sum is E_2-T_4.

Consider a nonintrinsic zero-sum six-tuple containing an opposite pair.
Deleting that pair leaves a nonintrinsic, opposite-free four-tuple. There
cannot be two disjoint opposite pairs: deleting them would leave a third
opposite pair and make the word intrinsic. Consequently all opposite
pairs in the six-tuple belong to one value pair {u,-u}, and deleting one
of them gives a unique four-element multiset.

Fix such a four-element multiset. Its ordered weight is
\(4!/\prod_v a_v!\). Extend it by one opposite value pair. If that pair
is absent, the ordered weight is multiplied by 30. If one of its values
already occurs a times, the multiplier is 30/(a+1). Therefore the total
multipliers, summing over all n/2 possible opposite value pairs, are

| Four-term multiplicities | Extension multiplier |
|---|---:|
| 1111 | 15n-60 |
| 211 | 15n-50 |
| 31 | 15n-75/2 |

The fractional last coefficient causes no integrality issue: A_31 is a
multiple of four. The uniqueness just proved prevents double counting
between base multisets. It follows that

\[
 E_3-T_6-R_6=(15n-60)A_{1111}+(15n-50)A_{211}
                         +(15n-75/2)A_{31}.              \tag{3}
\]

A pattern-31 relation is (u,u,u,-3u), in four placements, so A_31=4J.
For pattern 211, choose the six repeated-position pairs and then choose
v,w with v+w=-2u. From r_2(-2u)=r_2(2u), remove (-u,-u), and remove
( u,-3u ) and ( -3u,u ) when 3u belongs to S. Hence

\[
 A_{211}=6U-12J.
\]

Substitution into (3) proves (1), and subgroup invariance proves (2).

## 2. The next arithmetic input, with its exact payoff

The positive extension count also gives

\[
 \boxed{T_6(n)+R_6(S)\le E_3(S)
       \le T_6(n)+15n\bigl(E_2(S)-T_4(n)\bigr)+R_6(S).}   \tag{4}
\]

In particular, if E_2(S)=T_4(n), then

\[
 \boxed{E_3(S)=T_6(n)+R_6(S).}                            \tag{5}
\]

For a subgroup in the working window n^4/4<=p<=n^4, suppose the two
quantitative inputs

\[
 E_2(H)-T_4(n)\le A n^2 L,\qquad R_6(H)\le D n^3 L
                                                               \tag{P6}
\]

hold with L>=1. Formula (4) gives

\[
 E_3(H)\le[15+(15A+D)L]n^3.
\]

The self-contained Konyagin inequality proved in
[sigma-subgroup](sigma-subgroup-2026-09-05.md), at orders (3,3), then gives

\[
 \boxed{M(H)\le[15+(15A+D)L]^{1/9}n^{8/9}.}               \tag{6}
\]

Thus a uniform (P6), with L a fixed power of log n and A,D independent
of n and p, would extend the 8/9 exponent to every prime where those
hypotheses hold. On the primes with intrinsic fourth energy, its only
remaining input is the bound on R_6.

This is a reduction, not a proof of (P6). Moreover, 8/9 remains far from
the square-root target. The logarithmic-depth condition (SG) in
[subgroup-target](subgroup-target.md) is still needed for that target;
no assertion of equivalence to the official prize is made.

## 3. What caused the known finite exceptions

The sigma audit lists fourteen quartic primes at n=64 that fail its
particular bound E_3<=(15+log 64)64^3. Every one of those primes has

\[
 E_2-T_4=1536=24n.
\]

At thirteen of them, the exact part of the sixth-energy excess forced by
four-term relations is 1,382,400. At p=11127041 it is 1,397,760, with
r_2(2)=5. All fourteen have 3 outside H. These are exact counts of actual
subgroups, independently reconstructed in this lane.

| p | E_2-T_4 | Opposite-containing sixth excess | R_6 |
|---:|---:|---:|---:|
| 6878593 | 1536 | 1382400 | 276480 |
| 7041409 | 1536 | 1382400 | 276480 |
| 7177601 | 1536 | 1382400 | 276480 |
| 7204033 | 1536 | 1382400 | 353280 |
| 7884353 | 1536 | 1382400 | 276480 |
| 7987009 | 1536 | 1382400 | 276480 |
| 8019073 | 1536 | 1382400 | 276480 |
| 9190913 | 1536 | 1382400 | 345600 |
| 9877633 | 1536 | 1382400 | 276480 |
| 10219457 | 1536 | 1382400 | 322560 |
| 11127041 | 1536 | 1397760 | 122880 |
| 12942337 | 1536 | 1382400 | 276480 |
| 13640513 | 1536 | 1382400 | 276480 |
| 14721281 | 1536 | 1382400 | 276480 |

In particular R_6<2n^3 at every listed exception. The largest excess is
therefore principally the forced extension of four-term relations, not
an especially large opposite-free sixth count. This explains the finite
exception mechanism. It proves no uniform upper bound on either class
of relations and does not extrapolate that mechanism to all primes.

For the separate resonant test p=6700417, n=64, the decomposition is

\[
 E_3=T_6+698880+367680,
 \qquad E_2-T_4=768,\quad r_2(2)=3,\quad 3\notin H.
\]

## 4. Intrinsic fourth energy does not force a good sixth energy

The following obstruction uses genuine subsets of prime fields and the
same dyadic cardinality and quartic range. It deliberately drops the
multiplicative-subgroup hypothesis. This shows exactly which structure a
proof cannot discard.

Let n=2k with k>=8 a power of two. Choose an odd prime q with
k<=q<2k. For t=0,...,k-1 define integers

\[
 c_t=t+2q(t^2\bmod q),\qquad b_t=8q^2+c_t,
 \qquad B=\{b_t:0\le t<k\}.
\]

Let p be any prime in n^4/4<=p<=n^4 with p=1 modulo n, and regard

\[
 S=B\cup(-B)\subset\mathbb F_p.
\]

Then

\[
 \boxed{|S|=n,\quad S=-S,\quad
 E_2(S)=3n^2-3n,\quad E_3(S)\ge n^4/20.}                \tag{7}
\]

**Proof.** If c_a+c_b=c_c+c_d as integers, reduction modulo 2q gives

a+b=c+d

as integers, since both sides lie between 0 and 2q-2. The remaining
parts give a^2+b^2=c^2+d^2 modulo q. The sum and sum of squares recover
the product modulo the odd prime q, hence the two unordered residue
pairs are equal. Since all indices are in [0,q), the unordered integer
pairs are equal. Thus B is an integer Sidon set, with repetitions
allowed in the two-term sums.

We have 8q^2<=b_t<10q^2 and

\[
 6\max B<60q^2<240k^2<4k^4\le p,
\]

where k>=8 is used in the last strict inequality. So a signed sum of at
most six elements of B vanishes modulo p only if it vanishes as an
integer. A signed four-term zero sum must have two positive and two
negative entries, because 3min B>max B. The Sidon property then forces
opposite pairing, proving the fourth-energy assertion.

A signed six-term zero sum must have three positive and three negative
entries: 4min B>2max B excludes the nearest alternative. Consequently

\[
 E_3(S)=20E_3(B).
\]

The k^3 ordered triple sums from B take at most 6q^2+1<=25k^2 distinct
integer values. Cauchy--Schwarz therefore gives

\[
 E_3(B)\ge k^6/(25k^2)=k^4/25,
 \qquad E_3(S)\ge(4/5)k^4=n^4/20.
\]

This proves (7). Formula (5) also gives

\[
 R_6(S)\ge n^4/20-T_6(n).
\]

Thus E_3/n^3 is unbounded along these parameters, and cannot be bounded
by a fixed power of log n. The theorem itself is conditional only on
the displayed elementary prime parameters. To obtain unbounded dyadic sizes, take arbitrarily large odd primes q
and let k be the largest power of two at most q; then k<=q<2k
automatically. Nonemptiness of the quartic progression for all
sufficiently large dyadic n uses the already sourced prime-count
specialization in [pass 4](parallel4-subgroup-2026-09-04.md); no new
prime-gap assertion is required here.

These sets are not subgroups: in particular 1 is absent from S. They
satisfy the symmetry, field-size and fourth-energy constraints, but not
multiplicative closure. They are not counterexamples to (P6) for H, to
(SG), or to the prize.

Exact finite instances are:

| n | q | p | E_2(S) | E_3(S) |
|---:|---:|---:|---:|---:|
| 16 | 11 | 16417 | 720 | 54880 |
| 32 | 17 | 262337 | 2976 | 626000 |
| 64 | 37 | 4194433 | 12096 | 7597000 |
| 128 | 67 | 67109633 | 48768 | 119956280 |
| 256 | 131 | 1073748737 | 195840 | 1811245240 |

The asymptotic obstruction follows from the construction and inequalities,
not from these five examples.

## 5. Verification and limitations

Run [the exact verifier](../experiments/parallel21_subgroup_next_input_2026_09_05.py).
The [results](../results/parallel21_subgroup_next_input_2026_09_05.json)
record 44 subgroup fixtures, nine further symmetric-set fixtures, five
Sidon constructions, 34 multiplicity-classification checks, and 86,922
admissible weighted triple-pair terms in the direct R_6 enumeration.
All 53 general decomposition checks and all 44 subgroup specializations
pass exactly. Direct r_2/r_3 convolution is independent of the
opposite-free enumeration, which uses unordered triples and exact
factorial multiplicities. All energies, primality checks, weights and
bounds use integers; no floating-point tolerance is used. Source and
script hashes are recorded in the result.

The asymptotic arguments are ordinary mathematics and have not been
formalized in Lean or independently refereed. The conditional payoff
uses the already reconstructed Konyagin inequality. During this lane,
[Di Benedetto et al., Sections 4--5](https://arxiv.org/html/2003.06165v1#S4)
was freshly inspected to check whether a simple trilinear substitution
could improve 8/9; the earlier
[amplification ledger](analytic-bounds-and-amplification.md) already
records why that substitution is weaker. No new result from that source
is used in (1)--(7), and no source archive or central project document was
modified by this lane.

The next substantive task is an upper estimate for R_6(H) exploiting
multiplicative closure, while retaining the explicit four-term resonance
contribution in (2). The current argument supplies the reduction and a
counterexample to dropping that structure; it supplies no such upper
estimate. The uniform subgroup target, the exceptional-prime problem,
the classical two-set conjecture, and the official-prize bridge remain
unproved.
