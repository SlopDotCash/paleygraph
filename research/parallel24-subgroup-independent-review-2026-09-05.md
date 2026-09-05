# Independent review: the finite dyadic subgroup classification

**Verdict: no mathematical correction requested.** The reviewed proof and
verifier support the stated complete finite classification for subgroup
orders 4, 8, and 16 in their quartic prime windows. They do not establish
an estimate for growing subgroup order or an asymptotic saving toward
the full conjecture.

This review independently rederived the intrinsic product-ratio ledger,
the dyadic norm bound, the coverage argument, and the all-partitions
conclusion. Fresh computations used separately written code without
importing the production verifier. They did not rerun its complete
690-field calculation.

## Reviewed bytes

| Artifact | SHA-256 |
| --- | --- |
| [Proof](parallel24-subgroup-unbalanced-2026-09-05.md) | `ffe63245f28880ff53ace7210bf54f3baafca77c110b5aec0ac42cd3bf99bf1b` |
| [Verifier](../experiments/parallel24_subgroup_unbalanced_2026_09_05.py) | `3d0ffb8477da636b3466b97360e907ff554f54255beeb4951898b4464b96cca5` |
| [Results](../results/parallel24_subgroup_unbalanced_2026_09_05.json) | `c75846a83de374db55cdd41ef3f5a394a18af059314795f89cd2e605762e7336` |

Every hash in the result artifact's `input_sha256` mapping matched its
current file. Python syntax parsing succeeded. This review changes none
of the reviewed files.

## Intrinsic fibers and the zero orbit

For an intrinsic zero-sum six-word, negating the last three coordinates
gives either a permutation pair of triples or multisets

    x={a,-a,z},       y={b,-b,z}.

The latter representation has unique internal opposite value pairs when
the triples are not permutations. Squaring identifies the n/2 opposite
value pairs bijectively with H^(2). A prescribed nonidentity square ratio
therefore admits exactly n/2 ordered choices of the two different
opposite value pairs. There are n−4 external choices of z, each of
weight 36, and four overlapping choices of weight 18. This proves

    I_1=6n³−9n²+4n,
    I_ρ=18n(n−2) for ρ∈H^(2)\{1},
    I_ρ=0 otherwise.

The formula includes n=4: its two nonzero intrinsic fibers are
I_1=256 and I_(−1)=144, summing to T_6(4)=400. There is no absent
external-z case being inadvertently counted at this endpoint.

The three-way decomposition into intrinsic, opposite-containing
nonintrinsic, and opposite-free words is necessary in general. The
reviewed note retains that category and correctly specializes the
[pass-21 decomposition](parallel21-subgroup-next-input-2026-09-05.md).
When E_2=T_4, deleting any opposite pair leaves an intrinsic four-word,
so every opposite-containing six-word is intrinsic. This justifies
discarding the middle category in the finite classification.

The orbit identity also has the correct normalization. For fixed
x∈H_j, the pairs u,v∈H with s=1−u−v∈H_j are in bijection with
triples a,b,c∈H summing to x, by

    c=x/s,       a=−cu,       b=−cv.

Inversion of a,d in ad−a−d gives exactly the same pair count on
that coset, because ad∈H. Hence W_j=r_3(x). At zero,
w(0)=r_2(1) and r_3(0)=nw(0). The factors n in the nonzero
coset energy and n² in its zero term are therefore correct. This is
an energy identity, as the note explicitly acknowledges.

## Norm bound and the small-order exclusions

For n=2^k and d=n/2, an opposite-free multiset reduces modulo X^d+1
to coefficients whose nonzero absolute values are its multiplicities.
Removing their positive gcd is legitimate when p>s, where s is the
word length: that gcd is at most s and so is invertible modulo p.
Every prime in the finite windows exceeds six.

The Eisenstein argument after X↦X+1 proves the irreducibility of
X^d+1. Thus a nonzero reduced polynomial of degree below d has
nonzero integer norm. For 0<|i−j|<d, summing ζ^{u(i−j)} over odd
u modulo n gives zero. Consequently the mean squared complex absolute
value of the d conjugates is Q, and AM–GM gives |Norm(f)|≤Q^(d/2).
Reduction at the order-n element of the prime field implies that p
divides the norm. No cancellation or equidistribution assumption enters
this argument.

Primitive content removal is essential to the stated pattern bounds:
Q≤10 for four entries and Q≤26 for six entries. For example, (4,2)
reduces to primitive multiplicities (2,1), with Q=5; it must not be
assigned Q=20. At six entries the only patterns with primitive Q>10
are (5,1), (4,1,1), (3,2,1), and (3,1,1,1).

The numerical inequalities 10^(n/4)<n⁴/4 for n=4,8,16 and
26^(n/4)<n⁴/4 for n=4,8 hold strictly. They exclude all
opposite-free four-words in the three windows, and all opposite-free
six-words in the first two. A nonintrinsic zero-sum four-word is
automatically opposite-free, so the passage to E_2=T_4 is valid.

## Independent coverage certificate

The translation action on opposite-free exponent multisets is free for
dyadic n. Any nontrivial stabilizer subgroup contains translation by
n/2, which would force opposite exponents into the support. Therefore
the number of translation orbits of opposite-free six-multisets is

    (1/n) Σ_(s=1)^min(n/2,6) binom(n/2,s) 2^s binom(5,s−1).

Here one chooses s opposite value pairs, one sign from each, and a
positive composition of six over the chosen support. This independent
formula gives the following totals.

| n | All opposite-free six-multisets | Translation orbits | Eligible primes |
| ---: | ---: | ---: | ---: |
| 4 | 24 | 6 | 16 |
| 8 | 608 | 76 | 95 |
| 16 | 27008 | 1688 | 579 |

For all 1,770 stored representatives, separate code checked that the
word is opposite-free, is the least of **all n translations**, and is
distinct from every other stored representative. The preceding counts
therefore certify full orbit coverage; no sampling premise is needed.
The code also reconstructed every primitive coefficient vector and
content, every stored prime factorization product, and every eligible
prime-factor list. A separate sieve through 65,536 verified every listed
factor's primality and independently reproduced all 690 eligible prime
fields, including the quartic-window and divisibility restrictions.

The production norm recursion correctly pairs X and −X and reduces
e(Y)²−Yo(Y)² modulo Y^(d/2)+1. Its second norm implementation uses
exact rational elimination on the multiplication matrix. Both were
reviewed in source. Independently written integer Bareiss elimination
recomputed the three displayed exceptional norms:

| Exponent multiset | Determinant |
| --- | ---: |
| (0,0,0,0,1,10) | 67426 = 2·33713 |
| (0,0,0,0,1,14) | 74402 = 2·37201 |
| (0,0,0,0,1,4) | 83042 = 2·41521 |

This fresh audit recomputed those three determinants, rather than all
1,770 production determinants. The full determinant agreement is a
recorded production check whose implementations and coverage were
reviewed here.

## Independent field checks and every-partition conclusion

Separate code built H directly as {h: h^n=1} in each tested field,
enumerated all ordered triples, and joined triples with equal sums.
This counts product-ratio fibers without normalizing the first entry.
A second direct enumeration traversed all unordered six-multisets,
retained zero-sum opposite-free ones, and restored their exact
multinomial weights. Both routes agreed:

| p | n | E_3 | R_6 | Opposite-free multisets | Three-subset checks |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 73 | 4 | 400 | 0 | 0 | 0 |
| 33713 | 16 | 51040 | 480 | 16 | 320 |
| 37201 | 16 | 51040 | 480 | 16 | 320 |
| 41521 | 16 | 51040 | 480 | 16 | 320 |

Every intrinsic fiber in these four fields matched the ledger, and the
residual fibers in each exceptional field had masses 144,144,96,96.
The 960 fresh partition checks test all 20 three-subsets of every
opposite-free multiset. They independently confirm the stronger
all-partitions conclusion, not just the production verifier's chosen
witness at ten complementary splits.

The proof also establishes that conclusion without requiring this
additional enumeration. Every eligible norm representative has odd
exponent sum. Translation adds a multiple of six, and a new primitive
generator changes exponent parity by an odd multiplier, preserving it.
Thus the sixfold product is nonsquare within H. Since −1 is a square
within an order-16 group, each three-versus-three ratio

    −(product on I)/(product on Iᶜ)

has that same nonsquare class and cannot equal one.

Finally, E_2=T_4 ensures uniqueness of the unordered pair with any
nonzero prescribed sum: two different such pairs would give a
nonintrinsic zero-sum four-word. The displayed pair summing to −4
at each exceptional prime therefore yields exactly 16 scalar copies
of a (4,1,1) multiset, each of weight 6!/4!=30. This proves R_6=480
with no hidden multiplicity or stabilizer factor.

The verified result is a finite arithmetic classification. The norm
bound grows exponentially with n while the prime window grows
polynomially; no extension to all dyadic orders is supplied. The
general sixth-energy estimate and the full goal remain unproved.
