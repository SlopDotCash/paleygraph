# Compatible shifted-product equations: a weighted exceptional-prime bound

For each fixed dyadic order n>=4, this pass proves an unconditional
bound for the collision excess summed over eligible primes. It implies
an absolute bound on the number of primes in the quartic range that can
fail the triangle criterion from pass41. It does not bound every prime,
prove the full Paley conjecture, or settle the official prize.

This specializes the workspace's earlier
[cyclotomic prime-average argument](cyclotomic-prime-average.md) while
retaining the multiplicative compatibility of six-term relations. The
underlying norm characterization of circularity is published in
[Ke–Kiechle, 2023, Proposition 3.4 and Corollary 3.7](https://admjournal.luguniv.edu.ua/index.php/adm/article/download/2130/pdf).
No literature novelty is claimed. The root and incidence lanes independently
checked the orbit and valuation bookkeeping before the subagents reached
their usage limit. Root subsequently sharpened the numerical constant
by an exact second-moment calculation and checked finite certificates.

## 1. The weighted norm budget

For a prime p congruent to 1 modulo n, let H_p be the unique subgroup
of F_p^* of order n. Put

    R_p=(H_p-1) minus {0},
    X_p=E^times(R_p)-[2(n-1)^2-(n-1)].

Thus X_p counts ordered solutions of

    (a-1)(d-1)=(b-1)(c-1), a,b,c,d in H_p minus {1},

excluding the two trivial pair matchings. Define

    M=(n-1)^4-2(n-1)^2+(n-1),
    S=8n^2(n-1)^2-2n^4.

Then

    sum_(p prime, p=1 mod n) X_p log p
       <= (M/2) log(S/M)
       <= M log 6.                                             (1)

The sums have finite support for every fixed n. No prime-counting
theorem or density assumption enters (1).

Let zeta be a primitive complex n-th root and let A be the M ordered
exponent tuples (a,b,c,d), each entry in {1,...,n-1}, excluding the
two pair matchings. For t in A put

    F_t=(zeta^a-1)(zeta^d-1)-(zeta^b-1)(zeta^c-1).

Every F_t is nonzero. Indeed, for a complex unit u,
bar(u-1)=-(u-1)/u. If two nonzero products of shifted units are
equal, conjugation forces their unshifted products to be equal.
The original equality then forces their sums to be equal. Their
unordered pairs are therefore the roots of the same quadratic
polynomial, which is exactly an excluded trivial matching.

The product

    P_n=product_(t in A) F_t

is a nonzero algebraic integer fixed by every cyclotomic Galois
automorphism: each automorphism permutes the whole tuple set A.
Consequently P_n is a nonzero rational integer. Expanding F_t gives
six unit monomials, so |P_n|<=6^M.

For p=1 mod n, the cyclotomic polynomial splits into distinct roots
modulo p. Fix one primitive root and its p-adic lift. Each of the X_p
vanishing tuple factors is divisible by p in that embedding, giving

    v_p(P_n)>=X_p.                                             (2)

Equivalently, use integer multiplication determinants for every F_t.
Their total product is |P_n|^(phi(n)); summing the nullities over
all primitive roots gives phi(n)*X_p. Divide phi(n) exactly once.
This equivalent proof avoids using any unproved fixed-root converse
from a single norm divisibility statement.

Summing (2) over p proves the weaker budget in (1). To sharpen it,
write r_a=zeta^a-1. Orthogonality gives

    sum_a r_a=-n, sum_a |r_a|^2=2n.

Summing over every ordered tuple, including the trivial tuples whose
F_t is zero, yields

    sum_t |r_a r_d-r_b r_c|^2
       =2(n-1)^2 (sum_a |r_a|^2)^2-2|sum_a r_a|^4
       =S.

The arithmetic-geometric mean inequality on the M positive squared
absolute values gives |P_n|^2<=(S/M)^M. Together with (2), this is
the stronger budget in (1). In particular the right-hand side is
asymptotic to (log 6)n^4/2.

## 2. An absolute exceptional-prime count for the triangle criterion

For any L>1 and T>0, (1) immediately gives

    #{p>=L: p=1 mod n, X_p>T}
       <= M log(S/M)/(2T log L).                              (3)

Take L=c n^4, with fixed c>0 and n large enough that L>1, and
T=n^(61/31). The number of such exceptions in any quartic interval
[c n^4,C n^4] is

    O_c(n^(63/31)/log n).                                     (4)

For every other eligible prime in that interval, the purely
combinatorial implication from
[pass41](parallel41-level-refinement-2026-09-06.md) yields

    W(H_p)<<n^(17/3).                                         (5)

This is an unconditional absolute exception bound for the intermediate
triangle estimate. A proportion of eligible primes requires an additional
denominator estimate; none is assumed here. More importantly, (3)-(4)
leave individual exceptional primes possible. A single X_p of quadratic
size at a quartic prime costs only O(n^2 log n) of the O(n^4) budget.
The budget therefore does not exclude the obstruction needed for a
uniform theorem.

The improvement over applying the older general sixth-relation budget
comes from preserving the shifted-product compatibility. There are only
M=O(n^4) compatible exponent tuples here. A generic six-term relation
argument counts a larger family and loses this information. This is a
scoped specialization of the existing method, not a new principle that
controls logarithmic-depth moments or every exceptional prime.

## 3. Complete fixed-order certificates

The [exact checker](../experiments/parallel42_norm_budget.py) factors the
product of all nontrivial tuple norms for n=4,8,16. It compresses ordered
tuples by unordered shifted products and their exact multiplicities.
Dividing every aggregate norm valuation by n/2 recovers the factorization
of |P_n|. Every split prime dividing P_n is then checked independently
using direct ratios in (H_p-1) minus {0}.

The complete split-prime exceptional lists are:

| n | Primes with X_p>0 |
|---|---|
| 4 | 5 |
| 8 | 17, 41 |
| 16 | 17, 97, 113, 193, 241, 257, 337, 353, 401, 433, 449, 577, 593, 641, 881, 1217, 2113, 2129 |

Completeness follows from factoring every tuple norm, not from a search
up to a numerical cutoff. The program also checks extra split primes,
446 independent multiplication determinants, the exact sum S, and both
integer product inequalities. The [certificate](../results/parallel42_norm_budget_2026_09_06.json)
stores all prime valuations, including a strict valuation loss at n=16,
p=17: X_p=2730 while v_p(P_n)=2856. Thus the inequality in (2) must
not be silently replaced by equality.

These are circularity exceptions. They are different from the shorter
lists of primes with excess additive quadruples in the earlier note.

## 4. Off-diagonal collisions can survive a zero diagonal count

One resulting exact example is

    p=353, n=16, H=<304>.

It satisfies n^2<p. Its positive coset multiplicities a(C) consist of
one value 1 and seven values 2. Thus every diagonal rich-cell term is
zero. Nevertheless,

    X=X_dist=72, B=507,
    W=364800, W_(rho>=3)=82944.

There are twelve rich cells, each with rho=3 and T_b=144; all have
three distinct edge cosets. Therefore no inequality

    X_dist <= C*(X-X_dist)

can hold with any finite constant C for all subgroups in n<sqrt(p).
The example lies outside the narrower quartic window, so it does not
alone refute that window-restricted comparison. It is not a Paley
counterexample. Its role is to rule out a proposed general diagonal
domination step and to retain the off-diagonal arithmetic in future
arguments. Direct full-field counts and the shifted-ratio calculation
independently reproduce it.
