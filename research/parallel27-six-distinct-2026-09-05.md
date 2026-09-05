# Exact isolation of the six-distinct subgroup remainder

**No uniform bound on the six-distinct remainder or full conjecture is
proved.** This pass proves that the distinct, opposite-free balanced
classes are disjoint, gives an exact sparse formula for the remaining
count D6, and checks that formula against direct enumeration. It reaches
six new quartic-window cases at orders 512 and 1024 without enumerating
all triples in the primary algorithm. No literature novelty is claimed.

Use the definitions from [pass 25](parallel25-subgroup-growing-orders-2026-09-05.md):
H is a dyadic multiplicative subgroup of F_p*, n=|H|≥4, p is odd, and
D6 counts ordered zero-sum six-words whose entries are distinct, have
no opposite pair, and are product-unbalanced at all ten unordered
three-versus-three partitions. The exact formulas below require only
these subgroup conditions. The computational complexity estimate and
new samples use n^4/4≤p≤n^4.

## 1. A distinct opposite-free word has at most one balanced partition

Call I|Iᶜ balanced if the product of the entries in I equals minus the
product of the entries in Iᶜ. Suppose two different unordered partitions
are balanced. Choose the side of each so that their intersection has
two positions: different three-subsets of a six-set have intersection
one or two unless equal or complementary, and complementation exchanges
one and two. Relabel the six nonzero entries to write the equations as

    abc + def = 0,      abd + cef = 0.

Set A=ab and B=ef. Multiplying the first equation by c and the second
by d, then subtracting, gives

    A(c²−d²)=0.

Since A≠0 in a field, c=d or c=−d. Both contradict the hypotheses.
Therefore a distinct, opposite-free word is balanced at at most one
partition. Its sum is irrelevant to this argument.

The [Lean file](../prove2me/Check_paley_balanced_split_collision.lean)
checks this algebraic implication and its no-two-splits corollary over
an arbitrary field. The combinatorial relabeling of arbitrary partitions
is the ordinary argument just given; it is not encoded in that file.
All three printed statements have only standard axioms in the
[verification record](../results/parallel27_independent_check_2026_09_05.json).

Consequently, if B* counts ordered distinct opposite-free zero-sum words
balanced at one fixed partition, then

    R6,distinct,balanced = 10 B*.                         (1)

The earlier union bound becomes an equality on this distinct subset.
Repeated words are excluded: their balanced partition classes can overlap.

## 2. Exact removal of repeated words

For a six-word h, let ν(h)=Σ_a binom(m_a(h),2), where m_a is its
multiplicity. If h is repeated, ν(h)>0. Mark an equal positional pair
and normalize its common nonzero value to one. The remaining ordered
four entries a,b,c,d lie in H and sum to −2. Each word with repetitions
is marked ν(h) times, while there are 15 possible position pairs and
n possible repeated values. Giving a marked word weight 1/ν(h) yields

    R6,rep = 15n Σ 1/ν(1,1,a,b,c,d),                    (2)

where the sum is over ordered a,b,c,d in H with a+b+c+d=−2 and with
the full six-word opposite-free. This is an equality, unlike the
earlier upper estimate using ν≥1. The normalization preserves ν and
opposite-freeness. Grouping this sum by multiplicity pattern gives exact
counts for the individual repeated patterns as well.

To evaluate (2), store pairs according to their sum. Join the bucket
at s with the bucket at −2−s. Unordered pairs can be stored with weight
one for equal coordinates and two otherwise. The product of the two
weights gives exactly the number of ordered four-tuples represented.
Use rational arithmetic for 1/ν; the final total is an integer by (2).

The number of ordered candidates is exactly

    r4(2) = Σ_s r2(s)r2(−2−s) ≤ E2(H).                 (3)

Negation symmetry identifies r4(−2) with r4(2), and Cauchy proves the
inequality. No full six-word enumeration is used in this step.

## 3. Sparse evaluation of the fixed balanced count

The [pass-23 bijection](parallel23-subgroup-upper-2026-09-05.md) gives
the normalized word for a fixed balanced partition:

    (ad, b, c, −bc, −a, −d),
    ad−a−d = bc−b−c,         a,b,c,d∈H.                (4)

Each normalized quadruple corresponds to n ordered words, one for each
common scaling. Both equations of balance and zero sum hold in (4).
Filtering (4) for distinct entries and opposite-freeness gives exactly
B*/n. Store (a,d) by z=ad−a−d and join each bucket with itself.
The number of ordered candidates is

    Σ_z w(z)² = E×(H−1),
    w(z)=#{(a,d)∈H²:ad−a−d=z}.                         (5)

Again unordered pairs with their weights preserve the exact count.
Combining with (1) gives the whole distinct-balanced count.

## 4. The full energy and exact D6 formula

For each nonzero H-coset K let W_K=Σ_{z∈K}w(z). The earlier orbit
identity gives

    E3(H) = n²w(0)² + n Σ_K W_K².                     (6)

Cosets are identified exactly by z↦z^n: its kernel in F_p* is H.
Raising instead to (p−1)/n generally identifies the wrong cosets. An
initial implementation made that error; the finite direct checks rejected
it before any result was saved. The final code uses z^n, rechecks its
kernel, and passes the independent counts below.

Let T4=3n²−3n and T6=15n³−45n²+40n. The previous exact ledger for
nonintrinsic words containing an opposite pair is

    J6 = (15n−60)(E2−T4) + 60n(r2(2)−1) −30n·1_(3∈H).

The four types are disjoint and exhaustive, so (1), (2), and (6) yield

    D6 = n²w(0)² + nΣ_K W_K² − T6 − J6 − R6,rep −10B*.  (7)

This is an exact evaluation formula. The orbit term still contains the
unbounded sixth energy. Subtracting the exactly evaluated other types
does not supply a uniform upper bound for the remainder.

There is a further exact check. A six-element subset has a scaling
stabilizer in H whose order divides both six and n. Any nontrivial
such stabilizer has order two and contains −1, contradicting
opposite-freeness. Thus its scaling orbit has size n. Distinct words
have 6!=720 orderings, so both D6 and the distinct-balanced count are
divisible by 720n. D6/(720n) counts unordered six-set scaling orbits.

## 5. Cost and finite evidence

With H supplied, the arithmetic-operation cost of (7) is

    O(n² log n + r4(2) + E×(H−1)),

using O(n²) storage. Sorting supplies the bucket operations within the
n² log n term; the implementation uses hash tables. The same logarithm
also accounts for powering the coset labels. The bound counts field
and comparison operations and does not assert unit bit
cost for arbitrarily large primes or include finding a subgroup generator.
The imported [MRSS Theorem 3](https://arxiv.org/html/1712.00410v1#Thmtheorem3)
bound E2(H)≪n^(49/20)log(n)^(1/5) applies for n≤√p. The
[shifted-energy theorem](https://arxiv.org/html/1504.04522v1#Thmsatz6)
gives E×(H−1)≪n²(1+log n) for n²<p. These hypotheses hold in the
quartic window. They therefore give O(n^(49/20)polylog n) field
operations for the count. This is a computational saving, not an
improved upper bound on the answer D6. The relevant source statements
were rechecked in primary HTML; their earlier archives are preserved.

The [sparse verifier](../experiments/parallel27_six_distinct_2026_09_05.py)
checks the partition lemma on 36,800 distinct opposite-free six-sets,
without imposing zero sum. It compares the exact formula with full
six-multiset enumeration in seven small fields, including three dense
subgroups outside the quartic window, and with all four saved
pass-25 subgroup cases. Its [results](../results/parallel27_six_distinct_2026_09_05.json)
include six new cases, chosen as the first eligible prime at or above
each of n^4/4, n^4/2, and 3n^4/4 for n=512 and 1024:

| n | p | D6 | D6/(720n) |
|---:|---:|---:|---:|
| 512 | 17179869697 | 1843200 | 5 |
| 512 | 34359753217 | 737280 | 2 |
| 512 | 51539610113 | 0 | 0 |
| 1024 | 274877908993 | 5898240 | 8 |
| 1024 | 549755860993 | 2211840 | 3 |
| 1024 | 824633730049 | 1474560 | 2 |

In all six new cases E2=T4 and J6=R6,rep=10B*=0, so E3=T6+D6
exactly. These cases have E3<15n³. Earlier cases at orders 64 and 128
already violate that literal bound; it cannot be asserted uniformly.
The chosen samples do not locate the worst quartic prime at either new
order, justify extrapolation, or estimate a proportion of exceptional
primes.

A [separate C++ implementation](../experiments/parallel27_direct_six.cpp)
directly enumerates every normalized ordered-six-word class in all 17
cases, with exact multiplicity weights and actual checks of every split.
It does not use equations (1), (2), (6), or (7) to compute its counts.
The [comparison script](../experiments/parallel27_independent_check_2026_09_05.py)
also checks the Lean algebraic core, with results linked above. Both
implementations are by root; this is independent implementation, not a
separate-author or human review. The complete combinatorial formula is
an ordinary proof with exact finite checks, not a Lean theorem.

The missing mathematical step remains a uniform upper bound on (7),
together with sufficient control of the lower-order terms. Even a cubic
sixth-energy bound alone would give only a partial maximum-period saving;
the square-root target needs substantially more. No classical Paley,
full spectral, or official-prize implication is established by these
finite counts.
