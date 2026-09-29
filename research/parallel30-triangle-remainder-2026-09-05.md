# Zero-sum triples and a counterexample to D6 <= n^3

**The full Paley conjecture and Proximity Prize remain unproved.**
This note disproves the literal extrapolation D6<=n^3 from the complete
small-order census. It does not disprove D6=O(n^3), a bound with a
logarithmic factor, or either Paley target. The uniform upper-bound
problem remains open. The argument is ordinary mathematics with exact
finite checks, not a Lean theorem or a claim of literature novelty.

The explicit quartic-window example is

    p=215535361, n=128, H=<25525303> in F_p^*.

Two distinct counting implementations give

    D6=10967040=(5355/1024)n^3 > n^3.

A short triangle argument already supplies a lower bound exceeding n^3,
independently of the full six-word enumeration. The observed inequality
in [pass 28](parallel28-pass-summary-2026-09-05.md) remains true for all
28,774 cases actually checked at orders 4 through 64; the new example
shows why it cannot be promoted to an all-order theorem with constant one.

## 1. Distinct zero-sum triangles

Let p>3 be prime, H<=F_p^* have dyadic order n>=4, and define

    kappa=#{x in H : 1-x in H},
    tau=kappa-3*1_(2 in H).

The number of ordered zero-sum triples in H is n*kappa: normalize the
third entry to -1. Exactly 3n of these have repeated entries when 2 is
in H, and none do otherwise. Indeed their form is (a,a,-2a) in one
of three positions, and all three entries cannot coincide since p>3.
Thus the number of unordered, distinct zero-sum triples is

    |T|=n*tau/6.                                         (1)

Such a triple has no opposite pair. Its multiplicative scaling
stabilizer in H has order dividing both 3 and n: its action on the
three nonzero elements is free. Hence the stabilizer is trivial,
every scaling orbit has size n, and

    a=|T/H|=tau/6 is a nonnegative integer.                (2)

## 2. A six-set has at most one zero-triple partition

Suppose a set of six distinct elements has sum zero. If two different
unordered partitions into triples both have zero sums, choose the
side of each partition so that their intersection has two elements.
Subtracting the two three-term sums then equates their different third
elements, a contradiction. Therefore the zero-triple partition is unique.

For an opposite-free six-set of nonzero elements, every proper nonempty
zero-sum subset would have size three. Sizes one and five are excluded
by nonzero entries, and sizes two and four by the absence of opposite
pairs. Consequently the six-distinct, opposite-free, fully unbalanced
count D6 splits exactly into

    D6 = D6_tri + D6_primitive,                           (3)

where D6_tri contains one zero-triple partition, and D6_primitive has
no proper nonempty zero-sum subset. Both count ordered six-words.
This is a distinction within the actual remainder, not a new name for
the whole unknown energy.

## 3. Pairing triangle orbits

Choose one representative A for each of the a triangle scaling orbits.
For two representatives A,B and t in H, consider A union tB. At most
nine values of t cause an overlap and at most nine cause an opposite
pair: the respective forbidden ratios are x/y and -x/y for x in A,y in B.

There are ten unordered three-versus-three partitions of the six labeled
positions. For one partition, let j of its three positions come from tB.
Product balance has the form

    t^j alpha+t^(3-j) beta=0,
    t^(2j-3)=-beta/alpha,

with alpha,beta in H. The exponent is one of -3,-1,1,3 and is coprime
to the dyadic order n. Therefore exactly one t in H solves each such
equation. At most ten further relative scalings are product-balanced.
There are at least max(0,n-28) good relative scalings for each ordered
pair of triangle representatives.

Let G(A,B) denote their exact number. Every ordered pair of actual
triangles from the two orbits can be written uniquely as (sA,stB),
with s,t in H. The zero-triple partition uniqueness shows that every
good six-set is counted by exactly two ordered triangle pairs.
Multiplying by 6!=720 gives the exact identity

    D6_tri=360 n sum_(A,B) G(A,B).                       (4)

Together with (2), this proves

    10 n max(0,n-28) tau^2 <= D6_tri <= 10 n^2 tau^2.    (5)

In particular, for n>=56 a uniform D6=O(n^3) bound would imply
kappa=O(sqrt(n)). For the triangular contribution alone the implication
is an equivalence, by both sides of (5). More precisely,

    D6_tri=10 n^2 tau^2+O(n tau^2),

with absolute error at most 280 n tau^2. The six-term remainder thus
already contains the square-root-size intersection question for H and
1-H. That intersection estimate is not proved by this note.

## 4. A sharper bound from one triangle orbit

For A=B, define R=A/A={x/y:x,y in A}. It has at most seven elements.
Overlap or an opposite pair occurs exactly for t in R union (-R).
For a mixed product partition, choose the element x of A on the side
with only one A-element and the element y of the other copy on the side
with only one tA-element. Cancellation of the total product of A gives

    t=-(y/x)^2.

The unmixed partition requires t^3=-1, whose unique solution in H is
t=-1. It is included in -R^2, since 1 is in R. Thus all bad relative
scalings are exactly

    R union (-R) union (-R^2).                           (6)

Each set has at most seven elements, and -1 belongs to the last two.
The union therefore has size at most twenty. A single triangle orbit
already gives

    D6 >= 360 n max(0,n-20).                             (7)

This follows either directly from its good ordered triangle pairs, or
from the corresponding diagonal term of (4). No estimate for kappa or
for the remaining triangle orbits is needed to use this lower bound.

## 5. The explicit field and short certificate

The [certificate](../results/parallel30_verification_2026_09_05.json)
checks primality of p=215535361 by trial division through floor(sqrt(p)).
It also checks

    67108864=n^4/4 < p < n^4=268435456,
    (p-1)/n=1683870,
    g=25525303,
    g^64=-1 mod p,    g^128=1 mod p,
    g^19=190010057 mod p,
    1+g+g^19=p.

Hence g has exact order 128 and A={1,g,g^19} is a distinct zero-sum
triangle in H. Formula (7) alone proves

    D6>=360*128*108=4976640 >128^3=2097152.              (8)

This disproves the constant-one bound without enumerating all six-words.
The arithmetic discovery came from
Norm_(Q(zeta_128)/Q)(1+zeta_128+zeta_128^19)=215535361.
An independent 64-by-64 integer multiplication determinant verifies
that norm. The norm search was exploratory; the field and triangle
checks in this certificate establish the counterexample directly.
No completeness assertion about the search's factor candidates is needed.

There is also an exact arithmetic dependence among the relation
polynomials in this example. Put R=Z[X]/(X^64+1) and f=1+X+X^19.
Evaluation at g maps R onto F_p, with kernel of additive index p.
Since f(g)=0, the principal ideal fR is contained in that kernel.
Multiplication by f has integer determinant of absolute value p, so
fR also has additive index p. The two subgroups are therefore equal.
Every polynomial relation at this generator, including each primitive
six-term relation, is an integral multiple of f in R. In particular,
counting these different relations as independent prime divisibility
constraints would be unjustified. This does not give a bound on how
many short multiples fR contains.

For the exact triangular count, the exponents of R in this generator are

    D={0,1,18,19,109,110,127} mod 128.

The bad exponent set is D union (64+D) union (64+2D), which has exactly
twenty elements. All 108 remaining scalings are checked directly against
distinctness, opposite-freeness and every product partition. The exact
intersection count is kappa=6 and 2 is not in H, so tau=6. Thus there
is exactly one triangle orbit, and equality holds in (8) for D6_tri.

The full finite decomposition is:

| Ordered-word category | Count |
|---|---:|
| Intrinsic three-opposite-pair words, T6 | 30,725,120 |
| Nonintrinsic words containing an opposite pair, J6 | 22,855,680 |
| Opposite-free repeated words | 2,438,400 |
| Six-distinct opposite-free product-balanced words | 1,751,040 |
| Fully unbalanced words with a zero-triple partition | 4,976,640 |
| Fully unbalanced words with no proper zero-sum subset | 5,990,400 |
| Total E3 | 68,737,280 |

The last two remainder categories have 54 and 65 free scaling orbits
of unordered six-sets, respectively, for 119 in total. Their ratios are

    D6_tri/n^3=1215/512,
    D6_primitive/n^3=2925/1024,
    D6/n^3=5355/1024.

Therefore deleting the triangular contribution would still not restore
the literal constant-one bound in this example. The primitive uniform
upper-bound problem remains open as well.

## 6. Verification and scope

The existing sparse energy decomposition and the existing direct
normalized-six-word C++ enumeration first agreed on the full D6 count.
The new [direct program](../experiments/parallel30_direct_six.cpp) also
classifies the zero-triple partition by testing its actual sums, without
using the orbit formula. It agrees with a separate enumeration of all
triangle pairs in five fields. The [verifier](../experiments/parallel30_verify_2026_09_05.py)
checks the general counting identities and bounds in those fields, and
checks the explicit arithmetic and norm determinant exactly.

These are independent implementations by the same author, not an
independent-author review. No new Lean job or Prove2Me submission was
launched. The ordinary counting arguments remain available for review
in this note; their asymptotic scope is not certified by finite sampling.

The useful change to the frontier is a concrete rejection of the
constant-one extrapolation, plus a quantified necessary intersection
estimate and an exact primitive-relation remainder. No upper-bound
exponent improves. Neither this finite witness nor its lower bound
refutes an unspecified cubic constant or proves or disproves Paley.
