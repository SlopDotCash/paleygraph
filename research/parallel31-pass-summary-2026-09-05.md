# Thirty-first research assessment

**The full Paley conjecture and the Proximity Prize remain unproved.**
This pass determines the minimum number of triangle equations needed
to generate each of the 119 six-term remainder orbits in the pass30
field. It does not improve an upper bound or cancellation exponent.

The [ordinary proof](parallel31-short-multiples-2026-09-05.md) uses the
unique integer quotient F_S/f in Z[X]/(X^64+1), with f=1+X+X^19.
The absolute sum of its coefficients is exactly the minimum number of
scaled zero-sum triangles whose entries cancel down to S. This follows
from injectivity and coefficient counting, not from a search that failed
to find a shorter derivation.

The complete classification at p=215535361, n=128 is:

| Type | Minimum triangle equations | Scaling orbits |
|---|---:|---:|
| Has a zero-triple partition | 2 | 54 |
| Primitive | 4 | 20 |
| Primitive | 38 | 1 |
| Primitive | 40 | 9 |
| Primitive | 42 | 17 |
| Primitive | 44 | 18 |

For example, the primitive exponent set {0,8,16,45,62,100} needs at
least 38 triangle terms. Its quotient is provided explicitly. No use of
other scales, repetitions, or negative triangles shortens it, since the
quotient is unique. Thus 45 primitive orbits cannot be explained by a
few triangle identities in this formal cancellation model.

For a primitive six-set, every minimum cancellation graph is connected
and has cycle rank L/2-2, where L is the minimum triangle count. The
20 short primitive orbits have tree derivations, while the other 45
require cycle rank 17 through 20. This is a finite structural result,
not an infinite obstruction or a statement about the separate spectral
necklace graphs.

The [classifier](../results/parallel31_short_multiples_2026_09_05.json)
and [separate certificate checker](../results/parallel31_certificate_check_2026_09_05.json)
pass 2,174 checks in total. The latter uses direct trinomial row equations,
tests every proper subset for zero sum, and verifies the field arithmetic
and product partitions independently. All 119 orbits account for exactly
the full D6 count previously checked by sparse and C++ enumeration.
Both implementations are by root; separate-author review and Lean
verification remain outstanding. No new Lean or hosted proof job ran.

The initial higher-energy check found no stronger input through the
specific direct substitutions examined: the source's repeated-difference
energies have different definitions, while coset norm splitting and
one-step Cauchy interpolation give only the stated weaker estimates.
No claim is made that every possible higher-energy approach has been
exhausted.

The previous turn is progress. This pass adds exact minimum-length
certificates and shows why knowing the generating triangle alone does
not control the primitive six-term count by a few local patterns.
The next mathematical input must count short outputs of these larger
cancellation networks, or control the relevant spectral quantity by
another method. The uniform subgroup, classical Paley, and prize
bounds remain open in this project; the goal stays active.
