# Eighth pass: all length-six degree-two necklace words

**Progress: the remaining 42 fixed length-six words are now bounded.**
All 729 words have a uniform O(p^(7/2)) estimate for primes p≡1 mod 4.
This completes a finite-depth step; it does not prove the full Paley
conjecture, a growing-depth signed aggregate, or the prize reduction.

The [common contraction theorem](parallel8-necklace-2026-09-05.md)
uses a second quadratic middle convolution. The first produces a
rank-four irreducible sheaf. Either singleton-anchor twist followed
by another convolution produces rank six. Its pairing with an
irreducible rank-two Legendre sheaf cannot have an invariant.
Curve cohomology then supplies square-root cancellation in the
remaining summation variable.

The same lemma handles all three previously unresolved classes:

| Representative used in the proof | Catalog representative | Words | Uniform bound for the whole class |
|---|---|---|---|
| ABBACC | AABCCB | 18 | 48p^(7/2) |
| ABCBAC | ABACBC | 18 | 47p^(7/2) |
| ABCABC | ABCABC | 6 | 45p^(7/2) |

Two applications are direct matrix identities. The nested application
uses the previously proved adjacent-pair graph normalization and
retains its −p t₄+2 term. Both infinity corrections, the finite
spike, zero-valued character masks and all exceptional fibers are
included. The second infinity correction need not be evaluated:
its weight bound is sufficient.

The first proof draft was tested against the original matrices.
The final [verifier](../experiments/parallel8_necklace_2026_09_05.py)
completed on p=5,13,17,29,41 using integers and exact signs in Q(√p).
It records 4 monodromy cases, 1,900 finite-stalk checks, 1,140 compact
trace-vector checks, 11,180 pointwise bounds, 1,140 inner-core bounds,
3,420 formal coefficient checks, 1,140 generic and 120 exceptional
fiber checks, 15 direct trace identities and 210 orbit values.
At p=5 it also evaluates the original nested graph and all three
literal necklace sums.

The [results](../results/parallel8_necklace_2026_09_05.json) and
[integration audit](../results/parallel8_pass_audit_2026_09_05.json)
pin the exact files and source hashes. The finite checks verify
formulas and examples, not the source theorems or asymptotic claims.
The proof uses the same Katz middle-convolution and curve-cohomology
inputs audited in the preceding pass. Root developed and checked this
extension locally; the parallel workers remain interrupted by their
account usage limit.

The subgroup estimate remains the seventh-pass
M≤(17+log N)^(1/9)N^(8/9) on the existing density-one prime class.
No new subgroup or classical arbitrary-set bound is claimed here.

The next necklace obstacle is control as length grows, including
the weighted sum of all words required by the spectral argument.
Knowing every word at length six does not supply that estimate.
A useful next question is whether the rank and local-monodromy
calculation can be organized for longer chains while retaining
every punctual correction and tracking equal-rank invariant cases.
The stronger subgroup energy budget, exceptional primes, full
restricted operator and exact official prize bridge also remain open.
