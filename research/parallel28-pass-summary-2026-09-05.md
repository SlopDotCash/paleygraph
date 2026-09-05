# Twenty-eighth research assessment

**The full Paley conjecture and the Proximity Prize remain unproved.**
This pass improves the subgroup cancellation baseline recorded in this
project by combining existing analytic estimates. It also completes the
small-order census of the distinct six-term remainder.

For a multiplicative subgroup H of order n>=4 in a prime field, with
n^4/4<=p<=n^4, the [ordinary proof](parallel28-moment-recurrence-2026-09-05.md)
derives

    max_(a != 0) |sum_(h in H) exp(2 pi i a h/p)|
        << n^(71/72) (log n)^(1/18).

It combines the known third-energy bound, an invariant-set energy
estimate, and the mixed (3,6) Konyagin inequality. A nonnegative weighted
version of the energy estimate removes a logarithmic loss from the
recurrence; the proof explicitly retains the mass at zero. The typeset
published recurrence independently gives the same power with a weaker
logarithmic factor. The resulting exponent is smaller by 1/320 than
2849/2880, the prior project baseline. This is a derived consequence of
literature inputs. No claim of literature novelty or current-best
status is made. No separate-author review or Lean verification of this
analytic proof has been completed.

The [source ledger](../results/parallel28_source_scope_2026_09_05.json)
pins the published and preprint formulas and the two existing inputs.
The published formula numbers differ from the preprint. Complete relevant
typeset pages were visually inspected; the recurrence's size condition
applies to its iteration, not the one-step inequality used here.

The exponent ledger also checks that doubling and interpolating these
energy inputs, followed by the same mixed-moment inequality, supplies a
maximum saving of 1/72. This is a precisely limited statement about those
estimates. It does not rule out using additional information to do better.
The desired subgroup power remains 1/2+o(1), and neither a general
two-set Paley estimate nor the prize certificate follows from this bound.

## Complete finite remainder census

The [exact census](../results/parallel28_complete_remainder_2026_09_05.json)
enumerates every eligible prime in the quartic window at each listed
order, with an independent sieve and exact trial-division primality checks:

| n | Eligible primes | Maximum D6 | Maximum D6/n^3 |
|---|---:|---:|---:|
| 4 | 16 | 0 | 0 |
| 8 | 95 | 0 | 0 |
| 16 | 579 | 0 | 0 |
| 32 | 3,705 | 23,040 | 45/64 |
| 64 | 24,379 | 184,320 | 45/64 |

D6 counts the distinct, opposite-free zero-sum six-words that are
unbalanced at every three-versus-three product split. The 28,774 cases
all satisfy D6<=n^3; this does not prove that assertion at growing n.
All 28,753 overlapping E2/E3 values agree with the earlier sigma data.
The smaller old energy windows omitted 21 of the n=4 and n=8 cases;
the new census includes them.

At n=32 exactly eleven primes have nonzero D6, each with one scaling
orbit. At n=64 the counts of primes with 0,1,2,3,4 scaling orbits are
22,650, 1,536, 169, 20, and 4. The maximum is attained at
p=4,885,633; 6,356,033; 7,636,481; and 13,928,833. The previous high-E3
exceptions at p=6,700,417 and p=7,204,033 have D6=0. Thus the largest
full energies are not automatically the hardest cases for this remainder.

A separate C++ enumeration agrees in 75 selected fields, covering every
observed tuple of category counts and all positive-D6 maximizers.
It directly tests the six entries and ten product splits. It is an
independent implementation by the same author, not separate-author review.
The [verification record](../results/parallel28_verification_2026_09_05.json)
also checks exact weighted examples, all invariant sets in seven small
fields, the origin contribution, and the exponent arithmetic. These finite
tests do not prove the uniform analytic statements.

The previous goal turn is classified as progress. This turn changes the
quantitative subgroup baseline and replaces selected small-order examples
with a complete finite census. The full objective remains active and
unachieved. Prove2Me's five accepted elementary proofs remain recorded;
this analytic argument has no new hosted verification verdict. The three
workers last reported terminal usage-limit errors, so this pass was
completed by root. No model change, reset, credit purchase, public release,
mission launch, community message, or memory write occurred.

The next proof-search step is to scrutinize the weighted lifting and
mixed-moment combination independently, then seek arithmetic information
beyond the now quantified recurrence envelope. The exceptional D6 bound,
the classical signed correlations, the full spectral estimate, and the
scalar prize bounds remain open.
