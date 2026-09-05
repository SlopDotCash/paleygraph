# Thirtieth research assessment

**The full Paley conjecture and the Proximity Prize remain unproved.**
This pass finds an explicit failure of the tentative constant-one
six-term bound, and isolates how zero-sum triples contribute to it.
No cancellation exponent improves.

For p=215535361, n=128 and H=<25525303>, the
[proof and certificate](parallel30-triangle-remainder-2026-09-05.md) give

    D6=10967040=(5355/1024)n^3 > n^3.

This is a quartic-window prime subgroup. It refutes extrapolating
D6<=n^3 from the complete census through n=64. It does not refute an
unspecified O(n^3) bound, the earlier finite census, or Paley.

The counterexample also has a short proof independent of the full count:
g has order128 and 1+g+g^19=0 in the field. One zero-sum triangle orbit
produces at least 360n(n-20)=4976640 fully unbalanced ordered six-words,
already exceeding n^3. The exact remainder splits into 54 triangle-pair
scaling orbits and 65 primitive six-relation scaling orbits. The latter
alone still exceeds n^3, so deleting triangle pairs does not recover
the literal bound.

For general dyadic subgroups, write kappa=|H intersect (1-H)| and
tau=kappa-3*1_(2 in H). A unique-partition argument proves

    10n max(0,n-28) tau^2 <= D6_tri <= 10n^2 tau^2.

Consequently a uniform cubic upper bound for the full remainder would
already imply kappa=O(sqrt(n)). That intersection estimate and the
remaining primitive-relation upper bound are both unproved here.

The norm of 1+zeta_128+zeta_128^19 is exactly p. Its multiplication
matrix has determinant p, so all polynomial relations at this generator
belong to the same principal ideal. Different short relations cannot
be presumed to give independent prime divisibility constraints.

The [verification record](../results/parallel30_verification_2026_09_05.json)
contains 13,201 exact checks: trial-division primality, generator order,
the zero triangle, all good relative scalings, an independent integer
norm determinant, and agreement of triangle-pair enumeration with direct
six-word enumeration in five fields. The existing sparse count agrees
with the separate direct implementation on every full remainder category
of the witness. These are independent implementations by the same author;
separate-author review and Lean verification remain outstanding.

The prior turn is classified as progress. This pass adds an arithmetic
counterexample and a structural lower bound that change the next proof
search step. No new Lean process, hosted proof job, public submission,
or memory write occurred. Prove2Me's accepted elementary proofs are
unchanged. The full goal remains active and unachieved.

The next useful input must control the number of short relations within
these arithmetic dependencies, allowing the observed constants, or
bound the signed spectral quantity by another argument. Treating the
finite census as a uniform theorem or treating relation orbits as
independent prime factors would not close the gap.
