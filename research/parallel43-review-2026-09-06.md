# Pass 43 proof and certificate review

Root completed this review. The previous parallel lanes were unavailable
after an account usage-limit failure. Independent computational methods
below mean distinct calculations, not a separate reviewer. The ordinary
proofs are not Lean-certified or externally peer reviewed.

## Energy-based exceptional set

The exact normalization is W=n sum rho*T, sum T=A2^3, E_*=nA2.
Positivity therefore gives W<=R_max E_*^3/n^2. Root visually reread
Lemma2, equation18, on printed page197 of the published Shkredov2013
source. Its three-single-coset application satisfies the actual required
ratio-image cardinality n^2. Its left side is n*rho and its right side
is O(n^(5/3)). This proves R_max=O(n^(2/3)); it does not follow by
discarding repeated fibers in arbitrary unions of cosets.

Root reread the full normalized-tuple norm and AGM argument for the
fourth-energy prime average in
[the existing note](cyclotomic-prime-average.md). The intrinsic count is
3n^2-3n, with no principal-frequency subtraction beyond that definition.
All excess terms are nonnegative. Threshold n^(7/3) gives exponent
4-7/3=5/3 for the absolute prime count. Inserting the same energy scale
into W gives 3*(7/3)-4/3=17/3. The c-dependent large-n restriction ensures
the quartic interval lies in n^2<p. It does not prove any denominator
estimate for eligible primes. The previous uniform49/20 energy input is
insufficient for this criterion; the uniform W exponent does not improve.

## Quartic witness

The exact [checker](../experiments/parallel43_pointwise_inputs.py) verifies
primality by trial division and exact order by modular exponentiation.
It computes shifted energy and every rich cell through the established
bijective product-to-triple map. A separate full-field rho count and
direct difference calculation verify W, its low/high split, the symmetry
identities, and the off-diagonal excess. Thus the zero denominator in
X_dist<=C*(X-X_dist) is an actual subgroup property inside the quartic
window, with a strictly positive numerator. The same example has minimum
additive energy and does not violate the triangle estimate itself.

The single-coset incidence normalization is independently enumerated on
one rich cell: the ratio image has n^2 elements and the retained count
is 3n. This finite test checks the implementation of the source hypotheses;
the uniform justification is the single-coset group action in the proof.

## Valuation identity

The prime-power lift has derivative n*g^(n-1) invertible modulo p, so
each lift step is unique. Every nonidentity shifted root is a unit and
distinct labels stay distinct modulo p. The two trivial pair matchings
are therefore the same at every precision. Nontrivial product classes
refine as precision increases. Summing each nonzero tuple factor's finite
p-adic valuation gives sum_r X_r exactly, with no extra cyclotomic-degree
division: mathcal_P_n is already the rational tuple product.

The checksum-verified pass42 integer norm factorization supplies 21
independent valuations at n4,8,16; the new product-hashing lifts match
every one. At n16,p17 the sequence is read from the result file and its
sum equals2856 although X_1=2730. The three larger quartic sequences
are (72,0), (114,0), and (720,0) at orders32,64,128. Once a precision
has zero excess, refinement proves all later terms zero; the stopping
condition is mathematical, not a timeout or extrapolation.

The [results](../results/parallel43_pointwise_inputs_2026_09_06.json) retain
all precision counts and the four distinct exact energy/triangle checks.
The published source PDF/image and local dependencies are
[pinned](../results/parallel43_source_scope_2026_09_06.json). No asymptotic
prime-energy saving, full Paley bound, or official prize conclusion follows
from those finite certificates. The full goal remains active.
