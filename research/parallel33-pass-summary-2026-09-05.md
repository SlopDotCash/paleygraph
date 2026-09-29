# Thirty-third research assessment

**The full Paley conjecture and the Proximity Prize remain unproved.**
This pass makes two improvements to the treatment of repeated entries.

First, the [ordinary proof](parallel33-centered-distinct-moments-2026-09-05.md)
shows that the centered distinct-coordinate count approximates the full
centered subgroup moment with relative error O(s^2/sqrt(n)) at depth s.
The exact error is a finite product, and the estimate is uniform for
s=O(log n) in the quartic window. The two counts use their respective
principal terms: (n)_(2s)/p for distinct tuples and n^(2s)/p for all
tuples. Inclusion-exclusion subtracts the principal term of every
collision partition before applying Holder. This avoids leaving a large
uncentered collision term in the centered target.

Consequently a Gaussian upper bound for the centered distinct-coordinate
count would suffice for the existing subgroup moment target at logarithmic
depth. That upper bound is still missing; distinct tuples may contain
opposite pairs, and they are not being identified with the narrower D6
remainder.

Second, subgroup averaging gives R_6<=15 sqrt(E_2 E_3) for all repeated
zero-sum six-words. The previously used energy estimates yield

    R_6 << n^(129/40) (log n)^(3/5),

improving the earlier n^(69/20)(log n)^(1/5) bound on the opposite-free
repeated subset. This saves 9/40 in that component's exponent. It does
not improve the full sixth-energy or period exponent.

The [exact verifier](../results/parallel33_verification_2026_09_05.json)
passes 740 checks across 33 set/order cases in about0.63seconds. A
subset dynamic program independently verifies the partition-inversion
counts. The checks include a characteristic-three counterexample
outside the required 2s<p hypothesis and exact rational parameter
bounds at logarithmic depth. The constants can be large at modest
orders; no practical small-order approximation is asserted when the
explicit error exceeds one.

The [recent-source check](parallel33-source-check-2026-09-05.md)
found that dense relative uniformity fails directly on the original
sparse subgroup because of diagonal terms. A different Fourier-versus-
Bohr theorem applies in size but leaves an additional intersection
estimate. Neither source supplies a new upper bound in this pass.

The previous turn is progress through the uniform triangle-length
bound. This pass is progress through an improved repetition estimate
and an exact centered reduction valid at growing moment depth.
Both arguments have author review and finite checks; separate-author
review and Lean verification remain open. No new Lean or hosted proof
job ran. The remaining centered distinct-coordinate upper bound,
classical signed spectral problem, and scalar prize bounds are unproved.
The full goal remains active and unachieved.
