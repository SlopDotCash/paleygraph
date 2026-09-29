# Sixteenth pass: the sharp inputs still need a prime-field estimate

**Progress: a quantified field restriction for the current spectral route.**
We have not improved the prime-field clique bound or proved the full goal.
The new result explains why combining the sharp elliptic correlations
and ambient projection identity, as proposed after pass fifteen, is
insufficient without another input using prime field order.

The [derivation](parallel16-square-field-barrier-2026-09-05.md) uses actual
Paley graphs over F_q with q=ℓ². The subfield F_ℓ is a clique of size √q;
this is classical. Its centered indicator is an exact positive Paley
eigenvector, with squared-mass fraction 1−1/ℓ on the subfield. After
localizing at any fixed a anchors in that subfield, the centered norm
converges to 2^a/(2√(2^a−1)). For a≥2 this exceeds the desired edge 1.
At three anchors the limit is 4/√7. These outliers persist as both the
field size and characteristic tend to infinity.

These same examples retain the actual ±1 Paley entries, full projection
identity, complete elliptic kernel and exceptional border, original
four-puncture Kummer realization, and sharp all-order twisted-correlation
constants. Thus they retain the two separate constraints lost by the
abstract constructions of passes twelve and fifteen. No signs or
projections are modified. The elliptic character-sum argument works over
finite fields F_q and still needs independent mathematical review.

The obstruction also reaches the exact transfer matrix. Its spectral
radius tends to √(2^a−1). For a=3 and ℓ≥29, the full energy at depth j
exceeds (63/16)^j, and the even signed trace exceeds (63/16)^j−q.
Consequently neither supplies the required long-depth estimate on this
family. A theorem uniform over square fields cannot close this route.

This is **not a counterexample to the prime-field conjecture**. The
[primary spectral formulation](https://arxiv.org/html/2303.16475v1#S2)
explicitly takes its asymptotic parameter through primes. The elementary
subfield clique is established background, not a novelty claim. Our
contribution here is checking its consequences for the exact collection
of inputs and normalizations being considered in this project.

The [verifier](../experiments/parallel16_square_fields_2026_09_05.py)
constructs five quadratic fields, with orders 49,121,289,361,841, using
exact arithmetic. It checks 938,165 entries of each ambient identity,
25 localizations, 4,983 subfield-eigenvector entries and 12 exact power
traces. Three complete elliptic examples check 266,496 kernel entries,
16,656 squared-kernel entries, and 12 sharp certificates covering all
character twists for selected shift supports. Some sharp bounds attain
equality, so the certificates explicitly handle singular positive
semidefinite matrices. See [results](../results/parallel16_square_fields_2026_09_05.json)
and the [artifact audit](../results/parallel16_pass_audit_2026_09_05.json).

The next substantive step is a quantitative prime-field estimate that
excludes this concentration mechanism. Qualitative absence of a proper
subfield alone is not that estimate. The earlier logarithmic depth range,
classical arbitrary-set cancellation, uniform subgroup target and
exceptional primes, and exact official Reed–Solomon prize bridge remain
unresolved. The final live check confirms that all three parallel workers
are terminal at account usage limits; this pass ran locally.
