# Tenth pass: all-word quadratic control and its dimension limit

**Progress: there is now a simultaneous bound for arbitrary signed combinations
of the original word kernels, together with a proved limit on this approach.**
The [written argument](parallel10-word-aggregate-2026-09-05.md) is
source-dependent and checked locally. Independent mathematical review
remains outstanding. The full Paley conjecture and prize goal remain open.

The rank-growth argument from pass nine can be run backwards. Applying
inverse middle convolution reveals the last label through the signs of
local Jordan-block count differences. Removing that label and repeating
recovers the whole word. Different words therefore give geometrically
inequivalent principal sheaves. Arithmetic self-duality fixes the diagonal
correlation term; curve cohomology bounds every off-diagonal term.

This gives a uniform quadratic estimate over a complete word family,
valid for every complex coefficient vector. The estimate also controls
restriction to any further disjoint set of prescribed adjacency conditions.
The raw character paths are retained: their difference from the principal
trace functions is bounded using the preceding strict weight gap.

An exact finite recurrence now tracks generic raw ranks, including all
boundary contributions. For two anchors, the sum Q_m of squared raw ranks
over all length-m words begins 1,17,199,2001 and grows like
((9+√65)/2)^m. The general growth rate Λ_a is explicit. The whole
word family is an asymptotic isometry when Λ_a^m=o(√p). With s fresh
adjacency conditions, a sufficient relative-error condition is
2^s(a+s)Λ_a^m=o(√p), for fixed a. This extends the earlier bounds
for a span indexed by one kernel rank to all label sequences at once.

There is also an unconditional obstruction to extrapolating this
isometry to the depth needed for the spectral edge. The word map has
(2^a−1)^m columns but only p−a−1 available field points. Once there
are more columns than rows, a nonzero coefficient vector evaluates to
zero. No relative isometry error below one can then hold for every
coefficient vector. The checker supplies an exact integer null vector
for the original paths at p=13, with 27 words and 10 evaluation points.
This does not refute cancellation for the particular spectral coefficients.

The next step must use those particular coefficients, their algebraic
relations, or an appropriate quotient of word space. An unrestricted
all-coefficient isometry cannot itself settle the growing-depth spectral
criterion. The full cyclic trace also still needs its rank-one J terms
and exceptional anchor fibers handled together.

The [verifier](../experiments/parallel10_word_aggregate_2026_09_05.py)
and [results](../results/parallel10_word_aggregate_2026_09_05.json)
record 7,541 word recoveries and distinct local signatures, 773 full
numerical constituent-inventory comparisons, 114 integer convolutions,
573 literal path entries and 172 zero masks. Seventy exact
positive-definiteness certificates check 35 Gram or cell bounds;
six have relative error below one. The separate integer null vector
certifies the dimension obstruction on actual paths. Finite checks
do not prove the imported sheaf theory or the asymptotic statements.
The [artifact audit](../results/parallel10_pass_audit_2026_09_05.json)
pins final inputs and confirms preservation of pass nine.

The subgroup exponent remains 8/9 on the existing density-one quartic
prime class. The uniform subgroup target, exceptional primes, stronger
centered energies, classical arbitrary-two-set conjecture and exact
official prize bridge are unchanged. Root completed this pass locally;
the previously launched workers remain terminal at their account limit.
