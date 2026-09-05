# Twenty-second pass: parallel results and stronger obstruction checks

The three parallel research lanes completed concrete results, and the
root developed and checked a stronger weighted obstruction. **No full
Paley proof, improved all-vector spectral bound, new cancellation
exponent, or official prize implication follows.**

## Actual spectral progress: one specified direction is controlled

The [spectral note](parallel21-spectral-next-input-2026-09-05.md)
proves a uniform statement for actual two-anchor clique neighborhoods
in prime fields p≡1 mod4. For their normalized uniform vector, the
fraction of total Fourier mass on nonzero quadratic-residue frequencies
is 3/8+O(p^(−1/2)), and belongs to [1/4,1/2) for p≥73.

An exact Jacobi-sum evaluation gives the neighborhood quadratic form,
with both anchor corrections retained. This removes the exponentially
small determinant-certificate loss on this one direction. It does not
control the remaining vectors, the coupling to them, or the complete
operator spectrum. The note writes that coupling as an exact
variance of cubic character traces.

The [spectral verifier](../experiments/parallel21_spectral_flat_direction_2026_09_05.py)
and [results](../results/parallel21_spectral_flat_direction_2026_09_05.json)
check 79 actual prime neighborhoods, 1,742 affine edges, 1,418
quadratic-trace identities and 31,822 quadratic-map fibres.
A [separate-agent spectral review](parallel22-spectral-independent-review-2026-09-05.md)
verified the algebra and found an omitted normalization input: naming
the positive projection as residue frequencies requires the positive
quadratic Gauss evaluation. The review supplies it from
[NIST DLMF 20.11.2](https://dlmf.nist.gov/20.11.E2). Root verified
the formula and amended the spectral note's source scope; its formulas
are unchanged. The reviewed original bytes are
[preserved](../results/parallel22_spectral_reviewed_snapshot_2026_09_05.json).

## Classical target: an exact quartic-correlation energy interface

The [classical note](parallel21-positive-upper-review-2026-09-05.md)
uses the full field-translate identity to express the missing sixth
aggregate through squarefree quartic sums. For n≥6 and n⁴≤p, put

    U_C(t)=Σ_(Q⊂C, |Q|=3) K(Q⊔{t}),
    L_out(C)=Σ_(t∉C) U_C(t)².

Then

    L_out≤20pT_6+p²n³,       20pT_6≤L_out+3p²n³.

Thus a uniform bound L_out=O(p²n³) would supply the missing positive
T_6=O(pn³) bound. It has not been proved. Correlating two disjoint
quartic families returns the original six-root character sum, so the
lower degree of each individual summand does not remove the difficulty.

The [verifier](../experiments/parallel21_positive_upper_review_2026_09_05.py)
and [results](../results/parallel21_positive_upper_review_2026_09_05.json)
cover 216 actual field/set cases, five Sidon size-slice examples and
1,600 pairwise Gram identities. A
[separate-agent review](parallel22-quartic-star-independent-review-2026-09-05.md)
found no mathematical defect and confirmed both directions and the
absence of a circular upper-bound claim. The convenient Sidon transfer
threshold remains n≥25; the smaller numerical fixtures do not certify
that transfer. A single sixth-moment estimate would not prove full Paley.

## Subgroup target: isolate the remaining six-term relations

The [subgroup note](parallel21-subgroup-next-input-2026-09-05.md)
separates intrinsic opposite-pair relations, extensions of four-term
relations, and six-term relations with no opposite pair. It proves an
exact decomposition, including repeated-entry corrections.

If the fourth energy is intrinsic, the sixth energy is exactly its
intrinsic count plus the opposite-free six-term count R_6. A uniform
R_6=O(n³ polylog n), together with the stated fourth-energy control,
would extend the existing conditional 8/9 exponent to the primes
satisfying those inputs. This upper estimate remains unproved and
8/9 remains short of the square-root target.

All fourteen recorded order-64 exceptions have fourth-energy excess
24n. Their sixth-energy excess is mainly forced by extending those
four-term relations, and every one has R_6<2n³. These are finite
facts, not an extrapolation to every subgroup.

A separate construction has dyadic size in the same quartic prime
regime and exactly intrinsic fourth energy, but sixth energy at least
n⁴/20. These actual symmetric sets are not multiplicative subgroups.
They show why a proof cannot discard multiplicative closure.

The [verifier](../experiments/parallel21_subgroup_next_input_2026_09_05.py),
[results](../results/parallel21_subgroup_next_input_2026_09_05.json) and
[lane audit](../results/parallel21_subgroup_next_input_audit_2026_09_05.json)
cover 44 subgroup fixtures, nine other symmetric sets, five constructions
and 86,922 admissible terms in an independent opposite-free enumeration.
A [separate-agent subgroup review](parallel22-subgroup-independent-review-2026-09-05.md)
found no mathematical correction. It also bounds the constructed
non-subgroups' sixth energy above by (5/2)n⁴, confirming their
fourth-power scale. This is not an upper bound for subgroups.

## Root result: the weighted obstruction survives every individual order

The [new proof](parallel22-all-orders-obstruction-2026-09-05.md)
constructs a positive rational measure on sign rows at every fixed
moment order 2r≥6, with k=⌊p^(1/(r+1))⌋. Its cube component has
density between 1/2 and 3/2. It has the exact Gram matrix pI−J
and satisfies every individual distinct-product Weil-shaped bound,
simultaneously through degree k.

For every real coefficient vector it also has Gaussian-order bounds
at all lower even moments through 2r−2. Nevertheless the next
moment, for the sum of all columns, divided by pk^r is at least
k/2 and tends to infinity. The columns can have B_h labels throughout
the permitted SS-B* relation range.

This measure is rationally weighted and does not obey actual
field-translate or character-factorization identities. It therefore
rules out a bootstrap from just the retained relaxed constraints;
it does not refute a prime-field character estimate. Compared with
pass 21's model it retains all individual orders and all lower
coefficient-vector bounds, with different fourth-moment constants.

The [root verifier](../experiments/parallel22_all_orders_obstruction_2026_09_05.py)
and [results](../results/parallel22_all_orders_obstruction_2026_09_05.json)
check 2,300 explicit weighted rows, 508 all-subset correlations,
536 unsigned moment formulas, 332 coefficient-vector bounds and
17 compressed models on genuine prime-size slices at r=3,…,9.
They also check the binomial and label-threshold ledgers.

A [separate-agent proof review](parallel22-all-orders-independent-review-2026-09-05.md)
found no mathematical correction. Its
[independently written verifier](../experiments/parallel22_all_orders_independent_audit_2026_09_05.py)
and [results](../results/parallel22_all_orders_independent_audit_2026_09_05.json)
check 37,116 explicitly weighted rows across nine models. Its larger
fixtures use integer, nonprime parameters and are labeled accordingly;
the general prime quantifier is proved algebraically. The root's
separate large fixtures use primes. Neither evaluates a large Paley
character kernel.

## Review status and next step

Root checked the four derivations directly and the relevant primary
prime-counting input for the subgroup construction. The pinned proof
snapshots retain initial workflow-status annotations; the completed
verifier and review records linked here supersede those annotations.
Separate-agent review is not independent human certification or formal
proof verification.

The [final audit](../results/parallel22_pass_audit_2026_09_05.json)
records input hashes, review scope, exact counts, preserved pass-21
artifacts and all 26 preceding source-manifest entries. One NIST
metadata entry records the reviewed Gauss normalization; direct archival
retrieval returned 403, so no local source content hash is claimed. The late-finishing worker
notes retain their pass-21 filenames; this assessment is their first
combined root audit.

The next decisive step remains a uniform character-sensitive upper
estimate: the off-set quartic energy in the classical lane, opposite-free
sixth relations exploiting multiplicative closure in the subgroup lane,
or nonconstant vectors and their coupling in the spectral lane.
The full classical conjecture, uniform subgroup square-root target,
spectral edge and exact official Reed–Solomon prize bridge remain
unproved. The full goal is active and unachieved.

