# Nineteenth pass: reducing the classical target to relation-free sets

The full classical Paley conjecture can be reduced to sets with no
nontrivial equal h-term sums, for any fixed h. For h=2 these are
additive Sidon sets. This is a proved transfer, not a proof of the
required cancellation on that restricted class.

The [derivation](parallel19-relation-free-reduction-2026-09-05.md)
gives an exact greedy extension rule. It partitions any input set
into B_h blocks of a chosen size k with at most
L_h(k)=(k−1)+Σ_{j=1}^h(k−1)^(2h−j) leftover elements.
Keeping the leftover pairs explicitly proves that cancellation for
all polynomial-sized B_h pairs would imply cancellation for arbitrary
polynomial-sized pairs. The converse follows by restriction.

This narrows the existing sufficient moment criterion. It is enough
to establish, for an unbounded sequence of fixed r_j,

M_(2r_j)(C) ≤ C_j p^(1+β_j) k_j^r_j,

only for B_(h_j) sets C of size k_j=⌊p^(1/(r_j+1))⌋, provided
β_j→0 and h_j/r_j→0. The relation orders can remain at h_j=2
or grow, for example like √r_j. **This restricted moment estimate,
called SS-B, is unproved.** The transfer does not infer it merely
from minimal additive energies.

The cost of extraction is retained. At fixed r,h, the resulting
range is ε>max(1/(r+1)+β,(2h−1)/(r+1)). In particular a
fourth-moment estimate only for Sidon sets does not give the earlier
ε>1/3 result through this method. The earlier forced-row obstruction
at size p^(1/r) does not refute the smaller SS-B slice.

The simplest mixed Weil product was also checked. Its four ordinary
and exceptional cases reproduce the pass-18 Hilbert–Schmidt identity;
no additional character-sum saving follows from that product.

The [exact verifier](../experiments/parallel19_relation_free_2026_09_05.py)
and [results](../results/parallel19_relation_free_2026_09_05.json)
cover 2,406 complete forbidden sets and 38,006 candidate extensions;
172 partitions with 33,276 greedy extensions; 26,298 minimal-energy
checks; 74 two-sided transfers and 296 moment transfers at orders
2,4,6,8; 12 arbitrary bounded kernels; and 112,707 independently
multiplied matrix quotients. Interval-size restrictions and rational
exponent choices are also checked. The large moment orders in the
exponent ledger are arithmetic parameter checks, not evaluations of
the unproved high-moment hypothesis. The [final audit](../results/parallel19_pass_audit_2026_09_05.json)
pins these files and preserves pass 18.

No sharp extraction theorem is imported. Shkredov's primary arXiv v1
HTML was reviewed for Sidon terminology and context and archived as
a background source. No new PDF was rendered or visually reviewed.
Independent mathematical review of the new deductions is outstanding.

This gives a more focused sufficient input for the full classical
conjecture. It does not establish a new arbitrary-set exponent, the
subgroup square-root target, a stronger clique bound, or the official
Reed–Solomon prize bridge. The full goal remains active and unproved.
Root worked locally; the three workers were stopped at account limits
in the most recent live check, during pass 18.
