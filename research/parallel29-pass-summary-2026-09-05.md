# Twenty-ninth research assessment

**The full Paley conjecture and the Proximity Prize remain unproved.**
This pass proves an exact-centered subgroup energy recurrence in ordinary
mathematics and determines the limit of several ways of iterating it.
The subgroup cancellation exponent recorded in pass 28 does not improve.

For an odd prime p, a multiplicative subgroup H of order n, and its
ordered additive convolution energies, put T_s=E_s-n^(2s)/p. The
[proof](parallel29-centered-recurrence-2026-09-05.md) gives

    T_(2s) <= D n^(2s-1/2) T_s

with an absolute D independent of s. A point--plane encoding and an
elementary incidence variance identity remove the uniform contribution
exactly. Layer decomposition extends the estimate to weights; positive
and negative parts are handled separately, and the origin is retained.
The argument does not assume the origin-free signed-function statement
previously refuted by the project. Its only non-elementary new input is
the already archived published point--plane theorem.

This removes the principal-term cap from the previous recurrence. However,
its deficit increases by only one half each time the order doubles.
The best power from the centered mixed-moment gate remains 71/72 in the
quartic window. Feeding that amplitude bound back into higher moments and
using exact products of coset moments also fails to improve it: an explicit
piecewise linear envelope is closed under all these operations. This is
a limit of those particular upper-bound rules, not an obstruction to a
stronger theorem using additional information.

The proof also states why the prime-field invariant-set estimate cannot
be transferred automatically to extension fields: a prime subfield
supplies an unbounded ratio against the proposed extension-field bound.
This familiar obstruction prevents a direct application to the official
prize parameters. It does not refute the prize or either Paley target.

The [verification record](../results/parallel29_verification_2026_09_05.json)
separates exact incidence, variance, centering, origin and exponent checks
from approximate Fourier checks. The [source ledger](../results/parallel29_source_scope_2026_09_05.json)
reuses pinned primary sources. These are bounded checks of the bookkeeping,
not a formal proof of the analytic estimate. Separate-author review and
Lean verification of the new argument remain outstanding.

The previous turn is progress: the live performance investigation established
that the single-file check had finished, identified competing proof and
documentation builds, and located active kernel checking under heavy CPU
contention. This pass added no Lean job. The small standard-library verifier
runs in one process. Existing parallel agents were inspected and still
report terminal usage-limit errors; root completed this pass locally.
The earlier accepted Prove2Me theorem is terminal and was not resubmitted.

The next useful mathematical step is an input that exceeds the explicit
closed energy envelope, or direct control of the signed spectral quantities.
More iteration of the listed estimates is not sufficient. Growing-order
subgroup remainders, arbitrary-input classical correlations, the full
spectral operator, and scalar prize bounds are still unproved.
