# Pass 39: a finite constant obstruction and two structural restrictions

The Paley conjecture, the full Proximity Prize, and the uniform subgroup
period estimate remain unproved. The current energy exponent49/20 and
ordinary period exponent71/72 are unchanged.

The [principal-shell construction](parallel39-shell-structural-input-2026-09-06.md)
proves that a cyclotomic element of norm kp, with p not dividing k,
in a split evaluation kernel forces an actual nonzero coset with
V<=k^2 and eta>=n-4 pi^2 k^2. Here n>=4 is a power of two.
The known norm-prime trinomial at p=215535361,n=128 yields
a=38468180 and V=1 exactly. Rational comparisons establish
eta(a)>2 sqrt(n log(p/n)), refuting the literal finite constant C=2.
It does not refute a bound with an unspecified absolute constant or
an asymptotic bound allowing finite exceptions.

The construction has a useful converse restriction: every valid
M<=U(n,p) forces k^2>=max(0,(n-U)/(4 pi^2)). In any range with
U=o(n), including the quartic window under the existing scoped inputs,
k>=(1-o(1))sqrt(n)/(2 pi). Near-prime norm templates therefore cannot
scale in that window. No new upper bound follows from this restriction.
The [independent review](parallel39-independent-principal-shell-review-2026-09-06.md)
checks the norm-prime theorem and finite example with a separate
64-by-64 determinant computation. Root reran the
[exact certificate verifier](../experiments/parallel39_principal_shell.py).

The [edge-coset reduction](parallel39-edge-coset-reduction-2026-09-06.md)
bounds the coincident-edge part of the weighted triangle sum by
3 F3* F4*/n, hence O(n^(17/3) log n) under the existing inputs.
Three distinct edge cosets remain uncontrolled, and distinctness
alone does not repair the old incidence argument. The universal
bound received an [independent review](parallel39-independent-coset-coincidence-review-2026-09-06.md)
and direct finite checks.

Independent reviews also validate the
[single-degree moment comparison](parallel39-independent-single-degree-review-2026-09-06.md)
and [shell inversion sections1-4](parallel39-independent-shell-review-2026-09-06.md).
The latter required making the intended domain A=H,a!=0 explicit;
that clarification is now in the source. These are reviews by separate
agents, not external peer review or Lean certification. The positive
moment bound and uniform thin-annulus estimate remain open.

During this pass the user asked about a long Lean run. A fresh native
process sample confirms active kernel checking, together with severe
machine contention. The [timestamped diagnosis](lean-performance-refresh-2026-09-06.md)
separates that tooling issue from proof-search progress. No new Lean
build or Prove2Me submission was started for these mathematical results.

Next mathematical work should address the positive moment bound, the
distinct-edge weighted sum with its actual level-set structure, or
uniform shell deviations subject to the new arithmetic restriction.
The finite C=2 obstruction supplies no asymptotic disproof.
