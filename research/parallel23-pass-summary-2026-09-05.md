# Twenty-third pass: actual partial upper bounds and spectral coupling

**The full Paley conjecture and official prize implication remain unproved.**
This pass establishes an almost-all classical upper bound, bounds a specified
part of the subgroup relations, and evaluates the spectral coupling left open
in pass 22. It does not improve a worst-case cancellation exponent. Actual
character estimates now control some previously unestimated quantities, with
their exceptional cases explicit.

The user-authorized parallel work used three agents plus root. All four
research lanes and four separate-agent reviews are complete. Each lane has
a proof note and exact finite verifier. The combined audit below reconciles
the reviewed bytes, source inputs, corrections and verification scope.

## Classical: an upper bound for almost all Sidon sets

For each large odd prime, let n=⌊p^(1/4)⌋ and choose a Sidon n-set C uniformly.
The [classical proof](parallel23-classical-upper-2026-09-05.md) shows

    P[M_6(C)>15pn³ | C Sidon] ≤ (216/5+o(1))n²/p.

Every good set also satisfies the actual quartic-star estimate
L_out(C)≤(2/3)p²n³. This holds at each sufficiently large prime, without
averaging over primes. The proof uses an explicit variance bound, the full
character transform, a proved subset-slice variance inequality, and a Sidon
collision count.

This extends the local typical-set fourth-moment result to the sixth moment
and the requested Sidon slice. The pass-20 transfer needs the bound for
**every** input of the prescribed size. A vanishing exceptional proportion
can still contain an entire adversarially chosen structured family. No
worst-case Paley exponent follows.

The [verifier](../experiments/parallel23_classical_upper_2026_09_05.py) and
[results](../results/parallel23_classical_upper_2026_09_05.json) check 3,551
exhaustive small-slice sets, 15,042 insertion identities, 1,780 full-transform
identities, 38 slice-variance checks and 19 variance certificates. Larger
actual-field checks cover 287 sampled Sidon sets and five complete
quartic-star energy fixtures. Five further large-prime records evaluate
analytic probability certificates only. The proof establishes the almost-all
quantifier; the sampling success rate does not establish a universal bound.

## Subgroup: multiplicative closure controls balanced products

The [subgroup proof](parallel23-subgroup-upper-2026-09-05.md) separates
six-term relations according to a product constraint. A bijection turns
equal-sum, equal-product ordered triple pairs in H into shifted
multiplicative-energy solutions with an additional factor n=|H|. An existing
Shkredov theorem then bounds the opposite-free relations admitting a balanced
three-plus-three partition by O(n³(1+log n)) in the quartic window.

The unbalanced remainder is not controlled at that order. At the actual
subgroup p=1073748737, n=256, the relevant cyclotomic triple-collision count
is zero, but 368640 opposite-free six-term relations remain. The verifier
provides an explicit unbalanced witness. This shows that eliminating the
balanced portion does not eliminate the missing term. It is not a
counterexample to a Gaussian-order upper bound.

The [verifier](../experiments/parallel23_subgroup_upper_2026_09_05.py) and
[results](../results/parallel23_subgroup_upper_2026_09_05.json) cover nine
actual subgroups, 121,355 normalized weighted six-term entries and 560
product-ratio identities. Small groups also receive direct triple-bucket and
cyclotomic-matrix checks. Explicit bounds for the circular subgroup class
remain at the fourth-power scale. No improved total cancellation exponent
is claimed.

## Spectral: the controlled direction does not decouple

Let B be the compressed residue-frequency projection on an actual two-anchor
prime-field neighborhood, and u its unit uniform vector. Pass 22 showed
u^TBu→3/8. The [spectral calculation](parallel23-spectral-coupling-2026-09-05.md)
now establishes

    ||(I−uu^T)Bu||→1/8,       u^TB²u→5/32.

A known square-parameter Legendre-family second moment transfers to the
neighborhood through exact parameter symmetries and the complete character
transform. This evaluates the restricted trace variance left open previously.

Deleting the cross terms between u and its orthogonal complement has
operator error tending to 1/8. Thus u is not an approximate eigenvector.
This does not rule out an argument retaining the coupling or controlling
a larger subspace. No all-vector spectral upper bound or longer
complete-power range follows.

The [verifier](../experiments/parallel23_spectral_coupling_2026_09_05.py) and
[results](../results/parallel23_spectral_coupling_2026_09_05.json) check 79
actual prime fields, 36,465 parameter-symmetry and square-fibre cases, 9,057
neighborhood rows and 2,369 deleted-block entries. The asymptotic input is
Grove's published theorem. The verifier does not independently compute
Hecke traces or prove that source theorem.

## Prize: a restricted bridge and the exact missing counts

The [prize note](parallel23-prize-bridge-2026-09-05.md) identifies lists for a
monic received polynomial of degree k+a with fibres of the first a
elementary symmetric sums of its agreement roots. A multivariate character
bound estimates this restricted family. Linear subgroup periods cover a=1.

A pigeonhole argument shows that a uniform square-root bound for all
polynomial phases of degree proportional to the domain size cannot provide
that extension. At the pinned official profile, a nonzero phase of degree
at most 26215 exists with real subgroup sum above 180224 on the 262144-point
domain. This is an existence theorem, not a computed production-size
polynomial or a large-list assertion.

For arbitrary official words, the note states the exact missing certificate
as a sum of two maxima: distinct remainder explanations for 16 rows and
nontrivial remainder-ratio scalars for 8 rows. Their combined integer budget
is 274980728111395087. Neither maximum is bounded here. The source is the
archived official contract at commit b34c0131cfa36b51111521541d7d3e35c8791082.
No newer-contract or current-leaderboard claim is made.

The [verifier](../experiments/parallel23_prize_bridge_2026_09_05.py) and
[results](../results/parallel23_prize_bridge_2026_09_05.json) check 841
Newton/subset identities, ten monic-list classifications, 182
multivariate-count bounds, two extension-field centre lists, 24 general
list equalities and 72 MCA ratio equalities. Production-size MCA and list
maxima are not evaluated.

## Review and remaining work

Root checked the four arguments and their scope, including source
normalizations, classical variance constants, subgroup zero-product
corrections, and the exact-pinned definitions of Lambda and IsMCA.
The separate-agent reviews are complete for the
[classical](parallel23-classical-independent-review-2026-09-05.md),
[subgroup](parallel23-subgroup-independent-review-2026-09-05.md),
[spectral](parallel23-spectral-independent-review-2026-09-05.md), and
[prize](parallel23-prize-independent-review-2026-09-05.md) notes.
No mathematical correction was required. One report-only issue was fixed:
the classical verifier now labels the conditional Sidon certificate as
inapplicable outside its theorem range, including empty Sidon test slices.
Root reran that verifier; its reviewer reconstructed the old hashes after
reversing just the reporting changes, confirming unchanged mathematical
results and in-range certificates.

The prize reviewer separately checked 576 small cases of witness shrinking;
the subgroup reviewer separately recomputed the circular collision count,
checked its witness, and confirmed the strict opposite-free subset example.
The spectral reviewer checked the published normalization and optional
Hecke formula; it did not independently compute Hecke coefficients.
These are agent reviews, not human refereeing or formal verification.

The [final audit](../results/parallel23_pass_audit_2026_09_05.json) checks
23 recorded input hashes and 32 review input hashes. It preserves 24
previous proof/verification artifacts and all 27 earlier manifest entries.
The manifest now has 32 entries: Grove's primary HTML and four additional
ArkLib definition files at the already pinned contract dependency commit
were added. No PDF or newer official contract was substituted. The latest
live agent listing has all three children completed, with none running.

The next classical estimate must control the exceptional inputs or prove
a valid way to avoid them for every required input. The subgroup lane must
control the unbalanced relations, including exceptional primes. The
spectral lane must handle the remaining operator while retaining its
nonzero coupling. The prize lane needs bounds for arbitrary-word remainder
fibres and ratios with the correct field and radius quantifiers.

No estimate in this pass closes the full classical conjecture, uniform
subgroup square-root target, spectral edge, or official Reed–Solomon
certificate. The original goal remains active and unachieved.
