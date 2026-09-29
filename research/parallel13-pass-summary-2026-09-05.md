# Thirteenth pass: a larger logarithmic range for the full aggregate

**Progress: a source-dependent improvement of the complete spectral bound.**
The argument is checked locally and still needs independent mathematical
review. It does not yet improve the Paley clique bound or solve the full
conjecture or prize objective.

The [written proof](parallel13-principal-budget-2026-09-05.md) separates
principal constituents from their lower-weight boundary corrections.
The corrections have zero conductor at the moving endpoint y, whereas
every principal constituent has nontrivial inertia there. This excludes
principal/error tensor invariants and gives an extra factor of p^(-1/2)
in their averaged cross correlation.

Summing the actual principal local types gives an exact recurrence.
For a anchors, put h=2^(a−1), N=2^a−1. The sum of principal ranks
has exponential growth rate equal to the larger root λ_a of

\[
 \lambda^2-h(a+a^{-1})\lambda+N=0.
\]

A positive weighted count proves this rate at arbitrary depth, with
explicit upper and lower constants. It is strictly below the raw-rank
rate h(a+1)−1 used in the previous estimate. No boundary term has
been discarded to obtain the improvement.

Let Ψ_a=λ_a²/N and Θ_a=[h(a+1)−1]²/N. For fixed a≥2,
the complete normalized matrix energy and even trace are at most
(1+o(1))p throughout

\[
 j\le\min\left\{
   \frac{(1/2-\epsilon)\log p}{\log\Psi_a},
   \frac{(1-\epsilon)\log p}{\log\Theta_a}\right\}.
\]

For a=2 the controlling growth rate falls from 25/3≈8.33333 to
(19+5√13)/6≈6.17129. The leading coefficient of log p in the
available depth increases from approximately 0.23582 to 0.27474,
an increase of about 16.5%. This remains a bounded multiple of
log p; the needed j/log p→∞ is outside the proved range.

The [exact verifier](../experiments/parallel13_principal_budget_2026_09_05.py)
checks 8,516 enumerated words, 5,053 complete error inventories at y,
729 weighted recurrences and 1,458 coefficient identities in exact
quadratic fields. The new bound holds and is smaller in all 168 retained
finite energy cases. A p=10009 example is nonvacuous. These checks
support the written derivation; they do not themselves verify sheaf
weights or prove the asymptotic theorem. See
[results](../results/parallel13_principal_budget_2026_09_05.json) and
[artifact audit](../results/parallel13_pass_audit_2026_09_05.json).

The earlier projection obstruction remains valid. The next task is to
control the remaining principal Frobenius correlations at much longer
depths, or find a different argument for the extreme eigenvalues.
The subgroup exponent is unchanged at 8/9 on its existing density-one
quartic prime class; uniform subgroup cancellation, exceptional primes,
the arbitrary-two-set conjecture and the exact official prize bridge
remain open. The three workers remain unavailable at their previously
observed account limits; this pass was completed locally.
