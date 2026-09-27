# Pass 42 review and validation scope

The preceding completed pass made mathematical progress: it improved W
and quantified a sufficient collision-excess input. This pass adds a
weighted prime-exception budget and checks two specific obstructions.
The original goal is active and unproved.

## Review provenance

The incidence lane independently checked the compatible tuple count,
nonvanishing, and valuation bookkeeping for the coarse norm budget. The
source lane found and archived the primary circularity and cyclotomic
number sources. The classical lane returned the reflection projection
identities. All three reached the account usage limit before completing
their planned artifact work. Root wrote the final notes, checked the
source applications, and derived the sharper AGM constant and the
reflection recovery comparison. The latter steps have root review and
finite checks; no completed separate-agent review is claimed for them.

## Norm argument

Root checked that equality of two nonzero products of shifted complex
units forces both the same product and the same sum of the unshifted
units. It therefore forces the unordered pairs to coincide. This removes
exactly 2(n-1)^2-(n-1) trivial ordered tuples. The remaining set is Galois
invariant, so its product is a nonzero rational integer P_n.

At a split prime, every actual collision vanishes in a fixed primitive
root embedding. Each contributes at least one p-adic valuation, giving
X_p<=v_p(P_n). The alternate all-norm computation multiplies both the
total valuation and the counted roots by phi(n); division occurs once.
A norm divisor for one tuple implies existence of a vanishing primitive
root, not vanishing at an arbitrarily specified root for that same tuple.
The full Galois-invariant tuple product avoids that distinction.

Expanding the squared absolute values gives
S=8n^2(n-1)^2-2n^4. Trivial tuples contribute zero, so AGM applies to
the M positive remaining terms. This yields |P_n|^2<=(S/M)^M and the
weighted budget. Markov counting with threshold n^(61/31) and
log p bounded below by log(c n^4) gives the exponent 63/31. The pass41
pure combinatorial implication is used only at the nonexceptional primes.
No denominator estimate or conclusion about all primes is introduced.

The [norm checker](../experiments/parallel42_norm_budget.py) passed for
n=4,8,16. It compares recursive norms with 446 separate Bareiss
determinants, factors every tuple norm with exact multiplicities, checks
the exact second-moment sum and integer AGM inequality, and verifies all
split exceptional primes with a direct ratio computation. The 101 extra
split-prime checks are finite consistency checks. Completeness for the
three fixed orders follows from full norm factorization. The
[result file](../results/parallel42_norm_budget_2026_09_06.json) also
records the twelve rich cells at p=353, the zero diagonal histogram, and
the fact that this witness is outside the quartic window.

## Reflection calculation

Root independently derived the one-step sixth-moment recurrence,
including the nonzero M3^2 term. Self-adjoint idempotence applies because
the same center is used on both functions; it gives expected pairing
S/2^r. Zero mean of the character convolution is retained after every
projection. Interpolation uses the nonempty test set's mass m.

Every term in the sixth-moment recurrence is nonnegative, giving
B_r>=M6(F)/2^(5r). Substituting this lower bound in the explicit upper
certificate proves that this certificate is no better than direct
sixth-moment Holder. It is not a lower bound for the actual discrepancy.
The O(pn^3) smoothed-moment upper bound requires n^4<=p and imports
only the previously checked individual Weil estimate. The exact
projection and recovery identities themselves do not need n^4<=p.

The [reflection checker](../experiments/parallel42_classical_affine.py)
passed exact rational comparisons on 2,266 center tuples in six cases,
including cases with nonzero third moments. Its
[result file](../results/parallel42_classical_affine_2026_09_06.json)
distinguishes these identities from an original Paley bound.

## Remaining limits

The primary-source findings and version distinctions are in the
[source audit](parallel42-excess-source-audit-2026-09-06.md). No claim
is made that these results are novel, externally peer reviewed, or
formalized in Lean. Existing pass41 artifacts are preserved. The
weighted-triangle bound for every prime remains the pass41 bound,
and energy49/20, period71/72, full Paley, and official prize conclusions
remain unchanged. No process was stopped, reprioritized, or launched for
Lean; no Prove2Me mission or submission was started in this pass.
