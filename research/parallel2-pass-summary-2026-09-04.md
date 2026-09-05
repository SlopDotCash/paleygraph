# Second parallel pass: growing families and the aggregate gap

**Paley remains unproved.** This pass extends two three-label necklace
families to arbitrary length, derives an exact aggregate-to-clique
criterion, and rules out a specific classical fourth-moment approach.
The missing uniform cancellation estimates remain missing. No novelty,
new Lean certificate, prize submission, or full reduction to the prize
is claimed.

Three agents continued the subgroup, necklace, and classical moment
lanes. The primary agent derived the spectral transfer, reviewed their
proofs, and checked the analytic source used by the necklace lane.

| Lane | Proved or refuted in this pass | Remaining gap |
|---|---|---|
| [Growing necklaces](parallel2-necklace-2026-09-04.md) | For A={0}, B={1}, C={0,1}, the families A^rBC and A^rBAC have bounds of order p^((k+1)/2), uniformly in r. The large exceptional quadratic mode is removed exactly. | Arbitrary separating gaps and cancellation across the weighted sum of all words. |
| [Spectral transfer](parallel2-spectral-transfer-2026-09-04.md) | An exact Chebyshev trace formula retains the rank and anchor terms. It gives a quantitative condition implying clique bounds at growing localization degree, without a separate degree or localization-size estimate. | The required aggregate bound is unproved. Individual necklace bounds do not establish it, and a clique bound alone would not prove the full two-set conjecture. |
| [Subgroup mixed energy](parallel2-subgroup-2026-09-04.md) | A primitive-root polynomial encodes balanced mixed energy exactly. Its squarefree defect gives upper and lower bounds; no triple root implies B<=2k^2. A new order-512 quartic collision is certified. | The no-triple-root condition is itself a collision hypothesis, certified only in finite quartic cases. Logarithmic-depth centered moments remain uncontrolled. |
| [Classical fourth moments](parallel2-classical-2026-09-04.md) | An exact signed cross-ratio formula includes the random-set mean. Intervals in an unbounded prime family disprove a proposed sufficient L2 bound, even after this centering. | Directional cancellation against the elliptic trace vector. The failed L2 bound does not refute the logarithmic-moment hypothesis. |

The new necklace input is [Katz's Mellin-transform theorem](https://web.math.princeton.edu/~nmk/mellin186.pdf),
Corollary 4.2 and Theorems 15.1 and 15.3. Its hypotheses, local monodromy,
normalizations, and exceptional characters were checked. The primary PDF
is [archived locally](../sources/katz-finite-field-mellin.pdf); printed
pages 22, 52, and 53 were rendered and inspected. This is an explicit
dependency on an existing deep theorem, not a new proof of that theorem.

The spectral calculation makes the main distinction quantitative. Bounds
on separate words can still yield an exponentially large bound after
summing over all labels. The exact aggregate criterion identifies the
cancellation that would be needed to get beyond the present results.

All four lanes include exact verification scripts and results linked in
their notes. Computation checks the formulas and finite examples; it
does not replace the missing uniform proof. The original thin-subgroup
target, arbitrary-set Paley conjecture, and general prize connection are
unchanged.
