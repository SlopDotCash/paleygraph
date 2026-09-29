# Twenty-sixth research assessment

**The Paley conjecture and the Proximity Prize remain unproved.**
This pass makes the code-agreement projection proof machine-checkable and
corrects its attribution. It improves verification, not the unresolved
cancellation estimate. No percentage complete or time-to-proof estimate
is justified.

## What changed

Seven principal statements in the
[standalone proof](../prove2me/Check_paley_mca_projection.lean) pass Lean:
one nonzero row projection preserves every originally bad scalar; the
input-failure predicates are equivalent; and uniform bounds on bad-scalar
counts are equivalent for rows and scalars. These results include the
actual existence-of-agreeing-codeword witnesses on coordinate sets, with
an arbitrary linear code, finite field, and finite nonempty row set.
The [result and environment record](../results/parallel26_mca_projection_2026_09_05.json)
pins the proof and prints only standard axioms for all seven statements.

During the source comparison, root found that the already archived
[ArkLib Errors.lean](https://github.com/Verified-zkEVM/ArkLib/blob/e65197892890b8fd9b0dc05b8980273cf1d595cc/ArkLib/Data/CodingTheory/ProximityGap/Errors.lean#L1091)
already proves `mcaError_interleaved_eq`. Its proof also preserves all bad
scalars in one projection by avoiding at most |F| proper subspaces.
The pass-25 MCA conclusion is therefore a rediscovery. Its earlier
assessment did not identify this exact existing theorem; that attribution
is corrected here and in the central project records.

The [detailed note](parallel26-mca-formalization-2026-09-05.md) gives the
formal theorem, proof, interpretation of official definitions, and source
comparison. The Mathlib proof was compiled; the whole ArkLib dependency
closure was not rebuilt. No claim of a kernel-checked ArkLib import or
independent human review is made.

The list-size projection argument from pass 25 remains a written proof
with an author audit and finite verification. No matching finite
collision-count transfer was found in the few inspected source files;
this is not a novelty claim. The attempt to improve the prize obstruction
by coefficient pigeonholing also recovered an existing project result:
[list-to-winning-set.md](list-to-winning-set.md) already proves the
threshold at radius 122641/262144 using r=139503. No new obstruction or
stronger threshold is counted in this pass.

## Prove2Me and parallel execution

The [project setup](../PROVE2ME.md) is complete and reuses the existing
private draft and supported Lean environment. Four earlier elementary
proofs have ACCEPTED server verdicts in the setup readback. The newer
private incidence-lemma statement job
843b8456-f9ad-4158-809d-8378028cc958 remains PENDING in the
[latest server record](../results/prove2me_projection_server_2026_09_05.json).
It has no theorem ID, proof submission ID, or proof verdict yet. The
larger standalone projection proof has only local verification. No
duplicate request or public release was made.

All three workers still report terminal usage-limit errors. None is
currently researching or reviewing. Root continued with formalization
and the source audit. No model switch, usage reset, credit purchase,
automation, or memory write was performed. The full goal remains active:
the worker and queue conditions have not stopped meaningful local work.

## What remains mathematically necessary

The scalar reduction retains the full official extension field. It gives
no new upper bound on either scalar maximum and does not reduce that
field to its prime subfield. The combination-round budget and separate
spot check are not proved. No equivalence between the full Paley
conjecture and the official prize has been established.

The pass-25 signed, subgroup, and spectral bounds retain their earlier
limits. The signed moment still needs a uniform upper bound for actual
character inputs. The subgroup argument leaves six distinct unbalanced
relations and an error above the cubic target scale. The full spectral
operator still needs the elliptic saving in every relevant sector and
control of its nonzero border. This pass supplies none of those estimates.

## Evidence and preservation

The previous turn is classified as progress: Prove2Me setup and local
algebraic verification. This pass is also progress in verification and
source accuracy, without a stronger asymptotic conclusion. The
[prior-state snapshot](../results/parallel26_prior_state_2026_09_05.json)
and [final audit](../results/parallel26_pass_audit_2026_09_05.json) distinguish
immutable earlier proof artifacts from current documentation and the
deliberately refreshed server queue record. Earlier pass-25 proof and
result bytes are preserved. Two exact-commit source files are added to
the manifest; its 34 earlier entries remain intact.

The next productive research target remains one of the actual uniform
estimates above. Reproving interleaving invariance or the existing
coefficient-pigeonhole threshold does not advance that target further.
