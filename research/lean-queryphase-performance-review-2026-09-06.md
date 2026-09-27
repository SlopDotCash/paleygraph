# QueryPhase kernel-cost source review, September 6, 2026

This is an independent, read-only review of the v18/v19 proof sources and logs. It does not supply a live process measurement or a reliable ETA. No Lean process, existing source, cache, build setting, or generated output was changed, and no build was started.

## Evidence and limits

- [`/tmp/queryphase-preservation-v18.log`](/tmp/queryphase-preservation-v18.log:69) explicitly reports `(kernel) deterministic timeout` for `query_phase_step_preserves_fold`, declared at source line 502. This establishes that the v18 failure is a kernel-checking failure, rather than merely a slow tactic search.
- The reviewed [`v19`](/tmp/queryphase-preservation-v19.lean:502) has the same main statement and essentially the same proof. Its only changes from v18 wrap the local fiber-evaluation equalities in `as_aux_lemma`, at lines 682 and 889. The declaration spans lines 502-969, with the positive branch at 578-800 and zero branch at 801-969.
- The observed v19 log contains warnings, including one at source line 966, but no terminal success or timeout. A warning at a late source line is not a progress percentage: elaboration and kernel checking can be asynchronous. The running process's current state belongs to the root lane's live diagnosis.
- Source SHA256: v18 `6fad1a413df3a6385ee75dbb4e71dd939cf00d41a7c6f941d2fa51c8ba184885`; v19 `138b2b0eb454e2b0010b59fe92462ef41f21496a93d68c624a493182cf740abb`.

The relevant checkout guide, [`/private/tmp/proximity-duplex-doc/AGENTS.md`](/private/tmp/proximity-duplex-doc/AGENTS.md), requires serialized full builds/cache operations through `scripts/lake-locked.sh`. Its `lean-toolchain` selects `leanprover/lean4:v4.30.0-rc2`. This review inspected that installed version's source rather than inferring behavior from a different installed Lean release.

## Concrete expensive patterns in v19

The strongest source-level suspect is repeated conversion between propositionally equal dependent indices and functions inside one large declaration. This is a hypothesis consistent with the v18 kernel failure and the prior kernel-stack sample, not a measured attribution to one tactic.

| Source location | Pattern | Separately checkable replacement |
|---|---|---|
| [v19:510](/tmp/queryphase-preservation-v19.lean:510), statement through 558 | Large invariant and target types repeat `iterated_fold`, `getFoldingChallenges`, and suffix terms with freshly constructed dependent indices and arithmetic proofs. | Preserve the public statement, but give internal helpers concise statements using a canonical block source/destination. |
| [v19:577](/tmp/queryphase-preservation-v19.lean:577), branches through 969 | Both branches unpack the simulated query and rebuild its conversion to a fold. | First split the two branches into named theorems with explicit parameters; then replace common work with the logical API described below. |
| [v19:664](/tmp/queryphase-preservation-v19.lean:664), through 695 | `f_mid` transports the oracle function with a type equality, then unfolds `fiberEvaluations`, `getFiberPoint`, and suffix definitions to compare it to `fiber_vec.get`. | Reuse the existing logical query/fiber bridge rather than reopening these definitions in the preservation theorem. |
| [v19:728](/tmp/queryphase-preservation-v19.lean:728), through 800 | Fold transitivity is followed by independent rewrites of step count, destination, function, challenge vector, and suffix. Each has dependent arguments. | A standalone prefix-plus-block fold lemma with one canonical index representation; extract the challenge-append equality at 774-786 into a small function theorem. |
| [v19:888](/tmp/queryphase-preservation-v19.lean:888), through 927 | The zero branch expands casts and runs `congr! 8` twice across dependent subtypes, with `all_goals` and more casts. | Reuse the logical fiber bridge. If a residue remains, use typed congruence lemmas at the subtype-value level; the local `cast_val` and `suffix_val` proofs already identify the small reusable facts. |
| [v19:935](/tmp/queryphase-preservation-v19.lean:935), through 964 | The zero branch again reconciles `ϑ`, `(k_zero+1)*ϑ`, `k_zero*ϑ+ϑ`, and their suffix types. | Normalize these indices once at the helper boundary, rather than carrying them as separate expressions through matrix/fold rewrites. |

The linter's unused simp arguments are low-priority cleanup. Removing those warnings would not explain or demonstrably fix the kernel timeout. Likewise, raising the heartbeat limit only permits more work; it is not a source-level reduction in work.

## Existing imported bridge can replace repeated proof work

The strongest concrete reuse opportunity is already in the import closure. [`Soundness.lean:8`](/private/tmp/proximity-duplex-doc/ArkLib/ProofSystem/Binius/BinaryBasefold/Soundness.lean:8) imports the following modules, and v19 imports `Soundness` at line 7:

1. [`QueryPhasePrelims.lean:237`](/private/tmp/proximity-duplex-doc/ArkLib/ProofSystem/Binius/BinaryBasefold/Soundness/QueryPhasePrelims.lean:237) defines `logical_queryFiberPoints` as the function of pointwise oracle answers. This is exactly the conclusion of [`mem_support_queryFiberPoints`, v19:278](/tmp/queryphase-preservation-v19.lean:278), after function extensionality.
2. [`QueryPhaseFoldBridge.lean:45`](/private/tmp/proximity-duplex-doc/ArkLib/ProofSystem/Binius/BinaryBasefold/Soundness/QueryPhaseFoldBridge.lean:45) already proves `logical_queryFiberPoints_eq_fiberEvaluations` for arbitrary block index `k` using canonical source/destination expressions.
3. [`QueryPhaseFoldedValue.lean:44`](/private/tmp/proximity-duplex-doc/ArkLib/ProofSystem/Binius/BinaryBasefold/Soundness/QueryPhaseFoldedValue.lean:44) already proves `logical_computeFoldedValue_eq_iterated_fold`, again for arbitrary `k`. This encapsulates the matrix-form conversion that v19 repeats in both branches.

An internal preservation proof can therefore be organized as three small interfaces:

```text
successful simulated step output
  = logical_computeFoldedValue (... logical_queryFiberPoints ...)
  = iterated_fold of the current oracle over the block
  = iterated_fold of the first oracle over the full prefix.
```

The first equality needs a new, separately checked support/output lemma; the second has the existing theorem above; the third isolates the oracle consistency and prefix-transitivity reasoning from v19:706-800. This proposal preserves the mathematical argument while reducing how much dependent-index material is exposed in the outer proof. It has not been compiled or benchmarked in this review.

## Successful output already excludes a failing guard

There is also an avoidable logical dependency worth checking separately. [`checkSingleFoldingStep`, QueryPhasePrelims.lean:162](/private/tmp/proximity-duplex-doc/ArkLib/ProofSystem/Binius/BinaryBasefold/Soundness/QueryPhasePrelims.lean:162) queries the fiber, optionally executes a guard, and returns a folded value that does not depend on the old `c_cur`, except through that guard.

V19 already assumes `s'` belongs to the successful `OptionT` support, but at [619-635](/tmp/queryphase-preservation-v19.lean:619) it proves `h_guard_pass` afresh using the large incoming-fold invariant and `query_phase_consistency_guard_safe`. For a partial-correctness/output lemma, one can instead eliminate the successful guard witness from the support decomposition. If the guard fails there is no successful output; if it succeeds, the return expression is the same logical folded value.

This would let the output-identification helper omit `h_c_k_correct_of_k_pos` and `h_relIn` entirely. The outer theorem can retain its current signature initially, so callers do not need to change. This is a semantic refactor to validate in a small lemma, not a suggestion to remove the guard or weaken the verifier. It also matches the current docstring at v19:490, which says the output computation does not require the incoming fold value to be correct, although the current signature still includes that invariant.

Similarly, successful support of the whole bound computation forces its fiber-query result to be `some fiber_vec`. A focused helper using successful `OptionT` bind-support decomposition may therefore avoid the stronger zero-failure-probability detour at [v19:595](/tmp/queryphase-preservation-v19.lean:595) through 606 and [830](/tmp/queryphase-preservation-v19.lean:830) through 843. This would keep probability calculations out of a statement that only identifies an already-successful result. The exact simulator/`OptionT` rewrite sequence remains to be checked.

## What as_aux_lemma changes, and what it does not

Installed Lean's [`Init/Tactics.lean:14`](/Users/shawwalters/.elan/toolchains/leanprover--lean4---v4.30.0-rc2/src/lean/Init/Tactics.lean:14) documents the possible expression-size benefit: the proof term is not duplicated. Its implementation at [`Lean/Elab/Tactic/AsAuxLemma.lean:19`](/Users/shawwalters/.elan/toolchains/leanprover--lean4---v4.30.0-rc2/src/lean/Lean/Elab/Tactic/AsAuxLemma.lean:19) first runs the tactic and then calls `mkAuxTheorem`. [`Lean/Meta/Closure.lean:457`](/Users/shawwalters/.elan/toolchains/leanprover--lean4---v4.30.0-rc2/src/lean/Lean/Meta/Closure.lean:457) closes the term and returns a constant application. [`Lean/Meta/Tactic/AuxLemma.lean:58`](/Users/shawwalters/.elan/toolchains/leanprover--lean4---v4.30.0-rc2/src/lean/Lean/Meta/Tactic/AuxLemma.lean:58) creates a theorem declaration and calls `addDecl`.

Therefore v19's two wrappers can improve sharing, but they do not make either proof free: the auxiliary declaration itself still requires kernel checking. Nor do they isolate the remaining support, prefix-transitivity, challenge, and destination conversions in the outer theorem. Wrapping the entire existing proof in another auxiliary lemma would mainly move the work unless it also reduces duplication or exposes a smaller interface.

## Safest next measurement

The least invasive first experiment is to check named positive/zero branch lemmas separately, keeping their proofs unchanged, followed by the tiny case-split wrapper. This localizes the next terminal failure and gives component timings. The more promising reduction is then the already-imported logical bridge plus a successful-output helper.

Use the existing cached checkout and focused per-file checks, after the current attempt ends or as an explicitly coordinated replacement. Enable Lean's `profiler` and an appropriate `profiler.threshold`; the installed definitions are at [`Lean/Util/Profile.lean:16`](/Users/shawwalters/.elan/toolchains/leanprover--lean4---v4.30.0-rc2/src/lean/Lean/Util/Profile.lean:16). Record exact source hash, terminal exit status, wall/CPU time, and peak memory for each tested declaration. A failed or still-running attempt is not a speedup benchmark.

There is no defensible completion ETA from source size or the last warning. The heartbeat option counts small memory allocations in thousands, not seconds, as specified at [`Lean/CoreM.lean:28`](/Users/shawwalters/.elan/toolchains/leanprover--lean4---v4.30.0-rc2/src/lean/Lean/CoreM.lean:28). The current run's activity and resource pressure must be determined from fresh process samples; this note supplies concrete proof reductions to test, not a prediction that waiting a fixed number of minutes will succeed.
