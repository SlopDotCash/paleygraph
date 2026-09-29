# Lean performance diagnosis, September 6, 2026

Observed at approximately 01:20–01:23 EDT. This is a timestamped diagnosis, not a live monitor.

The original `/private/tmp/append-adapters-final.log` is terminal: `QueryPhase` took 5,478 seconds (91.3 minutes) and failed. `QueryPhaseFirstOracle` took 5,042 seconds (84 minutes). The failed build also lists `Steps.Relay` and `Steps.Commit`.

The newer process, PID 13687, checks `/tmp/queryphase-preservation-v19.lean` from `/private/tmp/proximity-duplex-doc`. At the saved snapshot it had run for 5 minutes 58 seconds, accumulated 6 minutes 5 seconds of CPU time, used approximately one full CPU core, and had 11.3 GiB resident memory. Earlier samples showed rising CPU time and memory. These observations establish activity, not useful mathematical progress or eventual termination.

A three-second native stack sample places an active worker in Lean's kernel type checker, particularly definitional equality (`is_def_eq_core`, `quick_is_def_eq`) and definition reduction. Other worker threads are waiting. The process is therefore CPU-bound during this sample, not waiting for a build lock. The machine has 16 logical CPUs and 128 GiB RAM. Swap-out counters did not increase across the observations; old occupied swap does not establish current memory thrashing.

The v18 repair log reports `(kernel) deterministic timeout` at line 502, the declaration of `query_phase_step_preserves_fold`. The current v19 file adds `as_aux_lemma` around two local proofs, at lines 682 and 889, but its main preservation lemma still spans almost half of the 977-line file. The same declaration is a strong candidate for the current bottleneck; the native sample itself does not expose its declaration name. The latest warning location is not a progress percentage.

There is no reliable remaining-time estimate. The current proof differs from the original 91-minute failed run. v18's source modification to last log write was about 29 minutes, but that is not an exact runtime benchmark and it ended with a timeout. No measured speedup for v19 is available yet.

Recommended next experiment: split preservation into separately checked lemmas for the zero and positive index branches; reduce dependent index casts through small, explicitly typed congruence lemmas; then benchmark only this dependency slice using the existing cache. Lean documents auxiliary lemmas as a way to prevent proof-term duplication, but this is a hypothesis to benchmark, not a guaranteed fix. Increasing parallel build counts will not directly accelerate the one busy kernel worker observed here.

For the next intentional run, enable component timings without doing a full build:

```sh
cd /private/tmp/proximity-duplex-doc
/usr/bin/time -l lake env lean -DautoImplicit=false -Dprofiler=true -Dprofiler.threshold=1000 /tmp/queryphase-preservation-v19.lean
```

Run that only after the current attempt finishes or as part of an intentional replacement experiment. It was not started during this diagnosis. The profiler options were checked against the installed Lean 4.30.0-rc2 source. Timings are reported by completed profiling scopes and are not a countdown. Full builds and cache operations must use the project's `scripts/lake-locked.sh` wrapper.

Evidence: [machine and process snapshot](../results/lean_performance_live_2026_09_06.json), [native stack sample](../results/lean_kernel_sample_2026_09_06.txt), [original build log](/private/tmp/append-adapters-final.log), [previous repair log](/tmp/queryphase-preservation-v18.log).

Primary references: [Lean auxiliary lemma tactic](https://lean-lang.org/doc/reference/latest/Tactic-Proofs/Tactic-Reference/), [Lean profiling options](https://lean-lang.org/doc/api/Lean/Util/Profile.html), [Mathlib guidance on large definitional-equality checks](https://leanprover-community.github.io/contribute/style.html).

Actions taken: inspected processes and logs, sampled the active process, and saved this report. No proof source, cache, running job, scheduler, or build settings were changed. This diagnostic provides no evidence that the Paley conjecture or Proximity Prize has been proved.
