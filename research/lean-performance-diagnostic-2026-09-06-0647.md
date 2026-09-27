# Lean diagnostic at 02:49 EDT, September 6, 2026

The current `queryphase-preservation-v19.lean` check is active, with no
evidence of a deadlock. It has no defensible completion ETA. This is a
point-in-time diagnosis, not a recurring monitor.

At 06:49:10 UTC, PID 13687 had run for 1h32m42s and accumulated
31m00.67s of CPU time. Resident memory was 23.1 GiB, down from 25.0 GiB
15 seconds earlier. CPU time increased by 1.65 seconds over that interval.
The machine has 16 logical CPUs and 128 GiB of RAM; its one-minute load
average was 93.33. Competing work is substantial. These observations do
not establish an exact uncontended runtime or slowdown factor.

A native stack sample at 06:47:38 UTC shows one worker executing Lean
kernel definitional-equality checks and reductions; the other sampled
workers wait. CPU activity establishes execution, not progress toward a
successful proof. The log has not changed since 05:17:10 UTC, and no exit
marker exists. Its final warning at line 966 is not a completion percentage.

The previous v18 log explicitly reports a kernel deterministic timeout
for `query_phase_step_preserves_fold`, whose declaration begins at line
502. The current source hash is
`138b2b0eb454e2b0010b59fe92462ef41f21496a93d68c624a493182cf740abb`.
It retains the large preservation proof and adds two `as_aux_lemma`
boundaries. There is no completed v19 benchmark from which to estimate
remaining work. Its source sets `maxHeartbeats` to 4,000,000; Lean's
installed option documentation defines heartbeats in terms of allocations,
not elapsed seconds.

The practical optimization is to shorten the proof term checked by the
kernel. The already imported
`logical_queryFiberPoints_eq_fiberEvaluations` and
`logical_computeFoldedValue_eq_iterated_fold` lemmas encapsulate conversions
that this proof repeats. A small successful-output lemma can identify the
returned folded value first, then apply these bridges. The incoming guard
does not need to be reproved safe merely to identify an already-successful
output. This refactor still needs Lean validation and timing; no speedup is
claimed. See the [source review](lean-queryphase-performance-review-2026-09-06.md)
for exact extraction points and correctness obligations.

A bounded independent reread identified a smaller first experiment: keep
the zero-branch support extraction through v19 line 855, then replace
lines 856–969 with a function-extensionality identification of
`fiber_vec.get` with `logical_queryFiberPoints`, followed by the existing
logical folded-value bridge. The remaining obligation is the zero-block
specialization and the identification of `getFirstOracle`. This would
remove two `congr! 8` blocks and repeated matrix/fiber expansion. The final
dependent-index normalization has not been compiled, so this is an exact
editing target rather than a validated replacement.

For the next deliberate attempt, use the warm dependency checkout and a
focused file check with `-Dprofiler=true -Dprofiler.threshold=1000`, recording
CPU time, wall time, peak memory, source hash, and terminal exit status.
These options were checked in the installed Lean 4.30.0-rc2 source. The
[official profiling talk](https://lean-lang.org/talks/community-meeting-may-2024.pdf)
also documents the structured trace profiler. Profiling helps localize
cost; it does not itself reduce it. Launching additional builds does not
parallelize the single busy kernel worker observed here.

No running job, priority, proof source, or cache was changed, and no new
Lean build was launched. The diagnosis is saved in the
[process observations](../results/lean_performance_diagnostic_2026_09_06_0647.json)
and [stack sample](../results/lean_kernel_diagnostic_2026_09_06_0647.txt).
The earlier [refresh](lean-performance-refresh-2026-09-06.md) remains
historical evidence. This diagnostic is not evidence of a Paley or
Proximity Prize proof.
