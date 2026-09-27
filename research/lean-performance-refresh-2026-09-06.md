# Lean performance refresh, September 6, 2026

Observed at 01:56 EDT (05:56 UTC). This is a saved diagnosis, not a
live monitor. The long-running process is PID13687, checking
`/tmp/queryphase-preservation-v19.lean` from
`/private/tmp/proximity-duplex-doc`.

At 05:56:58 UTC its elapsed time was40m30s, CPU time16m14.73s,
and resident memory19.6GiB. Over the preceding26.35seconds, CPU time
advanced2.99seconds: roughly11% of one core during that short interval.
The machine has16 logical CPUs and128GiB RAM. Its one-minute load
average was51.25; multiple other compilation and computation jobs
were active. Contention is materially extending wall time. This sample
does not establish an exact slowdown factor for an otherwise identical
uncontended run.

A three-second native stack sample again places a worker in Lean's
kernel `type_checker`, including `is_def_eq_core` and definition
reduction, while the main thread and other workers wait. CPU time
and memory both advanced. The process is active; the evidence does
not support calling it deadlocked. Activity does not establish useful
proof progress or eventual completion.

The log is3780bytes and last changed at05:17:10UTC; no exit marker
exists at the observation time. Its last warning at line966 is not a
completion percentage: Lean may elaborate later text while earlier
proof terms still await kernel checking. The prior v18 attempt reports
a kernel deterministic timeout at the declaration starting at line502.
The v19 source adds two auxiliary-proof boundaries but has no measured
successful runtime yet. The original full QueryPhase build ended in
failure after5478seconds (91.3minutes); that is not a prediction for v19.

There is no defensible remaining-time estimate from these observations.
The proof checker provides no remaining-work counter here, and both
proof structure and available CPU differ between attempts. Occupied
swap was13.7GB, but the system swap-out counter did not change between
the two snapshots; that short interval does not show sustained swapping.

The most concrete performance opportunity is to reuse the already
imported `logical_queryFiberPoints_eq_fiberEvaluations` and
`logical_computeFoldedValue_eq_iterated_fold` lemmas, replacing the
duplicated vector/fiber and matrix/fold conversions in the preservation
proof. A separately checked helper should identify the returned folded
value from successful support membership, then apply those existing
lemmas. The [source review](lean-queryphase-performance-review-2026-09-06.md)
gives the extraction points and limitations. This is a proposed refactor;
no speedup is claimed before Lean accepts it and a benchmark completes.

Splitting proofs can prevent duplicated proof terms, as explained in
the [Lean auxiliary-lemma documentation](https://lean-lang.org/doc/reference/latest/Tactic-Proofs/Tactic-Reference/).
Each extracted lemma still needs kernel checking. For the next deliberate
attempt, use the existing cached dependency slice and enable
`-Dprofiler=true -Dprofiler.threshold=1000`, recording both CPU and wall
time with `/usr/bin/time -l`. Do not add a full build just to test this
one declaration. Full builds and cache work must use the existing
`scripts/lake-locked.sh` wrapper. More simultaneous builds cannot
directly parallelize the one busy kernel worker observed here.

No existing proof source, cache, scheduler, process priority, or running
job was changed. No new Lean build was launched. Evidence is in the
[two process snapshots](../results/lean_performance_refresh_2026_09_06.json)
and [native stack sample](../results/lean_kernel_refresh_2026_09_06.txt).
The earlier [diagnosis](lean-performance-2026-09-06.md) is retained as
historical evidence. This compiler diagnosis gives no evidence that
the Paley conjecture or Proximity Prize has been proved.
