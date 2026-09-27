# Lean performance check, 08:42 UTC / 04:42 EDT

The existing v19 process is still CPU-active in Lean's kernel. It has no
terminal proof verdict and no defensible completion ETA. The machine also
has severe memory pressure. Waiting longer may still end in a timeout,
as happened to the preceding v18 attempt.

At 08:41:47 UTC, PID13687 had run for 3h25m19s and accumulated
50m44.54s of CPU time. Its source hash is unchanged from the earlier
diagnostic. Its 3,780-byte log was last modified at 05:17:10 UTC.
The fresh three-second stack sample again contains kernel
`type_checker::is_def_eq` calls in a worker, while the main thread and
most other workers wait. This establishes activity in the expensive
checking path; it does not establish a completion percentage.

The sample reports a 42.3 GB process footprint. The much smaller resident
set reported by `ps` is not a contradiction: resident memory and total
footprint are different measurements under compression and swapping.
The machine has 128 GiB RAM and 16 logical CPUs. At the memory check,
about 29 GB of swap was in use; load average was about115 at the process
snapshot. These are machine-wide observations. They support resource
contention as an additional concern but do not attribute all pressure to
this single process or supply a predicted speedup from freeing memory.

The [checkpoint](../results/lean_performance_checkpoint_2026_09_06_0842.json)
records the fresh evidence and hash of the raw stack sample. The prior
[source review](lean-queryphase-performance-review-2026-09-06.md) locates
the large dependent-type conversions and the existing logical fold API
that can replace repeated proof work.

## Prepared experiment

The [source transformer](../scripts/split_queryphase_benchmark.py) prepares
four standalone diagnostic candidates from the exact observed v19 file:

- Prefix only, to establish the cost of the imports and preceding lemmas.
- Positive branch, with its branch hypothesis explicit.
- Zero branch, with its branch hypothesis explicit.
- Both named branch lemmas followed by a small case-split wrapper with
  the original theorem signature.

The generator preserves both branch bodies apart from indentation, checks
the original source hash, enables the installed Lean profiler, and refuses
to overwrite an existing output directory. The generated candidates and
hashes are in the [manifest](../results/lean_queryphase_split_candidates_2026_09_06.json).
Their files are under `/tmp/paley-queryphase-benchmark-20260906-0845`.
Generation and source-preservation assertions passed. **The candidates
have not been Lean-validated or benchmarked.** Splitting a declaration
does not itself guarantee faster checking; this experiment isolates costs
before changing the mathematical proof.

Run these candidates sequentially once the existing attempt ends or is
deliberately replaced. Use the checkout's `scripts/lake-locked.sh env lean`
with `-DautoImplicit=false`, record wall/CPU time and peak memory, and impose
a wall-time cap on each diagnostic run. A timeout is a measured lower
bound on required time under those conditions, not a successful check.
If one branch remains expensive, replace its dependent-index conversions
with the existing `logical_queryFiberPoints_eq_fiberEvaluations` and
`logical_computeFoldedValue_eq_iterated_fold` interfaces identified in
the source review. That semantic refactor still requires validation.

No existing Lean process was stopped, reprioritized, or changed. No extra
Lean process, cache operation, full build, or Prove2Me submission was
started. This performance work does not change the Paley proof status.
