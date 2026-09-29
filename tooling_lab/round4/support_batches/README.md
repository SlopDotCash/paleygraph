# Certified symbolic interpolation batches

This prototype avoids visiting every interpolation basis when many bases reconstruct the same affine codeword track. It keeps a complete, checkable census of every qualifying scalar/codeword/support node, including symbolic whole-field tracks. The gain is structural: deliberately aligned piecewise inputs collapse to three interpolations, while generic or interleaved inputs can still require every basis.

Zero-suppressed decision diagrams, interpolation uniqueness and covering families are established methods. The work here is an implementation connecting those ingredients to the exact provenance census developed in round2. It is not a new list-decoding theorem, a worst-case complexity breakthrough, or a prize proof.

## Mathematical certificate

Let the evaluation points be distinct, and let the code consist of polynomials of degree less than k. For an input pencil `u0+z*u1`, interpolating u0 and u1 on a k-element basis gives a track `A+z*B`, where A and B are evaluation vectors of degree-less-than-k polynomials.

Define the track's common-coordinate set

`T(A,B)={i:u0[i]=A[i] and u1[i]=B[i]}`.

Every k-basis contained in T reconstructs **exactly the same pair (A,B)**. Therefore one interpolation plus its verified common-coordinate set certifies every such basis. Moreover, two distinct tracks have fewer than k common coordinates in common: otherwise both their intercept and slope polynomials agree on k distinct points and hence are identical. Their families of certified bases are consequently disjoint.

The round2 complete cover consists of all k-subsets of disjoint coordinate blocks. A track covers exactly

`Σ_blocks binomial(|T ∩ block|, k)`

selected bases. The verifier checks each track by interpolation, recomputes its maximal T, verifies distinctness, and sums these exact counts. If the sum is the cover's complete basis count, no selected basis is missing. This proof does not enumerate the skipped bases or trust the decision diagram.

The original cover remains complete because any s-coordinate support contains a selected k-basis. For a qualifying codeword at scalar z, its polynomial and the track polynomial at z agree on that basis, so interpolation uniqueness identifies them. Within each track, every coordinate imposes either no scalar constraint, an impossible constraint, or one field scalar. Exact scalar buckets recover all qualifying supports. Repeated nodes are deduplicated by scalar and full codeword; whole-field behavior follows the round2 convention.

## Symbolic search

Each block's remaining family of k-subsets is a reduced zero-suppressed decision diagram (ZDD). Its nodes decide whether a coordinate belongs to a basis; terminals represent the empty family and the singleton family containing the empty set. Shared subgraphs, exact family counts and union-of-support masks are cached.

At each step the algorithm chooses the first remaining basis, interpolates its track, and subtracts all represented subsets of T from every block family. Subtraction works on the diagram, so it can delete a huge basis family in one operation. The algorithm terminates only when every family is empty. A separate certificate readback verifies the resulting counts and output nodes without using ZDD operations.

The first implementation had substantial overhead on singleton batches. The measured iteration added a minimum-cardinality check: when a diagram branch contains fewer allowed coordinates than its smallest represented set, it cannot contain any subset of T and remains unchanged. This is a standard pruning condition, not a new mathematical invention.

## Exact tests and measured limits

All **16** saved round1 stacks and **7** extension-field controls match round2 in every scalar/codeword/support node and every qualifying-track certificate. Two additional GF9 cases force the symbolic and materialized whole-field branches. There are **4,976** exhaustive tiny family-subtraction and sequential-subtraction checks, including empty families and k=0. These tests pass for both the initial and pruned implementations.

The new inputs use three piecewise polynomial tracks that meet at the same codeword when z=1. Aligned cases assign whole cover blocks to a track; the interleaved control assigns coordinates cyclically among the three tracks. The generator receives only the resulting input words, not the planted polynomial labels. It uses its ordinary first-remaining-basis rule.

| Input | Selected bases | Tracks interpolated | Complete output |
|---|---:|---:|---|
| Generic n32, k4, s12 | 630 | 630 | No qualifying nodes |
| Aligned n48, k6, s28 | 546 | 3 | Exactly one node, z=1 |
| Interleaved n48, k6, s28 | 546 | 546 | Exactly one node, z=1 |
| Aligned n128, k8, s80 | 3,465 | 3 | Exactly one node, z=1 |
| Aligned n256, k16, s160 | 16,507,238 | 3 | Exactly one node, z=1 |
| Aligned n512, k32, s320 | 271,753,118,987,350 | 3 | Exactly one node, z=1 |
| Aligned n1024, k64, s410 | 357459200002524597110894827605335741280427934700 | 3 | Exactly one node, z=1 |

The first four cases were compared with complete round2 enumeration. The three largest baseline enumerations were **not run**: completeness is instead established by the exact disjoint-batch certificates and readback that does not use the ZDD. These large cases use multiplicative subgroups of F65537; their inputs are deliberately constructed from the selected cover blocks. The n1024 case allocates39,362 diagram nodes and performs18 subtraction-node visits after setup, certifying the listed basis count with three reconstructed tracks.

At n128 the saved full round2 census took9.13 seconds and the pruned batch census0.00244 seconds. At n1024 the batch census took0.492 seconds. Times vary substantially under concurrent workloads and exclude neither ordinary Python overhead nor the work within the census. They are observations about these fixtures, not general performance claims. Certificate readback and input generation are separate from the stored census timing.

The pruning iteration reduces exact subtraction-node visits from16,647 to6,210 on the generic case and from16,711 to7,770 on the interleaved case. Three alternating timing repetitions **did not show a speedup on those two cases**: initial/pruned medians were0.314/0.373 seconds and0.640/0.675 seconds. The additional cardinality checks have costs; fewer recursive visits do not by themselves prove lower runtime. The aligned n128 medians were0.00466/0.00199 seconds. See the complete timing samples in [iteration_results.json](iteration_results.json).

The interleaved failure is mathematically useful. No selected cover block contains k coordinates from any planted common set, so none of those three large same-track families is discoverable within the selected cover. Nevertheless every reconstructed track reaches the same global codeword at z=1. Exact same-track batching preserves too much distinction to exploit this scalar-only convergence. A future scalar-conditioned batch instrument must also control every other scalar before it can safely skip those tracks.

## Reproduce and inspect

From the repository root, using Python's standard library:

```sh
python3 tooling_lab/round4/support_batches/run_experiments.py
python3 tooling_lab/round4/support_batches/verify_exports.py
python3 tooling_lab/round4/support_batches/compare_backends.py
```

- [symbolic_batches.py](symbolic_batches.py): the pruned ZDD family engine, exact census and interpolation/count certificate verifier.
- [results.json](results.json): all25 legacy/extension/whole-field comparisons, full new-case batch certificates, planted inputs, exact costs and timings. Large outputs include every qualifying node and every batch witness. Legacy inputs remain in the round1 ledgers; their repeated full batch exports are summarized here with deterministic certificate hashes.
- [verification_results.json](verification_results.json): readback of every new export plus both whole-field branches, rejecting eight deliberate corruptions: wrong batch multiplicity, missing/duplicate track, false common coordinate, missing/duplicate node, mismatched problem parameters and a noncanonical field coordinate.
- [initial_results.json](initial_results.json), [symbolic_batches_initial.py](symbolic_batches_initial.py), [run_experiments_initial.py](run_experiments_initial.py): archived first implementation and measured failure. The archived runner is a historical source snapshot; use `compare_backends.py` to execute the archived engine against the final engine.
- [iteration_results.json](iteration_results.json): exact old/new output equality and operation/timing comparison.

The export verifier shares field and interpolation helpers with the generator; it is a readback verifier, not an independently implemented end-to-end decoder. The equality suite compares with round2, whose scalar/support census had a separate exhaustive independent review. All stated counts are finite exact computations, not sampling estimates.

An additional [independent round4 review](../novelty/support_batches_review.json) passed261 complete coefficient/scalar oracles,1,335 explicit selected-base Lagrange replays,162 whole-field cases including87 symbolic ledgers, and2,200 sequential ZDD operations against explicit set families. Its [review script](../novelty/review_support_batches.py) does not rely on the generator for the tiny codeword/scalar census. It found no engine correctness issue.

## Complexity, prior art and full-problem relevance

Write M for the complete selected-cover basis count and D for the number of distinct reconstructed tracks. The old instrument interpolates M bases. This instrument interpolates D, with roughly `O(D*(n*k+k²))` field work plus symbolic-family operations and output costs. Initial block-family diagrams require `O(Σ_block |block|*k)` nodes. Later diagrams can grow, and D can equal M. The certificate ledger also costs `O(D*n)` field entries. Nothing here removes exponential worst-case dependence when k scales with n.

The actual Proximity Prize profile uses a degree-six extension field and requires uniform control beyond any supplied input. Extension arithmetic is exercised by the inherited GF9/GF729 controls; the large measurements use F65537, not the full prize field. The known saved coset stacks over F1009 and the large prime each retain all1,287 distinct tracks; the F17 samples save only16 to64 interpolations. Thus the new batching mechanism has not compressed those hard obstruction inputs substantially. The tool can now exhaustively audit selected structured large pencils that were inaccessible by basis enumeration. It does not prove that arbitrary or obstructing prize inputs admit this compression.

Primary foundations include Bryant's [Graph-Based Algorithms for Boolean Function Manipulation (1986)](https://www.cs.cmu.edu/~bryant/pubdir/ieeetc86.pdf) and Minato's [Zero-suppressed BDDs and their applications (2001)](https://eprints.lib.hokudai.ac.jp/dspace/bitstream/2115/16895/1/IJSTTT3-2.pdf), which develops the ZDD method introduced in1993. Interpolation and the complete cover are inherited from the already reviewed round2 instrument.

A bounded local search found a prior repository note describing vectorized batches of k-subset interpolation in `docs/kb/deltastar-444-S1-list-growth-VERDICT-2026-06-15.md`; its named `probe_wfS1_engine.py` was not present in the current file inventory. The note is evidence that batching interpolation itself is existing local work. Targeted primary/web and local searches did not identify this exact same-track/ZDD/count-certificate implementation, but this is not proof of historical absence or a claim that its elementary mathematical ingredients are new.
