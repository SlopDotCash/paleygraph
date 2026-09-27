# Round 12: multiplicities, full-length certificates and actual pruning

This checkpoint completes three supplied-space tools. [Round 11](../round11/README.md) is preserved. No cluster discovery, general decoder, prize proof or historical novelty certificate is claimed.

The [capacity compiler](capacitated_pinning.py) retains zero columns and repeated projective directions. A nonnegative multiaffine basis polynomial is concave along each pairwise occupancy transfer, so some minimum fills or empties every class except possibly one. A subset transform evaluates these corners. This is a standard rounding mechanism; its precise argument and limits are in the [derivation](DERIVATION.md).

The complete [actual replay](actual_results.json) passes all **2,894** rank-three spaces spanned by quadruples of the saved 18 small candidate nodes, including the 2,839 spaces excluded by the preceding simple-matroid interface. It checks 1,620,640 raw triples and 7,971,600 agreement sets across 1,825 distinct dependency catalogues. The worst minimum is **50** informative triples, versus 120 in the 55 simple cases. The initial controls cover 596 patterns, 3,936 thresholds and 45,952 subsets. Complexity remains exponential in projective class count, capped at 22.

The [polynomial adapter](polynomial_certificate.py) exports self-contained certificates. The [standalone verifier](certificate_verifier.py) imports no constructor: it reconstructs classes by row reduction and checks small optima by searching all integer occupancies, without using corner reduction. For a large simple space it reconstructs all pair spans with separate grouping logic. The [interface review](certificate_verification.json) passes six certificates, rejects 12 corruptions and treats a verification budget limit as unfinished verification.

On the saved F65537, n=1024, degree<64, agreement>=410 space, all 523,776 pair spans certify

```
11,401,439 <= minimum informative triples in a 410-set <= 11,402,277.
```

The complete pair partition has 2,785 nontrivial lines and accounts structurally for all 178,433,024 triples. Those triples were not literally enumerated. The exact optimum is unresolved. The degree-only baseline is 7,023,974; the certified success probability is at least 1,628,777 / 25,490,432.

The [pruning sampler](pruning_sampler.py) turns that probability into an exact integer trial budget for a supplied affine three-dimensional space. It samples triples, solves their coordinate equations and checks each distinct candidate's complete agreement. This pruning mechanism and its union bound are known; see the [prior-art audit](prior_art.md).

The [actual runs](pruning_results.json) recover the complete lists for all 17 small received scalars and two large scalars. The small oracle enumerates all 4,913 space members. Each large received word is assembled from two degree<64 polynomials, both agreeing at 512 positions; any other such polynomial agrees at most 126 times, so the complete list is exactly those two. A [separate Gaussian replay](pruning_verification.json) checks all **9,234 queries** and **2,619 distinct candidate evaluations**. Eight altered traces are rejected, and a repeated-query control demonstrates why arbitrary replay queries do not inherit a uniform-sampling guarantee.

At the same per-run failure bound 2^-42, the large instance certificate requires **945 trials**, compared with **1,554** from the universal degree bound. Across all 19 runs the union of failure bounds is 49/2^46 < 2^-40, conditional on independent uniform draws. Recorded outputs are checked exactly on these particular words. This is a query reduction, not an end-to-end speedup: the full certificate took about 34 seconds to construct and 105 seconds to verify in the recorded run, while each large sampling run took about half a second after preparation. Timings depend on machine load and implementation.

Run from this directory with `/opt/miniconda3/bin/python3`:

```sh
python3 polynomial_certificate.py hard_rank3.input.json
python3 pruning_sampler.py hard_rank3.certificate.json received.json 40
python3 verify_round.py
```

`received.json` must be a canonical field-valued array of the certified length. Saved random traces can be replayed through the Python interface. The integration [manifest](manifest.json) binds all artifacts and checks that round 11 is unchanged. Earlier checkpoint files retain their historical, narrower scope statements; the standalone and pruning reviews supersede their pending-work notes.

The next design question is whether a much smaller deterministic family of informative triples can cover every large agreement set. A cheap certificate of that property could be more useful than an expensive estimate of its minimum density. This is a candidate tooling application of known covering and pruning ideas, not a claim of a new general derandomization theorem.
