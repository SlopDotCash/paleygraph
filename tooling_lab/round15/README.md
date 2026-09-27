# Round15: projective capacities, verified local repair, and exact decoding

This round turns a failed coordinate-query allocation into a certified optimum and a cheaper construction procedure. It then tests the result on received words at the agreement boundary. [Derivation](DERIVATION.md), [independent review](capacity_review.json), and [prior-art audit](prior_art.md) give the mathematical scope and evidence.

For a complete-arc partition with t queried blocks, u unqueried coordinates, projective class sizes c_j and z zero columns, necessarily

```
u >= z + sum_j max(c_j-t,0).
```

The [profile optimizer](projective_profile.py) combines this capacity constraint with the agreement condition t(d-1)+u<s and balances block sizes to minimize the query count. A cyclic assignment gives a valid allocation in the parallel-class relaxation. All higher-rank query independence remains an arithmetic obligation. If the lower bound is attained by an actual arc cover, it is optimal among complete-arc block partitions with the same size cap, not among arbitrary query families.

## A concrete obstruction on the actual domain

Use the saved1024-point domain over F65537 and the synthetic space span(1,P,XP), where P vanishes at62 domain points. The space has rank3, degree<64, no zero columns, one projective class of size62 and962 other classes of size one. This is a deliberately constructed polynomial control, not a discovered candidate cluster.

At agreement96, the unconstrained70070-query profile has47 queried blocks and one unqueried coordinate. It can accommodate at most48 members of any class. The class of size62 makes that profile impossible under every reshuffling.

The capacity-constrained optimum costs136155 queries with33 queried blocks and29 unqueried coordinates. The first allocation contains two dependent triples. [Bounded reshuffling](capacity_search.json) finds a valid allocation on attempt23. The [standalone certificate](root_class.certificate.json) passes both modular elimination and integer determinant checks. An independent enumeration of all feasible t,u checks173 profiles and establishes the same finite optimum.

## Repair instead of restarting the allocation

The [exchange prototype](exchange_repair.py) retains the actual dependent triples and considers swaps preserving block sizes and projective capacities. Only triples containing a moved coordinate can change. A swap is accepted when it strictly reduces the complete defect count; unqueried singletons can participate too. This is bounded local descent, with no completeness claim for the search.

Two swaps repair the initial allocation, using1682 additional rank tests after its136155-query audit. The observed restart run used3131565 search tests. Final verification and certificate checks are additional in both cases, and these counts do not establish a wall-clock speedup. The independent reviewer rescans both entire affected blocks, verifies all changes, and checks the [repaired certificate](repaired.certificate.json).

## Decoding evidence and limits

The completed [boundary tests](boundary_results.json) use the round14 pruner and independent Cramer/full-word reviewer on three deliberately planted received words. All408465 queries are independently replayed;408351 candidate words and418151424 coordinate values are checked.

| Received word | Known polynomial agreement | Returned list size | Explicit coordinate evaluations, including output | Full candidate scan |
|---|---:|---:|---:|---:|
| Boundary |96|1|923494|139373568|
| Late near miss |95|0|922753|139407360|
| True guard mismatch |96|1|138310|139370496|

The known candidate requires993 coordinate tests in each case, showing that the worst individual verification cost remains linear. In the first two cases it occurs in only one query; the true guard mismatch occurs four times. Frequency cutoffs or unconditional guard rejection would be unsound.

These are complete lists **inside the supplied space**. Agreement96 is below the ambient Johnson threshold for n1024,k64. The older two-piece ambient oracle has exclusion cap126 and cannot be used here. No full-ambient decoder or prize theorem is asserted.

The [capacity controls](parallel_controls.json) exhaust33841 coordinate partitions across47 capacity patterns and417 thresholds. The [search controls](search_controls.json) add32 realized polynomial configurations,1120 literal optimal-profile partitions, four rejected corrupt exchange histories, a successful singleton exchange and an explicit higher-flat obstruction. Six distinct points on one line violate a rank-two capacity even though projective-class counts alone permit the profile.

## Reproduce and continue

From the repository root, using the installed Python with NumPy:

```sh
/opt/miniconda3/bin/python3 tooling_lab/round15/capacity_search.py
/opt/miniconda3/bin/python3 tooling_lab/round15/exchange_repair.py
/opt/miniconda3/bin/python3 tooling_lab/round15/capacity_review.py
/opt/miniconda3/bin/python3 tooling_lab/round15/search_controls.py
/opt/miniconda3/bin/python3 tooling_lab/round15/boundary_runs.py
/opt/miniconda3/bin/python3 tooling_lab/round15/verify_round.py
```

The results are frozen by [the manifest](manifest.json); rerunning producers changes recorded timing fields and requires a fresh review and manifest. The earlier [preflight checkpoint](checkpoint.json) is a historical snapshot with an older README hash; it is not the final verification artifact.

The repeated class suggests a further tool: condition on its shared equation, count its exact agreements, and reduce the dimension before building queries. [Round16](../round16/README.md) implements the initial received-word-dependent compiler. Existing arc coloring, convexity, local search and affine conditioning are explicitly acknowledged. Global historical novelty is unestablished. Paley arithmetic tooling remains a separate pending lane, and neither prize is being proved.
