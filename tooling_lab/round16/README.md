# Round16: conditional decoding and profiles that permit early rejection

[Round15](../round15/README.md) shows that a62-coordinate projective class forces the minimum complete-arc partition cost up to136155 queries on a specific1024-point degree-bounded polynomial space. It also exposes many candidates that agree with all repeated coordinates but fail elsewhere.

The [conditional compiler](conditional_cover.py) uses the received word to split the parameter space. If the chosen projective direction is w and m_b coordinates ask for theta.w=b, choose a cutoff h. Values with m_b>h get individual affine branches of one lower dimension, each requiring s-m_b agreements outside the class. Every remaining candidate needs at least s-h agreements outside the class and is covered by one shared residual branch. Each branch receives a separately checked arc certificate and an explicit parameter lift.

For span(1,P,XP), conditioning fixes the constant parameter. The two remaining directions are P and XP; outside the62 roots, all pairs are independent because their determinant is P(x)P(y)(y-x). The construction aims to replace a geometry-only all-s-set cover with a received-word-dependent family of covers. This distinction is essential: a lower query count here would not contradict round15's restricted optimum.

The initial scope is dimension three, with residual thresholds at least the residual dimension. Input or arithmetic cases outside the implemented scope are reported. If puncturing leaves fewer coordinates than the original degree bound, basis polynomials are reduced modulo the punctured domain's vanishing polynomial; their evaluations and original parameter interpretation are preserved.

The [compiled results](conditional_results.json) and [constructor-free review](conditional_review.json) now establish the following query families. The reviewer evaluates all original and punctured polynomials directly, checks every residual query with integer determinants, reconstructs the received histogram and verifies the affine parameter lifts and exhaustive branch routing.

| Received word | Cutoff h | Repeated-class frequencies | Residual queries | Fixed-equation queries | Total |
|---|---:|---|---:|---:|---:|
| Boundary |0|62|57570|13543|71113|
| Late near miss |0|62|57570|13543|71113|
| True guard mismatch |1|61,1|57760|13132|70892|

The first two covers reduce the query count from136155 to71113 on the same received words. This is a comparison of query counts in two different certificate models, not a timing measurement or a general decoding complexity theorem. The verifier does not certify optimal cutoff selection or optimality over all possible conditional query families.

The [full-scale decoder](conditional_pruning.py) now retains every branch query and decision. Its [independent review](decoding_review.json) solves the queries with integer Cramer determinants, fully evaluates every branch candidate and checks all parameter lifts. All three final lists match the frozen complete supplied-space baselines. The lists have sizes1,0,1. This validation also found a substantial cost regression that the query-count estimate missed.

## The first scaled result changed the profile design

The query-minimizing profile omits singleton checks in some branches. A false candidate can then require19 or29 coordinate tests before its agreement upper bound drops. On the guard-mismatch word this produces1467372 coordinate evaluations, substantially more than the geometry-only baseline's138310, despite using fewer queries.

The [revised compiler](guard_profiles.py) reserves at least one unqueried coordinate in each branch, allowing an early test within the same exact transcript-bound filter. It does not reject a candidate solely for disagreeing with a singleton. The [full independent review](guard_review.json) checks all branch candidates and verifies that the final lists still match both earlier decoders.

| Word | Baseline queries | First conditional queries | Revised queries | Baseline coordinate evaluations | First conditional evaluations | Revised evaluations |
|---|---:|---:|---:|---:|---:|---:|
| Boundary |136155|71113|71520|923494|452949|74419|
| Late near miss |136155|71113|71520|922753|450835|72345|
| True guard mismatch |136155|70892|73554|138310|1467372|76276|

Coordinate counts include adaptive candidate tests, branch output materialization and final whole-word checks; query solving and certificate construction/validation are separate. The experiment shows a useful finite tradeoff, not a universal speedup. The late near miss still takes661 tests inside its fibre. The [derivation](DERIVATION.md) explains both completeness and the cost failure.

The [exhaustive toy controls](toy_controls.json) decode every3125 received word over F5 in a length5, degree<4, dimension3 space with a repeated class. Every result is compared against all125 polynomials in the supplied space. The optimum uses one query for2500 words and three queries for625 words. Four corrupted certificates are rejected. These controls exercise polynomial reduction after puncturing, but their optimum always uses a single branch.

The [basis controls](basis_controls.json) exhaust another3125 words in each of two affine presentations: a nonfirst pivot and a projective direction with three nonzero entries. Thus9375 received-word/presentation cases pass a complete125-polynomial oracle. The [routing controls](routing_controls.json) deliberately test nonminimal cutoffs on125 further structured words. Their250 covers include25 mixed residual/fibre covers and100 covers with two fibres; all decoded unions match the complete oracle. Three corrupt routing certificates are rejected.

The [production-decoder controls](decoding_controls.json) add750 complete transcript/lift reviews with full125-polynomial oracles, including affine origins and375 multibranch cases. Four corrupted decoding artifacts are rejected. A tiny rank-two input correctly reports that no singleton-reserving profile exists. These checks supplement the earlier exhaustive toy tests, whose helper decoder predates the production implementation.

The source audit in [round15](../round15/prior_art.md) already identifies prior coordinate conditioning and affine-space pruning. Historical novelty is unestablished. The present advance is a checked received-word-dependent certificate, a verified failure of query-count-only optimization, and a profile revision that improves both measured query and verification counts against the baseline. No cluster discovery or full-ambient list guarantee is provided. The [final manifest](manifest.json) binds this evidence; the [earlier checkpoint](checkpoint.json) is historical. [Round17](../round17/README.md) returns to Paley arithmetic realizability.

```sh
/opt/miniconda3/bin/python3 tooling_lab/round16/conditional_cover.py
/opt/miniconda3/bin/python3 tooling_lab/round16/conditional_verifier.py
/opt/miniconda3/bin/python3 tooling_lab/round16/toy_controls.py
/opt/miniconda3/bin/python3 tooling_lab/round16/basis_controls.py
/opt/miniconda3/bin/python3 tooling_lab/round16/routing_controls.py
/opt/miniconda3/bin/python3 tooling_lab/round16/conditional_pruning.py
/opt/miniconda3/bin/python3 tooling_lab/round16/decoding_review.py
/opt/miniconda3/bin/python3 tooling_lab/round16/guard_profiles.py
/opt/miniconda3/bin/python3 tooling_lab/round16/guard_review.py
/opt/miniconda3/bin/python3 tooling_lab/round16/decoding_controls.py
/opt/miniconda3/bin/python3 tooling_lab/round16/verify_round.py
```
