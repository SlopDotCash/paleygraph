# Complete pencil lists without a codeword anchor

This prototype removes the exact-codeword-anchor requirement from one structured certificate route. It blindly discovers marginal polynomial covers, combines them into a verified affine polynomial-track partition, and compiles every qualifying codeword and maximal agreement support across the whole field. Both n1024 demonstrations have **no codeword anywhere on the input pencil**.

The two-track example nevertheless has two qualifying codewords at every scalar, or 131074 scalar/codeword nodes. The three-track example has exactly three qualifying scalars, with one codeword at each. These are complete ordinary agreement lists for simplified inputs. The official Proximity Prize MCA event, which imposes additional nonjointness conditions, is not certified.

Interpolation, root counting, affine coefficient arrangements, finite-field factorization and greedy set cover are known methods. The implemented contribution is the blind-discovery interface, a compact exact whole-field list compiler, and measured failure diagnostics. No new mathematical foundation or historical novelty is claimed.

## Reproduce and inspect

```sh
/opt/miniconda3/bin/python tooling_lab/round6/pencil_tracks/run_experiments.py
/opt/miniconda3/bin/python tooling_lab/round6/pencil_tracks/verify_results.py
```

`pencil_tracks.py` provides:

- `compile_certificate(F,domain,k,s,u0,u1,tracks)`: verifies supplied tracks and compiles a complete symbolic list if the strict cap holds.
- `discover_pencil(...)`: receives no planted labels or coefficients; uses the frozen round5 discovery tool for `u0,u1`, or two distinct scalar slices; intersects the partitions and optionally reassigns coordinates using full verified joint supports.
- `materialize_scalar(F,certificate,z)`: returns every qualifying codeword, its coefficients and maximal agreement support at z, including exact deduplication at collisions.
- `verify_certificate` and `verify_discovery`: source-level readback of the full partition, scalar arrangement, symbolic list and marginal algebra witnesses.

`results.json` saves full discovery and certificate records. `verification.json` binds sources, frozen dependencies and result bytes. This readback shares production algebra helpers; independent review is maintained separately under `../pencil_tracks_review/`. Earlier rounds and the Proximity Prize repository are unchanged.

## Mathematical contract

Let the domain have n distinct field elements, and let every code polynomial have degree<k. Suppose verified disjoint regions R_j partition all coordinates, and polynomials a_j,b_j of degree<k satisfy

\[
u_0(x_i)=a_j(x_i),\quad u_1(x_i)=b_j(x_i)\qquad(i\in R_j).
\]

At any scalar z the received word equals the polynomial `a_j+z b_j` on R_j. If a code polynomial g differs from **every** track polynomial at that scalar, each difference is nonzero of degree<k, hence

\[
|\{i:g(x_i)=u_0(x_i)+z u_1(x_i)\}|
\le \sum_j\min(|R_j|,k-1)=C.
\]

Therefore `s>C` forces every qualifying codeword to equal at least one track polynomial at that scalar. No anchor is used. Track collisions cause no gap in the argument: identical polynomials at a particular scalar are subsequently deduplicated.

For each track j and coordinate i, agreement is the affine equation

\[
(a_j(x_i)-u_0(x_i))+z(b_j(x_i)-u_1(x_i))=0.
\]

Its solution set is exactly the whole field, empty, or one scalar. Each track exports these three cases as always-coordinates, never-coordinates, and one-scalar buckets. Two distinct affine tracks collide exactly when the coefficient vector equation `(a_j-a_l)+z(b_j-b_l)=0` holds. That has no solution or one solution; tracks identical for every scalar are merged before compilation. Since n≥k and the domain is distinct, equality of the codeword evaluations is equivalent to equality of the degree<k coefficients.

Let E be the union of every coordinate exception and every track-collision scalar. Outside E, all supports are fixed and all retained track polynomials are distinct. Inside E, the compiler groups qualifying tracks by their exact coefficient vector. Consequently

\[
\#\text{scalar/codeword nodes}
=(q-|E|)L_{\rm generic}+\sum_{z\in E}L_z.
\]

No field-scalar enumeration is needed for this identity. A finite-event node references one representative track and all coincident qualifying track IDs. Its full codeword is `a_j+z b_j` on the saved domain; its maximal support is `always_j ∪ bucket_j[z]`. This symbolic recipe retains every coefficient, codeword and support without duplicating an n-entry array at every event.

## Blind discovery and its limits

The frozen round5 routine reconstructs a small monic bivariate annihilator, recovers polynomial branches by bounded root splitting/Hensel lifting, checks complete polynomial identities, and discovers a complete coordinate partition. Its exact rank/kernel/inconsistency witnesses are included unchanged. Marginal list completeness is **not** required here: only the recovered piece cover is needed, because the final joint partition has its own root-count proof.

The default marginals are u0 and u1. For finite slices at distinct z0,z1, a pair of recovered slice polynomials h0,h1 gives

\[
b=(h_1-h_0)/(z_1-z_0),\qquad a=h_0-z_0b
\]

on the intersection of their assigned supports. Empty intersections disappear and coincident affine tracks merge. The optional refinement computes every candidate's maximal **joint** support and applies ordinary greedy set cover, reassigning coordinates only where both polynomial equations are checked. It can remove unnecessary fragments caused by ambiguous marginal matches. Neither partition intersections nor greedy reassignment are asserted to recover the minimum joint cover.

A failure means one of two different things: a marginal cover was not discovered within the fixed ansatz/budgets, or a joint cover was found but its root-count cap is too high. Neither outcome means the true list is empty or that a better certificate cannot exist.

## Exact results

All n1024 cases use F65537, k64, s410 and the multiplicative subgroup domain. The planting data is kept outside discovery and revealed only for post-discovery comparison. All cases below have no exact code anchor.

| Fixture | Retained tracks | Cap | Coordinate/collision exception scalars | Generic list size | Qualifying scalars | Scalar/codeword nodes |
|---|---:|---:|---:|---:|---:|---:|
| n192, F257, two skew tracks, k8,s80 |2|14|139|2|257|514|
| n192, F257, three pair-colliding tracks, k8,s80 |3|21|3|0|3|3|
| n1024, two skew tracks |2|126|1011|2|65537|131074|
| n1024, three pair-colliding tracks |3|189|3|0|3|3|

“Skew” here means the two coefficient-space lines never meet at a common scalar. Their differences are `g(X)+z`, where g is nonconstant, so no scalar can make the two full polynomials identical. Each region has at least k points, which precludes a global code anchor.

The three-track fixture has coefficients `a_j=f+alpha_j*g`, `b_j=h+beta_j*g`, with alpha=(0,1,3), beta=(0,1,2). Pairwise collisions occur at three distinct scalars, while all three tracks never coincide. At n1024 a pair collision joins two regions of size approximately 341, exceeding s410; away from collisions every individual region is below threshold. Exact support records handle any additional shared roots of g.

The saved large blind runs take 12.54s and 16.97s, including marginal algebra, factor discovery, joint-cover verification and symbolic compilation. The full experiment suite takes 73.63s. These are machine measurements, not asymptotic claims. Every scalar was separately evaluated against all recovered tracks on both n192 cases; large cases use the completeness proof plus exact count comparison to the withheld planting witness and selected scalar materializations. No large-field codebook enumeration is claimed.

The compiler passes 833 exhaustive codeword/scalar checks: all 729 F3 pencils of length3 with k1, 64 arbitrary F5 cases, and 40 GF9 cases. They include 1779 nodes, 604 pencils without an anchor, 334 whole-field qualifying sets and 697 cases with track collisions. Blind discovery is additionally tested on 40 hidden two-track inputs: 13/20 F5 and 20/20 GF9 examples certify successfully, with all codewords and scalars checked for every success (464 nodes). Seven F5 cases remain explicitly incomplete. Five corrupted full certificates are rejected.

An inherited n48 example that does have an anchor matches the frozen round4 scalar-fiber certificate at all 65537 field scalars. The anchor route remains a simpler certificate when its premise holds; the new route handles additional structured pencils rather than superseding that reduction.

## Two measured iterations

1. **Ambiguous coordinate assignment.** For an F17 length12 pencil, separate marginal partitions create three joint regions, including a singleton ambiguity fragment. Their cap is 3, so threshold s3 fails the strict test. Full joint-support reassignment removes the unnecessary fragment, leaves two tracks with cap 2, and certifies 34 nodes across all 17 scalars. Every reassigned coordinate satisfies both polynomial equations exactly.
2. **Choice of marginals.** A length192 pencil has a two-by-two grid of joint regions. The default intercept/direction covers each require two polynomials, so their four intersections give cap 28<s80 and certify exactly two nodes at one scalar. Using scalar slices 0 and 1 makes the second slice require at least four polynomials; its exact monic inconsistency witnesses rule out all degrees allowed by the max-three discovery budget. That route fails while the intercept/direction route succeeds on the identical pencil. This is a parameterization-sensitive discovery limit, not a difference in the underlying list.

A three-by-three grid has nine recovered joint regions and cap 63 at s50, so compilation correctly refuses a complete-list certificate. Three generic length64 controls fail to produce small marginal covers. The split-marginal strategy does not remove the need for stronger joint discovery or a sharper completeness argument.

## Cost and connection to the full problem

For t verified tracks, there are at most tn coordinate exceptions and t(t-1)/2 pair-collision events. Their scalar values may coincide, so the actual set E can be much smaller. Agreement/collision construction uses O(tnk+t²k) field operations with ordinary polynomial evaluation. Finite-event deduplication adds O(|E|tk) field operations; this simple implementation additionally reconstructs/sorts supports during verification, costing up to O(|E|tn log n) ordinary operations. The stored ledger has O(tn+|E|t+tk) field elements or coordinate references, and does not scale linearly with q.

Blind discovery inherits round5's dense interpolation cost and explicit incomplete factor-search status; at fixed rate its dense linear algebra grows cubically with n. Marginal intersections can multiply piece counts and destroy the strict cap. GF9 arithmetic is tested exactly, but the inherited extension-field constructor remains a toy implementation that precomputes every element. There is no scalable official F_(p^6) backend or official-size prize certificate.

The new interface makes the lack of an anchor an addressable premise for structured pencils. It leaves the actual hard questions visible: discover a sufficiently small joint polynomial-track cover on arbitrary stacks, control uncovered coordinates without assuming such a cover, and connect ordinary proximity lists to the official nonjointness event. None is declared solved here.

## Existing mathematics and local prior art

The local [fiber-Chebyshev probe](/Users/shawwalters/proximityprize/scripts/probes/probe_rate_quarter_p1_fiber_chebyshev.py) already evaluates ratios of codeword-pair differences, distinguishes common roots, and verifies the k−1 foreign-fiber bound away from a proportionality scalar. Its [accompanying note](/Users/shawwalters/proximityprize/docs/kb/deltastar-466-rate-quarter-fiber-chebyshev-2026-07-11.md) records formal lemmas for both the root bound and the at-most-one exceptional scalar. It also explains why the bound fails for arbitrary input-relative ratios; the present compiler applies root counting only to actual polynomial differences on verified regions.

The [coefficient arrangement note](/Users/shawwalters/proximityprize/docs/kb/deltastar-466-design-matrix-affine-clusters-2026-07-09.md) explicitly constructs interpolation pencils and their incidence with scalar/codeword points. Round2 of this tooling lab already grouped common-coordinate supports and one-scalar agreement buckets. The new work compiles those known structures from a discovered small partition and handles all polynomial collisions explicitly.

[Ben-Sasson, Carmon, Ishai, Kopparty and Saraf, Proximity Gaps for Reed–Solomon Codes (2020; current revision2021)](https://eccc.weizmann.ac.il/report/2020/083/) analyzes classical algebraic decoders on formal elements of affine spaces. The report page records a proof correction in revision3 without a change to the main statements. This is relevant prior art for polynomial dependence on a scalar; the certificate here proves no new general proximity-gap theorem. The interpolation/factorization sources and local Hensel implementations audited in [round5](../../round5/piece_discovery/README.md) remain directly applicable.
