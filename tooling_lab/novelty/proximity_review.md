# Independent proximity prototype review — 2026-09-05

Reviewed `proximity/deformation_microscope.py`, `stack_provenance.py`, `verify_artifacts.py`, README, source metadata and saved deformation/stack results. No source in that lane or the original proximity checkout was edited.

## Mathematical findings

The projected-kernel formula is correct. Setting the leading coefficient to zero identifies the kernel of the leading-coefficient projection on `ker(T_ell)` with `ker(T_(ell-1))`; rank-nullity therefore gives the reported difference of nullities. The leading solution spaces form a decreasing filtration. This is a dimension of the initial syzygies extendable along the **specified** matrix trajectory, not the dimension of the full jet solution space or a proof of integration along an arbitrary curved root trajectory.

The tangent constraint `L M_1(v) R=0` is correct for preservation of **all** initial syzygies to first order: `M_1(v)c_0` must lie in the image of `M_0` for every `c_0` in its kernel. Over a field, the left kernel detects exactly failure of image membership. The distinction from subgroup-admissible tangent motion is essential and is stated correctly: for fixed p not dividing n, the roots of `X^n-1` are simple, and each root has zero first-order motion inside that fixed finite subgroup scheme. Ambient tangent dimension 15 is not a subgroup deformation dimension.

The binomial tangent count `2d+1` and one-root fragility arguments are valid under the stated assumptions of three distinct cosets and d>=2. The simple-root differential maps velocities bijectively to degree-below-d coefficient perturbations. The constant-syzygy condition imposes d-1 independent constraints. Replacing one root adds a nonzero degree-(d-1) polynomial to one binomial, which cannot lie in the span of the two unchanged binomials, `span(1,X^d)`.

The subset parity census is exact for the configured distinct evaluation points and `k<s<=n`. Vanishing dual checks is equivalent to interpolation by a polynomial of degree below k on that subset. A nonzero syndrome direction supplies at most one finite affine scalar. Both syndromes zero correctly marks a whole-field common support and an incomplete explicit node ledger. All saved stacks are in the complete-ledger branch.

Deduplicating by full agreement mask is safe: any mask has at least s>=k points, and two degree-below-k polynomials agreeing with the same word at k distinct points coincide. Consequently a decoded support of size a contributes exactly `binom(a,s)` subset certificates. Different decoded codewords cannot share a size-s certificate. This justifies the saved multiplicity check.

Affine-track generation from all pairs of distinct scalar/codeword nodes is complete for tracks containing at least two scalar values. Pair aggregation includes every node lying on such a track. For each non-common coordinate, a nonzero affine function of the scalar has at most one root, yielding `m(s-t)<=n-t` when t<s. The README correctly declines to sum capacities over overlapping tracks. The omitted projective point at infinity remains a real scope limitation for any later budget comparison.

## Independent finite checks

[review_proximity.py](review_proximity.py) provides brute-force oracles independent of the prototypes' elimination and syndrome algorithms:

- All 256 pairs of 2-by-2 coefficient matrices over F2, with a fixed noncommuting nonzero third coefficient: explicit enumeration of leading vectors and lifts through order two matches `lift_survival`.
- All 729 received-word pencils for RS(3,1) over F3 at agreement threshold two: enumerating every scalar and every codeword matches the census, including common-support cases.
- All 16 saved stacks, 95 nodes and 383 tracks: independently reconstruct every pair-defined affine track, recover its full membership, verify every agreement mask and exact non-common coordinate allocation, and compare the complete track key set.

All pass. Hashes and counts are saved in [proximity_review.json](proximity_review.json). This lane did not independently repeat the large-field exhaustive agreement-subset censuses, certify the external Lean files, or verify a prize theorem. The saved artifact verifier also shares algebra helpers with discovery, so it should be described as readback verification, not an independently implemented end-to-end census.

## Interpretation and remaining limits

No substantive mathematical correction was found in the README's elementary derivations. They are ordinary linear algebra, interpolation and first-order deformation facts; the document correctly makes no invention claim for them. The quantitative rows are finite experiments, not uniform bounds.

The phrase that the experiment demonstrates why “the joint instrument is necessary” is stronger than its evidence: it shows that the tested triple signatures do not determine the observed pencil incidence. Many possible additional diagnostics could supply that missing information. Prefer that precise statement over necessity of this particular implementation. Likewise, first-order preservation followed by second-order failure refutes a general heuristic in the tested directions, not every possible notion of stability.

The useful next tooling question is whether exact support and track provenance can be compressed with a certified completeness guarantee. The existing binomial coefficient subset census is an oracle for testing such compression on small codes, not a production-scale algorithm.
