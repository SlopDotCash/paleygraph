# Final round4 documentation and scalar-fiber review

Reviewed2026-09-05. Earlier review scripts, evidence and documents are unchanged. The new [review_scalar_fibers.py](review_scalar_fibers.py) and [scalar_fibers_review.json](scalar_fibers_review.json) record this follow-up.

The [marked-moments README](../marked_moments/README.md) accurately distinguishes its finite identity, actual information-loss witnesses, and the very small correction envelope in the reported critical-scale examples. Its q>=17 guard, distinction between cell-only data and genuine-matrix provenance, factor in the one-contraction reduction, and limitation of n=d tests agree with the independent review. Its million-prime claim of23 canonical relation classes and149-bit aggregate coefficients matches the saved scale record. The total change from the unmarked average is correctly distinguished from the relation correction cQ. No mathematical overclaim was found. The two links to the ablation README were temporarily pending during inspection; that file now exists.

## Scalar-fiber proof

Let C be the evaluation code of polynomials of degree less than k at n distinct field elements, and suppose `f=u0+zstar*u1` belongs to C on every coordinate. At any z different from zstar, the map

```
h -> g=f+(z-zstar)*h
```

is a bijection from C to C. It preserves the complete coordinate agreement support between h and u1 and between g and `u0+z*u1`, since multiplying a nonzero difference by the nonzero field scalar `z-zstar` cannot create or remove zeros. Thus one complete static list for u1 gives every nonanchor scalar fiber, including all supports. This is the usual linear-code symmetry, not a new algebraic principle.

At zstar there is exactly one qualifying codeword when s>=k: a second codeword would differ from f by a degree-less-than-k polynomial with at least k distinct roots. The prototype correctly rejects s<k. At that lower threshold the anchor fiber can have many codewords; the generator's recorded21-codeword example over F5 is consistent with this obstruction.

The anchor discovery routine is complete. Interpolate u0 and u1 on the first k coordinates and compute their two residual vectors. The pencil belongs to the code exactly when their affine combination is zero. A nonzero direction residual determines the only possible scalar, which is then checked on every coordinate. Zero direction residual gives either no anchor or the entire pencil in the code. These cases require no scan over scalars.

For the static list, the supplied partition assigns each coordinate to a degree-less-than-k piece polynomial matching u1 there. Equal piece polynomials are merged. Any polynomial different from every piece agrees on piece j at at most `min(region_size_j,k-1)` coordinates. Summing these disjoint-region bounds proves completeness whenever

```
s > sum_j min(region_size_j,k-1).
```

The inequality must be strict. At equality an unlisted polynomial can qualify, so the correct result is an unavailable certificate. The routine checks that the regions partition all coordinates, validates the piece values, and evaluates each piece globally to recover its maximal support. Regions are a proof device; they are not substituted for the true support.

If the complete static list has L codewords, every nonanchor scalar has exactly L qualifying codewords. Distinct static codewords remain distinct under the transport bijection. The complete symbolic node count is therefore `1+(q-1)*L`; the qualifying-scalar set is the singleton anchor if L=0 and the whole field otherwise. In particular, zero direction and an entire pencil inside the code are handled correctly. A whole-field symbolic family need not have the same evaluation word at different scalars, and the node identity includes its scalar.

`verify_certificate` reconstructs the static partition and transport from their explicit evidence, binds the field and problem data, and compares all exported mathematical fields. This is a shared-helper readback; the independent tiny oracle below is a separate end-to-end check. The event is **plain agreement on at least s coordinates**. These outputs do not certify a stronger prize event with additional nonjointness conditions.

## Independent execution

The new review uses separately implemented prime/GF9 arithmetic from the frozen support-review oracle and enumerates all degree-less-than-k coefficient vectors and all field scalars. It does not use the candidate's anchor finder, root bound or transport to derive expected nodes.

- **1,700** complete polynomial/scalar comparisons; **4,400** exact scalar/codeword/support nodes.
- **729** exhaustive length-three F3 pencil checks of anchor discovery, covering unique, absent and entire-code intersections.
- **108** complete static-list checks and **54** strict-cap unavailability checks.
- **904** whole-field cases, including **486** cases with the entire pencil in the code; **40** independently evaluated GF9 pencils.
- Five altered exports rejected: missing fiber, incorrect excluded scalar, false node count, duplicated partition coordinate, and false support.
- Six saved certificate readbacks, including the n1024 cases with one node and131,073 symbolic nodes. These large readbacks use the candidate verifier and its polynomial evidence; they are not independent enumeration of the enormous codebook.

The frozen engine SHA256 is `42792a1315c9aafdd34de92bd52258fd62cad5495d3584b78bf2d2ac047a5353`. No correctness issue was found in that source. Run:

```
python3 tooling_lab/round4/novelty/review_scalar_fibers.py
```

## Contribution and limits

This prototype does address the preceding interleaved fixture's specific failure: many interpolation tracks meet at one codeword, so processing tracks individually repeats information. The anchor transport replaces that work by one static direction list and a piece-partition certificate. The algebra and root-count bound are elementary known tools. The supported local contribution is their explicit certificate interface, complete support-preserving export, and tested application to the failure fixture.

Piece polynomials and their partition are supplied witnesses. The current code verifies them; it does not discover such a partition in arbitrary input. The anchor must be a complete codeword, and a pencil may have no such anchor. A failed strict root-cap test means this certificate is unavailable, not that the true list is empty. Thus the method does not solve generic list discovery or establish a uniform complexity improvement. Historical originality of the exact implementation remains unestablished.
