# Independent review of compressed agreement supports

Review date: 2026-09-05. Reviewed implementation SHA256:
`c92ad219679c17d704425d583a310729c3cc63f0ebfcc14418cd77c56fdda88d`.
The review covers the final prefix/suffix interpolation implementation and
the corresponding README. No actionable correctness error was found.

## Executable independent checks

Run `python3 tooling_lab/round2/novelty/review_compressed.py` from the Paley
workspace. Results are in `compressed_review.json`. The expected decoder
enumerates every scalar and every degree-less-than-k coefficient vector;
it does not use the candidate interpolation or track routines. GF(9)
arithmetic is independently implemented as F3[a]/(a^2+1).

- All 495 cover plans for n <= 9, every admissible k and s, and all three
  modes were checked against all 12,291 size-s agreement supports.
- All 729 pencils over F3 with n=3, k=1, s=2 were checked, supplemented
  by 98 seeded higher-dimension pencils over F5 and GF(9).
- Across partition, anchor, and full-k covers, with symbolic and
  materialized whole-field branches, all 4,962 complete node-set
  comparisons passed. This includes 756 symbolic whole-field cases.
- Node equality includes the scalar, entire codeword, and full agreement
  support. Symbolic tracks are expanded independently for comparison;
  ignoring them would be an incorrect interpretation of the output.

## Why completeness holds

For disjoint blocks, a support that contains no selected k-subset can
contain at most k-1 coordinates in each block (or the block size, if
smaller), and at most the number of outside coordinates. `verify_cover`
checks that this capacity is below s. Thus every qualifying support
contains a selected base. Every selected block's k-subsets are enumerated.

At such a base, uniqueness of degree-less-than-k interpolation forces a
qualifying polynomial at scalar z to equal A+zB, where A and B interpolate
the two pencil endpoints. At each coordinate, the affine residual either
vanishes identically, never vanishes, or vanishes at one explicitly
computed scalar. The common coordinates and scalar buckets therefore
give every qualifying node and its maximal agreement support. Tracks may
be deduplicated by their full pair (A,B), and nodes by (z,codeword), without
losing output. Different qualifying codewords at one scalar remain
distinct.

If a track has at least s common coordinates, it contributes every field
scalar. For large fields, the complete output is the union of the explicit
nodes and symbolic all-field tracks. The `nodes` array alone then is not a
complete ledger, correctly signaled by `node_ledger_materialized=false`.
The scope is finite affine scalars; no projective point at infinity is
claimed.

## Evidence limits and novelty

The checked capacity certificate proves universal support coverage. The
full completeness argument also uses the implementation's exhaustive
enumeration and interpolation logic; the certificate alone is not a
succinct certificate of all negative track results. The large stored
censuses were not independently rerun by this review. Their alternate
rotated cover is useful agreement evidence, not an independent decoder.

The final README correctly distinguishes a complete census for an
individual selected input pencil from an exhaustive census over all
pencils. It correctly reports that worst-case cover size can remain
exponential and that the large extension example is not the production
field. Turan covers, Lagrange interpolation, and affine parameter
classification are known ingredients. This review establishes the
correctness of the specialized local tool on the stated checks and
supports its elementary completeness argument; it establishes neither
global novelty nor a prize theorem.
