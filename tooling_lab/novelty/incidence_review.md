# Final incidence refinement review — 2026-09-05

Reviewed `observables/incidence_refinement.py` and `results_incidence.json`. The implementation builds the correct weighted K6: edge ij is the complete quadratic-character correlation of the four elements remaining after omitting i and j. Its rooted signature is the multiset of the six sorted incident-weight lists. The triangle feature is `tr(W^3)` with zero diagonal. Both are permutation invariants; neither is a newly invented graph invariant.

The independent standard-library script [review_incidence.py](review_incidence.py), importing no discovery code, replays the remaining F61 pair:

| Quantity | C={0,1,2,4,38,52} | D={0,1,3,8,10,21} |
|---|---:|---:|
| M2 | 330 | 330 |
| M4 | 4998 | 4998 |
| Boundary moment 2 | 14 | 14 |
| Boundary moment 4 | 86 | 86 |
| Additive energy | 74 | 74 |
| M6 | 101790 | 124830 |
| tr(W^3) | -1992 | 4152 |

The entire unordered quartic deck and the rooted signatures agree. All 720 vertex permutations were checked; none identifies the two weighted graphs. The unequal triangle traces already give a shorter certificate of that non-isomorphism. Saved expected values all match. Evidence is [incidence_review.json](incidence_review.json).

The implementation's refinement logic is sound: only previously ambiguous base fibers need additional work, because splitting a constant-target fiber cannot create ambiguity. Hence the reported 63-to-1-to-0 ambiguous-fiber sequence at F61 is meaningful for the full normalized census even though the additional feature evaluation runs only on the 5,550 records in initially ambiguous fibers. The saved `after` and `after_triangle` fiber counts refer to that affected subcollection, not all 455,126 records. The F257 experiment is explicitly sampled. This reviewer did not rerun either full refinement census.

This result identifies information lost by an unordered deck and then by a one-step rooted signature. It does **not** exhibit two inputs with isomorphic fully labeled quartic structures and different M6: the remaining graphs are non-isomorphic. The triangle feature uses products of already-defined quartic values, so it requires no new kind of character sum, but its uniform mathematical control may still be as difficult as the target. Zero remaining finite ambiguity proves no all-prime determination theorem, no critical-scale bound, and no new cancellation exponent. Both demonstrated n=6 fields are below p=n^4.

Code and saved-result wording are appropriately cautious on these points. A final top-level `tooling_lab/README.md` was not yet present at the first review read; any synthesis written afterward must preserve these limits.
