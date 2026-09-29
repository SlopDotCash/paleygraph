# Status after the first requested parallel pass

We have proved additional intermediate results and narrowed some failed
routes. We do not have a proof of either Paley formulation, or a verified
reduction to the full Proximity Prize. There is no justified percentage
complete or estimate of the remaining time. Finite-length results are not
fractions of a proof of an asymptotic conjecture.

At the user's request, three independent agents worked on subgroup
moments, localized necklaces, and classical arbitrary-set moments. The
primary agent worked on the prize connection and reviewed their arguments.
Each lane used separate proof, experiment, and result files. This was one
bounded parallel pass; it did not create a scheduled background task.

| Lane | Result from this pass | What remains missing |
|---|---|---|
| [Localized necklaces](parallel-necklace-2026-09-04.md) | Inversion transfers known bounds among all three choices of a two-label alphabet. This expands the length-six result from 64 to 189 words. A wheel identity resolves the first three-label word and completes all degree-two words through length three. | General longer words, including words using all three labels, and the uniform estimates needed to control extreme eigenvalues. |
| [Thin subgroups](parallel-subgroup-2026-09-04.md) | Projection onto the index-two symmetry gives a sharper mixed-energy inequality and improves a conditional induction constant from 22 to 7+4 sqrt(3). A centered version applies to convolution counts at arbitrary depth. | The required uniform upper bound on balanced mixed energy and logarithmic-depth moments. The improved constant does not change the missing exponent. |
| [Classical Paley moments](parallel-classical-2026-09-04.md) | Biased row patterns strengthen the necessary logarithmic allowance in any proposed moment bound. The construction is quantified for every sequence with p/n tending to infinity. | An upper bound for worst-case sets. This is a stronger obstruction, not a proof or refutation of the candidate logarithmic-moment hypothesis. |
| [Reed–Solomon connection](parallel-prize-subfield-2026-09-04.md) | For a base-field evaluation domain and base-field received words, the MCA-bad challenge set descends exactly to the base field. A broader subspace ratio-set restriction is proved. | Worst-case extension-valued words and a quantitative Paley-to-prize implication. The protocol winning set remains distinct from the MCA-bad set. |

The most direct extension of previous proved cases is the necklace
transfer. The other results improve the precision of the research problem
but do not close its main analytic gap. No new result here is claimed to
be novel in the literature or formally certified in Lean.

The earlier [prize winning-set suffix](list-to-winning-set.md) remains an
ordinary mathematical result for the pinned benchmark, beginning at
122641/262144. This pass does not sharpen that radius or establish its
optimality. It explains why a family can have a large protocol winning
set while having no MCA-bad challenges.

The next substantive targets are the balanced subgroup moment upper
bound, general mixed-necklace cancellation, and an arbitrary-small-set
Paley estimate. More numerical checks of already established identities
would not discharge any of those obligations.
