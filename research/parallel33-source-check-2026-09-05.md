# Scope of the recent higher-energy inputs

Two recent sources were checked while looking for a stronger six-term
upper bound. Neither supplies the missing estimate through the direct
applications below. This is not a claim that all uses of their methods
have been exhausted.

## Dense higher-energy uniformity

[Shkredov, *Some new results on the higher energies*, JNT 281 (2026)](https://lims.ac.uk/documents/some-new-results-on-the-higher-energies.pdf)
defines relative uniformity in Definition 7, equation (18), printed
page115. The published normalization was visually checked. Its
Theorem 2 obtains uniformity on an affine subspace; it does not assert
uniformity of the original set. The published PDF and the older arXiv
v2 text are archived separately.

Here is an elementary obstruction to imposing that dense normalization
on the original sparse subgroup. Write delta=n/p and f=1_H-delta.
At k=l=2, the full diagonal x_1=x_2 contributes

    p (sum_z f(z)^2)^2 = p n^2(1-delta)^2

to the source's unnormalized energy. Relative epsilon-uniformity would
require this to be at most epsilon^4 delta^4 p^4=epsilon^4 n^4.
Therefore

    epsilon^4 >= (p-n)^2/(p n^2).

In the quartic window with n>=4, the right side is at least
n^2/4-2/n>=7/2>1. Thus epsilon<=1 is impossible for the original H,
even at the lowest pair of even indices. In ambient dimension one,
the affine-subspace conclusion is allowed to be a singleton, which
does not imply the required global upper bound.

This calculation motivates keeping the sparse collision terms and
their exact principal contributions. It does not identify the source's
complete-bipartite correlation energies with our sum-equality moments.
The centered distinct-coordinate reduction in the main note has its
own elementary proof.

## Fourier coefficients and small doubling

[Shkredov, *On Fourier coefficients of sets with small doubling*, Theorem 1](https://arxiv.org/html/2412.11368v1)
gives an alternative between a large Fourier coefficient and concentration
on a regular Bohr set. For H, writing K=|H-H|/n and delta=n/p, its
condition 100K^2 delta<=1 follows from K<=n and p>=n^4/4 once n>=400.
So the issue is the conclusion, rather than failure of this size condition.

The large-Fourier alternative cannot provide the desired upper bound.
The Bohr-set alternative would need another estimate on subgroup
intersections with those sets. Its stated dimension and size guarantees
alone do not supply that estimate. No new bound is imported from this
alternative in the present pass.

## Imported estimates and verification status

The improved repeated-word estimate uses the same MRSS E_2 and E_3
inputs as earlier passes, under their unchanged size hypotheses. The
published 2026 paper and the Fourier-small-doubling paper are scope
checks only. The [source ledger](../results/parallel33_source_scope_2026_09_05.json)
pins the archived files and identifies the uses. No literature novelty,
independent-author certification, full Paley proof, or prize claim follows.
