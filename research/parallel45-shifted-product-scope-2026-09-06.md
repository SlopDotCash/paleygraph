# Shifted-product inputs for the surviving level concentration

The checked results below do not eliminate the remaining concentration
pattern. This is a comparison of specified theorems with our parameters,
not an exhaustive claim about the literature or a counterexample in a
prime field.

Let D be a set of d quotient cosets with T<=a(C)<2T, let Q be its preimage
in F_p*, and suppose |DD|<=K d in the quotient. Then

    R=H intersect (Q+1),  Td<=|R|=sum_D a<2Td,
    |RR|<=n,  |(R-1)(R-1)|<=|QQ|<=Knd.                  (1)

At the functional obstruction's scale T=n^(3/5), d=n^(1/5), the lower
mass is n^(4/5). A geometric interval in the quotient has bounded K.
These are the parameters to test if such a level were realized by an
actual subgroup. The abstract construction does not prove realization.

[Warren, *On products of shifts in arbitrary fields*, Corollary4,
published page248](https://msp.org/cnt/2019/8-3/moscow-v8-n3-p04-p.pdf)
gives |AA|+|(A+1)(A+1)| >> |A|^(11/9)/(log|A|)^(7/6) when
|A|<p^(1/4). Set A=R-1. In the quartic range, the tested mass scale
n^(4/5) meets the size condition for sufficiently large n. The resulting
lower power is44/45, below even the n upper bound for |RR| in (1).
Thus this consequence is consistent with the tested concentration scale.
No best-current-exponent claim is needed for this comparison.

A precise stronger input that would extend the previous mass cap to
low-doubling D is

    |RR| |(R-1)(R-1)| >> |R|^3                          (2)

for these particular intersections R, with any stated logarithmic losses.
Combining (1)-(2) would give T^3d^2<<Kn^2. At the displayed scale its
left side is n^(11/5), so bounded or polylogarithmic K would be excluded.
**Inequality (2) is an unproved sufficient input here.** It is not inferred
from Warren's theorem, small doubling, or the matrix identities.

A recent stronger-looking source has a different field hypothesis.
[Harrison–Mudgal–Schmidt, arXiv:2603.06483v1, Theorem1.5](https://arxiv.org/html/2603.06483v1)
controls intersections of low-doubling sets with non-coset varieties in
one-dimensional algebraic groups over the complex numbers. That theorem
does not state a prime-field analogue. No transfer with uniform dependence
on p is established in this project, so it cannot be inserted into (1).

Both primary sources are archived with hashes in the
[scope record](../results/parallel45_source_scope_2026_09_06.json).
The exact rational comparisons44/45<1 and11/5-2=1/5 are checked by the
[matrix experiment](../experiments/parallel45_weighted_matrix_algebra.py).
The full goal and every existing uniform exponent remain unchanged.
