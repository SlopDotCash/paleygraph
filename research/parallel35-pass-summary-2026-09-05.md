# Pass 35: the positive estimate at one degree is sufficient

**The full Paley conjecture and the Proximity Prize remain unproved.**
The previous turn made progress through the exact opposite-pair
transform and the degree-ten principal-term lower bound. This pass
strengthens the reduction: the needed upper bound can be stated at
one degree, without assuming a hierarchy of smaller-degree estimates.

For a symmetric A subset F_p^* of size n=2N, let T_s be its
nonprincipal even Fourier moment and let

    B_s=O_(2s)-2^(2s)(N)_(2s)/p,

where O counts ordered distinct opposite-free zero-sum words. For
every s>=1, the ordinary argument establishes

    B_s >= -(256sn)^s,
    |B_s| <= 4^s[T_s+(64sn)^s],
    T_s <= 16^s B_s+2(4096sn)^s.

Thus the negative side is already bounded at Gaussian scale, and a
one-sided upper estimate B_s<=(Ksn)^s at a single degree would give
T_s<=[(16K+8192)sn]^s. For a subgroup, choosing s of order log(p/n)
would then give the desired square-root scale up to logarithms.
The positive upper estimate remains unproved, and no period exponent
or full conjecture bound improves.

The [proof](parallel35-single-degree-comparison-2026-09-05.md)
applies Ravichandran's Theorem4.4 to a shifted polynomial for each
bounded real Fourier row. It handles small sizes separately and
averages only after obtaining pointwise inequalities. The primary
theorem was inspected in HTML and visually on PDF page12.

The lower bound also gives O_10>=n^6-(90+1280^5)n^5 when n>=90 and
p<=n^4. This improves the asymptotic possible shortfall to O(1/n)
relative to the principal term; its large constant makes the earlier
pass34 finite threshold better. It gives no matching upper estimate.

An exact p1153,n8 subgroup example shows that the averaged generating
polynomial is 1153-(1+2w)^4, with two nonreal roots. Applying a
real-root theorem after averaging would therefore be invalid. The
Newton recurrence was checked as a group-algebra identity, but its
leading correlation remains unestimated at the target degree.

The [verifier](../experiments/parallel35_verify_2026_09_05.py) passes
1053 exact checks in about1.71seconds:114 bounded rational rows,
78 exact Sturm root enclosures,43 finite-field cases,32 power-sum
identities and32 Newton identities, plus the subgroup counterexample.
These are same-author finite checks; separate-author review and Lean
verification of the uniform argument remain outstanding.

The three prior agents were inspected and retain terminal usage-limit
errors. No new Lean build or hosted proof job was launched. The root
agent can continue, and the full goal remains active and unachieved.
