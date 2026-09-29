# Pass 36: a lattice reformulation and a fast finite certificate

**The full Paley conjecture and the Proximity Prize remain unproved.**
The previous goal turn made progress by confirming that the long Lean
process had terminated with errors. This pass resumes the mathematical
work without adding another Lean build to the machine.

For a symmetric nonzero set of size n=2N, write V(a) for the sum of
squared centered residues of one representative from each opposite
pair, divided by p^2. The exact nonzero-frequency mean is
n(p+1)/(24p). If D is the maximum deviation from that mean and M_f
is the maximum centered real character sum, the ordinary argument gives

    12 p^2/(p^2-1) D <= M_f <= 36 p^2/(p^2+1) D.

The left constant improves from12 to24 if multiplication by2 preserves
the set. These comparisons also hold for every normalized l^r norm,
r>=1. The [proof](parallel36-shell-inversion-2026-09-05.md) centers the
Bernoulli Fourier series before applying an explicit, absolutely
convergent Dirichlet inverse. The exact zero extensions are essential.
Thus a uniform bound on these lattice-distance deviations would imply
the desired period bound without an exponent loss. That bound is open.

The distances are shortest squared norms in the p cosets of Z^N in
the dual of the cyclotomic evaluation lattice. For dyadic n=2N and
N>=2, a rank-one reduction modulo p and the cyclotomic norm imply
V(a)>=p^(-2/N). This only approaches1 in the quartic regime, while
the mean grows like N/12, so it does not give the required annulus.
The p1153,n8 example has exact norm p^3 and nearly attains this coarse
norm inequality. The calculation must use signed, ordered powers;
arbitrary positive representatives preserve V but not the polynomial norm.

Truncating the inverse toward zero yields a rigorous rational period
certificate from integer residue sums. For p6700417,n64, L4096 gives
M<43.832 and all frequency cosets other than H below39.838 in absolute
value. Direct rational cosine bounds at frequency1 identify the maximum
as eta(1). The integer certificate took about0.63seconds; the full
same-author verifier took about25.4seconds and made41574 exact assertions.
It checks4096 inverse coefficients, six general symmetric sets, twelve
full small-field implementations and28 dyadic norm cases, as well as
large-field extremizers. No numerical floating-point test is used for
acceptance, and no new Lean or hosted proof job was launched.

The project [already knew the maximizing coset and a tighter interval](period-polynomial-certificate.md)
from signed moments through degree24. The new computation independently
corroborates that conclusion with a different method; it does not improve
the best recorded enclosure of M. The initial comparison with the older
sqrt(1970) bound was incomplete and has been corrected. No benchmark
against a fresh run of the old full method was performed.

This is progress in reformulation and finite verification, not a proof
of the missing positive moment upper bound, a new uniform period
exponent, a classical Paley reduction, or the official prize result.
Separate-author review and Lean formalization of the new analytic
argument remain outstanding. The full goal stays active.
