# Seventeenth pass: a prime-field gap, and why it is too small

**Progress: an explicit prime-field spectral gap and a quantified limitation
of the determinant method.** The conjecture and official prize target
remain unproved. No stronger clique bound follows from this pass.

The [derivation](parallel17-prime-uncertainty-2026-09-05.md) imports the
classical prime-order Fourier minor theorem from Tao's primary paper.
It implies that a Paley projection compressed to m≤(p−1)/2 coordinates
has no eigenvalue exactly 0 or 1. This supplies a prime-field input
that was unavailable for the square-field examples of pass sixteen.
The imported uncertainty theorem is established background.

The Paley entries give an additional exact arithmetic invariant.
If X and Y are the two compressed conjugate projections, then
p^{m+1}det(X)det(Y) is a positive integer. Combining this divisibility
with positivity gives the explicit endpoint gap

\[
\delta_{p,m}=\frac{4^{m-1}}{p^m(p-m)}.
\]

Retaining the exact integer gives a stronger determinant certificate.
However, the actual ±1 entries yield an exact second-moment formula
which shows that even this stronger certificate is at most
(1−m/p)exp(−m(m−1)/p)/2. At the fixed-anchor neighborhood sizes
m∼p/2^a, that tends to zero exponentially. The desired spectral edge
for a≥2 needs a positive constant gap. This particular determinant
argument therefore cannot establish it.

That exponential bound concerns the certificate, **not the actual
Paley eigenvalue gap**. It does not show that the conjectured edge
fails. Nor does it exclude arguments retaining more information about
the minors or their ratios.

A separate elementary interval-frequency example over prime cyclic
groups shows why Fourier nonvanishing alone cannot imply quantitative
separation: a binomial vector has exponentially small leakage outside
half the frequencies. These examples preserve real circulant projection
structure and the ambient square identity. Their transformed entries
fail the Paley ±1 condition for sufficiently large primes, and their
coordinate interval is not asserted to be an anchor neighborhood.
They are not Paley counterexamples or examples satisfying every
hypothesis from preceding passes.

The [verifier](../experiments/parallel17_prime_uncertainty_2026_09_05.py)
passed 4,132 actual Paley compressions. It checked 41,420 positive
leading minors across both real embeddings, 4,132 independent integer
determinant identities, 4,132 exact second-moment and rational decay
bounds, and 518 minors after subtracting the rational gap. Two cases
beyond the rank threshold were correctly singular. Cyclotomic checks
cover the complete convolution identities of three interval projections;
four binomial examples certify leakage below one millionth.
See [results](../results/parallel17_prime_uncertainty_2026_09_05.json)
and the [artifact audit](../results/parallel17_pass_audit_2026_09_05.json).

For a concrete comparison, at p=257 with anchors {0,1,62}, the
32-coordinate neighborhood has a determinant-certificate gap below
0.01; the desired three-anchor edge needs a gap greater than 1/6.
Both comparisons are checked exactly. The uniform argument, rather
than this finite example, establishes the asymptotic limitation.

Root reviewed complete pages 2–5 of Tao's arXiv v6 PDF and its
versioned HTML. The new deductions have not received independent
mathematical review or formal verification. Five pass-sixteen artifacts
are preserved; the source manifest adds two entries and preserves
all preceding entries. The three workers were terminal at usage limits
in the last live check, so this pass ran locally.

The next estimate must quantitatively control quadratic-residue Fourier
mass on the actual anchor-defined neighborhoods, or provide another
route to the required long-depth bound. Qualitative support uncertainty
has now been made quantitative, but not strongly enough. The earlier
logarithmic depth range, classical arbitrary-set cancellation, uniform
subgroup target and exceptional primes, and exact official prize bridge
remain open. The full goal stays active.
