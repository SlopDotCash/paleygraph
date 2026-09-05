# Fourteenth pass: an exact elliptic model for three anchors

**Progress: a different exact model of the remaining spectral problem.**
This pass gives no stronger asymptotic operator or clique bound. The
full Paley and official prize goals remain unproved. The scalar
Fourier estimate uses standard curve cohomology; independent review
of its application and the preceding research remains outstanding.

For a normalized three-vertex clique {0,1,r}, its common neighborhood
is covered by the elliptic curve y²=x(x−1)(x−r) under P↦x(2P).
Every ordinary fiber has eight points. Sixteen rational 4-torsion
points account for the four exceptional fibers.

The [written derivation](parallel14-elliptic-model-2026-09-05.md) proves
that the complete lifted character matrix has entries

\[
M_{P,Q}=f(P+Q)f(P-Q),\qquad f(x,y)=\chi(y),\quad f(O)=0.
\]

The exceptional points form an explicit coupled block. Their coupling
has square-root size, so deleting it from the full operator would
change the spectral question. The exact centered restriction to the
ordinary fibers reproduces the previous normalized localized matrix.

Four elliptic translations split the neighborhood into four character
blocks, with explicit dimensions even when the action has fixed
points. At p=257, r=62, the previous 32-coordinate outlier certificate
reduces to an 8×8 integer matrix and the vector
(−1,1,−3,3,−3,1,2,2). Its squared norm is 38 and its character
quadratic form is −406. It reproduces the same normalized Rayleigh
quotient below −1. This remains a finite outlier, not a refutation
of asymptotic convergence.

On the quotient elliptic group G/G[2], every Fourier coefficient F_η
of f has |F_η|≤√p. The proof uses a tame rank-one Kummer sheaf
with four punctures and an unramified Lang-character twist; its
compactly supported first cohomology has dimension four. The
full sum is four times F_η. The coefficient bound is uniform in
the character, including the trivial one.

The transformed matrix still has off-diagonal entries. Each is a
sum of four products F_γF_{γβ}/|G/G[2]|. Their correlations, the
punctured restriction and the centered term remain to be controlled.
Scalar square-root estimates do not supply the missing operator
norm bound. No additive-phase Weil-representation theorem has been
assumed applicable to these different multipliers.

The [exact verifier](../experiments/parallel14_elliptic_model_2026_09_05.py)
passes on 27 curves: 136,704 kernel entries, 976 inverse root triples,
5,632 periodicity checks, 108 symmetry projectors, 352 positive
leading minors certifying the finite Fourier bounds, and 4,768
cyclotomic Fourier-matrix entries. Three exact off-diagonal witnesses
prevent a mistaken diagonal interpretation. The full exceptional
block, centered compression and reduced outlier certificate are
included. See [results](../results/parallel14_elliptic_model_2026_09_05.json)
and [artifact audit](../results/parallel14_pass_audit_2026_09_05.json).

The next obligation is a sharp asymptotic norm estimate for the
explicit elliptic kernel, or an alternative long-depth correlation
argument. The prior logarithmic range and subgroup exponent are
unchanged. Uniform subgroup cancellation, exceptional primes,
classical arbitrary-set cancellation, and the precise official
Reed–Solomon prize bridge remain open. The existing three workers
were rechecked and remain stopped by account usage limits; this
pass was completed locally without a reset or purchase.
