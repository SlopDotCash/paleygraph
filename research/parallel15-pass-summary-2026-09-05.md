# Fifteenth pass: higher correlations, and the information they miss

**Progress: a uniform higher-correlation theorem and an abstract obstruction.**
Neither result improves the Paley clique bound or proves the full goal.
The sharp spectral edge, classical conjecture, uniform subgroup target,
and exact official prize bridge remain unproved.

The [derivation](parallel15-elliptic-correlations-2026-09-05.md) extends
the individual elliptic Fourier bound to every product of translated
sign functions, with arbitrary character twists. If s distinct shifts
occur an odd number of times, and e other distinct shifts occur a
positive even number of times, the complete sum has absolute value
at most s√p+e. Every correction at a repeated shift's zero is explicit.
When every multiplicity is even, the answer is exact character
orthogonality minus the missing points.

The proof uses 4s disjoint tame punctures on the elliptic curve;
the compactly supported first cohomology has dimension 4s. The
four-to-one quotient gives the coefficient s. This applies the same
primary cohomological inputs used in pass fourteen and still needs
independent mathematical review. It also gives exact joint Fourier
convolutions and classifies every entry of the squared kernel.

There is a limitation to using only these estimates' square-root
orders. For every eligible sufficiently large prime, changing only
O(√p) signs can retain the even sign-kernel form, all four symmetries,
the full exceptional border, and all translated-product bounds with
a larger uniform constant. The resulting normalized operator has
norm at least (4/√7)(1−2/√p), which stays above 1.

This construction **does not preserve the sharp constants or the
ambient Paley identity**. No original four-puncture Kummer realization
is asserted for the modified function. It is an abstract obstruction,
not a Paley counterexample. Unlike the earlier projection modification,
it retains zero diagonals and ±1 entries in the quotient; unlike that
modification, it cannot preserve the full Paley projection. These are
separate examples, not one construction satisfying every constraint.
The final audit verifies the loss of the sharp bound directly: the
constant Fourier coefficients change from 3 to 387 at p=10009 and
from −1 to 2827 at p=65537. Their new squares exceed p. Thus
the examples leave the original sharp-constant approach open.

The [exact verifier](../experiments/parallel15_elliptic_correlations_2026_09_05.py)
checks 1,047 parity patterns, 162 supports with all character twists,
6,864 positive leading minors, 5,168 squared-matrix entries, and
100 cyclotomic product identities. Two planted examples retain
75,552 exceptional entries and pass 48 perturbation certificates.
At p=10009 the quotient has factors 2 and 1252; at p=65537 it has
factors 128 and 128, testing both interval and rectangle constructions.
Exact rational lower bounds for their centered normalized Rayleigh
quotients are 9899/7575 and 2287/1542, both above 1.
See [results](../results/parallel15_elliptic_correlations_2026_09_05.json)
and the [artifact audit](../results/parallel15_pass_audit_2026_09_05.json).

The next attempt must combine more of the actual arithmetic information:
the sharp constants, correlation identities and the ambient projection
constraints. Generic square-root cancellation with unspecified constants
is insufficient. This is a more precise restriction on the next argument,
not a claim that the needed combined argument has been found. The
previous logarithmic range is unchanged. The three workers were terminal
at account usage limits in the preceding live check; this pass ran
locally without restarting them or consuming a reset.
