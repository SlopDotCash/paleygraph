# Twelfth pass: the short-depth estimate cannot bootstrap from these inputs

**The full Paley and prize goals remain unproved.** This pass rules out
a proposed extension of the current estimate and supplies a small exact
finite counterexample to an overly strong edge bound. It does not improve
an asymptotic clique or subgroup bound.

The [written argument](parallel12-bootstrap-obstruction-2026-09-05.md)
modifies the actual Paley projection by rank at most two. It preserves
projection rank, the constant null direction, every selected anchor column,
the common-neighborhood mask, and S²=pI−J. It plants a zero or one
localized eigenvalue. The modified matrix fails the original character-entry
condition away from the anchors; it is a countermodel to a proposed
projection argument, not a counterexample to Paley.

All fixed-word traces change by at most 4k p^(k/2). For the complete
normalized transfer matrix, the energy obeys

\[
 H'_j\le(\sqrt{H_j}+2\sqrt2\,jN^{j/2})^2,
 \qquad H'_j\ge N^j,\quad N=2^a-1.
\]

Thus the established bound H_j≤(1+o(1))p survives throughout its
short-depth range, while the required long-depth condition fails for the
modified matrices. Projection identities, exact anchor data, and those
short-depth estimates cannot alone exclude extreme outliers. Further
quadratic-character information away from the anchors is necessary for
this route. This conclusion is conditional only in its use of the prior
short-depth bound; the perturbation calculation is elementary.

Separately, for the **actual** p=257 Paley matrix and anchor clique
{0,1,62}, a 32-coordinate vector with entries in {−3,…,3} has normalized
Rayleigh quotient

\[
 -812/(19\sqrt{1799})<-1,
 \qquad812^2-19^2\cdot1799=9905>0.
\]

The certificate even gives spectral radius greater than 9/8 for its
transfer matrix and H_j>(81/64)^j at this fixed prime. An exact edge
interval or polynomial-in-depth energy bound at every prime is therefore
false. A vanishing asymptotic error remains possible; this finite example
does not refute spectral-edge convergence or the Paley conjecture.

The [exact verifier](../experiments/parallel12_bootstrap_obstruction_2026_09_05.py)
checks six modified projections, 366 preserved anchor entries, 78 word
trace perturbations, 44 power trace perturbations, and 22 each of power
norm bounds, planted power identities, and energy recurrences. It separately
checks the actual p=257 certificate using integers. See
[results](../results/parallel12_bootstrap_obstruction_2026_09_05.json) and
[artifact audit](../results/parallel12_pass_audit_2026_09_05.json).

The next obligation is a uniform long-depth estimate using the full
character structure, with asymptotically vanishing edge error. The
preceding necklace and aggregate proofs still need independent review.
The subgroup exponent remains 8/9 on the existing density-one quartic
prime class; its uniform square-root target, the arbitrary-two-set
conjecture, and the exact official prize bridge remain open. The existing
parallel workers are unavailable because of their previously observed
account usage limits; this pass was completed locally.
