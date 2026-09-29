# Pass46: a row restriction and its precise limit

The full Paley conjecture and Proximity Prize remain unproved. This pass
turns a matrix identity into a new necessary inequality for the quotient
autocorrelation. It excludes one earlier abstract obstruction, but a
modified obstruction survives the total-excess version of the inequality.
No uniform exponent improves.

For any symmetric nonnegative integer matrix C with row sums at most n,
put h_i=sum_j C_ij(C_ij-1). For a subset D of d rows, let H_D=sum_D h_i,
R_D=max_(outside D)h_i, and X_D=sum_(D x D)(C_ij)_3. The
[row-capacity argument](parallel46-row-capacity-2026-09-06.md) proves

    H_D<=nd max(1,sqrt(R_D))+6^(1/3)d^(2/3)X_D^(2/3).

For actual subgroup data h_i is exactly the autocorrelation of a, with
the exceptional zero-row correction. This links a large row to the
repeated incidences available in the columns it touches.

Applied to the pass44 interval with Sidon filler, the inequality forces
X>=c n^(11/5), contradicting the available O(n^2 log n) source bound
as the parameter grows. This is an exclusion of that abstract family as
actual data at sufficiently large size, not a Paley counterexample.

Replacing the Sidon filler by a long interval preserves mass, parity,
the compared coset caps, Q=n^2, K=O(n^3), and correlation power61/15.
The modified family satisfies every total-X row-capacity inequality with
a formal budget O(n^2). No local X_D allocation or actual C is supplied.
Consequently the new total-excess condition alone does not improve the
standalone correlation exponent. The local allocations and off-diagonal
multiplication identities are still available for the next step.

The checker verifies6,656 subset inequalities on five actual matrices,
four abstract examples, two symbolic large-parameter calculations, and
399 conservative threshold checks covering every subset for the two
modified finite examples. All passed. The one positive lower certificate
in the actual-matrix checks is in a dense example used only for the general
matrix lemma, outside the range of the sparse source theorem.

Root completed the ordinary derivations and exact checks. Live agent
status still reports all three earlier lanes errored at their usage limit;
they were not restarted. No separate-agent, Lean, or external review is
claimed. No Lean process was launched, polled, or changed, and nothing was
submitted to Prove2Me in this pass. The current uniform W86/15,
energy49/20, period71/72, quartic absolute exception exponent5/3, and all
full-conjecture targets remain unchanged. The full goal stays active.
