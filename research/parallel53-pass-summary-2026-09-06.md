# Pass 53: the ambient spectral shortcut fails asymptotically

The Paley conjecture and Proximity Prize remain unproved. Recent
collision refinements did not improve the uniform exponent, so this
pass returns to the direct two-anchor spectral operator.

The actual target uses an elliptic kernel restricted to the common
neighbors C, followed by its S_3 symmetry average. On the nontrivial
sectors, the required upper bound is (2/3+o(1))p. Enlarging the kernel
to all nonzero quadratic residues makes it diagonal in multiplicative
characters, but destroys the hoped-for saving.

Its exact eigenvalues are Re(J(psi,chi)^2). Applying the checked
Lu–Zheng–Zheng discrepancy theorem proves

    lambda_max(K_Q on inversion-odd vectors)
       >= p-32pi^2 sqrt(p).

In particular this maximum exceeds 2p/3 for every eligible prime
p>=2^20, and its ratio to p tends to one. This is an asymptotic
obstruction to the ambient bound, not a counterexample to the actual
compressed, averaged operator. The source input and proof appear in
the [note](parallel53-ambient-elliptic-2026-09-06.md).

At p89 an exact integer inversion-odd vector has eigenvalue73>2p/3;
it has nonzero entries outside C. The
[checker](../experiments/parallel53_ambient_elliptic.py) also verifies
4,096 integer convolution/correlation coefficients and3,198 entries
of the original restricted operator identity across the stated fields.
All checks pass. The
[results](../results/parallel53_ambient_elliptic_2026_09_06.json) and
[audit](../results/parallel53_pass_audit_2026_09_06.json) separate these
finite checks from the sourced asymptotic proof.

The next direct spectral task is to control the interaction of the
C-support restriction with the S_3 average. Equidistribution over all
characters alone does not control the weights of an arbitrary supported
vector. No uniform spectral, period, energy, triangle or full-goal bound
improves here. The two-anchor edge itself would still be an intermediate
result, with the broader Paley and prize obligations retained.

Root completed ordinary derivations and exact checks. No independent
agent, Lean, external review or novelty claim is made. Previously
usage-limited agents were not restarted. No Lean process was polled or
changed and no Prove2Me submission was made. The goal remains active.
