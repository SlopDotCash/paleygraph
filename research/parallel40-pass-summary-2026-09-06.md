# Pass 40: exact limits of three proposed proof steps

The full Paley conjecture and Proximity Prize remain unproved. This pass
proves an algebraic converse with an explicit loss, refutes two proposed
shortcuts on actual quartic-window subgroups, and checks a structural
obstruction for the classical arbitrary-set route. No uniform energy or
period exponent improves.

The [shell lane](parallel40-shell-converse-2026-09-06.md) shows that an
actual short coset need not have a principal evaluation kernel. The
existing p=6700417,n=64,a=1 example has V=6120237/6700417<1, yet the
kernel is nonprincipal. The proof transfers a hypothetical generator
through Norm(2-X)=641*p and contradicts an elementary Gauss bound at641.
For a general coset polynomial F, the reciprocal T=tau*p/F is integral
when Norm(F)=p^(N-1)*tau, and has norm p*tau^(N-1). Here tau=1217,
and the cofactor1217^31 is the exact minimum within that scalar family.
This is not a minimum across all kernel elements; norm641*p is already
available by a different construction. A useful thin-annulus converse
would need to control this loss or use additional arithmetic structure.

The [incidence lane](parallel40-distinct-edge-input-2026-09-06.md)
finds positive repeated fibers on the actual old dyadic levels at
p=215535361,n=128, even after all equal edge cosets are excluded.
The correct incidence52736 becomes51968 if multiplicity is discarded.
An exact identity expresses the normalization collisions through
multiplicative intersections of the actual levels, with a valid upper
correction. That correction still needs a uniform bound; distinctness
and the use of difference levels do not make it vanish.

The [classical lane](parallel40-classical-structure-2026-09-06.md)
proves exact costs for translation envelopes, subsets, containers and
relative approximation of Sidon character convolutions. These concern
the sparse Sidon tests in the existing arbitrary-set reduction. The
relative L2 obstruction does not invalidate standard Croot-Sisask
tolerance: all shifts can satisfy its larger input-norm error bound.
The positive moment estimate for every Sidon test set remains missing.

The [root review](parallel40-independent-review-2026-09-06.md)
checks the proofs and independently reproduces the finite evidence,
including a different determinant algorithm and direct field incidence
enumeration. Its [results](../results/parallel40_independent_verification_2026_09_06.json)
also cover538 Sidon sets,5916 translations and42628 smoothing/set-pair
checks. The review is by agents in this task, not external peer review
or Lean certification. No Lean build or Prove2Me submission was launched.

The preceding turn is classified as progress: pass39 produced checked
mathematical results and a live compiler diagnosis. This pass is also
progress because the two concrete counterexamples change the available
next steps. The full objective stays active. Next work must estimate
the actual weighted collision correction, control the arithmetic loss
for shell bounds, or prove the character-sensitive Sidon moment input;
the rejected shortcuts cannot be assumed in those arguments.
