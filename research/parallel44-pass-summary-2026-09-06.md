# Pass 44: a structured pointwise case and the remaining interval case

The full Paley conjecture and Proximity Prize remain unproved. This pass
adds a restriction on actual difference multiplicities and proves the
intermediate triangle target under a specified structural hypothesis on
the large levels. It also shows why the new restriction alone cannot
improve the general functional estimate. No uniform exponent improves.

For a quotient subgroup coset D of size d, the
[two-subgroup intersection argument](parallel44-structured-levels-2026-09-06.md)
proves sum_(C in D)a(C)<<n^(2/3)d^(1/3). Thus an actual level coset
with multiplicity at least T satisfies T^3d^2<<n^2, stronger than
the general T^3d<<n^2 tail. This excludes the exact subgroup-supported
model that was sharp for the earlier functional inequality.

For three level pieces that are cosets of nested quotient subgroups,
the retained field incidence gives a weighted contribution
O(n^(17/3)d_min^(-1/3)d_mid^(-1)d_max^(-1)). In particular, if every
level above n^(23/60) can be partitioned into polylogarithmically many
such cosets from one chain, then W<<n^(17/3) up to logarithms. The
small-edge contribution is already bounded at that power using the
existing49/20 energy input. The required partition is unproved and is
not inferred from small doubling or a large correlation.

The [general-case review](parallel44-interval-and-matrix-review-2026-09-06.md)
constructs an explicit weighted interval family in prime-order cyclic
groups. It has mass n-1, Q=n^2, K=O(n^3), A2=O(n^(7/5)), and satisfies
every subgroup-coset mass cap, but retains correlation norm of order
n^(61/15). These are abstract functions, not actual subgroup difference
data; no incidence matrix is supplied. The construction shows that the
new caps alone cannot supply a power saving in the standalone functional
bound. It does not refute an improved actual triangle estimate.

Actual incidence matrices also satisfy C^2=nI-ne0e0^T+L(a), where
L(w) sums simultaneous shifts of C with weights w. A separate exact
identity gives ||L(w)||_F^2=n(n-1)(sum w)^2+(p-2n)sum w^2. Root derived
these from the existing coset-convolution algebra and direct pair counts.
They retain additional constraints not imposed on the interval model,
but no stronger spectral or triangle estimate is deduced here.

Exact checks cover 24 actual mass identities, three nested incidence
normalizations, two interval models, 21,364 matrix-square entries,
178 shifted inner products, and nine weighted Frobenius identities.
Root checked the ordinary proofs and the primary statement of Mit'kin's
lemma quoted in Shkredov1504.04522v1. The source conditions and constants
are [pinned](../results/parallel44_source_scope_2026_09_06.json).
This is root review and finite verification, not separate-agent, Lean,
or external peer review. The parallel lanes remained unavailable after
their usage-limit failure.

The preceding turn made progress through the sharper absolute
exception count and actual quartic witness. This pass makes structural
progress: it handles a specified pointwise case and identifies what
the new constraint still fails to exclude. The general W exponent86/15,
energy49/20, period71/72, and absolute prime-exception exponent5/3 remain
unchanged. No Lean build, process change or Prove2Me submission was made.
The full goal is active. The next task is to use actual matrix/arithmetic
constraints against general large-level concentration, without replacing
those levels by subgroup cosets as an unproved assumption.
