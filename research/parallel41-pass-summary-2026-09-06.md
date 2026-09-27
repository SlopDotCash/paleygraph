# Pass 41: a sharper collision bound and a quantified remaining input

The full Paley conjecture and Proximity Prize remain unproved. This pass
improves a uniform weighted-triangle estimate in the subgroup route,
reduces an arithmetic recovery loss, and verifies a smoothed character
moment estimate together with its unresolved recovery cost.

For H=-H of size n<sqrt(p), the [final correlation argument](parallel41-correlation-optimization-2026-09-06.md)
proves

    W << n^(86/15)(1+log n)^(8/15),

where W sums f(x)^2 f(y)^2 f(x-y)^2 over nonzero edges and
f=r_(H-H). All repeated incidence fibers are retained. This improves
our earlier n^6 estimate, but its power remains 1/15 above 17/3.
The proof combines an exact quotient-energy identity with asymmetric
mixed-level bounds and two convergent geometric sums, removing the
initial three extra logarithms. It needs neither the invalid injectivity
shortcut nor the stronger MRSS energy bound as a premise.

The exact dependence on nontrivial shifted multiplicative energy X
provides a more precise next input. The diagonal incidence row implies
sum_C a(C)^3 <=4(n-1)+(9/2)X. Combining it with the new correlation
estimate proves that

    X << n^(61/31), up to logarithms,

would reach the triangle power 17/3. This is an unproved sufficient
condition, not a bound obtained here. It requires a saving of 1/31 from
the currently available quadratic power. A separate exact subtraction
isolates the excess with three distinct edge cosets. These reductions
identify where an additional estimate must enter; they do not establish
the full arbitrary-set Paley or prize conclusions.

The [real-subfield recovery](parallel41-shell-arithmetic-recovery-2026-09-06.md)
turns a coset norm defect tau into a kernel element with cofactor
lambda^2*tau when supplied a real certificate of norm lambda*p.
In the known order-64 example, lambda=1 replaces 1217^31 by 1217.
The previously known cofactor 641 remains smaller. Exact determinants
also verify that coefficient centering can increase the norm defect,
so the new lattice correspondence does not solve norm minimization.
The [unit review](parallel41-shell-unit-distortion-2026-09-06.md) excludes
making the full correspondence a relative 1+o(1) Euclidean isometry
by choosing real units. It leaves estimates on actual cosets open.

For the classical route, the [smoothing argument](parallel41-classical-smoothing-2026-09-06.md)
gives a uniform expected sixth moment <=(325/32)pn^3 after averaging
ceil(n^(1/5)) random translates, when n^4<=p. Exact error identities
show that this does not supply the original unsmoothed moment bound.
The missing cancellation has been located in the recovery error.

The [root review](parallel41-independent-review-2026-09-06.md) checks the
ordinary proofs and primary-source hypotheses. Independent exact checks
cover six finite fields, including positive high-incidence cases, the
shell determinant and unit certificates, and 4,167 further smoothing
shift tuples. The original agent verifiers also passed on rerun. This
is internal mathematical review and finite verification, not Lean or
external peer review. No novelty or best-current-literature claim is made.

The preceding turn is classified as a verified wait plus diagnostic
progress: it freshly inspected the still-running Lean PID and saved a
source-level performance diagnosis. This pass is mathematical progress:
the W bound and arithmetic recovery changed, and the sufficient missing
collision estimate is sharper. No new Lean build, process change, or
Prove2Me submission was performed. Energy49/20, period71/72, the positive
original-moment estimates, and the full-goal status remain unchanged.

The next mathematical task is to bound the actual nontrivial collision
excess, or exploit restrictions on its joint concentration with the
actual difference levels. The standalone correlation inequality is
sharp for an abstract subgroup-supported function, so another use of
the actual prime-field subgroup structure is needed. That abstract
witness is not a Paley counterexample. The full objective remains active.
