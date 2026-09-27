# Pass 43: improved prime exceptions and a quartic obstruction

The full Paley conjecture and Proximity Prize remain unproved. This pass
improves the absolute exceptional-prime bound for the intermediate
triangle target W<<n^(17/3), from O_c(n^(63/31)/log n) to
O_c(n^(5/3)/log n). The bound for every prime remains open. The new
estimate is not a density statement or a claimed literature novelty.

The [energy argument](parallel43-energy-prime-exceptions-2026-09-06.md)
retains the exact positive triangle sum and gives
W<=R_max*E_*^3/n^2. The published single-coset intersection bound
R_max<<n^(2/3) therefore makes E(H)<<n^(7/3) sufficient. The earlier
fourth-energy prime-average theorem bounds the number of primes that
can fail this energy condition. Root visually checked the published
lemma and verified its cardinality hypothesis for three single cosets.
The invalid union-of-cosets application from pass38 is not used.

The [quartic witness](parallel43-quartic-collisions-2026-09-06.md) has
p=278177,n=32,H=<160164>, minimum additive energy2976, and
X_dist=X=72 with zero diagonal rich mass. It refutes every constant-factor
domination of off-diagonal by diagonal rich mass inside the actual
quartic window. All excess in the shifted-product equation uses four
distinct entries. Its twelve rich cells have rho3 and T_b144. It does
not violate the desired triangle bound, and is not a Paley counterexample.

An exact local refinement gives v_p(mathcal_P_n)=sum_r X_r, where
X_r counts collisions modulo p^r among compatibly lifted roots of unity.
This explains higher valuation rather than assuming equality with X_1.
Hashing unordered products computes each precision using O(n^2) modular
arithmetic and storage. Twenty-one comparisons with independently
factored integer norms passed, along with three larger quartic examples,
the witness's full-field triangle audit, and four energy/triangle checks.

This pass was completed by root because the parallel lanes had reached
the account usage limit. The [review record](parallel43-review-2026-09-06.md)
separates ordinary proofs, imported hypotheses, and finite checks.
No separate-agent or Lean review, external peer review, uniform energy
improvement, or full-conjecture proof is claimed. Energy49/20 and
period71/72 remain unchanged. No Lean build, process change, or Prove2Me
submission was performed; prior research artifacts are preserved.

The previous turn is classified as progress: it changed the mathematical
state and also verified the active Lean process. The present turn is
mathematical progress. The next task remains control of exceptional
individual primes or a stronger use of actual joint incidences; the
diagonal-domination shortcut is now ruled out in the target range.
