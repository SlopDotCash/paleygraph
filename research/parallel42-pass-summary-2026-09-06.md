# Pass 42: an exceptional-prime budget with the pointwise gap retained

The full Paley conjecture and the Proximity Prize remain unproved. This
pass bounds the total nontrivial shifted-product collision excess across
primes at each fixed dyadic subgroup order. It also rules out a proposed
diagonal-domination shortcut and verifies the recovery cost of a concrete
affine smoothing method. None of these results improves the uniform
energy exponent 49/20 or period exponent 71/72 used in this project.

For n>=4 dyadic, put

    M=(n-1)^4-2(n-1)^2+(n-1),
    S=8n^2(n-1)^2-2n^4.

For each prime p=1 mod n let H_p be its order-n subgroup and let
X_p count nontrivial ordered multiplicative collisions in
(H_p-1) minus {0}. The [compatible norm argument](parallel42-excess-equations-2026-09-06.md)
proves the unconditional weighted budget

    sum_(p=1 mod n) X_p log p <= (M/2)log(S/M).

The product of all compatible nontrivial shifted-root differences is a
nonzero rational integer. Its p-adic valuation pays for every collision;
an exact second moment bounds its absolute value. Keeping the product
compatibility uses O(n^4) exponent tuples. This specializes an existing
cyclotomic norm method; no literature novelty is claimed.

Consequently at most O_c(n^(63/31)/log n) eligible primes in a fixed
quartic interval [c n^4,C n^4] can have X_p>n^(61/31). At all other
eligible primes the pass41 conditional estimate gives W<<n^(17/3).
This is an absolute count of exceptions. It does not assert their
proportion, exclude any particular exception, or bound every prime.
A quadratic excess at one quartic prime still fits comfortably inside
the budget. That pointwise arithmetic gap remains the next task.

Exact complete norm factorizations for n=4,8,16 find respectively
1,2,18 split exceptional primes. The checker includes 446 independent
integer determinants, 101 further split-prime checks, both archimedean
inequalities, and direct shifted-ratio counts at every exception.
Completeness at those three orders follows from factoring every relevant
norm, not from an ambient numerical search. At n=16,p=17 the valuation
strictly exceeds X_p, so it cannot be replaced by equality.

One resulting actual subgroup, p=353,n=16,H=<304>, has X_dist=X=72
and zero diagonal rich-cell count. Its twelve rich cells each have
rho=3 and T_b=144. This refutes domination of the off-diagonal count
by any finite constant times the diagonal count for the general range
n<sqrt(p). It lies outside the quartic window and is not a Paley
counterexample. The stronger finite witness prevents that general
shortcut from being mistaken for the missing arithmetic estimate.

The [classical affine calculation](parallel42-classical-affine-recovery-2026-09-06.md)
uses matched random reflection projections on the character sum and the
test set. It proves a positive expected sixth-moment bound, with exact
recovery of the original pairing in expectation. However, its explicit
Holder recovery certificate is always at least as large as direct
sixth-moment Holder, even with the exact smoothed moment. This closes
this particular proposed shortcut without excluding other affine methods.
Exact rational checks cover 2,266 reflection-center tuples in six cases.

The [primary-source audit](parallel42-excess-source-audit-2026-09-06.md)
checks the existential primitive-root quantifier, the exponential size
requirements of individual cyclotomic-number bounds, and the small-group
scope of the known shifted-energy estimate. None of the checked statements
supplies the required uniform X saving. Circularity cannot be assumed
as a premise for controlling its own collision excess.

Three parallel lanes returned their main findings, then failed at the
account usage limit. Root completed the proof notes, source application,
finite checks, and [review record](parallel42-independent-review-2026-09-06.md).
The coarse norm mechanism received a separate-agent check; the stronger
AGM constant and final certificates were completed by root afterward.
This is internal ordinary proof review and finite verification, not Lean
verification or external peer review. No new Lean build, process change,
or Prove2Me submission was performed in this pass. The full goal remains
active; the next advance needs control of individual exceptional primes,
or a genuinely stronger use of the actual incidence and difference levels.
