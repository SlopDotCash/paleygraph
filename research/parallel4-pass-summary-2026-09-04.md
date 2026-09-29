# Fourth parallel pass: an infinite subgroup class and larger kernel families

**The investigation has advanced on proved subproblems. The full Paley
targets and the general Proximity Prize bridge remain unproved.** Three
agents continued the necklace, subgroup, and classical-moment lanes in
parallel. The primary agent proved a kernel aggregate across ranks and
integrated the results. Independent review checked its singular stalk,
arithmetic normalization, and complex coefficients, and separately checked
the subgroup proportion and its prime-counting input. No claim of
literature novelty or formal verification is made.

## The strongest new subgroup conclusion

Let N tend to infinity through powers of two, and consider exactly

\[
 \mathcal P_N=\{p\text{ prime}:N^4/4\le p\le N^4,\ p\equiv1\pmod N\}.
\]

The [subgroup proof](parallel4-subgroup-2026-09-04.md) establishes that
the lower limiting proportion of these primes with

\[
 E_2(H_s)=3s^2-3s
 \quad\text{at every dyadic level }2\le s\mid N
\]

is at least

\[
 \boxed{1-\frac{\log2}{36}=0.9807459116\ldots .}
\]

Here E_2 counts ordered additive quadruples, and the displayed value is
the intrinsic contribution from opposite pairs. On this same class the
mixed energy between the two half-cosets is exactly (s/2)^2, and the
3+1 zero-sum count is zero at every level. Thus a condition previously
checked only at finite quartic examples now holds on a proved growing
class. The 98.07% is an asymptotic lower proportion for this particular
fourth-energy property. It is neither a guarantee at each finite N nor
a measure of how much of Paley has been proved.

The proof combines the existing cyclotomic norm budget with a new
boundary-product exclusion. Outside its negligible prime support,
every extra quadruple orbit has 24N ordered members. That quantum
turns the norm budget into the stated prime count. The exclusion matters:
the order-64 example has excess 12N because of repeated entries.
The asymptotic denominator follows from
[Thorner–Zaman, Corollary 3.1 and equation (3.2)](https://arxiv.org/html/2108.10878v2#S3.SS1),
specialized to modulus N, squarefree part 2, and the full quartic
annulus. Its hypotheses and constants were checked independently.

A separate density-one theorem gives
E_2(H_s)≤3s²−3s+N^(-1/2)s² simultaneously at all levels, outside a
proportion O(N^(-1/2)). More generally, the total primitive collision
mass is at most w_N/12 outside a proportion O(1/w_N), for every
w_N tending to infinity. These conclusions leave the worst-case primes
and the required high centered moments open. Minimal fourth energy alone
does not give the target maximum-period estimate.

## Two new analytic tools

The [necklace proof](parallel4-necklace-2026-09-04.md) allows a repeated
second exceptional label. For A={0}, B={1}, C={0,1}, r,h≥1,
k=r+h+1, and s=min(r,h),

\[
 \boxed{|N(A^rB^hC)|\le(s+1)p^{(k+1)/2}+|t_s-(-1)^s|.}
\]

The remainder is at most p^(s/2)+1. The bound is uniform in the ranks,
including ranks at least p, and decays after the usual necklace
normalization for k=o(√p). It covers 30 of the 90 length-five words
with multiplicities (2,2,1), after label transfers. The other 60 remain
unresolved. Their exact two-variable quartic reduction retains both
zero-coordinate corrections; a one-variable substitution and direct
planar duality do not apply to it.

The [kernel aggregate proof](parallel4-kernel-aggregate-2026-09-04.md)
first gives uniform pairwise Mellin bounds. At equal rank it subtracts
the scalar invariant and restores the missing singular-stalk term:

\[
 \left|\sum_{t\ne0}\rho(t)
 \bigl(k_r(t)^2+p^{r-1}1_{t=1}-p^{r-1}\bigr)\right|
 \le(2r-2)p^{r-1/2}.
\]

For h_r=k_r/p^((r-1)/2), every multiplicative character ρ, and every
complex coefficient vector a of length R, this yields

\[
 \left|\frac1p\sum_{t\ne0}\rho(t)
       \left|\sum_{r=1}^R a_rh_r(t)\right|^2
       -1_{\rho=1}\|a\|_2^2\right|
 \le\left(\frac{3R^2-R-2}{2\sqrt p}+\frac2p\right)\|a\|_2^2.
\]

This controls a signed quadratic aggregate across ranks, uniformly for
R=o(p^(1/4)), including logarithmic R. It does not control the sum over
arbitrary necklace label patterns. Inserting a further translated
anchor by the current Mellin triangle inequality loses the saving.

Both analytic arguments import the hypergeometric and cohomological
theorems in [Katz's Section 2](https://web.math.princeton.edu/~nmk/g2hyper62finalcorrected.pdf)
and [Katz's monograph](https://web.math.princeton.edu/~nmk/Katz-GKM.pdf).
The relevant pages were visually inspected; an independent audit
rechecked the conjugation bars in the duality formula. The new proofs
retain actual tensor stalks rather than substituting their generally
larger middle-extension stalks.

## Classical structured sets

The [classical lane](parallel4-classical-2026-09-04.md) transfers the
established Burgess bound to the exact signed elliptic pairing. For
fixed 0<β<1/15, p^(13/30+β)≤n≤√p, and a set B differing from an
arithmetic progression by at most np^(-1/30) elements, it proves

\[
 |\langle W_B,E\rangle|\ll_\beta p^{1/2-\beta}n^4,
 \qquad |S(A,B)|\ll_\beta|A|n p^{-1/30}\quad\text{for every }A.
\]

The same transfer has a proper rank-two progression corollary using
[Alsetri–Shao, Theorem 1.1](https://arxiv.org/html/2509.07765v1#S1.SS1).
These are consequences of imported character-sum theorems for structured
sets. They do not resolve arbitrary small sets or the proposed strong
moment bound. The exponent 13/30 is optimal only within the specified
Burgess-to-fourth-moment calculation.

## What still stands between this and the goal

1. The dyadic spectral target needs centered moments at logarithmic depth
   and control for every eligible prime. The new fourth-energy class
   supplies neither condition.
2. The clique route still needs the full signed aggregate in the
   [second-pass spectral transfer](parallel2-spectral-transfer-2026-09-04.md).
   The new quadratic aggregate and individual word families are smaller
   objects. The two remaining (2,2,1) patterns are a concrete next test.
3. The arbitrary-two-set Paley estimate remains stronger than the
   conditional clique conclusion, and the reduction to the precise
   official Reed–Solomon challenges remains unverified.

## Verification

| Lane | Exact checks | Artifacts |
|---|---|---|
| Subgroup | 16 fixed-prime certificates, 68 dyadic levels, norm and boundary-product audits; includes both the 12N repeated-entry and 24N distinct-entry witnesses | [Script](../experiments/parallel4_subgroup_2026_09_04.py), [results](../results/parallel4_subgroup_2026_09_04.json) |
| Necklace | 192 block identities/bounds/reversals, 200 Mellin factorizations, 144 literal hypergeometric fibers, 32 norm certificates with 608 positive leading minors, 89 three-gap identities, 9 literal necklaces | [Script](../experiments/parallel4_necklace_2026_09_04.py), [results](../results/parallel4_necklace_2026_09_04.json) |
| Kernel aggregate | 42208 kernel entries, 83 norm certificates with 1452 positive leading minors, 162 signed aggregate cases, local stalk ranks 1 through 8 | [Script](../experiments/parallel4_kernel_aggregate_2026_09_04.py), [results](../results/parallel4_kernel_aggregate_2026_09_04.json) |
| Classical | 3869 exhaustive remainder/directional cases, 20 progression cases including 8 perturbations, 99 Burgess parameters, 4 rational exponent choices | [Script](../experiments/parallel4_classical_2026_09_04.py), [results](../results/parallel4_classical_2026_09_04.json) |

The finite tests supplement the ordinary proofs; they do not establish
the analytic asymptotics or the imported sheaf theorems. All verifier
runs and independent reviews completed. The
[integration audit](../results/parallel4_pass_audit_2026_09_04.json)
records source/input hashes, script syntax, and local links. No new
Lean proof, external submission, or claim of a completed goal is made.
