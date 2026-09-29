# Fifth parallel pass: a stronger subgroup exponent and all length-five words

**This pass gives an actual stronger maximum-period bound on a growing
class, and closes the remaining length-five necklace patterns. The full
Paley targets and the general Proximity Prize bridge remain unproved.**
Three agents worked on subgroup moments, necklaces, and classical
moments; the primary agent worked on direct anchor weights. Independent
reviews checked the subgroup amplification and the anchor-weight proof.
The primary agent independently audited all three agent proofs.
No claim of literature novelty or formalization is made.

## A stronger maximum-period bound on a density-one class

For N a growing power of two, retain exactly the quartic prime window

\[
 \mathcal P_N=\{p\text{ prime}:N^4/4\le p\le N^4,\ p\equiv1\pmod N\}.
\]

The [subgroup proof](parallel5-subgroup-2026-09-04.md) constructs a
class J_N⊂P_N with |J_N|/|P_N|=1−O(1/log N), on which

\[
 \boxed{M(H_N):=\max_{b\ne0}
   \left|\sum_{h\in H_N}e^{2\pi i bh/p}\right|
   \le C N^{23/24}(\log N)^{7/72}.}
\]

The exponent 23/24≈0.9583 improves the previously recorded uniform
baseline 2849/2880≈0.9892 on this density-one class. The previous
bound continues to apply to its complement. Neither bound is claimed
to be best in the literature. The requested target has exponent
1/2+o(1) and applies uniformly; that gap remains substantial.

The new arithmetic input is a sixth-energy theorem at every dyadic
level: outside a proportion O(1/w), all levels satisfy
E_3(H_s)≤15s³−45s²+40s+ws³. The principal Fourier term s⁶/p is
retained explicitly when converting this to centered moments. Taking
w=log N and intersecting with the preceding fourth-energy class
supplies E_2(H_N)≤4N² and E_3(H_N)≤(15+log N)N³.

The spectral consequence uses the fixed three-selection argument from
[Di Benedetto et al., Section 5 and Lemma 4.1](https://arxiv.org/html/2003.06165#S5),
with its energy parameters retained. The proof verifies the support
sizes after removing zero and uses admissible phase weights in the
trilinear estimate. The selections give
p≫N⁷Δ⁷²/(sqrt(A)B²(1+log N)⁵), where Δ=M/N,
E_2≤AN², and E_3≤BN³. This establishes the displayed exponent;
it is more than a substitution into the earlier conditional ledger.
The original [Petridis–Shparlinski Theorem 1.1](https://arxiv.org/html/1604.08469v4#S1.SS3)
was also checked: its ordering condition is handled by permuting the
sets and their factorwise weights, and no additional size hypothesis
is imposed on this trilinear theorem.
The prime count uses the previously checked
[Thorner–Zaman theorem](https://arxiv.org/html/2108.10878v2#S3.SS1).

The new maximum also improves higher centered energies by interpolation.
For example, on J_N,
E_4(H_N)−N⁸/p≪N^(59/12)(log N)^(43/36).
This remains above the desired N⁴ scale. The actionable next arithmetic
target is an averaged centered eighth-moment budget
Σ_(p∈P_N)(E_4−N⁸/p)log p≪N⁷(log N)^A.
It is unproved. Exact fourth energy alone is insufficient: the fixed
quartic example (p,N)=(262657,32) has intrinsic fourth energy but
additional six-term relations.

## The full three-gap necklace family

The [necklace proof](parallel5-necklace-2026-09-04.md) now treats every
word B A^(a−1) B A^(b−1) C A^(c−1), with a,b,c≥1, A={0}, B={1},
C={0,1}, and k=a+b+c:

\[
 \boxed{|N(w)|\le(3a+1)\min(b,c)p^{(k+1)/2}+5p^{k/2}+2.}
\]

An inner sum is a correlation of three hypergeometric kernels, one
pulled back by a fractional-linear map. Its quadratic inertia at
infinity excludes every global invariant, including equal-rank cases.
The proof retains the actual stalk at one, both excluded coordinate
fibers, and the constant term in the quartic correlation.

The bound is uniform in all ranks, including ranks at least p. After
label transfers it has normalized decay for k=o(p^(1/4)), including
logarithmic length. Together with earlier families it proves

\[
 \boxed{|N(w)|\le15p^3\quad\text{for every one of the 243
 degree-two words of length five}.}
\]

At length six the proved catalog covers 639 of 729 words. The remaining
90 have multiplicities (2,2,2). This counts individual estimates at
fixed lengths; it does not measure a fraction of Paley proved. The full
signed spectral aggregate is still unbounded.

## The translated-anchor estimate is now proved on the kernel span

The [direct anchor proof](parallel5-anchor-aggregate-2026-09-04.md)
removes the particular Mellin L1 loss left open in the previous pass.
For a nonempty set T of m distinct nonzero anchors, every multiplicative
character ρ, and ranks r,s≥1, it proves

\[
 \left|\sum_{t\ne0}\rho(t)
   \prod_{a\in T}\chi(t-a)k_r(t)k_s(t)\right|
 \le d_T(r,s)p^{(r+s-1)/2},
\]

where d_T=mrs if 1∈T, and d_T=mrs+r+s−1 otherwise. A direct
Kummer twist rules out invariants, so the weight need not be expanded
in Mellin characters before applying the trace bound.

This gives a simultaneous quadratic bound for arbitrary complex
linear combinations of the normalized kernels h_r=k_r/p^((r−1)/2).
For R=O(log p), their R-dimensional span has asymptotically the expected
L2 mass in every prescribed adjacency-sign cell with at most
(1/2−ε)log₂p anchors. The proof removes the half-valued contributions
at the anchors exactly. It is a bound on this kernel span, not on all
vectors supported on the cell. No approximation theorem for extremal
Paley eigenvectors by this span has been proved.

The necklace and anchor proofs import the hypergeometric and curve
cohomology inputs in [Katz, Section 2](https://web.math.princeton.edu/~nmk/g2hyper62finalcorrected.pdf)
and [Katz-GKM](https://web.math.princeton.edu/~nmk/Katz-GKM.pdf).
The relevant primary pages are archived and visually inspected; the
new calculations do not infer those deep theorems from finite tests.

## The classical moment lane narrows a sufficient test and excludes a bridge

The [classical proof](parallel5-classical-2026-09-04.md) shows that
moment estimates at one cardinality per rank suffice for full Paley.
For an unbounded sequence of fixed r, it is enough to establish
M_(2r)(B)≤C_r p^(1+β_r)|B|^r only when
|B|=floor(p^(1/(r+1))), with β_r→0. The exact subset-sampling
inequality transfers such bounds to larger and unequal rectangles.
These single-size estimates remain unproved; no equivalence or
converse is asserted.

An unconditional construction shows that, for every sufficiently large
prime, polynomial-size sets may have exactly minimal additive energies
through h<r while M_(2r)/(p|B|^r) diverges. The construction uses
power-sum digits and logarithmically many fully biased character rows.
It rules out the proposed implication from minimal lower additive
energies to a Gaussian higher character moment. It does not refute
the moment hypothesis with its necessary logarithmic term, or Paley.

## Verification and remaining work

All four lanes have exact companion verifiers and recorded results:

- Subgroup: [script](../experiments/parallel5_subgroup_2026_09_04.py),
  [results](../results/parallel5_subgroup_2026_09_04.json); 12 fixed
  prime cases and 47 moment levels, six/eight-term norm certificates,
  rational selection checks, and exponent/deletion audits.
- Necklace: [script](../experiments/parallel5_necklace_2026_09_04.py),
  [results](../results/parallel5_necklace_2026_09_04.json); inner sums,
  1580 inner checks, 472 full boundary/master identities and bounds,
  1215 length-five necklace checks, and nine literal necklaces.
- Anchor aggregate: [script](../experiments/parallel5_anchor_aggregate_2026_09_04.py),
  [results](../results/parallel5_anchor_aggregate_2026_09_04.json);
  8848 kernel entries, 294 all-Mellin operator certificates, 648
  quadratic aggregate checks and 1512 exact adjacency-cell checks.
  Equality cases retain 90 zero pivots in the semidefinite certificates.
- Classical: [script](../experiments/parallel5_classical_2026_09_04.py),
  [results](../results/parallel5_classical_2026_09_04.json); 12384
  subset-average checks, 92410 rectangle checks, eight structured
  witnesses, 38 energy checks, 24 moment checks and 99 exponent audits.

Finite checks supplement the ordinary proofs and source audits; they
do not prove the asymptotics or the imported analytic theorems. The
[integration audit](../results/parallel5_pass_audit_2026_09_04.json)
records input hashes, syntax, links, and coverage. No new Lean proof
or external submission is claimed.

The next open steps are centered eighth moments on the subgroup side,
the balanced (2,2,2) necklace patterns and full signed aggregate on the
clique side, and the single-size upper moment estimates for arbitrary
sets. The full two-set conjecture, uniform thin-subgroup target, and
the exact reduction to the official prize remain separate unresolved
requirements of the investigation.
