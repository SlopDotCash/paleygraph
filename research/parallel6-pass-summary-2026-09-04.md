# Sixth parallel pass: a stronger subgroup bound and a repaired kernel family

**There is measurable progress on restricted subproblems. Neither the full
Paley conjecture nor the requested thin-subgroup bound is proved. The
general reduction to the Proximity Prize remains unverified.** Three
agents worked on subgroup moments, necklaces, and classical moments;
the primary agent worked on multiple kernel seeds and reviewed all
three proofs. The subgroup and seeded-kernel arguments also passed
independent review.
No literature-novelty, optimality, or formalization claim is made.

## The subgroup exponent improves from 23/24 to 17/18

Let N grow through powers of two, and retain the precise prime window

\[
 \mathcal P_N=\{p\text{ prime}:N^4/4\le p\le N^4,\ p\equiv1\pmod N\}.
\]

The [new subgroup proof](parallel6-subgroup-2026-09-04.md) gives

\[
 \boxed{M(H_N)\le2^{1/6}(15+\log N)^{1/18}N^{17/18}}
\]

on the sixth-energy class K_N from pass five, whose relative size is
1−O(1/log N). An extra fourth-energy condition is no longer required.
The exponent 17/18≈0.9444 is smaller than the preceding 23/24≈0.9583.
The target is still N^(1/2+o(1)) **at every eligible prime**, so both
the strength and the exceptional-prime gap remain substantial.

The new argument uses subgroup invariance, Jensen, and an elementary
weighted bilinear inequality. For normalized r-term sum measures μ_r,
write u_r=μ_r(0), m_r=1−u_r, and

\[
 V_r=\frac{E_r-n^{2r}/p}{n^{2r}}
       -\frac p{p-1}(u_r-1/p)^2.
\]

For Δ=|η(a)|/n it proves, for every pair r,s≥1,

\[
 (\Delta^{rs}-u_r-u_s+u_ru_s)_+^2
 \le \frac{p-n}{n(p-1)}m_r^2m_s^2
       +(1-1/n)m_rm_s\sqrt{pV_rV_s}.
\]

This retains the origin atoms and the constant contribution on the
nonzero coordinates. At r=s=3, p≤n⁴ and E₃≤Bn³, it gives
Δ¹⁸≤8B/n. The pass-five sixth-energy theorem supplies B=15+log N.
The prime-density assertion continues to use the archived
[Thorner–Zaman theorem](https://arxiv.org/html/2108.10878v2#S3.SS1);
the new amplification needs no trilinear theorem or dyadic selections.

The resulting centered eighth energy satisfies
E₄−N⁸/p≪N^(44/9)(log N)^(10/9), improving 59/12 to 44/9 in
the exponent, but still exceeding the desired N⁴ scale. Feeding these
interpolated moments back into the same coarse inequality cannot improve
the maximum exponent further; the proof optimizes all fixed r,s.

The next arithmetic target remains

\[
 \sum_{p\in\mathcal P_N}(E_4(H_N)-N^8/p)\log p
                         \ll N^7(\log N)^A.
\]

It is unproved. Its payoff is now explicit: a centered eighth-energy
bound E₄−n⁸/p≤Dn⁴ would give M≤3^(1/4)D^(1/16)n^(15/16)
through the same proved theorem at r=1,s=4.

## One more balanced necklace class is controlled

The [necklace proof](parallel6-necklace-2026-09-04.md) proves a uniform
operator estimate for Q=S D_C S and W_C=S∘Q:

\[
 \|W_C+S\|\le2p,\qquad\|W_C\|\le2p+\sqrt p.
\]

A projective change of variables turns W_C+S into a principal
compression of an explicitly augmented single-anchor kernel. The
proof keeps the projective pole and diagonalizes all three added
coordinate modes. It reuses the previously checked
[Katz finite-field Mellin bound](https://web.math.princeton.edu/~nmk/mellin186.pdf).

Renormalizing the adjacent C vertices gives the exact identity

\[
 N(AABBCC)=f_0^T W_C S W_C f_1-p t_4+2,
\]

and hence |N(AABBCC)|≤7p^(7/2). The +2 is the retained
coincident-anchor correction. Its entire twelve-word orbit has exactly
the same value. Coverage at length six therefore increases from
639/729 to 651/729, leaving 78 balanced words in four classes.
This count concerns individual finite-length families; it gives no
percentage of the full conjecture solved.

The next two contractions have crossing or nested pairings and do not
obey the same matrix-chain estimate. Their exact identities are
recorded without claiming a bound. The full signed spectral aggregate
remains uncontrolled.

## The old kernel span misses a largest eigenvector

The [seeded-kernel proof](parallel6-seeded-kernels-2026-09-04.md) first
identifies a structural limitation: every old kernel satisfies
k_r(t⁻¹)=χ(t)k_r(t). Adding powers never leaves the corresponding
inversion-even subspace. At p=13, the common neighborhood of 0 and 1
is {4,10}; its largest restricted eigenvector is (1,−1), and every
old kernel misses it. This is an exact finite obstruction to the
universal claim that the old span captures the top eigenvector.

Allowing seeds k_(b,r)(t)=k_r(t/b) captures that omitted vector at
b=3 and gives uniform translated correlation bounds. For a nonempty
anchor set T, disjoint from the seeds b,c,

\[
 \left|\sum_{t\ne0}\rho(t)
   \prod_{a\in T}\chi(t-a)k_{b,r}(t)k_{c,s}(t)\right|
 \le(|T|rs+r+s-1_{b=c})p^{(r+s-1)/2}.
\]

The proof checks the invariant exclusion and actual singular stalks,
using the archived [Katz hypergeometric](https://web.math.princeton.edu/~nmk/g2hyper62finalcorrected.pdf)
and [curve cohomology](https://web.math.princeton.edu/~nmk/Katz-GKM.pdf)
inputs. It also treats empty anchor sets and the corrected equal-seed
diagonal separately. Colliding seeds and anchors are explicitly excluded:
already at rank one their correlation can have order p.

For L≤p^σ seeds and R=O(log p) ranks, arbitrary linear combinations
have the expected mass in every sign cell with
a≤(1/2−σ−ε)log₂p nonzero anchors, provided σ+ε<1/2.
An optional condition on χ(t) is included. This enlarges the proved
range, but still does not control all vectors supported on the cell.
The full rank-one seed matrix is invertible, which also proves that
the cell isometry cannot extend to all seeds with relative error
less than one. A full restricted-operator estimate needs another step.

## The classical moment result applies to typical sets

The [classical proof](parallel6-classical-2026-09-04.md) gives exact
formulas for the mean and variance of M₄(B) over n-element sets.
In the range n→∞, n=o(p), the variance is asymptotic to 24pn⁴.
At exactly n=floor(p^(1/3)), it proves

\[
 M_4(B)\le3pn^2
\]

for a proportion at least 1−(6+o(1))p^(−1/3) of those sets. It
also bounds the centered signed quartic pairing by 7p^(4/3) outside
O(p^(−1/3)) of the same slice. Both are upper estimates on the
stated class, with no averaging over primes.

They do not control adversarial sets. A fixed p=1009 example on the
same cardinality slice violates coefficient 3; it does not refute an
unspecified SS constant. The note explains why mean and variance alone
permit rare large values, and compares the result with the strong
cancellation already available for random sets. No worst-case Paley
threshold improves in this lane.

## Verification and the next obligations

Each lane has an exact verifier and recorded results:

- Subgroup: [script](../experiments/parallel6_subgroup_2026_09_04.py),
  [results](../results/parallel6_subgroup_2026_09_04.json).
- Necklaces: [script](../experiments/parallel6_necklace_2026_09_04.py),
  [results](../results/parallel6_necklace_2026_09_04.json).
- Seeded kernels: [script](../experiments/parallel6_seeded_kernels_2026_09_04.py),
  [results](../results/parallel6_seeded_kernels_2026_09_04.json).
- Classical moments: [script](../experiments/parallel6_classical_2026_09_04.py),
  [results](../results/parallel6_classical_2026_09_04.json).

The [integration audit](../results/parallel6_pass_audit_2026_09_04.json)
records final hashes, finite-check counts, proof reviews and local links.
Finite computations supplement the mathematical proofs; they do not
establish the imported source theorems or the unresolved asymptotics.

The next decisive steps are the averaged centered eighth-energy budget,
the remaining balanced necklace classes and full weighted aggregate,
and uniform single-size character moments for arbitrary sets. The
uniform subgroup target, full two-set conjecture, and exact prize
reduction remain separate open obligations. The goal stays active.
