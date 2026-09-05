# Known bounds and the loss in one amplification argument

**Status: the uniform subgroup target and the Paley conjecture remain open.**
This note supplies an applicable literature baseline, checks the exponents
in a proposed extension of its method, and separates a centered-moment
formula from an issue in its broader function statement. No new subgroup
cancellation theorem or current-best claim is made.

Write `n=|H|`, `M=max_{b≠0}|η_b|`, `Δ=M/n`, and let `E_j` count ordered
equalities between two j-term sums from H. Work in the specified window
`n⁴/4≤p≤n⁴`, with dyadic n≥4. This is a concrete working window, not an
independently verified reformulation of the prize.

## An unconditional baseline

[Di Benedetto et al., arXiv:2003.06165v1](https://arxiv.org/html/2003.06165),
Theorem 3.1, gives `M≲n^(2689/2880)p^(1/72)` for `p^(1/4)<n<p^(1/2)`;
the symbol allows `p^ε` for every ε>0. Its §5 uses three dyadic selections,
the energy inputs `E₂≪n^(49/20)(log n)^(1/5)` and `E₃≪n⁴log n`, and a
weighted trilinear estimate. The quantities checked in that proof yield
`p≳n^(191/40)Δ^72` in its nontrivial branch. These are existing results.

The strict size hypotheses hold throughout our window: a prime is not the
fourth power n⁴, and `p≥n⁴/4>n²` for n>2. Substitution gives, for every
η>0, a constant C_η such that uniformly over these subgroups

\[
 M\le C_\eta n^{2849/2880+\eta}.
 \tag{B}
\]

Indeed choose the source's ε=η/4 and use p≤n⁴. The exact exponent is
`2689/2880+4/72=2849/2880`. This establishes a saving of `31/2880` from
the trivial exponent one, while the target is `1/2+o(1)`.
The implied constants do not give a usable finite bound without further
work; (B) must not replace any existing exact finite certificate.

## What improving the energy inputs can supply in this ledger

Here is an **abstract conditional calculation**, motivated by that proof.
For fixed positive integers r,s,ℓ, suppose dyadic selections with parameters
Δ₁,Δ₂,Δ₃ yield all the following inequalities, ignoring subpolynomial
factors:

\[
\begin{aligned}
 |X|&\gtrsim n^{2r-e_r}\Delta^{2r}/\Delta_1^2,\\
 |Y|&\gtrsim n^{2s-e_s}\Delta_1^{2s}/\Delta_2^2,\\
 |Z|&\gtrsim n^{4\ell-e_{2\ell}}\Delta_2^{4\ell}/\Delta_3^2,\\
 \Delta_1&\gtrsim\Delta^r,\quad
 \Delta_2\gtrsim\Delta_1^s,\quad
 \Delta_3\gtrsim\Delta_2^{2\ell},\\
 |X||Y||Z|^{1/2}\Delta_3^4&\lesssim p.
\end{aligned}
\tag{L}
\]

The e_j are exponents in full energy bounds `E_j≤n^(e_j+o(1))`.
We do **not** assert that these selections exist for arbitrary r,s,ℓ.
In particular their zero-sum deletions, thresholds, and weights require
proof. Centered energies cannot simply replace full energies in a
Cauchy–Schwarz bound for the support of a convolution.

Multiplying the first three lines and the factor Δ₃⁴ produces

\[
 p\gtrsim n^A
 \Delta^{2r}\Delta_1^{2s-2}\Delta_2^{2\ell-2}\Delta_3^3,
 \qquad
 A=2r+2s+2\ell-e_r-e_s-e_{2\ell}/2.
\]

All powers of the auxiliary Δ's are nonnegative. Successive substitutions
give `p≳n^A Δ^(8ℓrs)`. Thus the saving supplied by (L) is

\[
 B(r,s,\ell)=\frac{A-4}{8\ell rs}.
\tag{E}
\]

The actual source parameters r=s=3, ℓ=1, e₃=4, e₂=49/20 give
`A=191/40`, recovering `B=31/2880`. Keeping those selection orders but
granting the optimal energy powers e₃=3,e₂=2 only gives `B=1/24`.

There is also an exact bound on this formula even if all full energies
are granted their smallest possible powers. Diagonal solutions and
Cauchy–Schwarz give `E_j≥n^j` and `E_j≥n^(2j)/p`. Consequently in this
window the optimistic power is at least `max(j,2j−4)`. Inserting those
values into (E) gives

\[
 B(r,s,\ell)\le
 \frac{\min(r,4)+\min(s,4)+\min(\ell,2)-4}{8\ell rs}
 \le\frac1{16}.
\tag{F}
\]

For a positive numerator, replacing r,s by min(r,4),min(s,4) and ℓ by
min(ℓ,2) preserves the numerator and decreases the denominator. Thus
checking the 32 remaining triples is exhaustive. The maximum 1/16
occurs at (1,4,1) and (4,1,1). If r,s≥2, the maximum is 3/64, at
(2,4,1) and (4,2,1). Nonpositive numerators cannot exceed either maximum.
The attached exact rational audit checks every reduced case.

This is a limitation of the exponent **supplied by (L)**: even its
optimistic best amplitude exponent is 15/16. It is not a lower bound on
actual M, a proof that the optimistic energy bounds hold, or a limitation
on every amplification method. Additional information discarded by (L)
could still be useful.

## A published centered formula, with a source qualification

[Shkredov, arXiv:1802.09066v2](https://arxiv.org/html/1802.09066v2),
Theorem 3, states for r=2^k

\[
 E_r-n^{2r}/p\le
 A_k n^{2r-(k+7)/2}E_2,
 \qquad A_k=2^{3k^2}(C_*\log^4p)^{k-1}.
\tag{C}
\]

This is a relevant centered formula. However, the general-function
Theorem 25 as rendered in the inspected v2 HTML permits a mass at zero,
while equation (57) uses `||f||₂²≤||f||₁²/|H|`. We do not use that
unqualified statement as an independently checked proof dependency.
The following is our own test of those displayed hypotheses, not a
claim to refute the subgroup statement (C).

Take H=F_p* and `f=1_{0}−1/p`. It is invariant and has zero mean, and

\[
 \|f\|_1=2(1-1/p),\qquad
 \|f\|_2^2=1-1/p,\qquad f*f=f.
\]

Hence every positive-order convolution energy `T_j(f)=1−1/p`.
At k=2, the displayed general-function bound would require

\[
 1\le65536 C_*(\log p)^4(1-1/p)^4/\sqrt{p-1}.
\]

For fixed C_* the right side tends to zero, a contradiction. A support
restriction or a separate origin term is needed in that general statement.
The script independently checks the convolution idempotence and norms
using rational arithmetic in four small fields. This test does not show
that (C) fails, nor establish a repaired version of the whole proof.

Even **granting (C)**, its quantitative implication does not approach our
goal. Coset repetition gives

\[
 n M^{2r}\le p(E_r-n^{2r}/p).
\]

For a fixed k and `E₂≤n^(e₂+o(1))`, (C) would give

\[
 M\le n^{1-b_k+o(1)},\qquad
 b_k=\frac{k+1-2e_2}{4\,2^k}.
\tag{G}
\]

For e₂=49/20, the maximum is `11/1280` at k=5, smaller than the
verified baseline saving `31/2880`. Even granting e₂=2 only gives a
maximum `1/64`, at k=4 and k=5. These are global maxima: with
`a=2e₂−1`, the difference `b_(k+1)−b_k=(a+1−k)/(8·2^k)` changes sign
at the stated indices.

Taking r on the order of log n instead makes this expression's saving
`O(log log n/log n)`, not a fixed power near 1/2. The explicit factor
`A_k^(1/(2r))` has logarithm `O((log log n)²/log n)` there; hiding that
factor cannot improve the conclusion. Thus this direct use of (C)
would still not supply the required logarithmic-depth estimate (SG).

## Other scope checks

[Kowalski–Untrau, arXiv:2505.22059](https://arxiv.org/html/2505.22059),
Theorem 3.8, treats growing **prime** subgroup order d with
`d=o(log p/log log p)`, proving a Gaussian distribution over frequencies.
Those hypotheses exclude the dyadic quartic window. Distributional
convergence also does not itself bound the maximum. The inspected HTML
has a v2 arXiv header dated July 28, 2025 and an internal date of August
24, 2026; those dates are not treated as interchangeable.

For the separate clique formulation, a uniform `O(log p)` bound is
already impossible: the Graham–Ringrose lower bound is
`ω(P_p)≥c log p log log log p` infinitely often, as recorded in
[Magsino–Mixon–Parshall, introduction](https://arxiv.org/html/1907.05971).
The stronger candidate relative to a subpolynomial bound is
**polylogarithmic**. The original Graham–Ringrose proof was not audited
here; this statement is attributed through that primary research paper.

## Verification and remaining obligation

Run `python3 experiments/analytic_bound_ledger.py`. It checks the exponent
elimination, the complete reduced optimization, the centered-formula
optimization, and the exact origin example; results are saved in
`results/analytic_bound_ledger.json`. All four inspected HTML sources are
archived with hashes under `sources/analytic-bounds-2026-09-04/`.
These checks certify algebra and provenance, not all cited proofs, a
new analytic bound, or Lean formalization.

The next arithmetic input must control centered high moments or add
information beyond the full-energy ledger above. Neither treating the
unqualified function theorem as a black box nor optimizing its weak
exponent supplies square-root cancellation. The uniform conjecture and
the full reduction to the official prize remain unproved.
