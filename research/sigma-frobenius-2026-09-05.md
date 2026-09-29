# The non-square tuple sum as a signed sum of Frobenius traces

**Status: a synthesis. Two identities are PROVED and checked exactly; the
framing they support is HEURISTIC and proves nothing about the conjecture.
It unifies the tuple and dual reformulations of the sigma pass and states
the one input both need in the same language.**

Notation as in [the tuple note](sigma-tuple-2026-09-05.md): `A,B ⊆ F_p`,
`0 ∉ A`, ratio set `R = {−b/a : (a,b) ∈ A×B}`, and for `D ⊆ R`

    f_D(t) = Π_{r∈D} (t − r),     W_D = Σ_{t∈F_p} χ(f_D(t)).

The tuple note proves `T_ns(A,B;k) = Σ_{D⊆R} c(D)·W_D − Corr` with explicit
coefficients `c(D)` built from the signed ratio-class sums, and that the
conjecture is equivalent to `T_ns ≤ p^{−2kδ'}(mn)^{2k}` for a suitable `k`
(with the dilate `t = 0` excluded; the difference between `Σ_{t∈F_p}` and
`Σ_{t≠0}` is a bounded boundary term absorbed in `Corr`).

## 1. Point counts (PROVED)

For every `D`, `W_D = N_aff(C_D) − p`, where `N_aff(C_D)` is the number of
`(t,y) ∈ F_p²` on the hyperelliptic curve `C_D : y² = f_D(t)`. This is the
identity `#{y : y² = u} = 1 + χ(u)`; the verifier checks it for 528 random
branch sets at `p ≤ 200`. For `|D| = d ≥ 2` distinct roots the curve has
genus `⌊(d−1)/2⌋` and `|W_D| ≤ (d−1)√p` (Weil), checked on all 454 such
cases (maximum ratio to the bound `0.949`).

So `T_ns` is a signed combination of affine point-count deviations of the
`2^{|R|}` hyperelliptic curves branched on subsets of the ratio set, with
weights `c(D)` that carry the only dilation-breaking information about the
pair `(A,B)`.

## 2. Subgroup branch loci are Fermat quotients (PROVED)

Let `H ≤ F_p^*` have order `m`, `c ≠ 0`, and `D = cH`. Then
`f_D(t) = t^m − c^m`, and with `a = c^m`,

    W_{cH} = Σ_{ψ^m = 1, ψ ≠ 1} ψ(a) χ(−a) J(ψ, χ),   J(ψ,χ) = Σ_v ψ(v)χ(1−v),

with no boundary term: the trivial-character contribution `−χ(−a)` cancels
the `t = 0` term `+χ(−a)` exactly. Proof: `#{t : t^m = u} = Σ_{ψ^m=1} ψ(u)`
for `u ≠ 0`, then substitute `u = av`. The verifier checks this for every
subgroup order `m ≤ 39` of every prime `p ≤ 120` and five shifts each
(730 cases, error below `10^{−13}` in floating-point cyclotomic arithmetic).

Consequently, for the subgroup case `B = H` of the conjecture, the Weil
terms entering `T_ns` are exactly sums of Jacobi sums `J(ψ,χ)` over the
character group of order `m = |H|`, weighted by `ψ(a)`. This is the same
object as the "Statement J" of [the dual note](sigma-dual-2026-09-05.md):
both reformulations of the subgroup case reduce to cancellation in
`Σ_ψ ψ(u) J(ψ,χ)` over a subgroup `Ψ` of characters, which has more terms
than the original sum and known moduli `√p` but unknown phases.

## 3. What the framing says and does not say (HEURISTIC)

- For a "random" ratio set the Katz–Sarnak philosophy predicts that the
  normalized traces `W_D/√p` behave like independent centered variables of
  variance about `|D|−1` as `D` varies, which is exactly the pattern-level
  random-sign scale the tuple note observes on random rectangles
  (`|T_ns|/R_D` of order one, no growth in `p`). Nothing is proved by this;
  equidistribution theorems apply to algebraic families, not to the
  combinatorial family `{C_D : D ⊆ R}` with prescribed weights.
- For structured pairs the traces are correlated through the Jacobi sums
  of §2, and the complete-rectangle identity `T_ns = (#biased dilates)·(mn)^{2k}`
  shows there is no cancellation to find; the conjecture only needs those
  cases to have `mn ≤ p^{o(1)}`, which is again the open subgroup case.
- The one input every route in this workspace now needs, in one sentence:
  a bound on `Σ_{ψ∈Ψ} ψ(u) J(ψ,χ)` for subgroups `Ψ` of the character group
  of order `p^{1−ε}`, uniform in `u`, with a power saving over `|Ψ|√p`.
  No source examined in the pass addresses it; it is equivalent to
  Bourgain's Problem 5 in Chang's survey.

## Verification

`experiments/sigma_frobenius_2026_09_05.py` (standard library, under one
second) writes `results/sigma_frobenius_2026_09_05.json`.
