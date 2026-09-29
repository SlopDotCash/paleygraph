# Stepanov wave: shared brief (2026-09-26)

Read this first, then `research/sigma-brief-2026-09-05.md` (rules: honesty
labels, exact checks, precise citations, no shared-file edits, no git, no
sub-agents), then `HANDOFF.md` and `research/sigma-pass-summary-2026-09-05.md`
§6 and §9. Your file prefix is **`stepanov`** (not `sigma`, not `parallel`).

## Why this wave

Every averaging method (Fourier, Weil moments, spectral traces, dilation
moments) is blocked by dilation invariance: its bound is the same for all
`p−1` dilates `tA`, so it cannot single out `t = 1`
(`research/sigma-crux-2026-09-05.md`). Every field-agnostic method is blocked
by the `F_{p²}` test (`research/sigma-barriers-2026-09-05.md`).

Stepanov's method, as used by Hanson–Petridis (HP, arXiv:1905.09134, local
text `sources/sigma-hanson-petridis-1905.09134.txt`, Theorem 1.2 and its
one-page proof in §2), is the only known argument that is both **local**
(it bounds a single pair `(A,B)`) and **prime-sensitive** (its conclusion
fails over `F_{p²}`, where `A = B = F_p` gives `|A||B| = p² > (p²−1)/2 + p`).
Its limitation is that it needs a **complete** biclique.

## Notation

`p` odd prime, `χ` Legendre symbol, `χ(0) = 0`, `d = (p−1)/2`, `Q` the
nonzero squares. For `A, B ⊆ F_p`, `m = |A|`, `n = |B|`,
`S(A,B) = Σ_{a,b} χ(a+b)`, `r = |B ∩ (−A)|` (number of zero sums), and for
`b ∈ B` let `e_b = #{a ∈ A : χ(a+b) = −1}` (bad partners).

**HP (Theorem 1.2, `d = (p−1)/2`).** If `A + B ⊆ Q ∪ {0}` then
`mn ≤ d + r`. Proof: with `c_k` chosen (Vandermonde) so that
`G(x) = Σ_k c_k (x+a_k)^{M−1} ≡ 1`, the polynomial
`F(x) = −1 + Σ_k c_k (x+a_k)^{d+M−1}` has degree exactly `d`, a root of order
`M` at each `b ∈ B∖(−A)` and of order `M−1` at each `b ∈ B ∩ (−A)`.
At a point `b` with bad set `E(b)` one gets, for `0 ≤ j ≤ M−1`,
`F^{(j)}(b) = −2·(D)_j · Σ_{k∈E(b)} c_k (b+a_k)^{M−1−j}` with `D = d+M−1`
(check this yourself); so one bad partner already kills the root at `b`.

Explicitly `c_k = 1/Π_{l≠k}(a_k − a_l)` (all nonzero), by partial fractions.

## The target

**Robust Hanson–Petridis, RHP(f, g).** For all odd primes `p` and all
`A, B ⊆ F_p` with `m = |A| ≤ √p` such that `e_b ≤ e` for every `b ∈ B`:

    m·n − r  ≤  f(e/m)·(p−1)/2 + g(m),

where `f : [0, 1/2) → [1, ∞)` satisfies `f(η) → 1` as `η → 0` and
`g(m) = o(p)` uniformly in `m ≤ √p` (e.g. `g(m) = O(m)`; **erratum
2026-09-26:** the earlier example `O(m²)` is of order `p` at `m ≈ √p` and
is too weak, as the `algebra` worker pointed out).

Known: `f(0) = 1, g = 0` (HP, exact). The second moment gives, for
`b ∉ −A`, `F_A(b) ≥ m − 2e`, hence the count of such `b` is at most
`m(p−m)/(m−2e)²`, i.e. `f(η) ≤ 2/(1−2η)²` up to the `r` term: the second
moment is off from HP by exactly the factor 2 at `η = 0`. RHP asks that
the factor-2 gain of Stepanov's method survive a small fraction of bad sums.

**Consequence to be proved rigorously (target theorem, CONDITIONAL on RHP).**
For every `κ > 0` there is `η(κ) > 0` such that for all large `p` and all
`A, B` with `|A||B| ≥ (1/2 + κ)p` (and the balanced regime
`|A|, |B| ≤ p^{1/2+o(1)}`, the unbalanced one being Karatsuba's),
`|S(A,B)| ≤ (1 − η(κ))|A||B|`. Sketch: negative sums
`= (mn − S − r)/2`; Markov gives `(1−1/λ)n` values of `b` with
`e_b ≤ λ·η·m`; apply RHP to that sub-rectangle, using `r ≤ min(m,n)`.
Chung's bound gives a constant saving only for `mn ≥ (1+κ)p`; HP gives only
`S ≤ mn − 2` for `mn > d + ω`. Nothing is known in between; an unconditional
proof of this consequence would be a genuinely new result for arbitrary
sets (not the full conjecture). Verify this claim of novelty against the
literature before asserting it.

The `F_{p²}` test: RHP with `r` retained is false over `F_{p²}` already at
`η = 0` (take `A = B = F_p`), so any proof must use the prime field, as HP's
does (degree `D ≤ p−1`, binomials and factorials nonzero mod `p`).

## Deliverables (per worker)

- `research/stepanov-<slug>-2026-09-26.md` — Status line first
  (PROVED / CITED / CONDITIONAL / HEURISTIC / REFUTED / OPEN), then precise
  statements, complete proofs, witnesses, and open obligations.
- `experiments/stepanov_<slug>_2026_09_26.py` — exact verifier (standard
  library, numpy, scipy, sympy; C++ helpers allowed if the source is saved
  under `experiments/` and the verifier re-checks their stored witnesses),
  under ten minutes, writing `results/stepanov_<slug>_2026_09_26.json`.
- A final message of at most 400 words.

Write the note and verifier within your first few actions and update them
incrementally: workers in earlier passes were killed by account usage limits
and anything not on disk was lost. Never `pkill` by a pattern that could
match another process's verifier.
