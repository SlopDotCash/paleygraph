# Stepanov wave: a robust Hanson–Petridis theorem and constant cancellation above p/2

**Status: the two-set Paley conjecture remains open. This wave proves,
with exact verification, a robust form of the Hanson–Petridis inequality
and, as a consequence, the first constant-factor cancellation bound for
arbitrary sets in the window `p/2 < |A||B| < p`, where Chung's bound is
trivial. The threshold `1/2` is sharp. The argument has been read line by
line by the orchestrator, re-proved independently by a referee agent,
and checked by three independent programs. It has not been reviewed by a
human mathematician, and targeted literature
searches found no prior statement of it; novelty is not established.**

## 1. The result (PROVED; `research/stepanov-robust-2026-09-26.md`)

Notation: `p` odd prime, `χ` the Legendre symbol, `d = (p−1)/2`,
`S(A,B) = Σ_{a∈A,b∈B} χ(a+b)`, `m = |A|`, `n = |B|`, `r = |B ∩ (−A)|`,
`e_b = #{a ∈ A : χ(a+b) = −1}`, `N₋ = Σ_{b∈B} e_b`.

**Theorem A (Hankel-minor Stepanov inequality, Thm 2.1).** For every
`A` with `1 ≤ m ≤ (p+1)/2` and every integer `0 ≤ e ≤ (m−1)/2`,

    Σ_{b ∈ F_p, e_b ≤ e} (e+1−e_b)·(m − (3e+e_b)/2 − [b∈−A])  ≤  (e+1)(d−e).

For `e = 0` this is exactly Hanson–Petridis.

**Theorem B (supersaturation, Cor. 2.3).** For all `A, B` with
`|A| ≤ (p+1)/2` and `0 ≤ e ≤ (m−1)/2`: `N₋ ≥ (e+1)·[n − (d−e+r)/(m−2e)]`. Every element of
`B` beyond the Hanson–Petridis bound forces about `e+1` non-residue sums.

**Theorem C (robust Hanson–Petridis, Thm 2.4).** If every `b ∈ B` has at
most `ηm` bad partners, `η ≤ 1/8`, then
`|A||B| − r ≤ (1−√(2η))^{−2}·(p−1)/2 + |A|`. (The robust note assumes
`|A| ≤ (p+1)/2`; the referee proves the bound without it, §8.2 of
`research/stepanov-referee2-2026-09-27.md`.)

**Theorem D (constant cancellation, Thm 2.5).** For `p ≥ 11`,
`0 < κ ≤ 3/2`, `u = (1+2κ)^{−1/2}` and all `A, B ⊆ F_p` with
`|A||B| ≥ (1/2+κ)p`:

    |S(A,B)| ≤ [ 1 − (1−u)² + (√((1/2+κ)p)+1)/(2(p−1)) ]·|A||B|.

The saving is about `κ²` for small `κ`. It cannot exceed about
`2κ/(1+2κ)`. The threshold `1/2` is sharp over all set sizes: for
`p ≡ 1 mod 4`, `A = {0,1}` with its complete partner set has
`|A||B| = (p+3)/2` and bias `1 − O(1/p)`. Whether `1/2` is sharp when both
sets grow is OPEN. The bound is vacuous for small `p` because of the
`O(p^{−1/2})` term: at `κ = 1/2` it bites from `p ≥ 47`, at `κ = 0.1` from
`p ≥ 2741` (referee).

**Corollary E (Paley graph, `p ≡ 1 mod 4`).** Taking `B = −A`: every set of
at least `√((1/2+κ)p)` vertices spans an induced subgraph with edge
density between `c(κ)/2 − o(1)` and `1 − c(κ)/2 + o(1)`, where
`c(κ) = (1−(1+2κ)^{−1/2})²`. Previously Hanson–Petridis gave only "not a
clique" in the range `√(p/2) < |A| < √p`, and Chung's bound gave density
bounded away from 1 only for `|A| ≥ (1+κ)√p`. (The corollary is an
immediate consequence of Theorem D; it is not separately verified.)

**Where the idea is.** Near a point `b` the normalized derivatives
`F^{(s)}/(D)_s` of the Hanson–Petridis polynomial split as `2ρ_s + ε_s`.
The part `ρ_s = Σ_{k∈E(b)} c_k(x+a_k)^{D−s}` is a sum of `e_b` pure powers,
so its Hankel matrix has rank `≤ e_b` **identically in `x`**, like
Reed–Solomon syndromes with error locators `(x+a_k)^{−1}`. The part
`ε_s` vanishes to order `≥ m−s` at `b`. Expanding the `(e+1)×(e+1)` Hankel
determinant row by row, every surviving term keeps at least `e+1−e_b`
rows of `ε`s, so the determinant vanishes to order about `(e+1−e_b)m` at
`b`, while its degree is exactly `(e+1)(d−e)`. The exact degree is a
Krattenthaler determinant whose factorials all lie below `p`; that is
where the prime field enters.

## 2. Verification

| check | cases | failures |
|---|---|---|
| worker verifier (`experiments/stepanov_robust_2026_09_26.py`), every step | 1,208,814 | 0 |
| independent exhaustive C check of Theorem A, every `A ∋ 0`, `|A| ≤ (p+1)/2`, all `e`, `p ≤ 29` (`experiments/stepanov_star_exhaustive_2026_09_27.c`) | 1,074,401,579 | 0 |
| Theorem D on the 165 most biased `F_p` rectangles found by the adversary | 165 | 0 |
| orchestrator's line-by-line reading of Thm 2.1 steps 1–5, Lemma 2.2, Thms 2.4–2.5 | — | none found |
| independent referee (`research/stepanov-referee2-2026-09-27.md`): own proof of every step; full polynomials `H_{e+1}` in 2,509 cases; `Λ ≢ 0` for every `p ≤ 1500`, `m ≤ (p+1)/2`, `e ≤ 12`; exhaustive `p ≤ 23`; Corollary E at 40 primes | 257,753,636 | 0 |

Theorem A is tight at `e = 1` for `p ≤ 13` (ratio exactly 1 in the
exhaustive check). Literature: Hanson–Petridis follow-ups (Kalmynin,
Rudnev–Tyrrell, Yip–Yoo, Kim–Yip–Yoo, Elsholtz–Wurzinger) treat
complete containment only; no robust form or constant-bias bound in this
window was found (searches of 2026-09-26 and 2026-09-27).

## 3. The other directions of the wave

- **Algebra** (`stepanov-algebra`, 2,118,597 checks): exact lifting of
  `S` from its residue when `|A||B| < p/2`; the kernel
  `binom(d,j) ≡ (−1/4)^j binom(2j,j)`; Hanson–Petridis for every sign
  pattern; every Johnson-type list-decoding bound collapses to the second
  moment.
- **Pencil** (`stepanov-pencil`, root, 18,191 checks): the
  Hanson–Petridis polynomials of all subsets `A∖E` form one explicit
  linear system, and a Riemann–Hurwitz inequality
  `(2M−2)N₀ + (M−2)N₁ ≤ 2d−2`. It gives only constant 2 and is superseded
  for the robust question by Theorem C, whose weights `(e+1−e_b)` come from
  using the rank of the error part identically in `x` rather than its
  values at `b`.
- **Adversary** (`stepanov-adversary`, partial, 61,052 checks): exhaustive
  data for `p ≤ 61` over 774 million sets; the best robust ratios found at
  `η = 1/16` are `0.70` (`p = 503`) and `0.54` (`p = 1009`), far inside
  Theorem C.
- **Lean** (`stepanov-leanhp`, 102,524 checks): Hanson–Petridis
  Theorem 1.2 for every proper divisor `d`, the Paley clique bound
  `|A|(|A|−1) ≤ (p−1)/2`, and the sharp `|A| = 2` example, all proved in
  Lean 4 with only the standard axioms, compiled against the two Mathlib
  builds on this machine. Nothing was submitted: `~/prove2me_workspace`
  and its API key no longer exist.

## 4. What this does and does not do for the conjecture

It moves the constant-cancellation threshold for arbitrary sets from
`|A||B| ≈ p` (Chung) to `|A||B| ≈ p/2`, which is the best possible
threshold, and it does so with a local, prime-sensitive argument. It does
not give a power saving and does not reach below the square-root scale.
The degree of the Hankel minor, `(e+1)(d−e)`, is still proportional to `d`
times the order, so the method is capped at `|A||B| ≍ p`, like
Hanson–Petridis itself.

## 5. Open obligations

1. Human review of `research/stepanov-robust-2026-09-26.md` §2, in
   particular Step 3 (multilinear expansion) and Lemma 2.2 (the
   Krattenthaler evaluation).
2. Close the gap `κ² ≲ η(κ) ≲ 2κ` for the saving.
3. The bias statement for characters of order `k` (Proposition 2.9 gives
   the inequality, not the bias bound).
4. Formalize Theorems A–D in Lean, building on the Hanson–Petridis
   formalization; submit to prove2me once a workspace and key exist.
5. Whether any variant of the Hankel-minor idea can reduce the degree
   budget below `(e+1)d`; that would be needed to reach below `√p`.
