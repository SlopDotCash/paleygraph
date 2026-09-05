# Cancellation across tuples: the two-set conjecture as a second-level cancellation statement for Weil sums

**Status: PROVED — (a) the exact tuple decomposition of the dilation `2k`-th moment over `t ∈ F_p^*` with the value of every tuple term (sign, odd-root set, even-root correction), (b) the exact reformulation "Paley ⟺ the non-square tuple sum is dominated" with explicit constants and quantifiers (Theorems 1.4–1.5; the converse direction requires the degenerate dilate `t = 0` to be excluded, which the prompt's `Σ_{t∈F_p}` does not do), (c) the exact `D`-decomposition `T_ns = Σ_D c(D)W_D − Corr` with `c(D)` explicit as a coefficient of `Π_{D} sinh(σ_r x)·Π_{∉D} cosh(σ_r x)` (signed class sums, not multiplicities), the bound `|E₂ − Corr| ≤ (k+1)·C(2k,2)·R_ν·(mn)^{2k−2}` on the elementary part, and the full dilation invariance of the decomposition; (d) exact bookkeeping of the three scales (triangle, pattern-random-sign, target) showing the conjecture asks for a saving of `(2k−1)p^{1/2+2kδ}` over the triangle inequality while allowing a loss of `(mn)^k p^{−1/2−2kδ}` over square-root cancellation. REFUTED (witnesses) — the literal statements "`W(τ) = p − #roots` for square `τ`", "`W(τ)` depends only on the odd-root set", and "`c(D)` counts tuples / is a function of `ν`"; the ordered-tuple scale `R` as the random-sign null; the strong conjecture `|T_ns| ≤ C_k√p·R_ν^k` (flat-profile witness `A = {1}`, `B = Q`, ratio `≈ 3√p → ∞`); and the prompt's form `C_k(mn)^{2k} + C_k√p R_ν^k` as anything that could imply the conjecture (it is consistent with all data but is implied by the trivial bound in the Paley regime, §5.4). SUPPORTED, not proved — square-root cancellation at the pattern level for random rectangles: `T_ns/R_D = O(1)` with slope `0` in `p` (§3), and Conjecture T(k) of §5, which holds with `C_2 ≤ 7.8`, `C_3 ≤ 36` on every one of the 663 computed instances, implies `C(ε,δ)` for explicit `δ`, and implies the open subgroup case with an explicit exponent. OPEN — the conjecture itself, Conjecture T(k), and any non-trivial upper bound on `T_ns` that beats the triangle inequality by a power of `p`.**

Worker `tuple`, prefix `sigma`, 2026-09-05. Verifier `experiments/sigma_tuple_2026_09_05.py` (standard library + numpy; ≈ 3.5 minutes), output `results/sigma_tuple_2026_09_05.json`. Every identity and inequality below is checked exactly (integer arithmetic) by the verifier; the number of checks of each kind is recorded there. Stored rectangles are read from `results/sigma_crux_2026_09_05_search.json` (the crux pass's search witnesses; `results/sigma_crux_2026_09_05.json` itself stores only summaries). Notation of the crux note `research/sigma-crux-2026-09-05.md` is used where possible, with one change of convention stated in §0.

---

## 0. Setup, and one deviation from the prompt

`p` is an odd prime, `χ` the Legendre symbol with `χ(0) = 0`, `A, B ⊆ F_p` with `0 ∉ A`, `m = |A|`, `n = |B|`, `Π = A×B`, and

    g(t) = Σ_{a∈A} Σ_{b∈B} χ(ta + b)   (t ∈ F_p),    S(A,B) = g(1),    g(t) = S(tA, B).

**Class data.** For `(a,b) ∈ Π` put `r(a,b) = −b/a`. For `r ∈ F_p` let `C_r = {(a,b) ∈ Π : r(a,b) = r}`,
`ν_r = |C_r|`, `σ_r = Σ_{(a,b)∈C_r} χ(a)`, `𝓡 = {r : ν_r > 0}` (the ratio set, `L = |𝓡|`), `R_ν = Σ_r ν_r²`,
`R_χ = Σ_r σ_r²`. Elementary facts: `Σ_r ν_r = mn`; `ν_r ≤ min(m,n)` (in a class, `a` determines `b` and
`b` determines `a`), so `mn ≤ R_ν ≤ mn·min(m,n)`; `|σ_r| ≤ ν_r`; and, since `χ(ta+b) = χ(a)χ(t − r(a,b))`,

    g(t) = Σ_{r∈𝓡} σ_r χ(t − r).                                                        (0.1)

**Tuples.** For `τ = ((a_i,b_i))_{i=1}^{2k} ∈ Π^{2k}` put `P_τ(t) = Π_i (ta_i + b_i)`, `c_τ = Π_i a_i ≠ 0`,
`j_τ(r) = #{i : r(a_i,b_i) = r}` (the multiplicity pattern, a `2k`-multiset on `𝓡`),
`D(τ) = {r : j_τ(r) odd}` (the odd-root set; `|D(τ)|` is even because `Σ_r j_τ(r) = 2k`),
`E(τ) = {r : j_τ(r) even, > 0}`. `τ` is *square* if `D(τ) = ∅`.

**The deviation: `t` ranges over `F_p^*`.** The prompt defines `W(τ) = Σ_{t∈F_p} χ(P_τ(t))`. The term `t = 0`
is the degenerate dilate `0·A = {0}`, where `g(0) = m·Σ_{b∈B} χ(b)`; for `B ⊆ Q` (the nonzero residues) this
equals `mn` whatever the conjecture says, so the full-`F_p` moment `Σ_{t∈F_p} g(t)^{2k}` is `≥ (mn)^{2k}` for
such `B` and the converse direction of the reformulation (Theorem 1.4(ii)) would be false. Throughout this
note, therefore,

    W(τ) := Σ_{t∈F_p^*} χ(P_τ(t)) = Σ_{t∈F_p} χ(P_τ(t)) − Π_i χ(b_i),
    W_D  := Σ_{t∈F_p^*} Π_{r∈D} χ(t − r)   (D ⊆ 𝓡),      W_∅ = p − 1,

and all moments are `Σ_{t∈F_p^*}`. Every statement below transfers to the prompt's convention by adding the
explicit term `g(0)^{2k} = (m·Σ_b χ(b))^{2k}` to the moment and `Π_i χ(b_i)` to each `W(τ)` (the verifier
records `g(0)` and the square-tuple part `Q_k(0)` of that term for every instance). With this convention

    Σ_{t∈F_p^*} g(t)^{2k} = Σ_{τ∈Π^{2k}} W(τ) = T_sq + T_ns,
    T_sq = Σ_{τ square} W(τ),   T_ns = Σ_{τ not square} W(τ),
    T_tri = Σ_{τ not square} |W(τ)|,   R = (Σ_{τ not square} W(τ)²)^{1/2},   N_ns = #{τ not square}.

---

## 1. The tuple decomposition and the exact reformulation (Task 1)

**Lemma 1.1 (value of a tuple term; PROVED).** `P_τ(t) = c_τ · Π_{r∈𝓡} (t − r)^{j_τ(r)}`, and

    W(τ) = χ(c_τ) · W̃(τ),      W̃(τ) = Σ_{t∈F_p^*, t∉E(τ)} Π_{r∈D(τ)} χ(t − r).

If `τ` is square, `W̃(τ) = p − 1 − |E(τ)∖{0}|`. If `τ` is not square,

    W̃(τ) = W_{D(τ)} − Σ_{e∈E(τ)∖{0}} Π_{r∈D(τ)} χ(e − r),

and (Weil) `|W_D| ≤ (|D| − 1)√p + 1` for `D ≠ ∅`; moreover `W_{{r,r'}} = −1 − χ(rr')` for `r ≠ r'`.

*Proof.* `ta + b = a(t − r(a,b))`, which gives the factorisation. For fixed `t`, `χ((t−r)^j)` equals `1` if
`j` is even and `t ≠ r`, `0` if `t = r`, and `χ(t − r)` if `j` is odd. Hence `χ(P_τ(t)) = χ(c_τ)·1[t ∉ E(τ)]·
Π_{r∈D(τ)} χ(t − r)` (for `t ∈ D(τ)` both sides vanish). Summing over `t ∈ F_p^*` gives the formula; for a square
`τ` the product is empty and the count of admissible `t` is `p − 1 − |E(τ)∖{0}|`; for a non-square `τ`,
`Σ_{t≠0, t∉E} = Σ_{t≠0} − Σ_{t∈E∖{0}}`. The polynomial `Π_{r∈D}(t − r)` is squarefree of degree `|D| ≥ 2`, not a
constant times a square, so `|Σ_{t∈F_p} χ(Π_{r∈D}(t−r))| ≤ (|D|−1)√p` (Weil, brief); removing `t = 0` costs at most
`1`. For `|D| = 2`, `Σ_{t∈F_p} χ((t−r)(t−r')) = −1` (brief) and the `t = 0` term is `χ(rr')`. ∎

**Corrections to the prompt's description (REFUTED as stated, witnesses in the verifier's `brute_force`
section).** (i) For a square `τ` the value is `χ(c_τ)(p − 1 − #nonzero distinct roots)`, not `p − #roots`: at
`p = 7`, `A = {6,1,2}`, `B = {4,5}`, `τ = ((6,4),(6,4),(6,5),(2,4))` has ratios `4,4,5,5` (square, two nonzero roots)
and `W(τ) = χ(6·6·6·2)(7−1−2) = χ(5)·4 = −4`. The *sum* `T_sq` is nevertheless nonnegative (Proposition 1.2).
(ii) `W(τ)` is not a function of `D(τ)` alone: the even-root correction is nonzero (`p = 7`, `A = {1,2}`,
`B = {0,2,1}`, `k = 2`: `Σ_D c(D)W_D = −1104`, `Corr = −96`, `T_ns = −1008`), and the sign `χ(c_τ)` is not a function
of the pattern (Remark 2.3).

**Proposition 1.2 (the square part; PROVED, after crux Theorem 3.1/Proposition 3.2, re-proved here).**

    T_sq = Σ_{t∈F_p^*} Q_k(t),   Q_k(t) := Σ_{τ square} χ(P_τ(t)) = (2k)!·[x^{2k}] Π_{r∈𝓡, r≠t} cosh(σ_r x) ≥ 0,

hence `0 ≤ T_sq ≤ (p−1)·N_k`, where `N_k = #{square τ} = (2k)!·[x^{2k}] Π_r cosh(ν_r x) ≤ (2k−1)!!·R_ν^k`.

*Proof.* By Lemma 1.1, `χ(P_τ(t)) = χ(c_τ)·1[t ∉ E(τ)]` for square `τ`. Group square tuples by their (all-even)
pattern `j`: the positions carrying each class can be chosen in `(2k)!/Π_r j_r!` ways and then the sum of
`χ(c_τ) = Π_i χ(a_i)` over the class members factorises as `Π_r σ_r^{j_r}`. So `Q_k(t) = Σ_{j even, j_t = 0}
(2k)!/Π j_r! · Π σ_r^{j_r} = (2k)!·[x^{2k}] Π_{r≠t} Σ_{j even} (σ_r x)^j/j!`, the stated coefficient. All Taylor
coefficients of `cosh` are nonnegative, so `Q_k(t) ≥ 0`, and they increase when `|σ_r|` is replaced by `ν_r ≥ |σ_r|`,
so `Q_k(t) ≤ (2k)![x^{2k}]Π_r cosh(ν_r x) = N_k` (the same count with `χ ≡ 1`). Finally, every square tuple is the
image of a perfect matching of the `2k` positions (`(2k−1)!!` choices) and, for each matched pair, an ordered pair of
elements of `Π` in a common class (`R_ν` choices): match the positions of each class arbitrarily. Hence
`N_k ≤ (2k−1)!! R_ν^k`. ∎

**Theorem 1.3 (one direction, exact; PROVED).** Let `η ∈ (0,1]`, `k ≥ 1`. If `T_sq + T_ns ≤ (η mn)^{2k}` then
`|g(t)| ≤ η mn` for every `t ∈ F_p^*`, in particular `|S(A,B)| ≤ η mn`; if the hypothesis is strict, so is the
conclusion. In particular this holds if `T_ns ≤ ½(η mn)^{2k}` and `T_sq ≤ ½(η mn)^{2k}`.

*Proof.* `Σ_{t≠0} g(t)^{2k} = T_sq + T_ns` and every summand is `≥ 0`, so `g(t)^{2k} ≤ (η mn)^{2k}`. ∎

**Definition.** `C(ε,δ)`: there is `p₀` such that for all `p > p₀` and all `A,B ⊆ F_p` with `|A|,|B| > p^ε`,
`|S(A,B)| ≤ p^{−δ}|A||B|`. `C°(ε,δ)`: the same restricted to pairs with `0 ∉ A`. For `k ≥ 1`:
`M(ε,δ',k)`: there is `p₀` such that for all `p > p₀` and all `A,B` with `0 ∉ A`, `|A|,|B| > p^ε`,

    T_ns(A,B;k) ≤ p^{−2kδ'} (mn)^{2k}.

**Theorem 1.4 (exact converse and forward bookkeeping; PROVED).**

(i) If `k(ε − 2δ') > 1` and `M(ε,δ',k)` holds, then `C°(ε,δ)` holds for every `0 < δ < δ'`.

(ii) If `C(ε,δ)` holds, then for every `k ≥ 1` and all `p > p₀`, all `A,B` with `0 ∉ A`, `|A|,|B| > p^ε`:

    −(2k−1)!!·p^{1−kε}·(mn)^{2k}  ≤  −T_sq  ≤  T_ns  ≤  (p−1)·p^{−2kδ}·(mn)^{2k}.

In particular `M(ε, δ − 1/(2k), k)` holds for every `k > 1/(2δ)`, and `|T_ns| ≤ 2(2k−1)!!·p^{1−2k·min(δ,ε/2)}(mn)^{2k}`.

(iii) `C°(ε,δ)` implies `C(ε',δ')` for all `ε' > ε` and `δ' < min(δ, ε)`.

*Proof.* (i) By Proposition 1.2 and `R_ν ≤ mn·min(m,n)`,
`T_sq ≤ (2k−1)!!(p−1)(mn)^k min(m,n)^k = (2k−1)!!(p−1)(mn)^{2k}/max(m,n)^k < (2k−1)!!·p^{1−kε}(mn)^{2k}`.
With the hypothesis, `Σ_{t≠0} g^{2k} < [p^{−2kδ'} + (2k−1)!!p^{1−kε}](mn)^{2k}`. Since `1 − kε < −2kδ'`, for `p`
large the bracket is `≤ 2p^{−2kδ'} ≤ p^{−2kδ}` for any `δ < δ'`, and Theorem 1.3 gives `|S(A,B)| ≤ p^{−δ}mn`.
(ii) `g(t) = S(tA,B)` with `|tA| = |A|` for `t ≠ 0`, so `|g(t)| ≤ p^{−δ}mn` for all `t ∈ F_p^*` and
`T_sq + T_ns ≤ (p−1)p^{−2kδ}(mn)^{2k}`; the lower bound is `T_sq ≥ 0` and the bound on `T_sq` from (i).
`(p−1)p^{−2kδ} ≤ p^{−2k(δ − 1/(2k))}`. (iii) If `0 ∈ A`, `|S(A,B) − S(A∖{0},B)| = |Σ_b χ(b)| ≤ n ≤ p^{−ε}mn` and
`|A∖{0}| > p^{ε}` for `p` large (as `|A| > p^{ε'}`). ∎

**Theorem 1.5 (the equivalence, with quantifiers; PROVED).** For `0 < ε < 1` the following are equivalent:

(a) there is `δ₀ > 0` with `C°(ε,δ₀)`;

(b) there are an integer `k ≥ 1` and `δ' > 0` with `k(ε − 2δ') > 1` such that `M(ε,δ',k)` holds, i.e. *for one
`k = k(ε)`, the `2k`-th dilation moment over `F_p^*` is dominated by its square tuples up to `p^{−2kδ'}(mn)^{2k}`*.

Moreover in (b) the square tuples are automatically negligible: `T_sq ≤ (2k−1)!!p^{1−kε}(mn)^{2k} = o(p^{−2kδ'}(mn)^{2k})`.

*Proof.* (b)⟹(a) is Theorem 1.4(i). (a)⟹(b): replace `δ₀` by `δ₁ = min(δ₀, ε/3)` (still `C°(ε,δ₁)`), choose
`k > 1/(2δ₁)` and `δ' = δ₁ − 1/(2k) > 0`; then `k(ε − 2δ') = k(ε − 2δ₁) + 1 > 1`, and `M(ε,δ',k)` is
Theorem 1.4(ii). ∎

**Proposition 1.6 (the square-term threshold; PROVED).** `T_sq ≤ ½(η mn)^{2k}` holds whenever
`(mn)^{2k} ≥ 2(2k−1)!!(p−1)R_ν^k η^{−2k}`. Since `mn ≤ R_ν ≤ mn·min(m,n)`, a sufficient condition is
`max(m,n)^k ≥ 2(2k−1)!!(p−1)η^{−2k}`, and when all ratios are distinct (`R_ν = mn`) the condition is exactly
`(mn)^k ≥ 2(2k−1)!!(p−1)η^{−2k}`, i.e. `mn ≥ (2(2k−1)!!p)^{1/k}η^{−2}` — the prompt's `mn ≥ p^{1/k+o(1)}`. In the
regime of Theorem 1.5 (`k > 1/ε`) the condition holds for all large `p` with `η = p^{−δ'}`. (The verifier checks
`0 ≤ T_sq ≤ (p−1)(2k−1)!!R_ν^k` on every instance.)

**Remark 1.7 (what "second-level cancellation" means, exactly).** Write `N_ns = (mn)^{2k} − N_k ≤ (mn)^{2k}` for
the number of non-square tuples. Three scales for `T_ns`:

| scale | exact bound / value | size for distinct ratios, `mn = p^{2ε'}` |
|---|---|---|
| triangle `T_tri` | `≤ ((2k−1)√p + k + 1)·N_ns` (Lemma 1.1, `|E| ≤ k`) | `≈ (2k−1)√p·(mn)^{2k}` |
| pattern random-sign `R_D` (§2) | `(Σ_D c₊(D)² W_D²)^{1/2} ≤ ((2k−1)√p+k+1)·√((2k)!)·(mn)^{(3k−1)/2}·min(m,n)^{(k+1)/2}` (Lemma 1.7a) | `≈ 4.5·√p(mn)^k` (`k=2`), `61·√p(mn)^k` (`k=3`), observed |
| ordered-tuple `R` (prompt) | `(Σ_τ W(τ)²)^{1/2} ≤ ((2k−1)√p+k+1)√N_ns` | `≈ √p·(mn)^k` |
| target of Theorem 1.5 | `p^{−2kδ'}(mn)^{2k}` | `= p^{−2kδ'}(mn)^{2k}` |

**Lemma 1.7a (PROVED).** For every non-square `D`, `c₊(D) ≤ (2k)!/(2k−|D|)! · (2k−|D|−1)!! · min(m,n)^{|D|} · R_ν^{k−|D|/2}`,
hence `max_D c₊(D) ≤ (2k)!·(mn)^{k−1}·min(m,n)^{k+1}`, and

    R_D ≤ ((2k−1)√p + k + 1) · √((2k)!) · (mn)^{(3k−1)/2} · min(m,n)^{(k+1)/2}.

(When `max(m,n) ≥ 2k` the constant `(2k)!` can be replaced by `2k(2k−1)!!`, the value of the `|D| = 2` term; the
verifier tests this sharper form and it holds on every instance.)

*Proof.* A tuple with odd-root set `D` has at least one position in each class `r ∈ D`; fix the first such position
for each `r` (at most `(2k)!/(2k−|D|)!` position choices, `Π_{r∈D}ν_r ≤ min^{|D|}` element choices). Removing these
positions leaves a tuple whose multiplicities are all even, i.e. a square `(2k−|D|)`-tuple, of which there are at most
`(2k−|D|−1)!!R_ν^{(2k−|D|)/2}` (Proposition 1.2). This map is injective. With `R_ν ≤ mn·min` the bound is
`≤ C_{k,|D|}·(mn)^{k−|D|/2}min^{k+|D|/2} = C_{k,|D|}·(mn)^{k−1}min^{k+1}·(min/(mn))^{|D|/2−1} ≤ (2k)!·(mn)^{k−1}min^{k+1}`,
using `min ≤ mn`, `|D| ≥ 2` and `C_{k,d} = (2k)!/(2k−d)!·(2k−d−1)!! ≤ (2k)!`. (The ratio of consecutive bounds is
`B(d−2)/B(d) = max(m,n)/(2k−d+2)`, so for `max(m,n) ≥ 2k` the maximum is at `|D| = 2` with constant `2k(2k−1)!!`.)
Then `Σ_D c₊(D)² ≤ max_D c₊(D)·Σ_D c₊(D) ≤ max_D c₊(D)·(mn)^{2k}` and `|W_D| ≤ (2k−1)√p + 1`, while
`R_D² = Σ c₊(D)²W_D²`. ∎ (All three inequalities are checked exactly on every instance; a first version of this
lemma claimed the maximum at `|D| = 4` and was caught by that check.)

So the conjecture is equivalent to: *the `N_ns` Weil sums `W(τ)`, each of size up to `(2k−1)√p + k + 1`, cancel
among themselves by a factor `(2k−1)p^{1/2+2kδ'}` relative to the triangle inequality* — a saving of "more than
`√p`" — while it *allows* a loss relative to square-root cancellation at the pattern level: by Lemma 1.7a,
`R_D/(mn)^{2k} ≤ C_k√p/max(m,n)^{(k+1)/2} < C_k p^{1/2−ε(k+1)/2}` in the Paley regime, so `|T_ns| ≤ p^{c}·R_D` for
*any* `c < ε(k+1)/2 − 1/2` (and `k > 1/ε`) implies `C(ε,δ)` for every `δ < (ε(k+1)/2 − 1/2 − c)/(2k)`: the
conjecture sits strictly between "no cancellation" and "square-root cancellation" and is far weaker than the
latter — a loss of `p^{c}` with `c` up to `ε(k+1)/2 − 1/2 → ∞` (as `k → ∞`) over the pattern random-sign scale is
affordable. (This is not a proof strategy: no method is known that gives *any* power saving over the triangle
inequality for these sums with `k` fixed; §5–6.)

**Remark 1.8 (is this the "cancellation between words" of the spectral line?).** No, and no identity between the
two is claimed. The spectral line (`research/parallel9-all-degrees-2026-09-05.md`, `parallel2-spectral-transfer`)
bounds necklace traces `N = tr(D_{Z₁}S⋯D_{Z_k}S) = Σ_{x₁,…,x_k∈F_p} Π_i d_{Z_i}(x_i)χ(x_i − x_{i+1})` — complete sums in
`k` variables, one per word `(Z₁,…,Z_k)` of anchor subsets, with a Weil-type bound `p^{(k+1)/2}` and an exponential
constant `(2a+2)^{k−2}`; the aggregate criterion needs the *signed sum over all `(b−1)^k` words* to have
subexponential growth in `k` at fixed anchors `a`, and it targets clique bounds (`A = B = clique`). The present
statement concerns *one-variable* complete sums `W_D`, `D ⊆ 𝓡`, each bounded by `(|D|−1)√p`, weighted by the
combinatorial coefficients `c(D)` of §2, and needs a *power-of-`p` saving at fixed `k`*; it is *equivalent* to the
two-set conjecture (Theorem 1.5), whereas the necklace aggregate is a sufficient condition for a weaker statement.
Both are of the shape "a sum of complete character sums indexed by combinatorial data must cancel beyond the termwise
Weil bound", which is presumably why both are stuck; but the index sets (subsets of a product set `(−B)·A^{−1}`
versus words in anchors), the number of variables, and the regime (`p → ∞` at fixed `k` versus `k → ∞` at fixed `a`)
differ, and I found no change of variables making one a special case of the other.

---

## 2. Structure of the non-square sum (Task 3, proofs)

For `D ⊆ 𝓡` with `|D|` even define

    c(D)  := Σ_{j : D(j) = D} (2k)!/Π_r j_r! · Π_r σ_r^{j_r} = (2k)!·[x^{2k}] Π_{r∈D} sinh(σ_r x) · Π_{r∈𝓡∖D} cosh(σ_r x),
    c₊(D) := #{τ : D(τ) = D}                                 = (2k)!·[x^{2k}] Π_{r∈D} sinh(ν_r x) · Π_{r∈𝓡∖D} cosh(ν_r x),

where `j` runs over multiplicity patterns of total `2k` on `𝓡` (the generating function: odd multiplicities at
`r ∈ D` give `sinh`, even ones (including `0`) elsewhere give `cosh`). For `|D| = 2k`, `c(D) = (2k)!·Π_{r∈D} σ_r` and
`c₊(D) = (2k)!·Π_{r∈D} ν_r`. Note `|c(D)| ≤ c₊(D)`, and `c(D) = Σ_{τ: D(τ)=D} χ(c_τ)` is a *signed* count.

**Theorem 2.1 (exact `D`-decomposition; PROVED).**

    T_ns = Σ_{∅≠D⊆𝓡, |D| even ≤ 2k} c(D)·W_D − Corr,
    Corr = Σ_{j : D(j)≠∅} (2k)!/Π_r j_r! · Π_r σ_r^{j_r} · Σ_{e∈E(j)∖{0}} Π_{r∈D(j)} χ(e − r).

Splitting by `|D|`, `T_ns = T_W + E₂ − Corr` with the *Weil part* `T_W = Σ_{|D|≥4} c(D)W_D` and the *elementary part*
`E₂ = Σ_{|D|=2} c(D)W_D = −Σ_{r<r'} c({r,r'})(1 + χ(rr'))` (no character sum over `t` occurs in `E₂` or `Corr`).

*Proof.* Group the non-square tuples by pattern `j`. By Lemma 1.1 all tuples of pattern `j` have the same
`W̃(τ) = W_{D(j)} − Σ_{e∈E(j)∖{0}} Π_{r∈D(j)} χ(e−r)`, and `Σ_{τ of pattern j} χ(c_τ) = (2k)!/Πj_r!·Πσ_r^{j_r}` as in
Proposition 1.2. Summing over `j` with `D(j) = D` gives the coefficient `c(D)` of `W_D`, and the remaining terms are
`Corr`. `W_{{r,r'}} = −1 − χ(rr')` by Lemma 1.1. ∎

**Corollary 2.2 (the elementary part is small; PROVED).**

    |E₂ − Corr| ≤ (k+1)·C(2k,2)·R_ν·(mn)^{2k−2} ≤ (k+1)·C(2k,2)·(mn)^{2k}/max(m,n),

and `|T_W| ≤ ((2k−1)√p + 1)·N_W` with `N_W = #{τ : |D(τ)| ≥ 4} ≤ N_ns`.

*Proof.* A tuple with `|D(τ)| = 2 < 2k`, or with `E(τ) ≠ ∅`, has two positions `i < i'` in the same class; the number
of such tuples is at most `C(2k,2)·R_ν·(mn)^{2k−2}` (choose the positions, the ordered pair of same-class elements,
and the rest). Such a tuple contributes at most `|W_{{r,r'}}| ≤ 2` to `|E₂|` (only if `|D| = 2`, in which case
`|E| ≤ k−1`) and at most `|E(τ)∖{0}| ≤ k − 1` to `|Corr|`, hence at most `k + 1` in total; a tuple with `|D| ≥ 4`
contributes at most `|E| ≤ k − 2`. The second inequality is `R_ν ≤ mn·min(m,n)`. The bound on `T_W` is Lemma 1.1. ∎

**Remark 2.3 (the two naive forms are false; REFUTED with witnesses).** (i) The unsigned version
`Σ_D c₊(D)W_D` is not `T_ns` and not even `Σ_D c(D)W_D`: `p = 7`, `A = {6,1,2}`, `B = {4,5}`, `k = 2` gives
`Σ c(D)W_D = 64`, `T_ns = 16`, but `Σ c₊(D)W_D = −608`. The weights are the *signed class sums* `σ_r = Σ_{C_r} χ(a)`;
they equal the multiplicities `ν_r` only when `A ⊆ Q`. (ii) Dropping `Corr` is wrong: witness in §1. Both are checked
by full enumeration of all `(mn)^{2k}` tuples (37 cases, `p ≤ 13`, `k ≤ 3`; 362 312 tuple values compared with
Lemma 1.1), in agreement with the pattern-grouped computation for `T_sq, T_ns, T_tri, R², N_k, Corr` and
`Σ_D c(D)W_D`. (iii) The constant `k+1` of Corollary 2.2 cannot be lowered to `k`: the `gp` rectangle at `p = 1213`,
`(m,n) = (6,4)`, `k = 2`, `R_ν = 76` has `|E₂ − Corr| = 532 528 > 2·6·76·24² = 525 312` (the proved bound is `787 968`).

**Remark 2.4 (dilation invariance; PROVED).** Under `(A,B) ↦ (uA,B)`, `u ∈ F_p^*`: `r ↦ r/u`, `ν'_{r/u} = ν_r`,
`σ'_{r/u} = χ(u)σ_r`, `W'_{D/u} = Σ_{t≠0}Π_{r∈D}χ(t − r/u) = χ(u)^{|D|}Σ_{s≠0}Π_{r∈D}χ(s − r) = W_D` (`|D|` even), and
`c'(D/u) = χ(u)^{2k}c(D) = c(D)`; likewise `Corr`, `E₂`, `T_W`, `T_sq`, `T_tri`, `R`, `R_D` are unchanged. So the
entire decomposition is a function of the dilation orbit `{(uA,B)}` — it "sees" `t = 1` no more than any other
dilate, which is the exact content of crux Proposition 4.6 and of Theorem 1.5 (the moment statement is a statement
about the whole orbit). The prompt's sentence "`c(D)` depends on `ν`, the only place `A×B` enters beyond the ratio
set" should read: the decomposition depends on `(A,B)` exactly through the class data `(𝓡, ν, σ)`, i.e. through the
product set `𝓡 = (−B)·A^{−1}` with its representation function `ν = r_{(−B)·A^{−1}}` and its `χ(a)`-weighted version
`σ`. Nothing is broken by dilation; what the ratio set alone does not determine is `(ν,σ)` — and `min(m,n)`, which is
where Conjecture T(k) of §5 puts the saving.

---

## 3. Exact data (Task 2)

### 3.1 Instances and validation

The verifier computes, exactly, the decomposition of §§1–2 for 663 rectangles (700 pattern-grouped analyses
including the brute-force cases): on the 15 primes `p ∈ {101, 151, 211, 307, 401, 503, 701, 809, 1009, 1213, 1409,
1801, 2203, 2609, 2999}`, six constructions at `(m,n,k) ∈ {(6,6,2), (6,4,2), (4,4,2), (4,4,3)}` (all) and at
`(5,5,3)`, `(6,4,3)` (a subset), one instance per prime and construction (`random` twice):

* `random`: `A ⊆ F_p^*`, `B ⊆ F_p` uniform of sizes `m, n`;
* `interval`: `A = {1,…,m}`, `B = {u+1,…,u+n}`, `u` random;
* `gp`: `A = c·{r^i}_{i<m}`, `B = c'·{r^j}_{j<n}` with `r` a *quadratic residue* of order `> m+n` — the rectangle
  form of the crux note's biased-dilate mechanism: the ratio set collapses to `m+n−1` classes with triangular
  multiplicities and `σ_r = χ(c)ν_r` (all class sums of one sign);
* `gp_nr`: the same with `r` a non-residue: `χ(cr^i)` alternates, `σ_r ∈ {0, ±1}`, the profile is degenerate
  (`|S| ≤ 3`, `T_ns ≈ 0`); kept only to exhibit the classwise cancellation;
* `gp_biased`: `A = {r^i}_{i<m}` for a primitive root `r`, `B` = `n` random elements of `N(A) = {b : χ(a+b) = 1 ∀a∈A}`
  (bias `1` with the several biased dilates `r^j` of crux Theorem 6.2; `11–15` instances, where `|N(A)| ≥ n`);
* `greedy`: `A` random, `B` = the `n` shifts `b` maximising `Σ_a χ(a+b)` (bias `1` at every `p ≥ 307`; the
  "generic complete rectangle");

plus the 176 stored near-extremal rectangles of `results/sigma_crux_2026_09_05_search.json` with `mn ≤ 36`
(71 exhaustive maxima, 27 annealing `6×6`, 78 `k*` rectangles with `k ∈ {4,5,6}`; all have `|S| = mn`; one had
`0 ∈ A`, which was removed; their stored `S` is re-verified), all at `k = 2`; five lopsided witnesses `A = {1}`,
`B = Q` (`p ∈ {61,101,151}`, `k=2`; `p ∈ {41,61}`, `k=3`); and one subgroup rectangle (`A = H` of order `6` at
`p = 1009`, `B ⊆ N(H)`). For every instance the verifier checks the identities `T_sq + T_ns = Σ_{t≠0} g^{2k}`,
`N_sq + N_ns = (mn)^{2k}`, `T_ns = Σ_D c(D)W_D − Corr = T_W + E₂ − Corr`, `0 ≤ T_sq ≤ (p−1)(2k−1)!!R_ν^k`,
`0 ≤ Q_k(0) ≤ (2k−1)!!R_ν^k`, the termwise Weil bound `(|W_D|−1)² ≤ (|D|−1)²p`, `W_{{r,r'}} = −1 − χ(rr')`, the
triangle bounds of Remark 1.7, Corollary 2.2, and (Lemma 1.1 through brute force) 362 312 individual tuple values.
Lemma 1.7a is checked in its per-`D`, maximal and `R_D` forms (`3 × 695`). Total: 375 169 exact checks, 0 failures, 209 s.

Scales reported: `T_tri`, the prompt's `R`, the pattern scale `R_D = (Σ_D c₊(D)²W_D²)^{1/2}` and its Weil-part
version `R_{DW}` (`|D| ≥ 4`). `R` is *not* an admissible null scale: `W(τ)` depends only on the pattern, so the
`(2k)!` orderings of a distinct-ratio tuple carry one and the same Weil sum; for distinct ratios
`R_D/R ≈ √((2k)!)` at `k = 2` (observed `R_D/(√p(mn)^k)`: median `4.49`, `√24 = 4.90`; at `k = 3` the `|D| < 2k`
patterns, which have few `D`'s but large `c₊`, raise it to `60.8`). All ratios below use `R_D`.

### 3.2 The three regimes (medians over the 15 primes; `max` in brackets)

`k = 2`, `(m,n) = (6,6)`:

| family | bias | `T_W/(mn)^4` | `|T_ns|/R_D` | `|T_W|/R_{DW}` | `|T_ns|/T_tri` | `|E₂−Corr|/(mn)^4` | `T_ns<0` | `ρ_{T,W}` max |
|---|---|---|---|---|---|---|---|---|
| random (×2) | 0.08–0.11 | −0.01 | 0.73 / 1.03 [3.5] | 0.42 / 0.79 [2.4] | 0.0035 | 0.03 | 0.60 | 0.30 |
| interval | 0.22 | −0.02 | 0.73 [4.1] | 0.64 [4.3] | 0.0044 | 0.02 | 0.60 | 0.54 |
| gp (residue) | 0.22 | −0.44 | 2.17 [5.2] | 1.13 [4.6] | 0.185 | **1.00** | 0.67 | 6.02 |
| gp_nr | 0.06 | 0.000 | 0.00 | 0.00 | 0.0001 | 0.000 | — | 0.00 |
| gp_biased | 1.00 | **1.93** | 12.2 [21.6] | 12.3 [21.7] | 0.085 | 0.008 | 0.00 | 3.78 |
| greedy | 1.00 | **0.81** | 6.9 [12.8] | 7.2 [13.0] | 0.039 | 0.011 | 0.00 | 1.78 |
| stored anneal 6×6 (27) | 1.00 | 0.82 | 7.4 [18.4] | 8.7 [19.1] | 0.047 | 0.064 | 0.00 | 3.06 |
| stored k* 6×6 (25) | 1.00 | 0.85 | 8.9 [19.5] | 9.1 [20.8] | 0.056 | 0.081 | 0.00 | **9.17** |
| stored exhaustive 6×6 (10) | 1.00 | 0.80 | 6.2 [6.8] | 9.5 [10.5] | 0.074 | 0.27 | 0.00 | 3.22 |

`k = 3`, `(m,n) = (5,5)`:

| family | bias | `T_W/(mn)^6` | `|T_ns|/R_D` | `|T_W|/R_{DW}` | `|T_ns|/T_tri` | `|E₂−Corr|/(mn)^6` | `T_ns<0` | `ρ_{T,W}` max |
|---|---|---|---|---|---|---|---|---|
| random (×2) | 0.12 | −0.05 / −0.02 | 0.80 / 0.80 [1.7] | 0.80 / 0.58 [1.7] | 0.005 | 0.02 | 0.73 | 0.85 |
| gp (residue) | 0.28 | −0.94 | 1.93 [3.7] | 1.30 [3.3] | 0.199 | **0.75** | 0.73 | **33.5** |
| gp_biased (13) | 1.00 | **1.53** | 10.2 [16.1] | 10.5 [16.1] | 0.057 | 0.006 | 0.00 | 11.6 |
| greedy | 1.00 | **0.91** | 6.4 [19.9] | 6.5 [20.0] | 0.033 | 0.019 | 0.00 | 9.65 |

(`ρ_{T,W} = |T_W|·min(m,n)^k/(√p (mn)^{2k})` is the ratio of Conjecture T(k), §5. The `(6,4)` and `(4,4)` families
behave in the same way; all rows are in the results file, key `fits_by_type`.)

Three regimes, all forced or explained by the identities of §§1–2:

1. **Random and interval rectangles: square-root cancellation at the pattern level.** `|T_ns|/R_D` has median
   `0.6–1.0` and maximum `3.5` over 120 random instances at six `(m,n,k)`; the pooled log–log slope against `p`
   (60 instances per `k`, all sizes) is `+0.06 ± 0.11` (`k=2`) and `+0.05 ± 0.11` (`k=3`): `|T_ns| ≍ R_D ∝ √p(mn)^k`.
   Relative to the triangle bound, `|T_ns|/T_tri ≈ 0.004–0.02`, slope `−0.03 ± 0.14` / `−0.05 ± 0.13` (both scale
   as `√p`). The sign of `T_ns` is negative in `60–73 %` of the random instances: the moment is slightly *below* the
   independent-signs value `T_sq`, mostly because of the elementary part (HEURISTIC closed form, §6(3):
   `E₂ ≈ −6R_χ(g(0)² + g(∞)²) + 12R_χ²` at `k=2`, with `g(∞) := n·Σ_a χ(a)`; exactly, median
   `|E₂ − Corr| ≈ 0.02–0.06·(mn)^{2k}`, `p`-independent). The pure Weil part has
   `|T_W|/R_{DW}` median `0.4–0.8`, i.e. slightly *more* cancellation than independent signs, with a mild upward drift
   over this range (pooled slope `+0.30 ± 0.13`, `+0.38 ± 0.10`) that must flatten, since `|T_W| ≤ N_W((2k−1)√p+1)`
   is `O(√p)` at fixed `(m,n,k)`. The empirical mean square of `W_D` over top-size `D` is `0.99p` in every family
   (`k = 2` and `k = 3` alike), as the Katz–Sarnak `USp` equidistribution predicts for the trace (HEURISTIC, not used).

2. **Complete rectangles (bias 1): no cancellation, as the identity forces.** For greedy, gp_biased and all 176
   stored rectangles, `T_W/(mn)^{2k} ≈ 0.8–1.9` (median `0.78` over the stored set, max `4.7`), i.e.
   `T_ns ≈ Σ_{biased t} g(t)^{2k}`, the number of half-biased dilates being `1–25` (gp_biased: `15–25`, as in crux
   Theorem 6.2). Then `|T_ns|/R_D ≈ (mn)^k/(√((2k)!)√p)·#spikes ≈ 3–20` with fitted slope `−0.4` to `−0.6`
   against `p` (the spike is `p`-independent, `R_D ∝ √p`), and `|T_ns|/T_tri ≈ 0.03–0.09` with slope `≈ −0.5`. At
   `p ≤ 3000` and `mn ≤ 36` a single complete dilate exceeds the random-sign scale by `(mn)^k/√p ≈ 8–2800`; the
   crossover `(mn)^k = √p` is at `p ≈ (mn)^{2k}`, far beyond the range.

3. **Ratio-collapsed rectangles (`gp` with a residue ratio): little room for cancellation, and a large
   elementary part.** With only `L = m+n−1 ≤ 11` classes there are `n_D ≤ 385` distinct odd-root sets, each carrying
   a huge positive weight `c(D) = c₊(D)` (all `σ_r` have one sign). Here `|T_ns|/T_tri ≈ 0.2–0.5` (no more than a
   factor `5` below the triangle bound), `|T_W|/R_{DW} ≈ 1.1–1.8` (max `4.8`), the elementary part is of the full
   size `(mn)^{2k}` (median `0.75–1.2`), and `T_ns < 0` in `79 %` of the instances: the moment is up to `25 %` *below*
   the independent-signs value `T_sq ≈ (2k−1)!!(p−1)R_χ^k`, which is itself `≈ 30·(mn)^{2k}` at `(5,5,3)`,
   `p = 2999` (`T_sq = 1.86·10^{10}`, `T_ns = −3.76·10^9`, `(mn)^6 = 2.44·10^8`). This deficit is not asymptotic
   (at fixed sizes `|T_ns| = O(√p)` while `T_sq ∝ p`), but it is what dominates the constants of §5 at these sizes.

### 3.3 The prompt's three ratios, summarised

* `T_ns/T_tri`: random `≈ 0.004–0.02` (`∝ p^0`), complete `≈ 0.03–0.09` (`∝ p^{−1/2}`), collapsed `≈ 0.2–0.5`.
* `T_ns/R` (ordered-tuple scale): random `≈ 3–5` at `k=2`, `45–63` at `k=3` — the factor `√((2k)!)`-plus of §3.1;
  with the pattern scale `R_D` these become `0.6–1.0` at both `k`. Complete: `≈ 7–90` (`k=2`), `140–900` (`k=3`),
  `∝ p^{−1/2}`.
* `T_ns/(mn)^{2k}`: random `≈ ±0.01–0.15` (`∝ p^{+1/2}`), complete `≈ 0.8–1.9` (`∝ p^0`, the spike count),
  collapsed `≈ −1` to `−15` (dominated by the deficit below `T_sq`).

So `|T_ns|` is of order `R_D` (square-root cancellation across *patterns*) for unstructured rectangles, of order
`(#biased dilates)·(ηmn)^{2k}` — between `R_D` and `T_tri`, and equal to `T_tri·(#biased)η^{2k}/((2k−1)√p)` — for
biased ones, and within a factor `2–5` of `T_tri` for maximally collapsed ones. The difference between random and
biased rectangles is therefore exactly the spike, as Theorem 1.5 says it must be; no *additional* cancellation
phenomenon distinguishes them.

---

## 4. `D`-level statistics (Task 3, data)

**Alignment of `W_D` with the weights.** Let `cos(c,W) = Σ_D c(D)W_D / (‖c‖₂‖W‖₂)` over the odd-root sets `D ≠ ∅`
(so `T_W + E₂ = ‖c‖‖W‖·cos`), and `z = cos·√n_D` (the value `≈ ±1` expected for independent signs).

| family | `n` | median `z` | max `z` |
|---|---|---|---|
| random / random2 / interval / gp_nr | 90 / 30 / 60 / 60 | −0.34 / −0.50 / −0.41 / +0.07 | 1.98 / 1.15 / 2.44 / 1.03 |
| gp (residue) | 75 | −0.78 | 1.88 (min −2.68) |
| greedy | 90 | +2.32 | 15.1 |
| gp_biased | 76 | +4.15 | 12.6 |
| stored anneal / k* / exhaustive | 27 / 78 / 71 | +3.61 / +1.44 / +0.05 | 9.1 / 9.7 / 2.9 |

The unbiased families are consistent with independence (`|z| ≤ 2.7` throughout); the complete rectangles show a
coherent component at `2–15` standard deviations, i.e. the spike `(mn)^{2k}` is realised as a *weak alignment spread
over all `n_D ≈ 6·10^4` odd-root sets* (`cos ≈ 0.015` for a `6×6` complete rectangle at `p = 1009`, against
`0.004` for a random one), not as a few large `W_D`. The exhaustive witnesses (mostly `4×4`, `B ⊇ {0,1}`,
`A ⊂ ±Q` by the enumeration's normalisation) show no alignment: there `(mn)^{2k} = 65 536` is below `R_{DW}`.

**Pairwise correlation of overlapping Weil sums.** For `D ≠ D'` of size `2k` with `|D ∩ D'| = 2k−1` (up to `7.5·10^6`
ordered pairs per instance), the Pearson correlation of `(W_D, W_{D'})` over all such pairs, restricted to `L ≥ 16`:

| family | random | random2 | interval | greedy | gp_biased | stored anneal | stored k* | stored exhaustive |
|---|---|---|---|---|---|---|---|---|
| median | −0.0013 | −0.0016 | −0.0015 | −0.0015 | −0.0011 | −0.0027 | −0.0029 | −0.0160 |
| range | ±0.015 | [−0.011, 0] | ±0.012 | ±0.011 | ±0.011 | [−0.016, 0] | ±0.015 | [−0.022, −0.008] |

No family — biased or not — shows a detectable second-order dependence among the `W_D`; the `−0.016` of the
exhaustive witnesses is the `0 ∈ B`, `A ⊂ ±Q` normalisation artefact. The structure that distinguishes a complete
rectangle is *first-order* (alignment with `c(D)`, above), and `c(D)` is where `(ν,σ)` enter; the Weil sums
themselves look like independent signs in every family. (For `L ≤ 11`, i.e. the collapsed `gp` families, the pair
statistic is dominated by finite-size constraints and is not reported.)

---

## 5. One conjecture about `T_ns`, its consequences, and what refutes the alternatives (Task 4)

### 5.1 Statement (HEURISTIC / OPEN)

**Conjecture T(k).** For every integer `k ≥ 2` there is a constant `C_k` such that for every prime `p` and all
`A, B ⊆ F_p` with `0 ∉ A`:

    |T_ns(A,B;k)|  ≤  C_k · √p · (mn)^{2k} / min(m,n)^k .

By Corollary 2.2 this is equivalent (with a different `C_k`) to the same bound for the pure Weil part `T_W` plus the
elementary term `(k+1)C(2k,2)R_ν(mn)^{2k−2}`; the verifier tests both (`rho_T`, `rho_T_W`).

*Reading.* Since `R_ν ≤ mn·min(m,n)`, the right side is `√p·(mn)^k·max(m,n)^k`, i.e. `C_k·√p·R_ν^k·(max/min)^k`.
For a maximally collapsed rectangle (`R_ν ≈ mn·min`, e.g. `gp` or `A = H ⊆ Q`, `B ⊆ cH`) the pattern scale is
`R_{DW} ≈ c_k√p R_ν^k` and T(k) says *square-root cancellation across patterns up to a constant*; for a generic
rectangle (`R_ν ≈ mn`) it is weaker than that by `max(m,n)^k`, and that slack is exactly what allows biased dilates:
T(k) permits `#{t : |g(t)| ≥ ηmn}·η^{2k} ≤ C_k√p/min(m,n)^k`, nothing more. It is *not* a statement about the ratio
set alone — `min(m,n)` is the multiplicative-structure quantity of Remark 2.4.

### 5.2 What it implies (PROVED implications, exact exponents)

**Theorem 5.2.** Assume T(k) for one `k > 1/ε`. Then `C°(ε,δ)` holds for every `0 < δ < ε/2 − 1/(2k)`, with
`|S(A,B)| ≤ (2(2k−1)!!)^{1/(2k)} p^{1/(2k) − ε/2} mn` for `p ≥ p₀(k,ε)`.

*Proof.* For `|A|,|B| > p^ε`, `0 ∉ A`: `T_sq ≤ (2k−1)!!(p−1)(mn)^{2k}/max(m,n)^k < (2k−1)!!p^{1−kε}(mn)^{2k}` and
`T_ns ≤ C_k p^{1/2}(mn)^{2k}/min(m,n)^k < C_k p^{1/2−kε}(mn)^{2k}`. Hence `Σ_{t≠0} g^{2k} < 2(2k−1)!!p^{1−kε}(mn)^{2k}`
once `C_k p^{−1/2} ≤ (2k−1)!!`, and `|S(A,B)|^{2k} ≤ Σ_{t≠0} g^{2k}`. ∎

So T(k) for all `k` gives `C(ε,δ)` for every `δ < ε/2`: for `|A| = |B| = p^ε` this is `|S(A,B)| ≤ p^{o(1)}·mn/√min(m,n)`,
"square-root cancellation in the smaller set". For `m = n` this matches the size `n^{3/2}√(2 log p)` of the greedy
construction (`B` = the `n` best shifts of a random `A`; HEURISTIC), so T(k) is sharp up to `p^{o(1)}` in that
regime if it is true at all. Via crux Corollary 1.2, T(k) implies `n_p ≤ 2⌊p^{ε}⌋` for every `ε > 1/k`, and T(k) for
all `k` implies Vinogradov's conjecture `n_p ≤ p^{o(1)}`; T(k) is therefore at least as hard as Vinogradov.

**Theorem 5.3 (the subgroup case).** Assume T(k). Let `H ≤ F_p^*` be a subgroup of order `h` and `B ⊆ F_p`, `|B| = n`.
Then

    |S(H,B)| ≤ hn · [ (2k−1)!!(p−1)/(h·max(h,n)^k) + C_k√p/(h·min(h,n)^k) ]^{1/(2k)}.

*Proof.* `g(th) = S(thH,B) = S(tH,B) = g(t)` for `h ∈ H`, so `Σ_{t≠0} g^{2k} ≥ h·S(H,B)^{2k}`; bound the left side by
`T_sq + T_ns` as in 5.2 (with `m = h`). ∎

For `h = n = p^{1/3}` and `k = 3`: `|S(H,B)| ≤ (15p^{−1/3} + C_3p^{−5/6})^{1/6}·hn ≈ 1.6·p^{−1/18}·hn` — the open
one-set-a-coset case of the conjecture, with an explicit exponent. T(k) is thus strictly stronger than what the
workspace's subgroup line targets, and it is not implied by `C(ε,δ)` (which says nothing about `|A| ≤ p^{ε}` or
about the constants at small sizes).

### 5.3 Tests (663 instances; PROVED for these instances only)

| bound tested | `k` | instances | max ratio | worst instance | ratio `> 1` |
|---|---|---|---|---|---|
| T(k) on `T_ns` | 2 | 457 | **7.78** | stored `k*` `6×6`, `p = 337`: `S = 36`, 15 half-biased dilates, `R_ν = 96` | 118 |
| T(k) on `T_ns` | 3 | 206 | **35.2** | `gp` `5×5`, `p = 2999`: `T_ns = −0.20·T_sq` | 113 |
| T(k) on `T_W` | 2 / 3 | 457 / 206 | 9.18 / 33.5 | same two | 139 / 113 |
| strong form `C_k√p R_ν^k` on `T_ns` | 2 / 3 | 457 / 206 | 91.8 / 1209 | gp_biased `p=401`; greedy `5×5` `p=151` | 359 / 189 |
| prompt's form `C_k[(mn)^{2k} + √pR_ν^k]` | 2 / 3 | 457 / 206 | 5.76 / 13.5 | `gp` `6×4` and `5×5` at `p = 2999` | 80 / 63 |

For T(k), no family shows a ratio growing with `p`: the fitted slopes of `ρ_T` against `p` are `≤ +0.03` for every
family with `≥ 11` instances except `gp(5,5,3)` (`+0.25`) and `random2(6,6,2)` (`+0.30`), single-instance fits whose
siblings (`gp(4,4,3)`: `+0.03`, `random(6,6,2)`: `−0.13`) are flat; the complete rectangles have slope `−0.4` to `−1.1`
(spike fixed, bound `∝ √p`). Random rectangles sit below the bound by a factor `≈ max(m,n)^k` (`ρ ≤ 1.4`). The two
mechanisms that set the constants are (i) several biased dilates at small `p` (`ρ ≈ #spikes·min^k/√p`, decreasing in
`p`) and (ii) the ratio-collapsed `gp` rectangles, where `ρ ≈ (|T_W|/R_{DW})·√((2k)!)·(R_ν min/(mn)²)^k ≈ 3·27·0.3 ≈ 30`
at `(5,5,3)`; both are `p`-bounded for fixed sizes. HEURISTIC: the greedy construction forces
`C_k ≥ sup_p (2 log p)^k/√p ≈ (4k/e)^k` (`8.7`, `86` at `k = 2, 3`), consistent with the observed `7.8`, `35`.

**The critical untested regime.** T(k) is tight, up to `C_k`, exactly where the data are thinnest: collapsed
rectangles with `m = n = K → ∞` (`gp` with `K ≥ 8`, `A = B = H ⊆ Q`), and greedy rectangles with `n ≈ p^{1/k}`. A
ratio growing along any such family would refute T(k); §6.

### 5.4 The alternatives (REFUTED, or vacuous)

* **Strong form `|T_ns| ≤ C_k√p R_ν^k` (square-root cancellation across patterns for all rectangles): REFUTED.**
  (a) Every complete rectangle with `(mn)^k > C_k√p` violates it (spike); e.g. greedy `5×5` at `p = 151`, `k = 3`:
  ratio `1209`. (b) A *flat* profile violates it in the opposite direction: `A = {1}`, `B = Q` has `g(t) = −(1+χ(t))/2`
  for `t ≠ 0`, so `Σ_{t≠0} g^{2k} = (p−1)/2` and `T_ns = (p−1)/2 − T_sq` with `T_sq = 3(p−1)n² − …` (`k = 2`, `n = (p−1)/2`):
  exactly `T_ns = −153 120, −725 200, −2 475 300` at `p = 61, 101, 151`, ratios `21.8, 28.9, 35.8 ≈ (2k−1)!!√p → ∞`
  (`78.4, 102.4` at `k = 3`, `p = 41, 61`). The lower side `T_ns ≥ −T_sq ≈ −(2k−1)!!pR_χ^k` rules out any bound of size
  `√p·R_ν^k` uniformly; the `min(m,n)^{−k}` of T(k) is what absorbs this (`ρ_T ≤ 0.024` on these witnesses).
* **The prompt's form `|T_ns| ≤ C_k(mn)^{2k} + C_k√pR_ν^k`: consistent with all data (`C_2 ≤ 5.8`, `C_3 ≤ 13.5`;
  with `C_k = 1` it is violated by 80 + 63 instances, all `gp`), but it implies nothing.** In the regime of Theorem 1.5
  the term `(mn)^{2k}` is the trivial bound `|S| ≤ mn` itself: `Σ_{t≠0} g^{2k} ≤ T_sq + C_k(mn)^{2k} + …` only yields
  `|S(A,B)| ≤ (C_k + o(1))^{1/(2k)}mn`. Any conjecture that is to imply `C(ε,δ)` through Theorem 1.5 must have *every*
  term `o((mn)^{2k})` when `m, n > p^{ε}`; the `min(m,n)^{−k}` factor of T(k) is the minimal such correction consistent
  with the lopsided witnesses (`m = 1` forces a bound `≥ (2k−1)!!p·n^k`, which `√p·n^{2k}/1` provides).
* **`max(m,n)` in place of `min(m,n)`: REFUTED** by the same `A = {1}`, `B = Q` witnesses (it coincides with the strong
  form there).

---

## 6. Remaining obligations (OPEN)

1. **Any power saving over the triangle inequality.** By Theorem 1.5 the conjecture *is* the statement
   `T_ns ≤ p^{−2kδ'}(mn)^{2k}` for one `k`, i.e. a saving of `(2k−1)p^{1/2+2kδ'}` over `T_tri` for the signed sum of
   the one-variable Weil sums `W_D`, `D ⊆ (−B)·A^{−1}`, with the explicit weights `c(D)`. Nothing in this note or
   in the workspace proves a saving of `p^{θ}` for any `θ > 0` in any regime with `m, n → ∞`; the identities are
   exact, the inequalities are termwise Weil. The data say the saving is `√p`-type (`T_ns ≍ R_D`) for unstructured
   rectangles and *zero* for complete ones — which is the conjecture restated, not evidence for a mechanism.
2. **Test T(k) where it is tight**: collapsed rectangles with `m = n = K = 8,…,16` at `p ≈ 10^4–10^5` (cheap in `D`,
   `n_D = O(K^{2k})`; the cost is the `p×L` character table), the subgroup rectangles `A = H ⊆ Q`, `B = H` (maximal
   collapse *and* `|H|` biased dilates), and greedy rectangles with `n ≈ p^{1/k}`. Report `ρ_T` against `p` and `K`;
   growth along any family refutes T(k).
3. **The elementary part.** Prove the closed form `E₂ = −6R_χ(g(0)² + g(∞)²) + 12R_χ² + O(Σ_r σ_r⁴)` at `k = 2`
   (`g(∞) = n·Σ_a χ(a)`), its `k ≥ 3` analogue, and whether the bound of Corollary 2.2 is attained only for
   `A, B ⊆ ±Q` with `R_ν ≈ mn·min(m,n)`; a reformulation over the projective line (dilates `t = 0` and `t = ∞`)
   should absorb both degenerate dilates and the elementary part together.
4. **T(k) versus SI.** T(k) bounds the number of dilates of bias `≥ ½` by `C_k4^k√p/min(m,n)^k`, which is
   incomparable with crux Conjecture SI (`8 log₂ p`): weaker for `min(m,n) ≤ p^{1/(2k)}`, stronger beyond. Whether
   SI plus T(k) at small sizes is consistent with Theorem 6.2's `≥ ¼log₂p` biased dilates at `m ≈ ½log₂p`,
   `n ≈ √p log p` (it is: `m·min^k/√p → 0`) is checked above only heuristically.
5. **Exactness of the equivalence under the prompt's `F_p` convention.** Theorem 1.4(ii) fails with `Σ_{t∈F_p}`
   (`B ⊆ Q` gives `g(0) = mn`); if a full-`F_p` statement is wanted, `(m·Σ_bχ(b))^{2k} − Q_k(0)` must be carried as an
   explicit term (the verifier records both).
