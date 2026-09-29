# Sigma pass, direction `stress`: breaking Conjecture T(k) and Conjecture SI\* where they should be tight (2026-09-05)

**Status: REFUTED — (1) Conjecture T(k) of `research/sigma-tuple-2026-09-05.md` §5.1 (`|T_ns(A,B;k)| ≤ C_k √p (mn)^{2k}/min(m,n)^k`), for every `k ≥ 2` and every constant `C_k`, by two exact families: `A = B = Q` (the nonzero squares), where `T_ns = h^{2k}·Δ_k(Q)` with the closed form `Δ_k(Q) = h[1 − N_k(h) − N_k(h−1)]`, `h = (p−1)/2`, so that the ratio `ρ_T := |T_ns|·min^k/(√p (mn)^{2k})` equals `2(h−1)(3h−2)/(h√p) = 2.9968√p` at `p = 4001`, `k = 2` and tends to `(2k−1)!!·√p`; and `A = B = H` for a subgroup `H ⊆ Q` of order `h`, where `ρ_T` reaches `152.3` (`k = 2`, `p = 2833`, `h = 59`) and `2753` (`k = 3`, same pair) already at `p ≤ 3000`, `785.9` / `22 486` at `p = 13 183`, `h = 169`, and `334.8` / `7033` at `p ≈ 10^7`, `h = 724, 1024`, against the tuple note's observed `C_2 ≤ 7.8`, `C_3 ≤ 36`; the median of `ρ_T` over subgroups grows like `h^{0.65}` (`p ≤ 20 000`; `h^{0.52}` in the from-scratch `p ≤ 3000` scan) and the mean over the targeted family like `h^{0.54}` (`h = 16…1024`, `p ≤ 1.7·10^7`), flat in `p` at fixed `h`. (2) Conjecture SI\* of `research/sigma-si-2026-09-05.md` §4.2 (`|S| ≥ ½|A||B|` and `S² ≥ 8R ln p` imply `|D_{1/2}| ≤ 2 log₂ p`), by explicit "designed" `1×n` rectangles `A = {1}`, `B` = the `n` best shifts for a weighted random target set: at `p = 4001` (`n = 332`, `S = 166 = n/2`, `z² = 10.0 ln p`, `|D_{1/2}| = 32 > 23.9`), `p = 10 007` (`|D| = 47 = 3.54 log₂ p`), `p = 100 003` (`70 = 4.21 log₂ p`), `p = 1 000 003` (`124 = 6.22 log₂ p`), the ratio growing with `p`; the strongest reading (bias `1` at `t = 1`, `B ⊆ N(1)`, absolute threshold `|g(t)| ≥ ½mn`) also fails: `|D| = 42 = 2.53 log₂ p` (`p = 100 003`), `45 = 2.26 log₂ p` (`10^6`), `51 = 2.19 log₂ p` (`10^7`). PROVED — the moment route `T_ns = Σ_{t≠0} g(t)^{2k} − T_sq` with `T_sq` in closed form (which makes `T_ns` computable in `O(p·L)` instead of `O(p·L^{2k})`, validated against brute force and against the tuple verifier's `D`-decomposition), the subgroup reduction `T_ns(cH⁺, H; k) = |H⁺|^{2k}·Δ_k(H)`, `Δ_k(H) = Σ_{s≠0} T_H(s)^{2k} − (p−1−h)N_k(h) − hN_k(h−1)`, the closed forms for `Q` and for the elementary part `Δ_2(H) = 24Σ_{|D|=4}W_D − 2h(h−1)(6h−11) − 6hT_H(−1)²`, and the fact that the corrected statement T′(k) (exponent `k − ½` on `min(m,n)`) implies `C(ε,δ)` with exactly the exponents of tuple Theorems 5.2–5.3. SUPPORTED, not proved — T′(k) up to `p^{o(1)}`: `ρ′ := |T_ns|·min^{k−½}/(√p (mn)^{2k})` has rms `5.7` (`k = 2`) / `80` (`k = 3`) over the 6938 subgroup rectangles with `8 ≤ h ≤ √p` and maximum `60.5` / `1730` over all 8894 with `h ≥ 8`, the maximum being a single shifted-subgroup spike `M_H = 65 = 5.0√h` at `p = 13 183`, `h = 169`; the exponent `k − ½` is sharp (both refuting families have `ρ_T ∝ √min`). OPEN — T′(k) itself; any replacement of SI\* (the capacity heuristic of §2.4 says `|D_{1/2}| ≤ (32 ln 2 + o(1))·log₂ p ≈ 22 log₂ p` is the most that can be hoped for under the SNR hypothesis, and `≈ 8 ln 2·log₂ p ≈ 5.5 log₂ p` under the strongest reading).**

Worker `stress`, prefix `sigma`, 2026-09-05. Verifier `experiments/sigma_stress_2026_09_05.py` (standard library + numpy; imports the tuple verifier for the cross-check; 4626 exact checks, 0 failures, 192 s) → `results/sigma_stress_2026_09_05.json`. Heavy runs are stored in `results/sigma_stress_2026_09_05_search.json` (every stored witness is re-derived exactly by the verifier from `(p, A, B)` or `(p, h)`; the smallest SI\* witness also in pure Python). Notation is that of the tuple and SI notes: `χ` the Legendre symbol, `A, B ⊆ F_p`, `0 ∉ A`, `m = |A|`, `n = |B|`, `g(t) = Σ_{a,b} χ(ta + b)`, `S = g(1)`, `T_H(s) = Σ_{h∈H} χ(s + h)`, `M_H = max_{s≠0}|T_H(s)|`, `Q` the nonzero squares, `N(U) = {b ≠ 0 : χ(u + b) = 1 ∀u ∈ U}`, `D_{1/2} = {t ≠ 0 : 2|g(t)| ≥ |S|}`, `R` the ratio energy, `z² = S²/R`. All `t`-sums are over `F_p^*` (the tuple note's convention).

---

## 0. The computational route (PROVED, validated)

The tuple note computes `T_ns` through its `D`-decomposition, which costs `O(p·L^{2k})` (`L` = number of ratio classes) and confined its tests to `mn ≤ 36`. Its own Proposition 1.2 gives a route that costs `O(p·L)`:

**Lemma 0.1.** With `N_k(L′) := (2k)!·[x^{2k}] cosh(x)^{L′}` (the number of words of length `2k` over `L′` letters in which every letter occurs an even number of times; `N_2(L′) = 3L′² − 2L′`, `N_3(L′) = 15L′³ − 30L′² + 16L′`),

    T_ns(A,B;k) = Σ_{t∈F_p^*} g(t)^{2k} − T_sq,   T_sq = (p − 1 − |𝓡∖{0}|)·(2k)![x^{2k}]Π_{r∈𝓡} cosh(σ_r x) + Σ_{r∈𝓡, r≠0} (2k)![x^{2k}]Π_{r′≠r} cosh(σ_{r′} x).

*Proof.* `Σ_{t≠0} g^{2k} = T_sq + T_ns` by definition and `T_sq = Σ_{t≠0} Q_k(t)` with `Q_k(t) = (2k)![x^{2k}]Π_{r≠t}cosh(σ_r x)` (tuple Prop. 1.2); `Q_k(t)` takes the "full" value for `t ∉ 𝓡` and the value with the factor `r = t` removed for `t ∈ 𝓡∖{0}`. ∎

The verifier computes `Σ_t g^{2k}` from the value histogram of `g` and the cosh-coefficients with `Fraction`s, all in Python integers. Validation: (i) against a brute-force enumeration of all `(mn)^{2k}` tuples (own code, 15 cases at `p ≤ 13`); (ii) against `analyze()` of `experiments/sigma_tuple_2026_09_05.py` (the independent `D`-decomposition code path) on 98 rectangles at `p ≤ 401` including subgroup and gp rectangles, `k = 2, 3`; (iii) `g` from class data against the direct double sum. All agree exactly.

**Lemma 0.2 (subgroup reduction; PROVED).** Let `H ≤ F_p^*` have order `h`, `H⁺ = H ∩ Q`, `c ∈ F_p^*`, `A = cH⁺`, `B = H`, `m = |H⁺| ∈ {h, h/2}`. Then `g(t) = m·T_H(ct)`, and

    T_ns(cH⁺, H; k) = m^{2k}·Δ_k(H),    Δ_k(H) := Σ_{s∈F_p^*} T_H(s)^{2k} − (p − 1 − h)·N_k(h) − h·N_k(h − 1).

*Proof.* `χ(ta + b) = χ(a)χ(t + b/a)` and `b/a` runs over `c^{−1}H` as `b` runs over `H`, so `g(t) = (Σ_{a∈A}χ(a))·Σ_{u∈H}χ(t + u/c) = χ(c)m·χ(c)T_H(ct) = m·T_H(ct)`. The ratio classes are `r = −b/a ∈ −c^{−1}H`, each with `ν_r = m` and `σ_r = χ(c)m`, so Lemma 0.1 gives `T_sq = m^{2k}[(p−1−h)N_k(h) + hN_k(h−1)]` (`0 ∉ 𝓡`), while `Σ_{t≠0} g^{2k} = m^{2k}Σ_{s≠0}T_H(s)^{2k}`. ∎

In particular `T_ns(H,H;k) = T_ns(H,cH;k) = T_ns(cH,H;k)` for every `c` when `H ⊆ Q`: the SI note's `p = 97` witness `(2H, H)` and "`A = H`, `B = cH` at the extremal shift" are the same instance of `Δ_k(H)` as `A = B = H`. The verifier checks Lemma 0.2 on 110 `(p, H, c, k)` and the coset form `Σ_{s≠0}T_H^{2k} = h·Σ_{j<d}T_H(g^j)^{2k}` (`d = (p−1)/h`), which is what makes `p ≈ 10^7` feasible.

**Lemma 0.3 (elementary part for `A = B = H ⊆ Q`, `k = 2`; PROVED).** With `W_D = Σ_{s≠0}Π_{r∈D}χ(s + r)`,

    Δ_2(H) = 24·Σ_{D⊆H, |D|=4} W_D − 2h(h−1)(6h−11) − 6h·T_H(−1)².

*Proof.* For `s ∉ −H ∪ {0}` put `x_u = χ(s+u) ∈ {±1}` (`u ∈ H`); then `p_1 = p_3 = T_H(s)`, `p_2 = p_4 = h` for the power sums and Newton's identity gives `T⁴ = 24e_4 + (6h−8)T² − 3h² + 6h`; for `s ∈ −H` one `x_u` vanishes, `p_2 = p_4 = h − 1` and `T_H(s) = T_H(−1)` (`T_H(−u) = χ(u)T_H(−1)`). Sum over `s ≠ 0`, use `Σ_{s≠0}T_H(s)² = hp − 2h²` (Chung's identity minus `T_H(0)² = h²`) and subtract `T_sq` of Lemma 0.2 with `N_2(L′) = 3L′² − 2L′`. ∎ (Checked exactly for 20 subgroups with `4 ≤ h ≤ 24`, `p ≤ 1009`.)

So `|Δ_2(H) − 24Σ_{|D|=4}W_D| ≤ 18h³`, which is `≤ 36 h/√p · (√p h²)`: the elementary part is within the T(k) scale for `h ≤ 2√p` and the refutation below is carried by the Weil part.

---

## 1. Conjecture T(k) is false (REFUTED)

Recall (tuple note §5.1): **T(k)**: `|T_ns(A,B;k)| ≤ C_k·√p·(mn)^{2k}/min(m,n)^k` for all `p` and all `A, B ⊆ F_p`, `0 ∉ A`. Write `ρ_T = |T_ns|·min^k/(√p(mn)^{2k})`; the tuple note observed `ρ_T ≤ 7.78` (`k = 2`) and `≤ 35.2` (`k = 3`) on 663 instances with `mn ≤ 36`.

### 1.1 `A = B = Q`: the ratio is `(2k−1)!!·√p` (PROVED)

**Proposition 1.1.** Let `h = (p−1)/2` and `A = B = Q`. Then for every `k ≥ 1`

    T_ns(Q,Q;k) = h^{2k}·Δ_k(Q),   Δ_k(Q) = h − (p−1−h)N_k(h) − hN_k(h−1) = h[1 − N_k(h) − N_k(h−1)],

in particular `Δ_2(Q) = −2h(h−1)(3h−2)` and `Δ_3(Q) = −h(30h³ − 105h² + 137h − 62)`, and

    ρ_T(Q,Q;2) = 2(h−1)(3h−2)/(h√p) = 3√p·(1 + O(1/p)),    ρ_T(Q,Q;k) = (2k−1)!!·√p·(1 + O(1/p)).

*Proof.* `T_Q(s) = ½Σ_{x≠0}(1 + χ(x))χ(s + x) = ½(−χ(s) − 1)` for `s ≠ 0` (using `Σ_{x}χ(x)χ(x+s) = −1`), so by Lemma 0.2 `g(t) = hT_Q(t) = −h·1[t ∈ Q]` and `Σ_{t≠0}g^{2k} = h^{2k+1}`; `Δ_k` is Lemma 0.2 with `p − 1 − h = h`. `N_k(h) = (2k−1)!!h^k(1 + O(1/h))` gives the asymptotics; for `k = 2` the polynomial identity is elementary. ∎

Exact values (verifier `T_k_refutation_A_eq_B_eq_Q`; the exact inequality `|T_ns|²·min^{2k} > C²·p·(mn)^{4k}` is tested in integers with `C = 7.8`, `36`):

| `p` | `h` | `Δ_2(Q)` | `ρ_T` (`k=2`) | `ρ_T/√p` | `ρ_T` (`k=3`) | `ρ_T/√p` |
|---|---|---|---|---|---|---|
| 61 | 30 | −153 120 | 21.8 | 2.789 | 102.4 | 13.1 |
| 401 | 200 | −47 600 800 | 59.4 | 2.968 | 294.4 | 14.70 |
| 1009 | 504 | −765 606 240 | 94.9 | 2.987 | 472.7 | 14.88 |
| 4001 | 2000 | −47 960 008 000 | 189.6 | 2.9968 | 946.9 | 14.97 |
| 10007 | 5003 | −751 100 530 084 | 300.0 | 2.9987 | 1499.3 | 14.99 |

(`Δ_2(61) = −153 120` is the value the tuple note found for its `A = {1}`, `B = Q` witness — the same profile up to the factor `h^{2k}`; there `min = 1` absorbed it, here `min = h` does not.) The index-4 subgroup `Q_4` behaves the same way (`ρ_T = 173` at `p = 4001`, `k = 2`). So T(k) fails by `√p` in the large-set range; §1.2 shows it also fails by `√min(m,n)` in the range `m, n ≤ √p` that the conjecture is about.

### 1.2 `A = B = H`, `H ⊆ Q` a subgroup: the ratio grows like `√h` (REFUTED with exact witnesses)

By Lemma 0.2, `ρ_T(H,H;k) = |Δ_k(H)|/(√p·h^k)` and `ρ′ := ρ_T/√h = |Δ_k(H)|/(√p·h^{k+½})`. The verifier recomputes `Δ_k(H)` for **all** 1451 subgroups `H ⊆ Q` with `4 ≤ h ≤ 2√p` of all primes `p ≤ 3000` (full profile, 1 s), re-verifies a random sample of 250 of the 10 809 stored rows at `p ≤ 20 000` and all 136 stored rows of the targeted family (`h ∈ {16, 23, 32, 45, 64, 91, 128, 181, 256, 362, 512, 724, 1024}`, four primes `p ≡ 1 (mod 2h)` near each of `p ≈ h², h^{2.5}, h³`, `p ≤ 1.7·10^7`; coset method, cross-checked against the full profile for `p ≤ 2·10^5`).

**Witnesses (exact; each violates T(k) with the observed constants by the integer test):**

| `p` | `h` | `d` | `Δ_2(H)` | `ρ_T` (`k=2`) | `ρ′` | `Δ_3(H)` | `ρ_T` (`k=3`) | note |
|---|---|---|---|---|---|---|---|---|
| 2833 | 59 | 48 | 28 225 364 | 152.3 | 19.8 | 30 097 673 186 | 2753 | worst at `p ≤ 3000` |
| 13 183 | 169 | 78 | 2 577 125 616 | 785.9 | 60.5 | 12 461 603 622 936 | 22 486 | `M_H = 65 = 5.0√h` (one coset), worst `ρ′` |
| 38 011 | 181 | 210 | 2 128 060 440 | 333.2 | 24.8 | 5 486 248 729 260 | 4746 | targeted, `p ≈ h²` |
| 1 073 153 | 1024 | 1048 | 259 500 011 520 | 238.9 | 7.5 | 7 822 477 773 305 856 | 7033 | targeted, `p ≈ h²` |
| 14 109 313 | 724 | 19 488 | −659 287 950 384 | 334.8 | 12.4 | −6 662 500 737 718 104 | 4674 | targeted, `p ≈ h^{2.5}` |

**Growth.** Over the 6803 subgroup rectangles with `8 ≤ h ≤ √p`, `p ≤ 20 000` (`k = 2`): median `ρ_T = 9.9, 17.7, 28.0, 41.1` in the bands `h ∈ [8,15], [16,31], [32,63], [64,127]`, with band maxima `65.8, 188.8, 202.9, 304.3`; the log–log slope of the median of `ρ_T` against `h` is `+0.65` (109 values of `h`), while at fixed `h ∈ {8, 12, 16, 24, 32, 48}` the slope against `p` is between `−0.22` and `−0.02`. For `k = 3` the medians are `110, 220, 348, 509` (slope `+0.66`), maxima up to `6618`. In the targeted family the mean of `ρ_T` (`k = 2`) is `17.9, 25.4, 43.1, 52.2, 113.5, 136.2, 147.3` at `h = 16, 32, 64, 128, 256, 512, 724` (slope `+0.55`). The normalised ratio `ρ′ = ρ_T/√h` is flat: band medians `2.98, 3.81, 4.26, 4.50` (slope `+0.15`), rms `5.7` over the 6938 rows with `h ≤ √p`, and mean `3.8–7.7` at every `h` of the targeted family up to `1024`; for `k = 3` the `ρ′` medians are `34, 47, 53, 55`, rms `80`. In the range `√p < h ≤ 2√p` the same holds (medians of `ρ′`: `4.0–5.6`). The sign of `Δ_k` is negative in `58–81 %` of the cases (the elementary part `−12h³` of Lemma 0.3 and the variance deficit `Σ_{s≠0}T_H² = hp − 2h²`).

So: **for the subgroup family `ρ_T ≍ √h` up to bounded factors**, and no constant `C_k` works. Since `h` can be `p^ε` for any `ε` (Dirichlet), T(k) fails by a factor `p^{ε/2}` inside the Paley regime `|A| = |B| = p^ε` — exactly where Theorem 5.3 of the tuple note wanted to apply it.

**HEURISTIC explanation (coherence).** For `A = B = H` the profile is `g(t) = hT_H(t)`, constant on the `d = (p−1)/h` cosets of `H`, so `Σ_{t≠0}g^{2k} = h^{2k+1}Σ_{j<d}T_H(g^j)^{2k}` is `h` times a sum of only `d` terms. If `T_H(g^j)/√h` behaves like a standard Gaussian across cosets (the dual note's model, §3.4 there), the `2k`-th empirical moment fluctuates by `≈ √d·h^k·σ_{2k}` (`σ_4 = √96`, `σ_6 = √10170`), giving `|Δ_k| ≈ h·√d·h^k·σ_{2k} = σ_{2k}√p·h^{k+½}` and `ρ_T ≈ σ_{2k}√h`, `ρ′ ≈ σ_{2k}`. The observed `ρ′` (rms `5.7` and `80`; `p ≤ 3000` band medians `3.7–4.6` and `34–42`) are smaller than `σ_{2k}` (`9.8`, `101`) — the coset values are not independent (`Σ_j T_H(g^j)²` is fixed exactly) — but the `√h` law is what the data show. The tuple note's T(k) is "square-root cancellation across patterns" for collapsed rectangles; the `H`-symmetry forces the `W_D` to be constant on `H`-orbits of odd-root sets, which removes a factor `√h` of cancellation. The same mechanism, in weaker form, acts on the gp rectangles (§1.4: the profile is smooth along the orbit of the ratio with correlation length `K`).

### 1.3 The two regimes are one statement

Both refutations are the failure of the "independent-signs" value `T_sq ≈ (2k−1)!!·p·R_χ^k` to predict the dilation moment when the profile has a symmetry: for `Q` the profile takes two values and the moment is `h^{2k+1}` instead of `≈ 3ph^{3k}`; for `H` the moment is a `d`-term sum. In both cases `|T_ns| ≈ √(p·h)·R_χ^k·(bounded)` with `R_χ = h³` (`H ⊆ Q`), against T(k)'s `√p·R_χ^k` — the extra `√h = √min(m,n)`.

### 1.4 The other families (T(k) holds with the observed constants, or fails mildly)

All rows are exact (moment route; the verifier re-verifies every stored row with `p ≤ 40 009` and the `Q` rows). Maximum of `ρ_T` and of `ρ′` per family:

| family | sizes, primes | `k` | instances | max `ρ_T` (witness) | max `ρ′` | `ρ_T > C_obs` |
|---|---|---|---|---|---|---|
| gp, same residue ratio, `m = n = K` | `K ∈ {4,…,40}`, `p ∈ [1009, 10^6]` | 2 | 110 | 33.3 (`p = 20011`, `K = 24`, `T_ns = 899 437 307 200`) | 6.8 | 28 |
| | | 3 | 110 | 365 (same) | 74.6 | 52 |
| gp, `(K, 2K)` | `K ≤ 16`, `p ≤ 4·10^5` | 2 / 3 | 75 / 75 | 11.2 / 69.8 (`p = 400 009`, `12×24`) | 3.2 / 20.2 | 3 / 4 |
| greedy, `n ≈ p^{1/k}`, `m ∈ {n/2, n, 2n}` | `p ∈ [1009, 10^5]` | 2 | 8 | 0.78 | 0.13 | 0 |
| | | 3 | 21 | 9.5 (`p = 1009`, `10×10`) | 3.0 | 0 |
| unions of two cosets `(H∪uH, H)`, `(H, H∪uH)`, `(H∪uH, H∪uH)` | `6 ≤ h ≤ √p`, `p ≤ 6000`, 314 each | 2 | 942 | 14.8 / 14.8 / 58.9 (`p = 4001`, `h = 25`) | 3.0 / 3.0 / 8.3 | 43 / 43 / 84 |
| | | 3 | 942 | 132 / 132 / 695 (`p = 4463`, `h = 23`; `p = 4001`, `h = 25`) | 27.6 / 27.6 / 98.3 | 66 / 66 / 119 |
| random `K×K` | `K ∈ {8,16,32}`, `p ≤ 40009` | 2 / 3 | 11 / 11 | 0.14 / 0.21 | 0.05 / 0.07 | 0 |
| `(cH⁺, H)`, `p = 97` witness and relatives | `h ∈ {8,16,32}`, `c ∈ {2,3,5}`, `p ≤ 7681` | 2 / 3 | 87 / 87 | 34.1 / 342 (`p = 3137` resp. `7681`, `h = 32`) | 7.6 / 67.5 | 60 / 66 |
| `A = B = Q` | `p ≤ 10007` | 2 / 3 | 9 / 9 | 300.0 / 1499.3 | 4.24 / 21.2 | all |
| `A = B = H` (§1.2) | `p ≤ 3000` from scratch (`h ≥ 4`) | 2 / 3 | 1451 / 1451 | 152.3 / 2753 | 19.8 / 358 | 926 / 1136 |
| `A = B = H` (§1.2) | stored `p ≤ 20 000` + targeted | 2 / 3 | 10 945 / 10 945 | 785.9 / 22 486 | 60.5 / 1730 | most |

The gp family with equal ratio shows the same `√K` drift (the rms of `ρ_T` over the two instances per `(p, K)`, pooled over `p`, is `2.0, 2.8, 5.3, 7.7, 7.0, 13.4, 14.3, 12.2` at `K = 4, 6, 8, 12, 16, 24, 32, 40`, log–log slope `0.85`, noisy), reaching `ρ_T = 33 > 7.8`; the greedy family is harmless at these sizes (`ρ_T ≈ #spikes·(2 ln(p/n))^k/√p → 0`, as in the tuple note); unions of two cosets behave like single cosets with `min = 2h` or `h`.

### 1.5 What survives, and what it still implies (CONDITIONAL, exact)

The exponent of `min(m,n)` is pinned from both sides: for `A = B = H`, `ρ_T·min^{−θ}` is flat in `h` only for `θ = ½` (slope `+0.15` at `θ = ½`, `+0.65` at `θ = 0`); for `A = B = Q`, `|T_ns|·min^{k−θ}/(√p(mn)^{2k}) = (2k−1)!!√p·h^{−θ}(1+o(1)) → ∞` for every `θ < ½` and `→ (2k−1)!!√2` at `θ = ½`. The surviving candidate is therefore

> **Conjecture T′(k).** For every `k ≥ 2` and `η > 0` there is `C_{k,η}` such that for all `p` and all `A, B ⊆ F_p`, `0 ∉ A`: `|T_ns(A,B;k)| ≤ C_{k,η}·p^{½+η}·(mn)^{2k}/min(m,n)^{k−½}`.

The `p^η` is not decorative: `|Δ_k(H)| ≥ h·M_H^{2k} − …`, so a single coset with `M_H ≈ 5√h` (the `p = 13 183` witness) gives `ρ′ = 60`, and if `M_H ≈ √(2h log d)` (dual note's model) then `ρ′ ≲ (2 log p)^k√(h/p)`, bounded for `h ≤ p^{1−ε}` but not by an absolute constant; the Gaussian-fluctuation model likewise makes `sup ρ′ = ∞` at the rate `√log`. **T′(k) is SUPPORTED**: `ρ′ ≤ 60.5` (`k = 2`) and `≤ 1730` (`k = 3`) on all 10 945 subgroup rectangles (`h ≥ 4`, `p ≤ 1.7·10^7`; the maximum at `p = 13 183`, `h = 169 = 1.47√p`), `≤ 6.8 / 75` on gp, `≤ 3.0` on greedy, `≤ 8.3` on unions of two cosets, `4.24 / 21.2` on `Q`; nothing tested grows faster than `p^{o(1)}` at `θ = ½`.

**Theorem 1.2 (T′(k) has the same consequences as T(k); PROVED).** Assume T′(k) for one `k > 1/ε` and some `η ≤ (1−ε)/2`. Then `C°(ε,δ)` holds for every `0 < δ < ε/2 − 1/(2k)`, with `|S(A,B)| ≤ (2(2k−1)!!)^{1/(2k)}p^{1/(2k)−ε/2}mn` for `p ≥ p₀(k,ε,η)`; and for a subgroup `H` of order `h` and `|B| = n`,

    |S(H,B)| ≤ hn·[ (2k−1)!!(p−1)/(h·max(h,n)^k) + C_{k,η}p^{½+η}/(h·min(h,n)^{k−½}) ]^{1/(2k)}.

*Proof.* As in tuple Theorems 5.2–5.3: `T_sq ≤ (2k−1)!!p^{1−kε}(mn)^{2k}` and now `T_ns ≤ C p^{½+η−(k−½)ε}(mn)^{2k}`. The second exponent is at most the first iff `½ + η − (k−½)ε ≤ 1 − kε`, i.e. `η ≤ (1−ε)/2`, and the first is negative since `k > 1/ε`; so `Σ_{t≠0}g^{2k} ≤ 2(2k−1)!!p^{1−kε}(mn)^{2k}` for large `p`, and `|S|^{2k} ≤ Σ_{t≠0}g^{2k}`. The subgroup bound is the same computation with `Σ_{t≠0}g^{2k} ≥ h·S(H,B)^{2k}`. ∎

For `h = n = p^{1/3}`, `k = 3`, `η = 0`: `|S(H,B)| ≤ (15p^{−1/3} + Cp^{−2/3})^{1/6}hn` — the same `1.6·p^{−1/18}·hn` as before. The loss of `√min(m,n)` costs nothing in the exponents of the tuple note because there the square part `T_sq ∝ p` dominates the non-square part `∝ √p` by `p^{½−ε/2}`; the sharpness remark of tuple §5.2 (`|S| ≤ p^{o(1)}mn/√min(m,n)` for `k → ∞`, matching the greedy size `n^{3/2}`) also survives, because it comes from the square part `T_sq ∝ p(mn)^{2k}/max^k`, which is untouched.

---

## 2. Conjecture SI\* is false (REFUTED)

Recall (SI note §4.2): **SI\***: for all `p` and `A, B ⊆ F_p^*` with `|S(A,B)| ≥ ½|A||B|` and `S(A,B)² ≥ 8·R(A,B)·ln p`: `|D_{1/2}(A,B)| ≤ 2 log₂ p`. The SI note supports it with `≈ 371 000` high-SNR rectangles from six structured sources, maximum `1.07 log₂ p`.

### 2.1 Designed `1×n` rectangles (PROVED by exact computation)

Take `A = {1}` (then `R = n` since all ratios `b/1` are distinct, `S = T_B(1)`, `g(t) = T_B(t)` and `D_{1/2} = {x ≠ 0 : 2|T_B(x)| ≥ |T_B(1)|}`), `n = 32 ln p` rounded up to a multiple of 4 (so that `S = n/2` gives `z² = n/4 = 8 ln p` exactly at the threshold), and choose `B` to make many `x` biased: pick a random target set `X′ ⊆ F_p^*∖{1}` of size `M′`, score every `b ∉ −X′ ∪ {−1, 0}` by `s(b) = w·χ(1 + b) + Σ_{x∈X′}χ(x + b)` with anchor weight `w = 2`, let `B` be the `n` largest scores, and then improve `B` by swaps chosen on the near-threshold set (each swap re-evaluated exactly, the hypothesis enforced at every step). The rectangle satisfies the SI\* hypothesis by construction and `D_{1/2} ⊇` most of `X′`:

| `p` | `n` | `S` | bias | `z²/ln p` | `|D_{1/2}|` | `2 log₂ p` | `|D|/log₂ p` |
|---|---|---|---|---|---|---|---|
| 1009 | 224 | 112 | 0.500 | 8.68 | 15 | 20.0 | 1.50 (no violation) |
| 2003 | 244 | 122 | 0.500 | 8.02 | **23** | 21.9 | **2.10** |
| 4001 | 332 | 166 | 0.500 | 10.0 | **32** | 23.9 | **2.67** |
| 10 007 | 296 | 152 | 0.514 | 8.47 | **47** | 26.6 | **3.54** |
| 100 003 | 368 | 204 | 0.554 | 9.82 | **70** | 33.2 | **4.21** |
| 1 000 003 | 444 | 226 | 0.509 | 8.33 | **124** | 39.9 | **6.22** |

(Verifier `SI_star_refutation_designed`: every `B` is stored, `T_B`, `S`, `R = n`, the hypothesis and `|D|` are recomputed exactly; the `p = 2003` witness also in pure Python with `pow(·,(p−1)/2,p)`.) The ratio `|D_{1/2}|/log₂ p` **grows** with `p`; no constant `K` in `|D_{1/2}| ≤ K log₂ p` is safe from this construction unless `K ≳ 22` (§2.4). The tuple/crux "spike isolation" picture — one dilate at bias `≥ ½` and `O(log p)` others — is wrong at the level of constants: at `p = 10^6` the rectangle has `124` dilates at bias `≥ ¼` while the anchor sits at bias `0.51`.

### 2.2 The strongest reading also fails

One might object that `D_{1/2}` is relative (`|g(t)| ≥ ½|S|`, here bias `¼`) and that the anchor is only at bias `½`. Take instead bias `1` at the anchor (`B ⊆ N(1)`, so `S = n`, `z² = n`, and `n = ⌈8 ln p⌉` is the SNR threshold) and the absolute threshold `|g(t)| ≥ ½mn` (which coincides with `D_{1/2}` since `S = mn`). The same design (target set inside `N(1)`, no anchor weight needed) gives, after the targeted local search:

| `p` | `n = |B|` | `z²/ln p` | `|{x ≠ 0 : |T_B(x)| ≥ n/2}|` | `2 log₂ p` | ratio |
|---|---|---|---|---|---|
| 10 007 | 74 | 8.03 | 18 | 26.6 | 1.35 |
| 100 003 | 94 | 8.16 | **42** | 33.2 | **2.53** |
| 1 000 003 | 112 | 8.11 | **45** | 39.9 | **2.26** |
| 10 000 019 | 130 | 8.07 | **51** | 46.5 | **2.19** |

(Verifier `SI_star_strong_form`.) So even "a complete `1×n` rectangle at the SNR threshold has at most `2 log₂ p` half-biased dilates" is false from `p ≈ 10^5` on.

### 2.3 Families that do **not** break SI\* (the SI note's families, pushed further)

* **GP-type `(U_k, N(U_L))` at large `p` (task (a)).** For `p ∈ {10 007, 100 003, 1 000 003}`, 40/40/24 random ratios `r`, all `L ∈ [⌊log₂p⌋ − 9, ⌊log₂p⌋ + 3]` with `N(U_L) ≠ ∅` and `k ∈ {2, 3, 4, 6, 8, L/2, L}`: 1546 / 2283 / 1477 rectangles. Maximum of `|D_{1/2}|/log₂ p` among those with `z² ≥ 8 ln p`: `0.98` (`p = 10 007`, `r = 477`, `L = k = 10`, `|B| = 24`, `|D| = 13`), `0.90` (`p = 100 003`, `L = k = 14`, `|B| = 10`, `|D| = 15`), `1.05` (`p = 1 000 003`, `r = 48 998`, `L = k = 17`, `|B| = 11`, `z² = 8.04 ln p`, `|D| = 21`). This confirms the SI note's `≈ 1·log₂ p` ceiling for this family; the family is the wrong one — its `|D|` is the length of a complete bipartite design (`≤ log₂ p + O(1)`), not a half-biased one.
* **Subgroups (task (b)).** Over all 2258 primes `11 ≤ p ≤ 20 000` and all subgroups with `M_H ≥ |H|/2` (12 804 biased rectangles `(cH⁺, H)`): the largest `z²/ln p` is `2.84` (`p = 19 441`, `h = 30`, `M_H = 29`, `|D| = 270 = 19.0 log₂ p`); the largest `h/ln p` with `M_H ≥ h/2` is `7.36` (`p = 10 337`, `h = 68`, `M_H = 36`, `z² = 2.06 ln p`, `|D| = 340 = 25.5 log₂ p`); among the 264 rectangles with `h ≥ 4 ln p` the largest `z²/ln p` is `2.17`. So the SI\* hypothesis is never met by a biased subgroup rectangle (`z² = M_H²/h ≤ 2.84 ln p < 8 ln p`), exactly as the SI note's Theorem 2.2 predicts, and the hypothesis constant `8` could be lowered to `3` without a subgroup counterexample at `p ≤ 20 000` — but §2.1 shows lowering it does not matter.
* **Small rectangles by annealing (task (c)).** `p ∈ {1009, 2003, 4001, 10 007}`, `(m,n) ∈ {12×12, 8×12, 6×12, 12×8, 4×12, 3×12}`, three restarts of 6000–12 000 steps each, maximising `|D_{1/2}|` under the SI\* hypothesis (admissible rectangles need bias `≥ √(8 ln p/mn) ≈ 0.62–0.72`, so they are near-complete): best admissible `|D_{1/2}| = 5, 5, 6, 3` (`≤ 0.50 log₂ p`); `4×12` and `3×12` admit no rectangle at all. Small rectangles cannot refute SI\*: with `n` rows only, the per-dilate bias fluctuates by `1/√n ≈ 0.3` and the design capacity (§2.4) is `≈ 32 ln(p/n)/m ≤ 3 ln p` for `m = 12` before fluctuations.

### 2.4 What the constants have to be (HEURISTIC capacity bound, matching the data)

Choosing `B` as the top `n` of `p` rows of the `±1` matrix `(χ(x + b))_{b, x∈X}` by row sum, the mean row of `B` has per-coordinate bias `≈ √(2 ln(p/n)/|X|)` (the `(n/p)`-quantile of `N(0, |X|)` divided by `|X|`). A coordinate is in `D_{1/2}` iff its bias is `≥ ½·(anchor bias)`. Hence the number of designed half-biased dilates is

    |X| ≈ 2 ln(p/n) / β²,   β = threshold bias:   β = ¼ (SI\* as stated, anchor at ½) → |X| ≈ 32 ln(p/n);   β = ½ (strong form) → |X| ≈ 8 ln(p/n),

minus the loss to per-coordinate fluctuations (`±1/√n`), which the local search partly recovers. With `n ≍ ln p` this is `32 ln p (1 − o(1)) = 22.2 log₂ p` resp. `8 ln p = 5.55 log₂ p`; the observed `6.22 log₂ p` at `p = 10^6` (`32 ln(p/n)/log₂ p = 12.4` before losses) and `2.19–2.53 log₂ p` in the strong form (`8 ln(p/2n)/log₂ p = 3.1–3.6`) are `50–70 %` of the capacity. Conversely nothing in this construction gives more than `Θ(log p)` dilates, and a rectangle with `|D_{1/2}| ≫ log p` at high SNR would need `|X| ≫ ln(p/n)/β²`, i.e. row sums beyond the Gaussian tail of the `p` available rows — impossible for random-like rows and not achieved by any structured family tried. So the statement that *can* survive is `|D_{1/2}(A,B)| ≤ K·log₂ p` under the SNR hypothesis with `K ≥ 32 ln 2 ≈ 22.2` (relative threshold) or `K ≥ 8 ln 2 ≈ 5.5` (absolute threshold, bias-`1` spike); the hypothesis constant `8` is irrelevant to this (only `ln(p/n)` depends on it, logarithmically), and lowering it to `3` (the subgroup edge) changes nothing. Neither statement is tested beyond `p = 10^7`, and neither has any bearing on `t = 1` (SI note Prop. 2.3(4)).

---

## 3. Remaining obligations (OPEN)

1. **T′(k).** Prove or refute `|T_ns| ≤ p^{½+o(1)}(mn)^{2k}/min^{k−½}`. The subgroup case is the statement `|Σ_{s≠0}T_H(s)^{2k} − (p−1−h)N_k(h) − hN_k(h−1)| ≤ p^{½+o(1)}h^{k+½}`, i.e. that the `2k`-th moment of the shifted-subgroup sum over the `d` cosets fluctuates by no more than `√d·h^{k+1}·p^{o(1)}` — a "Gaussian moments across cosets" statement, strictly stronger than anything known about `M_H` at `|H| ≤ √p` (dual note §3.5). It is also exactly the input Theorem 1.2 needs; nothing here makes it more provable than T(k) was.
2. **Beyond subgroups.** Families with a large approximate symmetry group but `R_ν ≪ mn·min` (e.g. `A = H·A₀` with `A₀` random, `B` random) were not tested; the coherence heuristic predicts `|T_ns| ≈ √(p·|H|)R_χ^k`, which is *below* T′(k) when `R_ν ≈ mn`, but this should be checked.
3. **SI\*.** The constant in any log-type bound under the SNR hypothesis is at least `≈ 6.2` (observed) and heuristically `≈ 22`; whether `|D_{1/2}| = O(log p)` holds at all under the hypothesis is open, and the `1×n` construction should be pushed to `p ≈ 10^8` (cost `O(pn)` per evaluation) to see the ratio saturate. The SI note's Prop. 3.4 route (Gaussian `2k`-th dilation moments at `k ≈ ln p`) would give `|D_{1/2}| ≤ 4Γ`, so for the designed rectangles the dilation moments at `k ≈ ln p` exceed the Gaussian value by a factor `Γ ≳ 30`: a second, independent, exact witness that the "independent signs" model of the dilation profile fails at logarithmic order, complementing §1.
4. **The p = 13 183 subgroup.** `H` of order `169 = 13²` (`(p−1)/2 = 3·13³`) has `M_H = 65 = 5.0√h = 1.69·√(2h log d)`, far above the dual note's maximum `1.35`; whether prime-power-order subgroups with `13 | d` systematically have large `M_H` was not investigated.
