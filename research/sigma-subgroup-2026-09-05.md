# Sigma pass, direction `subgroup`: audit of the dyadic Gauss-period bound (2026-09-05)

**Status:** PROVED — the classical moment identities (Lemma 1.1–1.2) and the Karatsuba–Konyagin inequality
`M ≤ (pE_kE_l)^{1/(2kl)} n^{1−1/k−1/l}` for all `k,l ≥ 1` (Theorem 2.1, self-contained proof; CITED as Shkredov arXiv:1311.5726
Lemma 7 (15), attributed to Konyagin 2002); consequently pass 7's bound `M ≤ (B+2)^{1/9}n^{8/9}` under `E_3 ≤ Bn^3`, `p ≤ n^4`
is a *weaker* form of the 2002 inequality at `k=l=3` (`B^{1/9}n^{8/9}`), and every exponent in passes 5–7 (23/24, 17/18, 8/9,
conditional 7/8) is implied by Theorem 2.1 with the same energy inputs. AUDIT — the only class-dependent input of pass 7 is the
single-prime hypothesis `E_3(H_N) ≤ (15+log N)N^3`; no lemma fails outside the class, the hypothesis is simply unavailable there.
REFUTED (witnesses) — the class restriction is not vacuous: 14 of the 24,379 quartic primes at `N = 64` (e.g. `p = 7204033`,
`E_3 = 20.93·64^3 > 19.16·64^3`) lie outside it; the recurrence/level-descent route (i) cannot give any unconditional `c > 0`
(exact ledger); Konyagin's inequality with every cited Stepanov-type energy bound is trivial at `p = n^4` (exact ledger, `E_3 ≪
n^4 log n` sits exactly on the boundary). CITED — the best unconditional exponent for all primes at `p ≈ n^4` is Di Benedetto et
al. (2020) Theorem 3.1: `M ≤ c(ε)N^{2849/2880+ε}`, saving `31/2880 ≈ 0.0108`; no later improvement at `|H| ≈ p^{1/4}` was found.
OPEN — whether `E_3(H_N) = O(N^3 log N)` holds at *every* quartic prime (which would make 8/9 unconditional), and the target
`M ≤ C√(n log(p/n))` itself. DATA — 307,645 pointwise checks pass at 12,126 primes with exact periods; over all 51,780 scanned
primes `E_3/N^3 ≤ 20.93`; the ratio `M/√(n log(p/n))` has mean 1.29 and maximum 1.70 at `N = 64`, both below an iid-Gaussian
extreme-value null model (1.34 and 1.74), i.e. no sign of a growing constant through `N = 64`.
Verifier: `experiments/sigma_subgroup_2026_09_05.py` → `results/sigma_subgroup_2026_09_05.json` (≈175 s).

## 0. Plan of the note

1. What pass 7 claims, the prime class, and where the argument fails outside it (audit, Task 1).
2. Classical moment baseline written out exactly (Task 2).
3. Unconditional bounds for all primes at p ~ n^4 (Task 3).
4. Exact M(p,n) tables, ratio M/sqrt(n log(p/n)) (Task 4).
5. Verifier and results (Task 5).
6. Remaining obligations.

(Sections filled in below as work completes.)

## 1. Notation and the exact classical baseline (PROVED)

Throughout, `p` is prime, `H = μ_n ⊂ F_p^*` the subgroup of order `n | p−1`, `m = (p−1)/n`,
`e(t) = exp(2πit/p)`, `η(b) = Σ_{h∈H} e(bh)`, `M = M(p,n) = max_{b≠0}|η(b)|`, `Δ = M/n`.
For `k ≥ 1`, `r_k(x) = #{(h_1,…,h_k) ∈ H^k : h_1+⋯+h_k = x}` and
`E_k(H) = Σ_x r_k(x)^2 = #{h_1+⋯+h_k = h_{k+1}+⋯+h_{2k}}`. Since `n` is even, `−1 ∈ H`.

**Lemma 1.1 (moment identity, PROVED).** `Σ_{b∈F_p} |η(b)|^{2k} = p·E_k(H)`.
*Proof.* `η(b)^k = Σ_x r_k(x)e(bx)` and `|η(b)|^{2k} = Σ_{x,y} r_k(x)r_k(y)e(b(x−y))`; summing over `b`
and using `Σ_b e(bt) = p·1_{t=0}` gives `p Σ_x r_k(x)^2`. ∎

**Lemma 1.2 (coset constancy, PROVED).** `η(bh) = η(b)` for `h ∈ H`, so `b ↦ |η(b)|` is constant on
the `m` cosets of `H` in `F_p^*`, and `η(0) = n`. Hence, with `Q_k := (pE_k − n^{2k})/n = Σ_{cosets}|η_c|^{2k}`,

    M^{2k} ≤ Q_k = (p·E_k(H) − n^{2k})/n ≤ m·M^{2k}.                                   (1.1)

*Proof.* `η(bh) = Σ_{h'} e(bhh') = η(b)` because `h' ↦ hh'` permutes `H`. The `p−1 = nm` nonzero `b`
split into `m` cosets, each contributing `n|η_c|^{2k}` to Lemma 1.1; subtract the `b=0` term `n^{2k}`. ∎

**Corollary 1.3 (depth-`k` moment bound, PROVED).** If `p ≤ n^4` and `E_k(H) ≤ n^{2k−3−2kc}` then
`M ≤ n^{1−c}`. *Proof.* `M^{2k} ≤ pE_k/n ≤ n^{3}·n^{2k−3−2kc} = n^{2k(1−c)}`. ∎

What this asks of `E_k` at `p ≈ n^4`: the diagonal solutions alone give `E_k ≥ k!·n^k − O(n^{k−1})`
(all `k!` orderings of a `k`-tuple, distinct entries), and Cauchy–Schwarz gives `E_k ≥ n^{2k}/p ≥ n^{2k−4}`.
So `E_k ≤ n^{2k−3−2kc}` is compatible with the diagonal only if `c ≤ (k−3)/(2k)` and with the principal
term only if `c ≤ 1/(2k)`. Hence **no fixed depth `k ≤ 3` can give any saving by (1.1) at `p ≈ n^4`, and
the best a fixed depth can give is `c = 1/8` at `k = 4`, requiring `E_4(H) ≤ n^4` (i.e. `Ẽ_4 = E_4 − n^8/p = O(n^4)`)**.
Depth `k=2`: the Heath-Brown–Konyagin bound `E_2 ≪ n^{5/2}` gives `M^4 ≤ p·n^{5/2}/n`, i.e. `M ≪ p^{1/4}n^{3/8}`,
which is `≤ n` only when `p ≤ n^{5/2}`; at `p ≈ n^4` it reads `M ≪ n^{11/8}`, trivial. (Exact-rational
checks of every exponent in this paragraph: `ledger()` in the verifier, items (a),(b).)

## 2. The Karatsuba–Konyagin inequality and what pass 7 actually proved (PROVED / audit)

**Theorem 2.1 (PROVED here; classically Konyagin 2002, see §2.3 for the citation status).**
For every subgroup `H ≤ F_p^*` of order `n`, all integers `k, l ≥ 1`, and every `a ≠ 0`,

    |η(a)|^{kl} ≤ n^{kl−k−l} · sqrt( p · E_k(H) · E_l(H) ),   i.e.   M ≤ (p E_k E_l)^{1/(2kl)} n^{1−1/k−1/l}.   (2.1)

*Proof.* (i) `Σ_y r_l(y) η(ay) = Σ_{h_1..h_l} Σ_h e(ah(h_1+⋯+h_l)) = Σ_{h∈H} η(ah)^l = n·η(a)^l` by Lemma 1.2.
(ii) Hölder with total mass `Σ_y r_l(y) = n^l`: `n^k|η(a)|^{kl} = |Σ_y r_l(y)η(ay)|^k ≤ n^{l(k−1)} Σ_y r_l(y)|η(ay)|^k`.
(iii) Let `A(x) = p^{−1}Σ_t |η(t)|^k e(−tx)` (a real, possibly signed function). Fourier inversion gives
`|η(t)|^k = Σ_x A(x)e(tx)` and Parseval gives `Σ_x A(x)^2 = p^{−1}Σ_t|η(t)|^{2k} = E_k` (Lemma 1.1).
(iv) Bilinear orthogonality: for any `U, V : F_p → C` and `a ≠ 0`,
`|Σ_{x,y}U(x)V(y)e(axy)|^2 ≤ ‖U‖_2^2 Σ_x |Σ_y V(y)e(axy)|^2 = p‖U‖_2^2‖V‖_2^2` (Cauchy–Schwarz, then
Parseval in `x` since `x ↦ ax` is a bijection). With `U = A`, `V = r_l`:
`Σ_y r_l(y)|η(ay)|^k = Σ_{x,y}A(x)r_l(y)e(axy) ≤ sqrt(p E_k E_l)`. Combine with (ii). ∎

No positivity of `A` is used, so (2.1) holds for odd `k` as well; the "signed Fourier multiplier" of
pass 7 (`α_r = A/n^r`) is exactly step (iii). Pass 7's gate (7) is the same argument with the origin terms
of `A` and `r_l` split off and the remaining parts centred; it reads `Δ^{rs} ≤ 2/n + sqrt(p Ẽ_r Ẽ_s)/n^{r+s}`.
Since `Ẽ_r = E_r − n^{2r}/p ≤ E_r`, and at `p ≈ n^4`, `r = s = 3`, the subtraction `n^6/p ≤ 4n^2` is negligible
against `E_3 ≥ T_6(n) ≈ 15n^3`, the two bounds are numerically almost identical, and the classical one is
*better* by the additive term:

    pass 7 (13):  M ≤ (B+2)^{1/9} n^{8/9}          Theorem 2.1 (k=l=3):  M ≤ B^{1/9} n^{8/9}      (E_3 ≤ Bn^3, p ≤ n^4).

The verifier records the margin of (2.1) at every scanned prime for all `(k,l)` with `k,l ≤ 4` (`konyagin_margin`)
and the pass-7 gate margin (`gate_margin`); the pass-7 class bound `(17+log N)^{1/9}N^{8/9}` is checked
against the true `M` at every scanned prime, class or not (`pass7_class_bound_ok`).

**Audit conclusion for Task 1 (the class and the failing lemma).** Pass 7 claims: for dyadic `N`, on the class
`K_N = {p ∈ P_N : Σ_{2≤t|N, t dyadic} (E_3(H_t) − T_6(t))/t^3 ≤ log N}`, `P_N = {p ≡ 1 (N), N^4/4 ≤ p ≤ N^4}`,
`T_6(t) = 15t^3 − 45t^2 + 40t`, one has `M(p,N) ≤ (17+log N)^{1/9}N^{8/9}`, with `|K_N|/|P_N| = 1 − O(1/log N)`.
The only class-dependent input is the single-prime energy bound `E_3(H_N) ≤ (15 + log N)N^3` (pass 5, Theorem A,
via Markov applied to the prime-averaged norm estimate (N')). Everything else — Lemma 1.1, 1.2, Theorem 2.1 —
is unconditional for every subgroup of every prime. So the "lemma that fails outside the class" is not a lemma
at all: it is the hypothesis `E_3(H_N) = O(N^3 log N)`, which the norm estimate only controls *on average over p*.
Outside the class the unconditional input is Lemma 4.3 of Di Benedetto et al. (CITED below): `E_3 ≪ n^4 log n`
for `n < p^{1/2}`, and then (2.1) with `k=l=3`, `p ≤ n^4` gives `M ≪ (n^4·n^8 log^2 n)^{1/18} n^{1/3} = n (log n)^{1/9}`,
i.e. *nothing* — this is exactly Konyagin's 2002 threshold `n > p^{1/4+ε}` seen from the wrong side.

**Density-one is also cheaper than pass 5 states.** For the endpoint bound only the level `t = N` matters:
(N') at `k=6` gives `Σ_{p≡1(N)} (E_3(H_N) − T_6(N)) log p ≤ ½N^6 log 6`, so at most `N^3 log 6/(2w log(N^4/4))`
primes of `P_N` have `E_3(H_N) > T_6(N) + wN^3`; the all-levels sum in `K_N` is not needed for (2.1).

## 3. Exact scan results (final run, 175 s, 12 workers; all numbers from `results/sigma_subgroup_2026_09_05.json`)

Scan design (all exact integer energies; periods in double precision with recorded error ≤ 2·10⁻¹⁴):

| N | energy window (all p ≡ 1 mod N) | primes | quartic primes | periods | E₄ |
|---|---|---|---|---|---|
| 2,4,8,16 | [N^3.5, N^4.5] | 3, 31, 276, 2617 | 1, 10, 80, 579 | all | all |
| 32 | [185364, 5931641] | 24474 | 3705 | 7166 (all quartic + every 6th outside) | 2657 |
| 64 | [N⁴/4, N⁴] = [4194304, 16777216] | 24379 | 24379 | 2033 (every 12th, plus 6700417) | 255 |

**Witnesses: primes outside the pass-5/7 class in the quartic window.** For `N = 64` exactly 14 of the 24379
quartic primes violate `E_3(H_64) ≤ (15 + log 64)·64^3 = 5022494` (and hence the class condition):
6878593, 7041409, 7177601, 7204033, 7884353, 7987009, 8019073, 9190913, 9877633, 10219457, 11127041, 12942337,
13640513, 14721281. The largest is `p = 7204033` with `E_3 = 5486080 = 20.93·64^3` and `E_2 = 13632 = 3.33·64^2`.
Pass 5's Markov bound allows up to 4232 such primes; 14 occur. For `N ≤ 32` no quartic prime leaves the class
(for `N = 32`, `p = 194977 < 32^4/4` has `E_3 = 19.43·32^3` and class sum 5.80 > log 32, but it is outside the window).
So the class restriction in pass 7's theorem is **not vacuous** at `N = 64`; at the two exceptional primes where
periods were computed (11127041, 13640513) the true `M` (41.28, 47.55) still lies below the class value
`(17+log 64)^{1/9}·64^{8/9} = 56.59` — that is data, not a theorem.

**No used inequality fails.** 307,645 pointwise checks (moment identity to relative error ≤ 9·10⁻¹⁶; coset
constancy; `M^{2k} ≤ Q_k` for `k ≤ 4`; Theorem 2.1 for all ten pairs `k,l ≤ 4`; pass-7 gate (7) for eight `(r,s)`;
Jensen step; `a_3, u_3 ≤ 1/n`; pass-7 (13); level descent) pass at every one of the 12,126 primes with periods.
Minimum relative margins at `N=64`: Theorem 2.1 (3,3): 0.148; pass-7 gate (3,3): 0.786; Jensen: 0.226.

**Energies over ALL quartic primes (not only class primes).** Maxima of `E_k/N^k`:

| N | max E₂/N² | max E₃/N³ | max E₄/N⁴ (sampled) | k!: 2, 6, 24 | n^{2k}/p at p=N⁴: 1, N², N⁴ |
|---|---|---|---|---|---|
| 16 | 2.81 | 12.46 | 73.5 | | 1, 256, 65536 |
| 32 | 2.91 | 14.34 | 104.3 | | 1, 1024, 1048576 |
| 64 | 3.33 | 20.93 | 171.0 | | 1, 4096, 16777216 |

`E_2` is exactly intrinsic (`3N²−3N`) at 24359/24379 quartic primes for `N=64`; `E_3` is exactly `T_6(N)` at 21124.
`E_4/N^4` is dominated by the principal term `N^8/p ∈ [N^4, 4N^4]` plus the intrinsic count `T_8(N)/N^4 ≈ 105`:
the centred `Ẽ_4 = E_4 − N^8/p` stays `O(N^4)` in the sample (needed for exponent 7/8, §2).

**Konyagin exponents achieved by the true energies** (`log_N` of the bound (2.1), max over sampled primes, `N=64`):
(3,3): 0.967, (2,4): 0.962, (3,4): 0.949, (4,4): 0.945; true `max log_N M = 0.929`. Every pair with `min(k,l) ≤ 1`
or `(2,2),(2,3)` is trivial (`≥ 1`). So with the *actual* energies the best fixed-order Konyagin bound at `N=64`
is `M ≤ N^{0.945}`, against the truth `N^{0.929}` and the target `≈ N^{0.5+o(1)}`.

**Ratio `M/√(n log(p/n))`** (all sampled primes; quartic-window max in parentheses):

| N | count | max | mean | median | argmax p (class?) |
|---|---|---|---|---|---|
| 2 | 3 | 1.004 | 0.961 | 0.950 | 13 |
| 4 | 31 | 1.026 | 0.951 | 0.943 | 137 |
| 8 | 276 | 1.122 | 1.049 | 1.043 | 1601 |
| 16 | 2617 | 1.260 | 1.155 | 1.152 | 17377 |
| 32 | 7166 | 1.446 (1.446) | 1.251 | 1.246 | 376577 (in class) |
| 64 | 2033 | 1.697 (1.697) | 1.285 | 1.275 | 13640513 (**not** in class, E₃ = 20.6 N³) |

The maximum grows with `N`; the mean and median grow slowly (0.96 → 1.29). Interpretation in §4.

## 4. Unconditional bounds for ALL primes at p ≈ n⁴ (Task 3)

### 4.1 Citations (fetched sources archived under `sources/sigma-subgroup-2026-09-05/` with sha256 manifest)

- **[Sh14] CITED.** I. D. Shkredov, *On exponential sums over multiplicative subgroups of medium size*, arXiv:1311.5726v1,
  Lemma 7, formula (15): "for any positive integers l and m one has
  `M(Γ) ≤ p^{1/2lm} T_l^{1/2lm}(Γ) T_m^{1/2lm}(Γ)·|Γ|^{1−1/l−1/m}`", attributed there to [13] = S. V. Konyagin,
  *Estimates for trigonometric sums and for Gaussian sums*, IV Intern. Conf. Modern Problems of Number Theory (Tula 2001),
  Moscow 2002, and [18]. `T_k(Γ)` is the `2k`-variable count `#{a_1+⋯+a_k = a'_1+⋯+a'_k}` (his p. 3), i.e. our `E_k`.
  This is Theorem 2.1 verbatim. Same source, Corollary 5 (11): for `|Γ| ≤ p^{2/3}`, `E(Γ) ≪ |Γ|^{5/2}` ("It was proved in
  [11] (see also [13]) that E(Γ) = O(|Γ|^{5/2}), provided by |Γ| ≤ p^{2/3}", [11] = Heath-Brown–Konyagin, Quart. J. Math.
  51 (2000)); Theorem 6 (13), attributed to [13]: for `|Γ| < √p` and `d ≥ 2`, `T_d(Γ) ≪_d |Γ|^{2d−2+2^{1−d}}`.
  (His `E_3(Γ) ≪ |Γ|^3 log|Γ|` in (11) is the *cubic* energy `Σ_x r_{Γ−Γ}(x)^3`, not our `E_3`; not used.)
- **[MRSS] CITED.** Murphy–Rudnev–Shkredov–Shteinikov, *On the few products, many sums problem*, arXiv:1712.00410,
  Theorem 3: "`E(A) ≪ M^{8/5}|A|^{49/20} log^{1/5}|A|` … The same estimate, with M = 1 holds for a multiplicative
  subgroup `Γ ⊂ F_p^×`, with `|Γ| ≤ √p`"; Corollary 7: "`T_3(A) ≪ M^{12}|A|^4 log|A|`, the same holds with M = 1 if A is
  replaced by Γ, a multiplicative subgroup in `F_p^×` with `O(√p)` elements." (Their `T_k` = our `E_k`, their (20).)
- **[DGGGST] CITED.** Di Benedetto, Garaev, Garcia, Gonzalez-Sanchez, Shparlinski, Trujillo, arXiv:2003.06165v1
  (J. Number Theory 2020), Theorem 3.1: "Let `H` be a multiplicative subgroup of `F_p^*` of order `H` with `p^{1/2} > H > p^{1/4}`.
  Then `max_{(a,p)=1}|S_a(H)| ≲ H^{2689/2880} p^{1/72}`. In particular, when `H > p^{1/4}`, … `≲ H^{1−31/2880}`", where
  (their §2) "`A ≲ B` means `|A| < B p^{o(1)}`, or equivalently, for any `ε>0` there is a constant `c(ε)` … such that
  `|A| ≤ c(ε) B p^ε`". Its inputs are Lemma 4.2 (= [MRSS] Thm 3), Lemma 4.3 (= [MRSS] Cor 7) and the Petridis–Shparlinski
  trilinear bound (their Lemma 4.1). Their (1.2) records Bourgain–Garaev (2009): `≤ H^{1−175/9437184+o(1)}` at `H ∼ p^{1/4}`.
- **[BGK] CITED via Kowalski.** E. Kowalski, *Exponential sums over small subgroups, revisited*, arXiv:2401.04756,
  Theorem 1.1 (Bourgain–Glibichuk–Konyagin, J. London Math. Soc. 73 (2006)): "Let γ > 0 … There exists ν > 0, depending
  only on γ, such that for any prime p and any subgroup `H ⊂ F_p^×` with `|H| ≥ p^γ`, `Σ_{x∈H} e(ax/p) ≪ |H|p^{−ν}` for any
  `a ∈ F_p^×`". No explicit ν. Kowalski, Remark 1.2(3): "The dependency of the exponent ν on γ can be made explicit …
  currently the sharpest result … is due to Shkredov [11, Cor. 16]", [11] = Shkredov, Moscow J. Comb. Number Th. 8 (2019).
- **[Sh19] CITED, numerical reading flagged.** Shkredov, *Some remarks on the asymmetric sum–product phenomenon*,
  arXiv:1705.09703, Corollary 16: for `|Γ| ≥ p^δ`, `max_{ξ≠0}|Γ̂(ξ)| ≪ |Γ|·p^{−δ/2^{7+2δ^{−1}}}` (text extraction renders the
  exponent as "27+2δ−1"; the proof's last line `ρ ≪ |Γ|^{1−1/2^{k+2}}` with `k = ⌈2 log p/ log|Γ| + 4⌉ ≤ 2/δ + 5` fixes the
  reading `2^{7+2/δ}`). Its input is his Theorem 12, a bound for `T_{2^k}(Γ)` of the same shape as the centred formula (C)
  already discussed in `research/analytic-bounds-and-amplification.md`. At `|Γ| = p^{1/4}`: `k = 12`, saving `1/2^{14}`
  (proof line) or `1/2^{15}` (stated form) — explicit, but ~200× smaller than `31/2880`.

### 4.2 What these give at p ≈ n⁴ (PROVED consequences of the citations; exact-rational checks in `ledger()`)

Write `E_k ≤ n^{e_k}` (log factors dropped) and `p ≤ n^4`. Theorem 2.1 gives `M ≤ n^{κ(k,l)}`,
`κ(k,l) = (4 + e_k + e_l)/(2kl) + 1 − 1/k − 1/l`. With every unconditional input available for `n ≈ p^{1/4}`:

| inputs | pair (k,l) | κ | nontrivial? |
|---|---|---|---|
| `e_2 = 5/2` (HBK) | (2,1) | 11/8 | no |
| `e_2 = 5/2` | (2,2) | 9/8 | no |
| `e_2 = 49/20` (MRSS) | (2,2) | 89/80 | no |
| `e_2 = 49/20, e_3 = 4` (MRSS) | (2,3) | 83/80 | no |
| `e_3 = 4` (MRSS Cor 7) | (3,3) | 1 exactly (`×(log n)^{1/9}`) | no — this is Konyagin's 2002 threshold `n > p^{1/4+ε}` seen from below |
| `e_3 = 4, e_4 = 49/8` (Sh14 Thm 6) | (3,4) | 193/192 | no |
| `e_4 = 49/8` | (4,4) | 129/128 | no |

So **no Stepanov-type energy bound in the literature makes Theorem 2.1 nontrivial at `p = n^4`** — the known `E_3 ≪ n^4 log n`
sits exactly on the boundary. The tower route (i) is worse: the exact recurrence `E_2(H_{2k}) = 2E_2(K) + 6B + 8T` with the
only unconditional inputs `B ≤ E_2(K)`, `T ≤ E_2(K)` gives `E_2(H_{2k}) ≤ 16E_2(H_k)`, hence `E_2(H_N) ≤ 6·16^{log_2 N − 1} = (3/8)N^4`,
weaker than the trivial `E_2 ≤ n^3`; and the level-descent inequality `M_N ≤ 2^{j}M_{N/2^j}` combined with (2.1) at the lower
level gives `M_N ≤ min_j 2^{j/3}(E_3(H_{N/2^j})/(N/2^j)^3)^{1/9} N^{8/9}` (PROVED; the `2^{j/3}` bookkeeping is item (g) of
the ledger), which is never below the `j = 0` term unless a lower level has energy ratio smaller by `2^{3j}`; in the data the
minimum over `j` is attained at `j = 0` at every prime with `N ≥ 16`. **Route (i): REFUTED as a route to any `c > 0`.**

**Best unconditional exponent at `p ≈ n^4`, all primes (CITED, [DGGGST] Thm 3.1):** in the window `N^4/4 ≤ p ≤ N^4` one has
`N > p^{1/4}` (a prime is not a fourth power) and `N < p^{1/2}` for `N ≥ 4`, so for every `ε > 0`

    M(p, N) ≤ c(ε) · N^{2689/2880} · (N^4)^{1/72} · N^{ε} = c(ε) · N^{2849/2880 + ε},      2849/2880 = 1 − 31/2880.

`c = 31/2880 − ε ≈ 0.0108` is the largest explicit unconditional saving I can source (search for post-2020 improvements at
`|H| ≈ p^{1/4}` found none; [Sh19] Cor 16 is explicit but far smaller there). **The class restriction of pass 7 cannot be removed
by any cited input**: `M ≤ B^{1/9}n^{8/9}` needs `E_3(H_N) ≤ B N^3`, and the gap between the unconditional `N^4 log N` and the
needed `N^{3+o(1)}` is the whole content. Whether `E_3(H_N) = O(N^3 log N)` holds for *every* quartic prime is **OPEN**; the
prime-averaged norm estimate (N') bounds only the *sum* over primes and tolerates single primes with `E_3` as large as
`≈ N^6/(8 log N)`. Empirically (§3) `max E_3/N^3 = 20.93` over all 24,379 quartic primes at `N = 64`, i.e. `E_3 = O(N^3)`
with a small constant even at the 14 non-class primes — data, not a theorem.

## 5. The ratio M/√(n log(p/n)) and what the data say (Task 4; HEURISTIC interpretation, exact data)

Exact `M(p,n)` was computed at 12,126 primes (all `p ≡ 1 (N)` in `[N^{3.5}, N^{4.5}]` for `N ≤ 16`; all 3,705 quartic primes
and every 6th prime outside the window for `N = 32`; every 12th quartic prime, plus `p = 6700417`, for `N = 64`). The maximum of
the ratio over the sample grows with `N` (1.00, 1.03, 1.12, 1.26, 1.45, 1.70 for `N = 2,…,64`). Whether that growth means a
growing constant is exactly the question, so compare with the null model in which the `m = (p−1)/n` nonprincipal periods are
iid real `N(0, n)` (their exact mean square is `n(p−n)/(p−1)`): then `M/√(n log(p/n)) ≈ E[max_m |Z|]/√(log(p/n))`, and the
maximum over `S` sampled primes is the maximum of `mS` iid `|Z|`. Evaluated at the median-`p` prime (exact integral
`∫_0^∞ (1 − erf(x/√2)^m) dx`, verifier `gaussian_expected_max_abs`):

| N | m (median p) | primes S | model mean ratio | observed mean | model max over S | observed max |
|---|---|---|---|---|---|---|
| 8 | 782 | 276 | 1.305 | 1.049 | 1.818 | 1.122 |
| 16 | 8401 | 2617 | 1.323 | 1.155 | 1.851 | 1.260 |
| 32 | 31496 | 7166 | 1.331 | 1.251 | 1.852 | 1.446 |
| 64 | 162283 | 2033 | 1.338 | 1.285 | 1.738 | 1.697 |

Reading: (i) the observed mean sits *below* the iid-Gaussian value and rises toward it (1.05 → 1.29 against 1.30 → 1.34); the
limit of the model mean is `√2`; (ii) the observed *maximum* over the sample is below the null prediction at every `N`
(1.12, 1.26, 1.45, 1.70 observed vs 1.82, 1.85, 1.85, 1.74 predicted), and the gap closes as `N` grows because the trivial cap
`M ≤ n` bounds the ratio by `√(n/log(p/n))` (≈ 1.10 at `N = 8`, ≈ 2.3 at `N = 64`), which the Gaussian model ignores. So the
growth of the observed maximum with `N` is what the null model predicts for a growing sample whose cap is loosening, not
evidence of a growing constant: through `N = 64` the data are consistent with `M ≤ C√(n log(p/n))` with `C` of order `√2`.
The one prime that reaches 1.70 (`p = 13640513`) is a non-class prime with one extra quadruple orbit, i.e. exactly the
arithmetic resonance the target must tolerate. None of this is a proof: the sample
at `N = 64` covers 1/12 of the window, the null model is a heuristic, and `N = 64` is far from the sponsor's `n ≈ 2^{30}`.

Two further exact statistics recorded per prime: the smallest depth `k` at which the *centred* moment bound
`M ≤ Q_k^{1/(2k)}` beats the trivial `n` (6–9 at `N = 64`, 6–12 at `N = 32`, 8–17 at `N = 16` — always well above the fixed
depths 2, 3 that the literature's energy bounds address), and the (SG) constant `K_min = (Q_r/m)^{1/r}/(rn)` at depth
`r = ⌈log m⌉` (max over the sample 0.80, 0.73, 0.78, 0.94, 1.11, 1.50 for `N = 2,…,64`): (SG) with `K = 2` holds at every
sampled prime, but `K_min` is not monotone-bounded in this range, so the finite data neither support nor refute a uniform `K`.

## 6. Verifier

`experiments/sigma_subgroup_2026_09_05.py` (standard library + numpy; 12 worker processes; ≈175 s on this machine;
`SIGMA_QUICK=1` for a 2 s smoke run) writes `results/sigma_subgroup_2026_09_05.json` containing: the per-order summaries of
§3/§5, the full list of exceptional primes with their level-by-level `E_3/t^3`, all ten Konyagin pairs' margins, the pass-7 gate
margins, the Jensen/`a_3`/`u_3` checks, the moment identities, the coset-constancy check, the five exact norm-estimate (N')
certificates at `(n,k) ∈ {(4,4),(8,4),(16,4),(4,6),(8,6)}` (complete lists of primes with extra 4- and 6-term relations, by
Bareiss determinants over `Z[X]/(X^{n/2}+1)`, cross-checked by a direct scan), the exact-rational exponent ledger of §1, §2, §4
(items (a)–(l)), a compact table `[p, E_2, E_3, E_4, M, ratio, in_class, class_sum]` for every scanned prime, and sha256 hashes
of every input note. Every floating-point inequality is tested with relative tolerance `10^{−9}` and its margin recorded; energies
are exact integers (`≤ 2^{42}`).

## 7. Remaining obligations (OPEN)

1. **Single-prime sixth energy.** Prove `E_3(H_N) ≤ N^{3+o(1)}` for every prime `p ≡ 1 (N)` in `[N^4/4, N^4]`, or exhibit a
   prime with `E_3(H_N) ≥ N^{3+c}`. Either resolves whether the 8/9 exponent (Theorem 2.1 with `k=l=3`) is unconditional.
   The prime-averaged norm budget cannot decide this. Numerically `E_3/N^3 ≤ 20.93` through `N = 64`.
2. **Beyond 7/8.** Even `Ẽ_4 = O(N^4)` (OPEN) only gives `M ≤ N^{7/8}` through (1.1) or (2.1); the target `√(n log(p/n))`
   needs the centred moment bound (SG) at depth `≈ log m`, i.e. control of `E_k − n^{2k}/p` up to `(Kkn)^k` uniformly in `k`,
   which no fixed-order energy input can supply (§1, ledger items (a),(e)).
3. **Exceptional-prime structure.** All 14 exceptional primes at `N = 64` have intrinsic energies at every level `t ≤ 32` and a
   single extra four-distinct quadruple orbit at level 64 (`E_2 = 3N^2 − 3N + 24N`), i.e. they are exactly the window primes
   dividing the kernel discriminant `A_64` of `research/kernel-discriminant.md`; their excess `E_3 − T_6(N) ≈ 6.3N^3` is the
   combinatorial shadow of that one orbit. A bound `v_p(A_N) = O(N log N)`-type statement would control such primes' `E_2`
   but not directly `E_3`; no such bound is proved.
4. **Reduction to the prize** and the two-set Paley conjecture are untouched by this note.
