# Sigma pass, direction `extractor`: the two-set Paley conjecture as a two-source-extractor statement (2026-09-05)

**Status:** PROVED (complete elementary arguments, every step checked exactly): Theorem 1, the exact equivalence, for every rate `ε`, between the two-set conjecture (Satake's Conjecture 7 / property `𝒫(α,β)`) and each of (i) "`Ext_p(x,y) = ½(1+χ(x+y))` is an `(ε log₂p, p^{−δ})` two-source extractor for flat sources", (ii) the same for arbitrary sources of min-entropy `≥ ε log₂p` (convexity, with a constructive flat decomposition), (iii) the strong-extractor form; **this is a reformulation and contains no progress on the conjecture** — the direction (a)⇒(ii) is Chor–Goldreich 1988 as restated in Satake 2024, Proposition 15. PROVED: Theorem 3 (weighted Chung: `|E χ(X+Y)| ≤ √p·2^{−(H₂(X)+H₂(Y))/2}`, the Chor–Goldreich regime); Theorem 4 (`S(A,B)² ≤ p·E⁺(A,B)`, never below Chung, so additive-energy inputs cannot move Paley below rate 1/2 through the Fourier side); Theorem 5 (the Cauchy–Schwarz squaring trick that drives Bourgain's extractor is entropy-neutral for the Paley kernel: exact identity); Proposition 6 (no single-source condenser); Proposition 7 (the BIW product condenser reaches Paley only through the dilation family). CITED with fetched sources: Chor–Goldreich regime, Karatsuba regime, Barak–Impagliazzo–Wigderson (Theorem 1.1, Lemmas 3.1–3.2, §3.3), Bourgain 2005 via Rao's exposition (Lemmas 3.3, 3.4, 3.7, Corollary 3.5, Theorem 3.8) and Lewko (Theorems 2–3, Lemma 6, Theorem 9), Raz 2005 (via arXiv:2506.15547), Chattopadhyay–Zuckerman, Cohen, Ben-Aroya–Doron–Ta-Shma, Li 2023, Satake 2024 (Definitions 13–14, Proposition 15, Conjecture 9, Theorem 18, Corollaries 21 and 31, Lemma 28). REFUTED with witnesses: `H_avg` ("`|S(A,B)| ≤ 2·mean_{c∈F_p^*}|S(A,cB)|`"), the exact form in which a product condenser would have to control the two-set sum — at `p = 8009`, `|A| = |B| = 57`, the two-set bias is `0.376` while the mean over dilates is `0.0139`. OPEN: the conjecture; every published route below rate 1/2 changes the function (Bourgain/Lewko: Hadamard of a paraboloid encoding; BIW: iterated `ab+c` on ≥ 3 sources; the CZ line: non-algebraic), and in every one of them the `F_{p²}` test locates the prime-field input in an incidence or sum-product theorem, whereas Chung and Theorems 4–5 hold verbatim over `F_{p²}`.

Verifier: `experiments/sigma_extractor_2026_09_05.py` → `results/sigma_extractor_2026_09_05.json` (168,410 exact checks, 0 failures, 5.3 s). Sources fetched this pass are recorded in §2 with URLs; two PDFs (Rao's exposition, BIW) and one survey were text-extracted into the scratchpad. No file of any other pass was edited.

## 0. Conventions

`p` an odd prime, `χ` the Legendre symbol with `χ(0) = 0`, `S(A,B) = Σ_{a∈A,b∈B} χ(a+b)`, `f_A(b) = Σ_{a∈A} χ(a+b)`. A *source* is a probability distribution `P` on `F_p`; `H_∞(P) = −log₂ max_x P(x)`; `cp(P) = Σ_x P(x)²` and `H₂(P) = −log₂ cp(P)`; `U_A` is the flat source on `A ⊆ F_p`. "Rate `ε`" means `H_∞ ≥ ε log₂ p`, i.e. `max P ≤ p^{−ε}`; for `U_A` this is `|A| ≥ p^ε`. Statistical distance `SD(P,Q) = ½Σ|P−Q|`; for a bit `b`, `SD(b, U₁) = ½|E(−1)^b|`. `Ext_p(x,y) := 1` if `x + y = 0`, else `½(1 + χ(x+y))`; Satake's Definition 13 uses `χ(x−y)` with value `1` on the diagonal, which is the same object after `B ↦ −B` (checked exactly, §4/S2). Two-source `(k,η)`-extractor (Satake, Definition 7; Chattopadhyay, Definition 1.4): `SD(Ext(X,Y), U₁) ≤ η` for all independent `X, Y` with `H_∞ ≥ k`. *Strong* (Rao 2007, §1): `Pr_{y∼Y}[SD(Ext(X,y),U₁) > η] ≤ η`, and symmetrically in `x`. Satake's property `𝒫(α,β)` (his Definition 14): for every `S,T ⊆ F_p` with `|S|,|T| > p^α`, `|Σ_{s∈S,t∈T}χ(s−t)| ≤ C p^{−β}|S||T|` for a constant `C`; the workspace's `C(ε,δ)` is `𝒫(ε,δ)` with `C = 1` for `p > p₀`, and the two are interchangeable up to shrinking `δ`.

## 1. The reformulation (PROVED)

**Lemma 1.1 (bias versus statistical distance).** For `A, B ⊆ F_p` nonempty, with `Z(A,B) := #{(a,b) ∈ A×B : a+b = 0} ≤ min(|A|,|B|)`,

    E[(−1)^{Ext_p(U_A,U_B)}] = −(S(A,B) + Z(A,B))/(|A||B|),   SD(Ext_p(U_A,U_B), U₁) = |S(A,B) + Z(A,B)| / (2|A||B|).

*Proof.* `(−1)^{Ext_p(a,b)} = −χ(a+b)` when `a+b ≠ 0` and `= −1 = −(χ(a+b)+1)` when `a+b = 0`; sum over `A×B`. `Z ≤ min(|A|,|B|)` since `b = −a` is determined by `a`. ∎ (900 exact checks, S2, both sign conventions.)

**Lemma 1.2 (constructive flat decomposition).** Let `K ≥ 1` be an integer and `P` a distribution on a finite set with `P(x) ≤ 1/K` for all `x`. Then `P = Σ_{i=1}^{r} λ_i U_{S_i}` with `|S_i| = K`, `λ_i > 0`, `Σλ_i = 1`, and `r ≤ |supp P|`.

*Proof.* Maintain a nonnegative measure `P` of total mass `m` (initially `m = 1`) with the invariant `P(x) ≤ m/K` for all `x`. If `m > 0`, sort `P(x₁) ≥ P(x₂) ≥ …`; the invariant forces `|supp P| ≥ K`. Let `S = {x₁,…,x_K}` and `t := min(K·P(x_K), m − K·P(x_{K+1}))` (with `P(x_{K+1}) := 0` if the support has exactly `K` points). Then `t > 0`: `P(x_K) > 0`, and `K·P(x_{K+1}) = m` would force `P(x₁) = … = P(x_{K+1}) = m/K`, total `> m`. Replace `P` by `P − t·U_S` and `m` by `m − t`. Nonnegativity on `S` is `t ≤ K·P(x_K)`; the invariant off `S` is `t ≤ m − K·P(x_{K+1})`; the invariant on `S` is `P(x) − t/K ≤ (m−t)/K ⇔ P(x) ≤ m/K`, already true. Each step either annihilates `x_K` or raises `x_{K+1}` to the new maximal level `(m−t)/K`; atoms at the maximal level stay at the maximal level (they lose exactly `t/K`) and atoms at `0` stay at `0`; at most `K` atoms can sit at the maximal level (else the mass exceeds `m`), and if all `K` top atoms are maximal then `m = K·(m/K)` is exhausted by `S`, so `t = m` and the algorithm stops. Hence the potential "(#zero atoms) + (#maximal atoms)" strictly increases, and the algorithm terminates within `|supp P|` steps with `m = 0`, which is the claimed decomposition. ∎ (320 exact decompositions in rationals, S3, including boundary cases with atoms exactly at `1/K`.)

**Lemma 1.3 (bilinearity).** `β(P,Q) := Σ_{x,y} P(x)Q(y)(−1)^{Ext_p(x,y)}` is bilinear; hence if `P = Σλ_iU_{S_i}` and `Q = Σμ_jU_{T_j}` then `β(P,Q) = Σ_{i,j}λ_iμ_j β(U_{S_i},U_{T_j})` and `|β(P,Q)| ≤ max_{i,j}|β(U_{S_i},U_{T_j})|`. ∎ (320 exact checks, S3.)

**Lemma 1.4 (level sets).** For `A ⊆ F_p`, `t ≥ 0`, put `B⁺_t = {b : f_A(b) > t}`, `B⁻_t = {b : f_A(b) < −t}`. Then `S(A,B⁺_t) > t|B⁺_t|` and `S(A,B⁻_t) < −t|B⁻_t|` whenever the sets are nonempty, so `|B⁺_t| + |B⁻_t| ≤ (S(A,B⁺_t) − S(A,B⁻_t))/t` for `t > 0`. ∎ (736 checks, S4.)

**Theorem 1 (exact equivalence).** For `0 < ε < 1` consider the statements

- `C(ε)`: there are `δ > 0`, `p₀` such that for all `p > p₀` and all `A,B ⊆ F_p` with `|A|,|B| > p^ε`, `|S(A,B)| ≤ p^{−δ}|A||B|` (the two-set conjecture at rate `ε`);
- `F(ε)`: there are `δ > 0`, `p₀` such that for all `p > p₀`, `Ext_p` is an `(ε log₂p, p^{−δ})`-extractor for pairs of *flat* sources;
- `G(ε)`: the same for all pairs of independent sources of min-entropy `≥ ε log₂p`;
- `T(ε)`: the same, in the strong sense.

Then (1) `F(ε) ⇒ C(ε)`, and `C(ε′) ⇒ F(ε)` for every `ε′ < ε`; (2) `G(ε) ⇒ F(ε)`, and `F(ε′)` for some `ε′ < ε` implies `G(ε)`; (3) `T(ε) ⇒ G(ε)`, and `C(ε₁)` for some `ε₁ < ε` implies `T(ε)`. Consequently the following are equivalent: the conjecture for all `ε`; `Ext_p` is a flat-source extractor for every linear rate with polynomial error; `Ext_p` is an extractor for every linear rate with polynomial error; `Ext_p` is a strong extractor for every linear rate with polynomial error. (One may also take rate `ε` to mean `H₂ ≥ ε log₂p`; everything below goes through with `cp` in place of `max`, and Theorem 3 is stated in that form.)

*Proof.* (1) `F(ε) ⇒ C(ε)`: if `|A|,|B| > p^ε` then `U_A, U_B` have min-entropy `> ε log₂p`, so by Lemma 1.1 `|S(A,B)| ≤ 2|A||B|·p^{−δ} + Z ≤ (2p^{−δ} + p^{−ε})|A||B| ≤ p^{−δ″}|A||B|` for any `δ″ < min(δ,ε)` and `p` large. `C(ε′) ⇒ F(ε)`: flat sources of min-entropy `≥ ε log₂p` have `|A|,|B| ≥ p^ε > p^{ε′}` for `p > p₀`, and Lemma 1.1 gives `SD ≤ (p^{−δ}|A||B| + min(|A|,|B|))/(2|A||B|) ≤ ½(p^{−δ} + p^{−ε}) ≤ p^{−min(δ,ε)}`. (2) `G ⇒ F` is trivial. `F(ε′) ⇒ G(ε)`: let `K = ⌊p^ε⌋`; a source of min-entropy `≥ ε log₂p` has `max P ≤ p^{−ε} ≤ 1/K`, so Lemma 1.2 writes it as a convex combination of flat sources on sets of size `K > p^ε − 1 > p^{ε′}` for `p` large; Lemma 1.3 bounds `|E(−1)^{Ext}|` by the flat maximum, which is `≤ 2p^{−δ}` by `F(ε′)`; so `SD ≤ p^{−δ}`. (3) `T(ε) ⇒ G(ε)`: `SD(Ext(X,Y),U₁) ≤ E_{y∼Y} SD(Ext(X,y),U₁) ≤ η + η` by convexity of `SD`. `C(ε₁) ⇒ T(ε)`: let `X, Y` have min-entropy `≥ ε log₂p`, decompose `X = Σλ_iU_{A_i}` with `|A_i| = K = ⌊p^ε⌋ > p^{ε₁}` (Lemma 1.2), and put `t = p^{−δ}K`. If `|B⁺_t(A_i)| > p^{ε₁}` then `C(ε₁)` gives `|S(A_i,B⁺_t)| ≤ p^{−δ}K|B⁺_t| = t|B⁺_t|`, contradicting Lemma 1.4; so `|B^±_t(A_i)| ≤ p^{ε₁}`, and `Pr_Y[B^±_t(A_i)] ≤ 2p^{ε₁}·p^{−ε}`. Off these sets `|f_{A_i}(y)| ≤ t`, so by Lemma 1.1 `SD(Ext(U_{A_i},y),U₁) ≤ ½(p^{−δ} + p^{−ε})`. Hence `E_y SD(Ext(X,y),U₁) ≤ Σλ_i E_y SD(Ext(U_{A_i},y),U₁) ≤ η₀ := ½(p^{−δ}+p^{−ε}) + 2p^{ε₁−ε}`, and Markov gives `Pr_y[SD(Ext(X,y),U₁) > √η₀] ≤ √η₀`, i.e. strongness with error `p^{−δ‴}`, `δ‴ < ½min(δ, ε−ε₁)`; the same argument with the roles of `X` and `Y` exchanged. ∎

**Remark 1.5 (what this is and is not).** The implication `C ⇒ G` is the content of Satake 2024, Proposition 15 (attributed there to Chor–Goldreich 1988, Lemma 5 and Corollary 11; the original SICOMP paper was not fetched, so those internal numbers are UNVERIFIED at source), and Satake's Conjecture 9 (attributed to Chor–Goldreich) is exactly `G(ε)` for every `ε`. Theorem 1 adds the converse, the strong form, and explicit constants. It changes nothing about the truth of the conjecture; its only use is that it opens the two-source-extractor literature on the "half barrier" (Satake's Problem 8) to the `F_{p²}` test of the barriers note. Rao 2007, Theorem 5.1 (an argument of Barak) is the general form of (3): any two-source extractor with error `η` at entropy `k` is strong at entropy `k′` with error `2^m(η + 2^{k−k′})`.

## 2. Literature (CITED; every statement below was read in the fetched source named)

Sources fetched: Satake, *On the Paley RIP and Paley graph extractor*, arXiv:2405.08608 (HTML v1); Rao, *An exposition of Bourgain's 2-source extractor*, ECCC TR07-034, PDF from the author's page (text-extracted); Lewko, *An explicit two-source extractor with min-entropy rate near 4/9*, arXiv:1804.05451 (ar5iv), Mathematika 65 (2019) 950–957; Barak–Impagliazzo–Wigderson, *Extracting randomness using few independent sources*, FOCS 2004 / SICOMP 36 (2006) 1095–1118, PDF from Wigderson's page (text-extracted); Chattopadhyay, *A recipe for constructing two-source extractors*, SIGACT News survey PDF (text-extracted); ECCC TR15-119 (Chattopadhyay–Zuckerman), TR16-088 (Ben-Aroya–Doron–Ta-Shma), TR15-095 (Cohen) abstract pages; arXiv:2303.06802 (Li) abstract; arXiv:2506.15547 (restating Raz 2005) HTML. Not fetched: Chor–Goldreich SICOMP 17 (1988) 230–261; Bourgain, IJNT 1 (2005) 1–32; Raz STOC 2005 — for these the statements are taken from the secondary sources named and labelled accordingly.

| work | function analysed | min-entropy / error | Paley itself? | prime-field-only input | `F_{p²}` test |
|---|---|---|---|---|---|
| Chor–Goldreich 1988 (via Satake Prop. 15, Conj. 9; Rao Thm 3.1; survey p. 2) | inner product over `F_2^n` (Lindsey) and the Paley graph extractor `Ext_p` (introduced there) | `k₁ + k₂ > n`: rate `> 1/2` each; error `2^{−Ω(n)}` | yes | none: Lindsey/Chung, field-agnostic | holds over `F_{p²}` (consistent with being stuck at 1/2) |
| Karatsuba regime (survey p. 2–3; BIW footnote 7) | `Ext_p` | one source `(1/2+δ)n`, the other `C log n` (survey) / `|A| > p^{1/2+ε}, |B| > p^ε` for power-saving error (brief) | yes | none: Hölder + Weil | holds over `F_{p²}` |
| Barak–Impagliazzo–Wigderson 2004/06 | iterated `(a,b,c) ↦ ab + c` on `ℓ = (1/δ)^{O(1)}` sources; final step Lemma 3.2 (9 sources at `0.9 log|F|`) or, per footnote 17, any rate-`>1/2` two-source extractor or the leftover-hash lemma | every source rate `δ`; error `2^{−Ω(n)}` (Theorem 1.1) | no | Lemma 3.1 (`A·B + C` is `2^{−εm}`-close to min-entropy `min{(1+ε)m, 0.9 log|F|}`) uses Konyagin's Theorem 1.5 at exactly one place (§3.3, so stated by the authors); Theorem 1.4 (BKT) needs no subfield of size in `[|F|^{δ/2}, |F|)` | fails: `A = B = C = F_p ⊂ F_{p²}` gives `A·B + C = F_p`, no growth |
| Bourgain 2005 (via Rao Thm 3.8, Lemmas 3.3, 3.4, 3.7, Cor. 3.5; Lewko §1) | `Had(Enc x, Enc y)` with `Enc(x) = (x, x²) ∈ F_p²`, i.e. `xy + x²y²`, then a bit via `σ` (Rao Rem. 3.2, Lemma 4.1); Lewko: `ρ(x·y + (x·x)(y·y))` on `F²`, `ρ(x) = sign sin(2πσ(x))`, and `p ≡ 3 (mod 4)` (Lewko: needed because otherwise the zero set contains the line `{(t,it)}`) | two `(n,(1/2−γ)n)` sources, `γ` a universal constant; error `2^{−Ω(n)}`, `Ω(n)` bits | no | Rao Thm 2.17 (BKT/Konyagin line–point incidences in `F_p`, `K ≤ p^{2−β₀}` points and lines give `O(K^{3/2−α})` incidences), through Cor. 3.5 and Lemma 3.7 (`3X` is `|F|^{−Ω(1)}`-close to min-entropy `(1/2+γ)log|F²|` when `X` has `(1/2−γ)log|F|`); Thm 2.18 is the weaker `F_{2^p}` version | fails: for `X = U_{F_p}` in `F_{p²}`, `Enc(X)` and all `kX` live in `F_p² ⊂ F_{p²}²`, so `H_∞(3X) ≤ 2log₂p = ½log₂|F_{p²}²|`; also `−1` is a square in `F_{p²}` |
| Raz 2005 (via arXiv:2506.15547 Lemma 3; survey p. 3; Rao fn. 1) | `Ext(x,y)_i = G(x)_{(i,y)}`, `G` the small-bias generator `G(β,ν)_α = ν·Σ_{i<p′}(αβ)^i` (a fixed algebraic function: the `y`-th coordinate of an `ε`-biased encoding of `x`) | `k₁ ≥ (1/2+δ′)n₁ + 3log n₁ + log n₂`, `k₂ ≥ 4 log(n₁−k₁)` (as restated there; original UNVERIFIED), output `Ω(δ′·min(n₁/8, k₂/16))`, error `2^{−3m/2}` | no | none (small-bias + XOR lemma); the long source still needs rate `> 1/2` | not applicable; the `> 1/2` requirement is the same barrier |
| Lewko 2019 | Bourgain's `F²` map (Theorem 2) and its `F³` variant `ρ(x·y + (x·x)(y·y))` (Theorem 3) | rate near `21/44` resp. `4/9`; error `N^{−δ}` (`N = |F|^d`) | no | Lemma 6: `max_{λ≠0}|Σ a(x)b(y)e(λx·y)| ≤ |A|^{1/2}|B|^{1/2}|F|^{n/8}(Λ(A)Λ(B))^{1/8}` (additive energies `Λ`), fed by Theorem 9 (Rudnev–Shkredov): `A ⊂ P₃`, `|A| ≤ |F|^{26/21}`, `−1` non-square ⇒ `Λ(A) ≲ |A|^{17/7}`; `B ⊂ P₄`, `p^{4/3} ≤ |B| ≤ p²` ⇒ `Λ(B) ≲ |B|^{5/2}` (Rudnev's point–plane theorem, prime fields); Prop. 8: energy exponent `α` on `P_{d+1}` gives rate `(d+1)/(d(8−2α))`, and indeed `3/(2(8−34/7)) = 21/44`, `4/(3(8−5)) = 4/9` (arithmetic checked) | fails: the energy theorems are prime-field statements; sums of subfield points stay in the subfield |
| Chattopadhyay–Zuckerman 2016/19 (ECCC TR15-119 abstract; Annals 189) | resilient function composed with a non-malleable-extractor reduction | `log^C n` each; error `n^{−Ω(1)}`; one bit | no (not algebraic) | none; no field | not applicable |
| Cohen 2016/19 (ECCC TR15-095 abstract) | two-source *dispersers* via correlation breakers / challenge–response | polylog entropy | no | none | not applicable |
| Ben-Aroya–Doron–Ta-Shma 2016/17 (ECCC TR16-088 abstract) | Cohen's non-malleable extractor + a somewhere-random condenser with small entropy gap | `(log n)^{1+o(1)}` each; constant error (survey p. 3) | no | none | not applicable |
| Li 2023 (arXiv:2303.06802 abstract) | asymptotically optimal seeded non-malleable extractors ⇒ two-source extractors | `O(log n)` each | no | none | not applicable |
| Dodis–Li–Wooley–Zuckerman 2014 (via survey §3.1) | `Ext_p` as a *non-malleable* extractor | `k ≥ (1/2+δ)n`, error `2^{−Ω(n)}` | yes | Weil | rate `> 1/2` only |
| Satake 2024 | `Ext_p` (Definition 13), conditional on RIP of the Paley ETF `Φ_p`, `p ≡ 1 (mod 4)` | Theorem 18: `(p^{1/2+ε}, p^{−τ})`-RIP for some `τ > 0` ⇒ `𝒫(1/2−τ+γ, β)` for `0 < γ < τ`, some `β = β(ε,γ) > 0`; Corollary 21: then `Ext_p` is a `(αn, p^{−β})`-extractor for some `α < 1/2`; Corollary 31: Conjecture 29 (`δ_K = O(√(K/p)·log K·log p)`) ⇒ `Ext_p` is a `(αn,p^{−β})`-extractor for every `α < 1/2` | yes | none in the implication: Lemma 28 (`(K,δ)`-RIP ⇒ `|Σ_{u,v∈U}χ(u−v)| ≤ δ√p|U|` for `|U| ≤ K`) is linear algebra; the difficulty sits in the RIP hypothesis | passed in the right way: over `F_{p²}` the *hypothesis* is false, since `U = F_p` gives `Σ_{u,v∈F_p}χ(u−v) = p(p−1)` (every nonzero element of `F_p` is a square in `F_{p²}`), exceeding `δ√q·p` for `q = p²` |

Two remarks. (i) Satake 2020 (arXiv:2011.02907, Theorem 10; barriers note B7) proves the converse direction, `𝒫 ⇒` RIP beyond `√p` with exponent loss; together with Theorem 18 the RIP formulation and `𝒫` are equivalent up to exponents, so the RIP language is a third reformulation, not an independent input. (ii) Hanson's three-source result (`Σχ(a+b+c) = o(|A||B||C|)` at sizes `≥ δ√p`; sum-product note §1.3, arXiv:1509.04354 Theorem 1, not re-fetched here) is the additive analogue of §3.6 below and also sits at rate exactly `1/2`.

**Consolidated finding.** Every published two-source extractor for rate `< 1/2` either (a) applies the Hadamard/inner-product function to an *encoded* source on a paraboloid `P_{d+1} ⊂ F^{d+1}` (Bourgain, Lewko), (b) consumes at least three independent sources through `ab + c` (BIW and the multi-source line), or (c) is not an algebraic function at all (CZ and successors). None analyses `χ(x+y)` below rate 1/2, and the only paper that addresses `Ext_p` itself below 1/2 (Satake 2024) is conditional on an RIP statement equivalent up to exponents to the conjecture.

## 3. Chung's bound as an extractor bound, and the composition question (Task 3)

### 3.1 Theorem 3 (weighted Chung; PROVED)

Let `M = (χ(x+y))_{x,y∈F_p}`. Then `M` is symmetric and `M² = MMᵀ = pI − J` (`J` the all-ones matrix): `Σ_y χ(x+y)χ(x′+y) = p − 1` if `x = x′` and `−1` otherwise (3,345 entries checked). Hence for all real weights `α, β` on `F_p`,

    (Σ_{x,y} α(x)β(y)χ(x+y))² ≤ ‖α‖₂²·(p‖β‖₂² − (Σβ)²) ≤ p‖α‖₂²‖β‖₂²,

by Cauchy–Schwarz and `‖Mβ‖² = βᵀ(pI−J)β` (432 checks with signed and nonnegative integer weights). For sources `X, Y` this reads `|E χ(X+Y)| ≤ √p·√(cp(X)cp(Y)) = √p·2^{−(H₂(X)+H₂(Y))/2}`, non-trivial exactly when `H₂(X) + H₂(Y) > log₂p`. With Lemma 1.1, `Ext_p` is a `((1/2+δ)log₂p, ½(p^{−δ}+p^{−1/2−δ}))`-extractor for arbitrary (not only flat) sources. This is the classical Chor–Goldreich regime (Rao's Theorem 3.1 is the same computation for the Hadamard kernel over `F^l`), and it is the whole content of Chung's bound in extractor language.

### 3.2 The composition template, and why none of it touches Paley

A *condenser statement* has the shape: for independent sources `X₁,…,X_m` on `F_p` of min-entropy `≥ m₀` and a fixed map `f`, `f(X₁,…,X_m)` is close to a source of min-entropy `(1+c)m₀` (BIW Lemma 3.1 with `f = x₁x₂ + x₃`; Rao Lemma 3.7 with `f =` the threefold sum of `(x,x²)` into `F_p²`). Composed with Theorem 3 it bounds the bias of `χ(f(X) + g(Y))` by `√p·2^{−(H₂(f(X)) + H₂(g(Y)))/2}` plus the closeness error. For this to be a statement about `Ext_p` itself one needs `m = 1` and `f = g = id`.

**Proposition 6 (no single-source condenser; PROVED).** For every map `f : F_p → F_p` and source `X`: `max_z Pr[f(X) = z] ≥ max_x Pr[X = x]` and `cp(f(X)) ≥ cp(X)`. *Proof.* Merging atoms cannot decrease the largest atom or the sum of squares. ∎ (320 checks.)

So with two sources a condenser can enter a bound on `E χ(X+Y)` only through an inequality relating it to Paley-type expectations at *derived* sources. The literature uses three mechanisms: (i) Cauchy–Schwarz/Hölder squaring (Rao Lemma 3.3; Vinogradov, Karatsuba); (ii) Fourier expansion with additive energy (Lewko Lemma 6); (iii) exact identities from multiplicativity (dilations). Theorems 4 and 5 and Proposition 7 show that for the kernel `χ(x+y)` each mechanism is entropy-neutral or dominated by Chung.

### 3.3 Theorem 4 (the Fourier side is capped at Chung; PROVED)

Let `G = Σ_t χ(t)ζ^t` (`ζ = e^{2πi/p}`), `1̂_A(t) = Σ_{a∈A}ζ^{ta}`, and `E⁺(A,B) = #{(a,b,a′,b′) : a+b = a′+b′}`, `E⁺(A) = E⁺(A,A)`. Then

    G·S(A,B) = Σ_{t≠0} χ(t) 1̂_A(t) 1̂_B(t),   |G|² = p,   Σ_t |1̂_A(t)|²|1̂_B(t)|² = p·E⁺(A,B),

and therefore

    S(A,B)² ≤ p·E⁺(A,B) ≤ p·(E⁺(A)E⁺(B))^{1/2},   with   p·E⁺(A,B) ≥ p|A||B| ≥ |A||B|(p−|B|).

*Proof.* `Σ_t χ(t)ζ^{tu} = χ(u)G` for all `u` (substitute `t ↦ t/u` for `u ≠ 0`; both sides vanish at `u = 0`), which gives the first identity on expanding `1̂_A1̂_B`; `|G|² = p` is standard; the third identity is Parseval for the convolution `1_A*1_B`. Then `|S| ≤ p^{−1/2}Σ_t|1̂_A(t)||1̂_B(t)| ≤ p^{−1/2}·√p·(Σ_t|1̂_A|²|1̂_B|²)^{1/2}` by Cauchy–Schwarz with the constant vector, and `Σ_t|1̂_A|²|1̂_B|² ≤ (Σ|1̂_A|⁴)^{1/2}(Σ|1̂_B|⁴)^{1/2} = p(E⁺(A)E⁺(B))^{1/2}`. The lower bound is the diagonal `(a′,b′) = (a,b)`. ∎ The identities were checked exactly in `Z[ζ_p]` (integer vectors modulo `1+ζ+…+ζ^{p−1}`; 30 Gauss-sum, 60 Fourier and 120 Parseval identities for `p ≤ 61`), the inequalities (5,784 checks) on 1,446 pairs `(A,B)` from intervals, random sets, geometric progressions, greedy Sidon sets and subgroups at `p ≤ 8009`.

**Reading.** The Fourier route is the one through which "additive energy is small" (the output of every sumset condenser and of every incidence theorem) could enter; for Paley it delivers at best Chung, because the multiplier `χ(t)` is unimodular and `E⁺(A,B) ≥ |A||B|`. Numerically, the minimum over all tested pairs of `pE⁺(A,B)/(|A||B|(p−|B|))` is `8009/8005` (a pair of 4-element sets at `p = 8009` with `E⁺(A,B) = |A||B|`, i.e. no coincidences at all): even at its best the Fourier bound only reproduces Chung, and it exceeds `S²` by factors up to `3.5·10⁷`. Contrast Lewko's Lemma 6 for the Hadamard kernel: the same two applications of Cauchy–Schwarz give `|Σ_{A×B}e(λx·y)|⁸ ≤ |A|⁴|B|⁴|F|^nΛ(A)Λ(B)`, which *does* go below the barrier when the sets lie on a paraboloid `P_{d+1} ⊂ F^{d+1}` and entropy is measured in `F^d` (Prop. 8 of Lewko). Paley has no ambient dimension to trade: the source lives in `F_p` itself.

### 3.4 Theorem 5 (the squaring trick is entropy-neutral for Paley; PROVED)

For the Hadamard kernel Rao's Lemma 3.3 is `bias(X,Y)² ≤ bias(X−X, Y)`: one Cauchy–Schwarz replaces `X` by the *difference source* `X−X`, whose entropy grows for sources that "grow with addition"; this is the only place where BIW/BKT-type growth enters Bourgain's analysis (Rao §3.2.1–3.2.2). The verifier certifies this inequality exactly: with `h_y = Σ_x X(x)ζ^{xy}` and `w = Σ_y Y(y)|h_y|²`, the identity `W_Y·w − |Σ_{x,y}X(x)Y(y)ζ^{xy}|² = Σ_{y<y′}Y(y)Y(y′)|h_y − h_{y′}|²` holds in `Z[ζ_p]` (96 checks), and the right side is a sum of squared moduli.

For the Paley kernel the same step produces nothing new. For `x₁ ∈ F_p` and `d = x₂ − x₁ ≠ 0`, `x₁ + y ≠ 0`:

    χ(x₁+y)χ(x₂+y) = χ(d)·χ(d^{−1} + (x₁+y)^{−1}),

because `d^{−1} + (x₁+y)^{−1} = (x₂+y)/(d(x₁+y))` and `χ(1/z) = χ(z)`. Hence, for sources `X, Y`,

    |E χ(X+Y)|² ≤ Σ_{x₁} X(x₁)·Q(x₁),   Q(x₁) := Σ_{x₂,y} X(x₂)Y(y)χ(x₁+y)χ(x₂+y) = X(x₁)·Y(F_p∖{−x₁}) + T(x₁),
    T(x₁) = Σ_{u,v≠0} α_{x₁}(u)β_{x₁}(v)χ(u+v),   α_{x₁}(u) = χ(u)X(x₁+u^{−1}),   β_{x₁}(v) = Y(v^{−1}−x₁),

so `|E χ(X+Y)|² ≤ cp(X) + E_{x₁∼X}|T(x₁)|`, where `T(x₁)` is a signed Paley bilinear form whose weights are bijective relabellings of `X` and `Y`: `‖α_{x₁}‖₂² = cp(X) − X(x₁)²`, `‖β_{x₁}‖₂² = cp(Y) − Y(−x₁)²`, `‖α_{x₁}‖_∞ ≤ ‖X‖_∞`, `‖β_{x₁}‖_∞ ≤ ‖Y‖_∞`. (Exact checks: 6,300 identity checks over all `x₁` for `p ≤ 31`, random and structured weights; 48 checks of the inequality in rationals.) Applying Theorem 3 to `T` gives `|E χ(X+Y)|² ≤ cp(X) + √p·√(cp(X)cp(Y))`, weaker than Theorem 3 applied directly. Iterating the step only relabels again. Thus the derived sources have the entropies of the inputs (minus one atom), and every condenser statement about `X−X`, `X·X`, `XY+Z`, or `kEnc(X)` is inert for Paley through mechanism (i).

### 3.5 Proposition 7 (product condensers reach Paley only through dilations; PROVED, and `H_avg` REFUTED)

(i) For `C ⊆ F_p^*`: `Σ_{c∈C}Σ_{a∈A,b∈B}χ(ca + b) = Σ_{c∈C}χ(c)·S(A, c^{−1}B)` (multiplicativity of `χ`). (ii) `Σ_{a,b,c}χ(ab + c) = Σ_u r_{AB}(u)f_C(u)` with `r_{AB}(u) = #{ab = u}`, so by Theorem 3 `(Σ_{a,b,c}χ(ab+c))² ≤ E^×(A,B)·|C|(p−|C|)`; likewise `(Σ_{a,t,b}χ(a+t+b))² ≤ E⁺(A,T)·|B|(p−|B|)`. (iii) For a subgroup `H ≤ F_p^*`: `Σ_{a,b∈H,c∈C}χ(ab+c) = |H|·S(H,C)` and `E^×(H,H) = |H|³`: the product condenser does nothing for multiplicatively structured sources, and the three-source function `χ(ab+c)` contains the open shifted-subgroup problem (pass-summary §1.4). (iv) BIW's `A·B + C` does grow for generic sources (recorded: `log_p(1/cp(AB+C)) = 0.95–0.99` for random triples at rate `0.42–0.45`, `0.77–0.88` for triples of intervals or geometric progressions), and the six-source function `χ((ab+c) + (a′b′+c′))` obeys `bias ≤ √(E(AB+C)·(pE(A′B′+C′) − N⁶))/N⁶` by Theorem 3 (32 exact checks) — this is BIW's Theorem 1.1 with Paley as the final step, a statement about a six-variable function. (144 + 144 + 92 checks for (i)–(iii).)

(v) By (i), every quantity a product condenser can deliver is a signed or unsigned average over the dilation orbit `{S(A,cB)}_{c∈F_p^*}`; to reach `S(A,B)` one would need the specific dilate `c = 1` to be controlled by the orbit average. The hypothesis

    H_avg:  for all A, B ⊆ F_p^*,  |S(A,B)| ≤ 2·mean_{c∈F_p^*}|S(A,cB)|

is REFUTED. Witnesses (`N = ⌊p^{0.45}⌋`, `A` random, `B` the `N` elements of `F_p^*` maximising `f_A`; sets listed in the JSON):

| `p` | `N` | `S(A,B)` | `|S|/N²` | `mean_c|S(A,cB)|/N²` | `max_{c≠1}` | signed mean `|Σ_cχ(c)S(A,cB)|/((p−1)N²)` | random `C`, `|Σ_{c∈C,a,b}χ(ca+b)|/N³` |
|---|---|---|---|---|---|---|---|
| 503 | 16 | 141 | 0.551 | 0.0484 | 0.184 | 0.00125 | 0.0139 |
| 1009 | 22 | 244 | 0.504 | 0.0367 | 0.169 | 0.00027 | 0.00413 |
| 2003 | 30 | 408 | 0.453 | 0.0264 | 0.116 | 0.00013 | 0.00207 |
| 4001 | 41 | 691 | 0.411 | 0.0194 | 0.0886 | 0.0000061 | 0.00030 |
| 8009 | 57 | 1221 | 0.376 | 0.0139 | 0.0628 | 0.0000153 | 0.00285 |

The ratio `|S(A,B)|/mean` grows (`11, 14, 17, 21, 27`), so no fixed constant in `H_avg` survives, and the three-source sums of (i)–(ii) are two orders of magnitude below the two-source bias they would have to dominate. This is the dilation invariance of the crux note (Proposition 4.6 there: every termwise-Weil bound assigns the same value to all `p−1` dilates) in extractor language: the product condenser produces orbit averages, and `c = 1` is invisible to them.

### 3.6 The additive condenser and Hanson's regime

`Σ_{t∈T}S(A+t,B) = Σ_{a,t,b}χ(a+t+b)` is the three-source function of the sumset condenser `A ↦ A+T`. By (ii) its bias is `≤ √(E⁺(A,T)|B|(p−|B|))/(|A||T||B|)`, and since `E⁺(A,T) ≥ |A||T|` this is non-trivial only when `|A||T||B| > p` — total rate `1/2`, which is exactly Hanson's regime for `χ(a+b+c)` (no power saving). Below it, finite witnesses show the three-source function can be constant: at `p = 31391` (`n_p = 31`) with `A = T = B = [1,10]` (rate `0.222`), `Σχ(a+t+b) = 1000 = N³`; at `p = 18191`, `N = 9` (rate `0.224`), the sum is `729`; and for two intervals `S([1,N],[1,N]) = N²` at `p = 479, N = 6` (rate `0.290`), `p = 18191, N = 14`, `p = 31391, N = 15` (148,931 least-non-residue computations to `2·10⁶`, 10 witnesses verified exactly). These are finite witnesses only; their asymptotic version is Vinogradov's conjecture (crux note, Proposition 1.1).

### 3.7 Answer to Task 3

No published condenser argument can be composed with Chung's bound to give `Ext_p` itself a non-trivial bias bound at any rate `1/2 − c`, for the following exact reasons. (1) A single source cannot be condensed (Proposition 6), so a condenser must consume additional independent sources, and the composed function is then `χ(f(x₁,…,x_m) + g(y₁,…,y_n))`, not `χ(x+y)`. (2) The three mechanisms through which a multi-source condenser could re-enter a bound on `E χ(X+Y)` are inert for the Paley kernel: the Fourier/energy route is capped at Chung (Theorem 4); the Cauchy–Schwarz squaring route relabels the sources bijectively (Theorem 5); the multiplicative route yields only dilation-orbit averages, and the orbit average does not control the orbit member `c = 1` (Proposition 7, `H_avg` refuted). (3) The compositions that *do* work are written out exactly in §3.5(iv) (six-source Paley composition of BIW, rate `< 1/2` per source via Lemma 3.1) and in §2 (Bourgain/Lewko: Hadamard of the paraboloid encoding, where the entropy gain comes from the dimension shift `F^d → F^{d+1}` plus an incidence theorem); each changes the function. (4) The `F_{p²}` test agrees with (2): Chung, Theorem 4 and Theorem 5 hold verbatim over `F_{p²}`, while every working composition invokes a prime-field statement (Rao Theorem 2.17; BIW Theorem 1.5, used once in §3.3 of BIW; Rudnev–Shkredov Theorem 9 in Lewko).

## 4. Verification (`results/sigma_extractor_2026_09_05.json`; 168,410 checks, 0 failures, 5.3 s, seed 20260905)

| section | what is checked exactly | count |
|---|---|---|
| S1 | `M² = pI − J` entrywise, `p ≤ 31`; weighted Chung with signed/nonnegative integer weights; flat Chung | 3,345 + 432 + 216 |
| S2 | Lemma 1.1 in both conventions (`χ(a+b)`, Satake's `χ(a−b)`), `Z ≤ min`, `a−b ↔ a+b` | 900 |
| S3 | Lemma 1.2 greedy decompositions (`K ≤ 5`, `p ≤ 37`, atoms at the boundary `1/K`), Lemma 1.3 bilinearity and maximum, Proposition 6 | 320 + 320 + 320 |
| S4 | Lemma 1.4 for all thresholds `t ≤ |A|`, `p ≤ 53` | 736 |
| S5 | `GḠ = p`, `G² = χ(−1)p`, `G·S = Σχ(t)1̂_A1̂_B`, Parseval quartic and mixed identities in `Z[ζ_p]` (`p ≤ 61`); Theorem 4 inequalities on 1,446 pairs from five families at `p ≤ 8009` | 30 + 60 + 120 + 5,784 |
| S6 | Rao Lemma 3.3 via the exact Lagrange identity in `Z[ζ_p]`; Theorem 5 identities for every `x₁` and the inequality in rationals (`p ≤ 31`) | 96 + 6,300 + 48 |
| S7 | `E⁺(Enc A) = 2|A|² − |A|` exactly for `Enc(a) = (a,a²)`; threefold sums recorded | 20 |
| S8 | Proposition 7 (i)–(iv), subgroup collapse, six-source composition, `H_avg` witnesses | 144 + 144 + 92 + 5 + 32 + 5 |
| S9 | least non-residues `n_p` for all primes `≤ 2·10⁶`; interval witnesses verified by `pow` | 148,931 + 10 |

Recorded illustrations (HEURISTIC, exact numbers): for `Enc(A) ⊂ F_p²` with `|A| = p^{0.44–0.45}`, `H₂(3·Enc A)/log₂p = 1.01–1.11` against `H₂(2·Enc A)/log₂p = 0.77–0.81` for intervals, random sets, geometric progressions and Sidon sets — the growth that Rao's Lemma 3.7 asserts asymptotically, and which has no Paley counterpart because Paley has no encoding step.

## 5. Obligations (OPEN)

1. Fetch Chor–Goldreich 1988 and Bourgain 2005 and pin the internal statement numbers (Satake's "Lemma 5, Corollary 11"; Bourgain's theorem for `xy + x²y²`). The statements themselves are secured by three independent secondary sources each.
2. Theorems 4–5 and Proposition 7 cover the three mechanisms found in the literature (Fourier/energy, Cauchy–Schwarz/Hölder, multiplicative identities). A mechanism outside these — almost-periodicity (Croot–Sisask), or a "flattening" statement for the multiplicative convolution — could in principle relate `E χ(X+Y)` to a condensed source without passing through them; whether one exists is OPEN, and until a definition of "composition" wide enough to exclude it is fixed, §3.7 is a case analysis, not a theorem about all proofs.
3. Whether the Legendre-symbol output `χ(xy + x²y²)` of Bourgain's function is an extractor below 1/2: the Gauss-sum average over `t` loses a factor `√p` against Lewko-type bias bounds, so it is not implied by Theorems 2–3 of Lewko. OPEN.
4. Satake's Conjecture 29 (RIP constant of the Paley ETF) is unverified and, by Satake 2020 Theorem 10 and Satake 2024 Theorem 18, equivalent up to exponents to `𝒫`; it is not an independent route.
5. Nothing here improves any exponent: the conjecture, the subgroup case, and every route below rate 1/2 for two arbitrary sources remain OPEN exactly as in the pass summary.
