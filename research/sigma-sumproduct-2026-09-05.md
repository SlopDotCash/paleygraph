# Sigma pass, direction `sumproduct`: the additive-combinatorial (Burgess–Chang / sum-product) route to S(A,B) for arbitrary sets (2026-09-05)

**Status:** PROVED: Theorem E, a fully explicit form of the Burgess–Chang inequality chain with the expansion factor λ = |B − Z·I|/|B| and the ratio-energy deviation Φ = p²Σν²/(Σν)² tracked (Lemmas A–D, G, G′; Weil CITED); Corollary 3.1, the exponent bookkeeping θ = 1/(8−a−b) locating Chang's 4/9 exactly at her energy exponent 11/4 (and 2/5 at the AMRS exponent 5/2, 1/3 at the optimal exponent 2); Proposition F, a no-go: for every B, every Z, I, r, and *any* bound on the ratio energy, the chain's bound is ≥ |A||B|·|B|/E⁺(B)^{1/2}, so for additively unstructured B (E⁺(B) ≤ C|B|²) the chain cannot save more than √C. REFUTED with witnesses (p = 2003, |B| = 31, exhaustive over all singleton shifts): the hypothesis that some shift structure and a strong enough incidence/energy bound make the chain non-trivial for two arbitrary sets of size p^{α}, α < 1/2. CITED and settled: Chang's 4/9 theorem has the hypothesis |B+B| < K|B|; no published result gives θ < 1/2 for two arbitrary sets (Fouvry–Shparlinski–Xi 2024: Vinogradov and Karatsuba "still stand"). Numerically (737 exact checks, all pass): the energy deviation Φ is 1–26 for every tested family, so energy is never the bottleneck; the expansion factor λ is. OPEN: the conjecture itself; any route below θ = 1/2 for arbitrary sets.

## 0. Object and conventions

`p` an odd prime, `χ` the Legendre symbol mod `p` with `χ(0)=0`, so `χ(uv)=χ(u)χ(v)` for all `u,v ∈ F_p` and `χ(z)² = 1` for `z ≠ 0`.
For `A,B ⊆ F_p`: `S(A,B) = Σ_{a∈A}Σ_{b∈B} χ(a+b)`, `f_A(y) := Σ_{x∈A} χ(x+y)`.
`E⁺(X) = #{(x₁,x₂,x₃,x₄)∈X⁴ : x₁+x₂=x₃+x₄}`, `E^×(X) = #{x₁x₂=x₃x₄}`, `E^×(X,Y) = #{(x₁,x₂,y₁,y₂)∈X²×Y² : x₁y₁=x₂y₂}`.
An *interval* is `I = {1,…,H} ⊂ F_p`, `1 ≤ H < p`. `(2r−1)!! = 1·3·5⋯(2r−1)`.
Verifier: `experiments/sigma_sumproduct_2026_09_05.py` → `results/sigma_sumproduct_2026_09_05.json`. Every inequality below labelled PROVED was checked by exact integer computation there (see §5 for what was checked and the witnesses).

## 1. Literature pin (all CITED from local text extractions of the papers; statement numbers are those of the sources)

**1.1 Chang, Duke Math. J. 145 (2008), "On a question of Davenport and Lewis and new character sum bounds in finite fields", Theorem 6 (p. 5 and §3).** Exact statement: *Assume A, B ⊂ F_p such that (a) |A| > p^{4/9+ε}, |B| > p^{4/9+ε}; (b) |B+B| < K|B|. Then |Σ_{x∈A,y∈B} χ(x+y)| < p^{−τ}|A||B|, where τ = τ(ε,K) > 0, p > p(ε,K), χ non-principal.* The proof (her (3.10)–(3.24)) takes τ(ε) = ε²/(128K)-type and uses Freiman's theorem to place B in a proper GAP B₁ of dimension d ≤ K, a sub-progression B₀ and an interval I = [1,p^{δ}], δ = ε/(4d), with B − B₀I ⊂ B₂ := B₁ ∪ (−B₁). The energy input is her Proposition 1: for a box B (translated intervals allowed, max H_j < (√p−1)/2), E^×(B,B) < C^n (log p)|B|^{11/4}. **Settled: the small-doubling hypothesis (b) is part of the theorem; there is no arbitrary-two-set 4/9 theorem in the paper.** Her own introduction (p. 4–5): "Presently, such a result is only known (with no further assumptions) provided |A| > p^{1/2+δ} and |B| > p^{δ}. The problem is open even for the case |A| ∼ p^{1/2} ∼ |B|."

**1.2 Shkredov–Volostnov, arXiv:1606.00358 (2016/2017), Theorem 2 (main result).** *|A| > p^{12/31+δ}, |B| > p^{12/31+δ}, |A+A| < K|A|, |A+B| < L|B| ⇒ |Σχ(a+b)| ≪ (L log 2K/(δ log p))^{1/2}|A||B| for p > p(δ,K,L).* Two structural hypotheses (note |A+B| ≤ L|B| with |A| ≈ |B| forces |B−B| ≤ L²|B|²/|A| by Ruzsa's triangle inequality, so B is structured too); the saving is only logarithmic. Theorem 4 (three sets, |A|,|B|,|C| > p^{12/31+δ}, |A+A| < K|A|, |B+C| < L|B|) gives a power saving p^{−τ}, τ = δ²(log 2K)^{−3+o(1)}. Their Theorem 3 (Hanson's) is the three-arbitrary-set o(·) result. Their Theorem 10 is the AMRS collision bound (1.6 below) and Lemma 12 derives E^×(A) ≪ K^{3/2}|A|^{5/2} for |A±A| ≤ K|A|, |A|³K = O(p²).

**1.3 Hanson, arXiv:1509.04354, "Estimates for character sums with various convolutions".** Theorem 1: *A,B,C ⊂ F_p each of size ≥ δ√p ⇒ Σχ(a+b+c) = o_δ(|A||B||C|)* — three arbitrary sets, threshold exactly 1/2, no power saving (his proof uses Chang's theorem after BSG, hence needs E⁺ ≥ c|·|³). Theorem 2 / Corollary 1: the four-variable sum Σχ(a+b+cd) with |A|,|B|,|C|,|D| > p^{δ}, δ > 1/2 − 1/176, gets p^{−ε}. His statement of Chang's theorem (p. 4) carries the hypothesis |A+A| ≤ K|A|, and he writes (p. 5): "the proof of Theorem 1 relies on Chang's Theorem, which only allows one to estimate S_χ(A,B) past the square-root barrier under the hypothesis that |A+A| ≤ K|A|". Nothing for two arbitrary sets.

**1.4 Volostnov, Math. Notes 104 (2018) 197–203, arXiv:1712.09355** — as quoted in Schoen–Shkredov (1.5), Theorem 2: *|A|,|B| > p^{1/3+δ}, |A+A| < K|A|, |B+B| < L|B| ⇒ |Σχ(a+b)| ≪ p^{−δ²/C(K)}|A||B|.* Both sets structured.

**1.5 Schoen–Shkredov, arXiv:2004.01885 (2020), Theorem 3.** *|A| > p^{δ}, |B| > p^{1/3+δ}, |A||B|² > p^{1+δ}, |A+A| < K|A|, |B+B| < L|B| ≤ p^{δ/2}|B| ⇒ |Σχ(a+b)| ≪ exp(−c(δ⁴ log p/log² K)^{1/3})|A||B|.* Both sets structured. Their introduction: "Currently there are very few results regarding the above conjecture. The only affirmative answer was obtained in the following case |A| > p^{1/2+δ}, |B| > p^{δ}". Their Theorem 7 is Rudnev's point–plane theorem in the form |I(P,Π) − |P||Π|/p| ≪ |P|^{1/2}|Π| + k|Π| (|P| ≤ |Π|, k = max collinear points), Corollary 8 the derived ratio-energy bounds.

**1.6 Aksoy Yazici–Murphy–Rudnev–Shkredov, arXiv:1512.06613, Theorem 19 (image set theorem) with P = {(b,bc) : b∈B, c∈C}** (their (14) and Corollary 20's proof): *for A,B,C ⊆ F_p^* with |A||B||C| = O(p²), #{(b,b′,a,a′,c,c′) ∈ B²×A²×C² : b(a+c) = b′(a′+c′)} ≪ (|A||B||C|)^{3/2} + max(|A|,|B|,|C|)·|A||B||C|.* Built on Rudnev's Theorem 3 (arXiv:1407.0426): m points, n planes in P³, m ≥ n, n = O(p²), k = max collinear: I = O(m√n + km). Stevens–de Zeeuw (arXiv:1609.06284) Theorem 4: I(A×B, L) ≪ a^{3/4}b^{1/2}n^{3/4} + n for a ≤ b, ab² ≤ n³, an ≪ p².

**1.7 Petridis–Shparlinski, arXiv:1604.08469** — trilinear/quadrilinear sums with the additive character e_p(xyz) (Theorem 1.1 nontrivial for X = Y = Z ≥ p^{2/5}); not the multiplicative bilinear sum. **Shkredov–Shparlinski, arXiv:1803.08699, Theorem 1.1**: interval I = [1,X] × arbitrary S, bound in terms of E₃⁺(S,S,I); with the trivial energy bound it is non-trivial for S = X = p^{α}, α > 1/3 — one set an interval. **Alsetri–Shao, arXiv:2509.07765, Theorem 1.2**: for a GAP A ⊂ F_p of rank 2, E^×(A) ≪ (|A|² + |A|⁴/p) log p (the conjecturally optimal energy), and Theorem 1.1 (single sums over rank-2 GAPs, |A| ≥ p^{1/4+ε}). **Fouvry–Shparlinski–Xi, arXiv:2404.09295 (2024)**: trilinear/quadrilinear multiplicative character sums only; on the bilinear sum over arbitrary sets they state that Vinogradov's bound ‖α‖₂‖β‖₂p^{1/2} "has never been improved in full generality" and that the o(MN) conjecture "is far from proven, and the classical inequalities of Vinogradov and Karatsuba still stand".

**Conclusion of the pin (CITED).** For two *arbitrary* sets the best published θ is θ = 1/2: Vinogradov (|A||B| > p^{1+η}) and Karatsuba (|A| > p^{1/2+η}, |B| > p^{η}). Every sub-1/2 exponent (4/9 Chang; 12/31 Shkredov–Volostnov; 1/3 Volostnov, Schoen–Shkredov; 1/3 for interval × arbitrary) requires small additive doubling of at least one set, and 12/31 and 1/3 require it of both. The mission draft's suspicion is confirmed: Chang's 4/9 theorem is one-sided-structured, and no arbitrary-two-set version exists anywhere.

## 2. The chain, made fully explicit (PROVED)

Fix `A, B ⊆ F_p`, `Z ⊆ F_p^*` (nonempty), an interval `I = {1,…,H}`, an integer `r ≥ 1`. Define the *expanded domain* and *expansion factor*

    Y := B − Z·I = {b − zt : b∈B, z∈Z, t∈I},      λ := |Y|/|B| ≥ 1,

and the *ratio-energy weights* on `F_p²`

    ν(u₁,u₂) := #{(x₁,x₂,y,z) ∈ A×A×Y×Z : (x₁+y)/z = u₁, (x₂+y)/z = u₂},   Σν = |A|²|Y||Z|,
    Φ := p² Σν² / (Σν)²   (≥ 1 by Cauchy–Schwarz, since ν is supported on ≤ p² pairs).

**Lemma A (shift).** `Σ_{b∈B}|f_A(b)| ≤ (|Z|H)^{−1} Σ_{y∈Y} Σ_{z∈Z} Σ_{t∈I} |f_A(y+zt)|.`
*Proof.* For fixed (z,t) the map b ↦ b − zt is a bijection of B onto B − zt ⊆ Y, so Σ_{b∈B}|f_A(b)| = Σ_{y∈B−zt}|f_A(y+zt)| ≤ Σ_{y∈Y}|f_A(y+zt)|. Average over the |Z|H pairs (z,t). ∎

**Lemma B (Cauchy–Schwarz and the ratio identity).** `Σ_{y,z,t}|f_A(y+zt)| ≤ (|Y||Z|H)^{1/2}(Σ_{y,z,t}|f_A(y+zt)|²)^{1/2}` and, with `g(u₁,u₂) := Σ_{t∈I} χ(u₁+t)χ(u₂+t)` (real),

    Σ_{y∈Y,z∈Z,t∈I} |f_A(y+zt)|² = Σ_{(u₁,u₂)∈F_p²} ν(u₁,u₂) g(u₁,u₂).

*Proof.* The first is Cauchy–Schwarz. For the second, |f_A(y+zt)|² = Σ_{x₁,x₂∈A} χ(x₁+y+zt)χ(x₂+y+zt), and for z ≠ 0, χ(x+y+zt) = χ(z)χ((x+y)/z + t), so the summand equals χ(z)²χ(u₁+t)χ(u₂+t) = χ(u₁+t)χ(u₂+t) with u_i = (x_i+y)/z. Summing over t ∈ I first and then grouping the quadruples (x₁,x₂,y,z) by (u₁,u₂) gives the claim. ∎

**Lemma C (Hölder).** For ν ≥ 0 and any real g: `|Σ ν g| ≤ (Σν)^{1−1/r} (Σν²)^{1/(2r)} (Σ|g|^{2r})^{1/(2r)}.`
*Proof.* r = 1 is Cauchy–Schwarz. For r ≥ 2: Hölder with exponents 2r/(2r−1), 2r gives Σν|g| ≤ (Σν^{2r/(2r−1)})^{(2r−1)/2r}(Σ|g|^{2r})^{1/2r}; then ν^{2r/(2r−1)} = ν^{(2r−2)/(2r−1)}·ν^{2/(2r−1)} and Hölder with exponents (2r−1)/(2r−2), 2r−1 gives Σν^{2r/(2r−1)} ≤ (Σν)^{(2r−2)/(2r−1)}(Σν²)^{1/(2r−1)}; raise to (2r−1)/2r. ∎

**Lemma D (Weil moment).** `W_r := Σ_{(u₁,u₂)∈F_p²} g(u₁,u₂)^{2r} ≤ (2r−1)!! H^r p² + (2r−1)² H^{2r} p.`
*Proof.* Expanding, W_r = Σ_{t⃗∈I^{2r}} (Σ_{u∈F_p} χ(F_{t⃗}(u)))² with F_{t⃗}(u) = Π_{i=1}^{2r}(u+t_i), because the sums over u₁ and u₂ factor and χ is real. If every value occurs an even number of times in t⃗ then F_{t⃗} is a square and |Σ_u χ(F)| ≤ p; every such t⃗ arises from a perfect matching of the 2r positions (≤ (2r−1)!! choices) and a value on each pair (≤ H^r), so there are ≤ (2r−1)!!H^r of them. Otherwise F_{t⃗} is not a constant times a square, it has m ≤ 2r distinct roots, and Weil's bound (CITED: Iwaniec–Kowalski Thm 11.23/Cor 11.24; Chang's Theorem W) gives |Σ_u χ(F)| ≤ (m−1)√p ≤ (2r−1)√p. ∎
(For r = 1 this is exact: W₁ = H(p−1)² + H(H−1), verified in §5.)

**Theorem E (explicit Burgess–Chang chain).** With the notation above, for every choice of A, B, Z, I, r:

    |S(A,B)|²  ≤  (|A|²|Y|²/H) · (Σν²/(Σν)²)^{1/(2r)} · ((2r−1)!! H^r p² + (2r−1)² H^{2r} p)^{1/(2r)},

equivalently

    |S(A,B)|  ≤  λ |A||B| · Φ^{1/(4r)} · [ ((2r−1)!!)^{1/(2r)} H^{−1/2} + (2r−1)^{1/r} p^{−1/(2r)} ]^{1/2}.       (E)

*Proof.* |S(A,B)| ≤ Σ_{b∈B}|f_A(b)|. Lemma A, then Lemma B, then Lemma C applied to Σνg with Lemma D for Σ|g|^{2r}:
|S|² ≤ (|Z|H)^{−2}·|Y||Z|H·(Σν)^{1−1/r}(Σν²)^{1/2r}W_r^{1/2r}. Insert Σν = |A|²|Y||Z| to get the first display. For the second, write Σν² = Φ(Σν)²/p², so (Σν²/(Σν)²)^{1/2r} = Φ^{1/2r}p^{−1/r}, and use (X+Y)^{1/2r} ≤ X^{1/2r}+Y^{1/2r}: W_r^{1/2r} ≤ ((2r−1)!!)^{1/2r}H^{1/2}p^{1/r} + (2r−1)^{1/r}Hp^{1/2r}; multiply out (the p^{−1/r} cancels p^{1/r} in the first term and leaves p^{−1/2r} in the second), divide by H, take square roots, and use |Y| = λ|B|. ∎

*Remarks.* (i) Nothing is assumed about A. (ii) Chang's proof of Theorem 6 is exactly (E) with Z = B₀ (a sub-progression), Y ⊆ B₂ = B₁ ∪ (−B₁) (so λ ≤ 2e^{C(K)}), H = p^{δ}, r = ⌊10/δ⌋, and her (3.20) as the bound on Σν². (iii) The saving in (E) is at most the bracket, i.e. at most max(H^{−1/4}, p^{−1/(4r)}) ≥ p^{−1/(4r)} up to constants; the two obstacles to making (E) non-trivial are therefore
    (Obs-λ)  λ^{4r} Φ ≤ p^{1−η}   (needed for the p^{−1/2r} term), and   H ≥ λ⁴ Φ^{1/r}·(const)·p^{4δ} (for the H^{−1/2} term).

**Lemma G (Chang's (3.20), the bound on Σν² through multiplicative energies; PROVED, elementary).** For any A, Y ⊆ F_p and Z ⊆ F_p^* (ν defined with these Y, Z):

    Σν² ≤ |A|³ · max_{x,x′∈A} #{(y,y′,z,z′)∈Y²×Z² : (x+y)/z = (x′+y′)/z′} ≤ |A|³ · E^×(Z)^{1/2} · max_{x∈A} E^×(x+Y)^{1/2}.

*Proof.* Σν² counts (x₁,x₂,y,z,x₁′,x₂′,y′,z′) with (x₁+y)/z = (x₁′+y′)/z′ and (x₂+y)/z = (x₂′+y′)/z′. Choose (x₁,x₁′) (≤ |A|² ways), then (y,y′,z,z′) satisfying the first equation (≤ the max), then x₂ ∈ A freely (|A| ways); x₂′ = z′(x₂+y)/z − y′ is determined. For the second inequality, #{(x+y)/z = (x′+y′)/z′} = Σ_w r_{(x+Y)/Z}(w) r_{(x′+Y)/Z}(w) ≤ E^×(x+Y, Z)^{1/2}E^×(x′+Y, Z)^{1/2} by Cauchy–Schwarz, and E^×(X,Z) ≤ E^×(X)^{1/2}E^×(Z)^{1/2} (Cauchy–Schwarz again; Chang's Fact 1, [TV] Cor. 2.10). ∎

**Lemma G′ (direct collision count).** `Σν² ≤ |A| · T(A,Y,Z)`, `T(A,Y,Z) := #{(x,x′,y,y′,z,z′)∈A²×Y²×Z² : (x+y)z′ = (x′+y′)z}` (same proof, without fixing x₁′ separately). By 1.6 (AMRS Theorem 19 with the grid P = {(z,zy)}), if |A||Y||Z| = O(p²) and 0 ∉ Y, then T ≪ (|A||Y||Z|)^{3/2} + max(|A|,|Y|,|Z|)|A||Y||Z|  (CITED).

## 3. Exponent bookkeeping (PROVED as arithmetic; inputs CITED)

Take |A| = |B| = p^{α}, and suppose the shift structure is *free*: λ = p^{o(1)}, |Z| = |B|^{1−o(1)}, |Y| = |B|^{1+o(1)} — this is exactly what Freiman + Chang's (3.13) deliver when |B+B| ≤ K|B|, and what an interval delivers explicitly (Z = {1,…,M}, Y ⊆ {1−MH,…,N}, λ ≤ 1 + MH/N). Suppose the ratio energy satisfies Σν² ≤ p^{o(1)} |A|^{a}|B|^{b}. Then Φ = p^{2+o(1)}|A|^{a−4}|B|^{b−4}, and (Obs-λ) becomes |A|^{4−a}|B|^{4−b} ≥ p^{1+η}, i.e.

    **Corollary 3.1.**  θ(a,b) = 1/(8−a−b)   (equal sizes);  in general the condition is |A|^{4−a}|B|^{4−b} > p^{1+η}.

| input for Σν² (through Lemma G or G′) | (a,b) | condition | θ |
|---|---|---|---|
| Chang Prop. 1: E^×(GAP) ≪ log p·|GAP|^{11/4} ⇒ Σν² ≪ |A|³|B|^{11/4} log p | (3, 11/4) | |A||B|^{5/4} > p^{1+η} | **4/9** |
| Lemma 12 of 1.2 (from AMRS): E^×(P) ≪ K^{3/2}|P|^{5/2} for |P±P| ≤ K|P|, |P|³K = O(p²) ⇒ Σν² ≪ |A|³|B|^{5/2} | (3, 5/2) | |A||B|^{3/2} > p^{1+η} | 2/5 |
| Lemma G′ + AMRS Thm 19 (needs |A||B|² = O(p²)) ⇒ Σν² ≪ |A|^{5/2}|B|³ + |A|²|B|²max(|A|,|B|) | (5/2, 3) | |A|^{3/2}|B| > p^{1+η} and |A|²|B|² > p·max | 2/5 |
| Alsetri–Shao Thm 1.2 (rank-≤2 GAPs, incl. translated APs): E^× ≪ (|P|²+|P|⁴/p) log p ⇒ Σν² ≪ |A|³|B|² log p for |B| ≤ √p | (3, 2) | |A||B|² > p^{1+η} | 1/3 |
| the diagonal floor (see §5, structured sets): Σν² ≈ |A|²|B|³ | (2, 3) | |A|²|B| > p^{1+η} | 1/3 |
| absolute floor Σν² ≥ Σν = |A|²|B|² and Σν² ≥ (Σν)²/p² | (2,2) | |A||B| > p^{1/2+η} | 1/4 |

So: **4/9 is exactly the value forced by the exponent 11/4 in Chang's multiplicative-energy bound for boxes, inserted through the Cauchy–Schwarz splitting of Lemma G; the AMRS point–plane bound moves it to 2/5, and the optimal (Szemerédi–Trotter-strength) energy bound — proved for rank ≤ 2 GAPs by Alsetri–Shao — moves it to 1/3, which is also where the literature's best structured results sit.** All of these rows presuppose the free shift structure, i.e. small additive doubling of B; none of them says anything about arbitrary B. The bottleneck for *structured* B is the multiplicative-energy bound for GAPs of rank ≥ 3 (open beyond 32/13, cf. Alsetri–Shao's discussion); the bottleneck for *arbitrary* B is (Obs-λ), which is a theorem-level obstruction (§4), not an energy estimate.

Note on one-sided structure (CONDITIONAL on the CITED inputs, not claimed as new): rows 2–4 with Chang's Freiman step verbatim give, for arbitrary A and |B+B| ≤ K|B|, |A|,|B| ≤ √p, a power saving as soon as |A||B|^{3/2} > p^{1+η} or |A|^{3/2}|B| > p^{1+η} (θ = 2/5); for B a translated interval (no Freiman needed, λ ≤ 1 + MH/N explicitly) and Alsetri–Shao's energy bound, as soon as |A||B|² > p^{1+η} (θ = 1/3, which for intervals is already in Chang Thm 8/9 and Shkredov–Shparlinski Thm 1.1). I have not found the 2/5 one-sided statement in print; it is an assembly of cited results, and I do not claim it beyond that.

## 4. Why no energy bound can help for arbitrary sets (PROVED)

**Proposition F (no-go for the chain).** For every A, B ⊆ F_p, Z ⊆ F_p^*, interval I, r ≥ 1, the right-hand side of (E) satisfies

    RHS(E)  ≥  |A||B| · |B| / E⁺(B)^{1/2}.

In particular if E⁺(B) ≤ C|B|² then RHS(E) ≥ |A||B|/√C: the chain cannot save more than a constant factor, whatever Z, I, r are and whatever bound (even the exact value) is used for Σν².

*Proof.* Let H′ := Z·I (as a set), so Y = B − H′ and |H′| ≥ H (one z gives H distinct products). Cauchy–Schwarz on the representation function of B − H′: |B||H′| = Σ_{d∈Y} r_{B−H′}(d) ≤ |Y|^{1/2}(Σ_d r_{B−H′}(d)²)^{1/2}, and Σ_d r_{B−H′}(d)² = #{b−h = b′−h′} = Σ_x r_{B−B}(x) r_{H′−H′}(x) ≤ E⁺(B)^{1/2}E⁺(H′)^{1/2} ≤ E⁺(B)^{1/2}|H′|^{3/2}. Hence |Y| ≥ |B|²|H′|^{1/2}/E⁺(B)^{1/2}, i.e. λ ≥ |B||H′|^{1/2}/E⁺(B)^{1/2}. In (E), Φ ≥ 1 and the bracket is ≥ H^{−1/2}, so RHS(E) ≥ λ|A||B|H^{−1/4} ≥ |A||B|·|B|·|H′|^{1/2}H^{−1/4}/E⁺(B)^{1/2} ≥ |A||B|·|B|/E⁺(B)^{1/2} since |H′| ≥ H ≥ 1. ∎

Since E⁺(B) ≥ 2|B|² − |B| for every B, √C ≥ √2 − o(1) is the best case, and sets with E⁺(B) ≤ 3|B|² (random sets of size ≤ √p, Sidon-type sets) are the generic case. The same argument applies verbatim to the mirror chain that shifts A instead of B (interchange A and B in Lemma A). The chain is therefore *structurally* blind to the incidence input for unstructured sets: the loss is incurred in Lemma A before any energy appears.

**Refuted hypothesis (with witness, §5.4).** *H_chain(α): for every A, B ⊆ F_p with |A|,|B| ≥ p^{α} there exist Z, I, r for which the bound (E) — with the exact Σν² — is ≤ p^{−δ}|A||B| for some δ = δ(α) > 0 and all large p.* REFUTED for every α < 1 by Proposition F with any B having E⁺(B) ≤ 3|B|²; explicit witnesses (p, B, E⁺(B)) are recorded in the JSON. This is precisely the "E* would give 4/9 − c for arbitrary sets" scenario: the answer is that no E* exists inside this chain, because the obstruction is not the energy.

**Consequence for Task 3.** No published incidence/energy bound gives θ < 4/9 for two arbitrary sets, and none can through this chain (Proposition F). For sets with a free shift structure (small doubling, one side), published bounds already give θ = 2/5 (AMRS) and, for rank-≤2 GAPs, θ = 1/3 (Alsetri–Shao); for arbitrary sets the published state is θ = 1/2 and the additive-combinatorial route as it stands does not reach below it.

## 5. Numerical verification (exact integer computation; `results/sigma_sumproduct_2026_09_05.json`, 737 checks, 0 failures, 10.7 s)

**5.1 Lemma D.** For p ∈ {101, 199, 293}, H ∈ {3,5,8}, r ∈ {1,2,3}: W_r computed exactly from the p×p table of g and compared with (2r−1)!!H^r p² + (2r−1)²H^{2r}p (27 cases, all hold; e.g. p=101, H=8: W_r/bound = 0.909, 0.340, 0.099 for r = 1,2,3). For r = 1 the exact formula W₁ = H(p−1)² + H(H−1) holds in all 9 cases; for r = 2, H ≤ 5 the t-vector expansion Σ_{t⃗}(Σ_uχ(F_{t⃗}(u)))² reproduces W₂ exactly (6 cases).

**5.2 Theorem E and Lemmas A, B, C, G, G′, and the Proposition F inequalities**, on explicit configurations at p ∈ {101, 197, 293} (intervals × intervals, random × interval, random × random, random × subgroup coset; r ∈ {1,2}), on 36 random configurations (p ≤ 293, |A|,|B| ≤ 25, |Z| ≤ 4, H ≤ 5, r ≤ 3, B random/interval/AP), and on the larger configurations of 5.4–5.5 (p ∈ {1009, 2003}, |A|,|B| up to 176). Checked exactly in every instance: Lemma A; the identity Σ_{y,z,t}|f_A(y+zt)|² = Σνg (integer equality); Σν = |A|²|Y||Z|; Lemma C on the actual (ν,g); Theorem E in both forms; Φ ≥ 1; Σν² ≤ |A|³·max_{x,x′}#{…} ≤ |A|³E^×(Z)^{1/2}max_x E^×(x+Y)^{1/2}; Σν² ≤ |A|·T; |Y| ≥ |B|²|H′|²/(E⁺(B)E⁺(H′))^{1/2}; RHS(E) ≥ |A||B|²/E⁺(B)^{1/2}. Over the 36-instance sweep, max |S(A,B)|/RHS(E) = 0.066 and min RHS(E)/(|A||B|) = 1.91: at these primes the chain never beats the trivial bound (as it must not: its saving is at most p^{−1/4r}·λ^{-1}·Φ^{1/4r} and λ ≥ 1.6 in every instance).

**5.3 Exponent arithmetic.** θ(3,11/4) = 4/9, θ(3,5/2) = 2/5, θ(5/2,3) = 2/5, θ(3,2) = 1/3, θ(2,3) = 1/3 (exact rationals); (4/9)(9/8) = 1/2 (Chang's (3.24)); (2r−1)!! ≤ r^r for r ≤ 12 (my Lemma D constant is never worse than Chang's/SV's r^{2r}).

**5.4 Proposition F witnesses (REFUTED hypothesis H_chain).** Random B (seed 1) at (p,|B|) = (1009,16), (1009,22), (2003,21), (2003,31) — the sets are listed in the JSON, e.g. p = 2003, B = {177, 221, 335, 751, 758, 767, 805, 898, 1002, 1041, 1066, 1123, 1359, 1381, 1439, 1500, 1511, 1588, 1594, 1720, 1807}: E⁺(B) = 933 = 2.116|B|², so RHS(E) ≥ 0.688·|A||B| for every Z, I, r and every A. Explicit trials (|Z|,H,r) = (1,2,1), (3,3,2), (6,8,3), (12,16,4) give λ = 2.0, 8.6, 38.6, 83.8 (lower bounds from the proof: 1.12, 4.34, 17.6, 29.2), Φ_true = 221, 21.9, 7.2 (energy near its floor), and RHS(E)/(|A||B|) = 6.6, 12.8, 46.3 — while the actual |S(A,B)|/(|A||B|) = 0.032. **Exhaustive check (5.9):** for p = 2003, B = {26, 106, 147, 208, 231, 284, 321, 326, 372, 443, 504, 510, 572, 761, 779, 797, 835, 960, 1113, 1175, 1265, 1329, 1497, 1560, 1568, 1632, 1776, 1778, 1845, 1876, 1931} (|B| = 31 = p^{0.45}), the minimum of λ = |B − z·{1,2}|/|B| over **all** z ∈ F_p^* is 1.903 (attained at z = 153; the trivial maximum is 2), so even the smallest conceivable shift structure expands B by a factor ≥ 1.9 and RHS(E) ≥ 1.60·|A||B|. Same at p = 1009, |B| = 16: min λ = 1.8125 (z = 281).

**5.5 λ and Φ across set types (the "E*" test).** p ∈ {1009, 2003}, |A| = |B| = n ∈ {63, 110, 96, 176} (n = p^{0.60}, p^{0.68}), |Z| = H = n/10, r = 2; B ∈ {interval, random, subgroup coset (order ≈ n), interval ∪ random}, A ∈ {random, interval, AP of step 7, coset}. Findings (all values in the JSON):
- Interval B with its natural shift structure Z = {1,…,n/10}: λ = 1.56–2.64 (= 1 + MH/N up to overlap, exactly as predicted), Φ_true = 1.07–1.85 for random/coset/AP A, and 3.65–26 for A also an interval. So for one-sided structure the ratio energy is at (or within a small constant of) its Cauchy–Schwarz floor; the requirement λ^{4r}Φ ≤ p^{1−η} fails only because H = n/10 is far too short for p ≈ 2000 (H ≥ p^{1/r} is needed), i.e. only because p is small — consistent with Theorem E being non-trivial for such B when p → ∞.
- Random B, coset B, mixed B: λ = 7.4–21.5 with |Z| = H = n/10 (and 9–12 at n = 176), Φ_true = 1.07–1.31. The energy is essentially minimal; the expansion factor is what kills the chain, exactly as Proposition F says.
- Size scaling at p = 2003 (5.7), Y := B, |Z| = n/10: local exponents of Σν² in n are 5.46, 5.14, 5.23 for interval×interval (the structured diagonal z=z′, x₁−x₁′ = x₂−x₂′ = y′−y, i.e. row (2,3) of §3 with |Z| ∝ n), and 4.6 → 5.7 → 7.4 for random×random and interval×random (crossing over from the trivial diagonal Σν ∝ n⁴ to the uniform term (Σν)²/p² ∝ n⁸ as n grows past p^{2/3}). No configuration approaches the exponents 11/4 or 5/2 of the published energy bounds; those bounds are far from sharp on every tested family. HEURISTIC: this is consistent with the conjectural optimal energy (row (3,2)/(2,3) of §3, θ = 1/3 for one-sided structure), but exponents measured at p ≤ 2003 cannot discriminate 1/3 from 2/5.
- Multiplicative energy of translated intervals x + {1,…,N}, N ≈ √p/2 (5.6): max over 62 translates x of E^×(x+I) is 1025 (p = 1009, N = 16) and 2305 (p = 2003, N = 22), i.e. 0.50·N^{11/4} and 0.47·N^{11/4} (Chang's Prop. 1 holds with room), 0.58·N² log p and 0.63·N² log p, and 3.2·(N²+N⁴/p), 3.8·(N²+N⁴/p) (Alsetri–Shao's form holds with constant ≤ 4 on this sample). The 11/4 exponent is visibly not sharp.

**Verdict on E*.** Within the chain, the only energy quantity that matters is Φ; in every tested family Φ_true ≤ 26 ≪ p^{1−η}, so an energy hypothesis "E*" strong enough for 4/9 − c is *satisfied* by the data for every set type — and is irrelevant for arbitrary sets, because Proposition F shows the loss occurs before the energy enters. There is no E* for arbitrary sets to refute or confirm: the chain is refuted for them outright (5.4, 5.9).

## 6. What is established, and open obligations

PROVED here: Lemmas A–D, G, G′, Theorem E, Corollary 3.1, Proposition F (Lemma D and Theorem E rely on Weil's bound, CITED). CITED: §1 statements; AMRS Theorem 19; Alsetri–Shao Theorem 1.2; Chang Proposition 1. CONDITIONAL (not claimed as new): the one-sided-structure exponents 2/5 (via AMRS) and 1/3 (via Alsetri–Shao, rank ≤ 2) obtained by inserting the cited energy bounds into Theorem E + Lemma G with Chang's Freiman step. REFUTED: H_chain(α) for every α < 1 (Proposition F; witnesses in 5.4/5.9). OPEN: the Paley graph conjecture; any θ < 1/2 for two arbitrary sets.

Obligations remaining:
1. The one-sided 2/5 statement should be checked against the literature once more (I found no printed version; if new, it needs Chang's Freiman step written out with the K-dependence tracked, which I have not done) — CONDITIONAL until then.
2. Proposition F is a no-go for Theorem E's chain and its mirror (shifting A). It does not cover chains that shift a *convolution* (Croot–Sisask almost-periodicity, as in Shkredov–Volostnov), but those need |A+S| ≤ K|A| in Lemma 9 of 1.2 to make the almost-period set large; whether an almost-periodicity argument can be made to work for a set with E⁺ ≈ 2|B|² is OPEN (I expect not: the almost-period set has size |A|(2K)^{−O(ε^{−2}q)} with K ≈ |A| for unstructured A).
3. For arbitrary sets the additive-combinatorial route offers, as far as I can see, only the Hanson-type dichotomy (small E⁺ ⇒ ? ; large E⁺ ⇒ BSG ⇒ structured case). The missing half — cancellation in S(A,B) from *small* additive energy of B, with two sets — has no known mechanism; the natural candidate inequality |S(A,B)|² ≤ |A|·Σ_{a∈F_p}|f_B(a)|² = |A||B|(p−|B|) contains no energy at all. This is the precise obstruction to any energy-based E* for two arbitrary sets.
4. The measured exponents in 5.7 are HEURISTIC evidence only; larger p (≥ 10⁴, with |A||B|² ≫ p) would be needed to fit exponents to ±0.1.
