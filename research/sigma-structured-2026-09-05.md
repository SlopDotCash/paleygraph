# Sigma pass, direction `structured`: one multiplicatively structured set (2026-09-05)

**Status:** PROVED: (i) the two-set bound for pairs (A, H) with H a multiplicative subgroup is *equivalent* (Theorem 2.2) to the uniform shifted-subgroup bound max_{c≠0}|Σ_{h∈H} χ(h+c)| ≤ p^{−δ}|H|, and that bound transfers by the triangle inequality to cosets and unions of k cosets with loss k (Task 2); (ii) an elementary symmetric-set bound M(B) ≤ 2t + ((p|B|−|B|²)/|Sym_t(B)|)^{1/2}, giving M_H ≤ (p−|H|)^{1/2} for subgroups; (iii) an exact multiplicative amplification inequality whose only gain factor is the multiplicative energy E^×(A,S), which is maximal (E^× = |A||S|²) on exactly the H-invariant sets A that are extremal in (i); (iv) the coset-decomposition identity of Task 4 and, with exact witnesses, why it cannot be closed by any hypothesis on full subgroups. REFUTED (with witnesses): "a shifted-subgroup bound for H controls dense subsets of H" (p = 433, |H| = 24, M_H = 8, but the index-3 subgroup H' ⊂ H has M_{H'} = |H'| = 8); "coset decomposition by a subgroup of size ≈ p^{ε/2} reduces a generic B" (p = 1000003, |B| = 251, |H| = 6: all 251 cosets are singletons). OPEN: the uniform shifted-subgroup bound for |H| = p^ε (it is Bourgain's Problem 5 in Chang's 2010 survey, stated as open even at |H| ~ √p; no source examined, including the 2024–2026 ones, contains it), hence the entire Task-3 class |BB| ≤ K|B| is OPEN already at K = 1, and it also contains, at K = 2, the open consecutive-powers problem (Shparlinski's Problem 1). No new cancellation exponent is claimed.

Verifier: `experiments/sigma_structured_2026_09_05.py` (220,917 exact integer checks, 0 failures, 16 s) → `results/sigma_structured_2026_09_05.json`.

## 0. Notation

p an odd prime, χ the Legendre symbol with χ(0) = 0. For B ⊆ F_p and r ∈ F_p,

    T_B(r) := Σ_{b∈B} χ(r+b),    S(A,B) := Σ_{a∈A} Σ_{b∈B} χ(a+b) = Σ_{a∈A} T_B(a),

    M(B) := max_{r≠0} |T_B(r)|,   M_{2k}(B) := Σ_{x∈F_p} T_B(x)^{2k},   M_H := M(H) for a subgroup H ≤ F_p^*.

For B ⊆ F_p^* and s ∈ F_p^*: sB = {sb}. H⁺ := {h ∈ H : χ(h) = 1} (H itself if H ⊆ QR, else a subgroup of index 2 in H).

Background facts used freely (brief): Σ_x T_B(x)² = p|B| − |B|² (verified exactly, §6); Weil: |Σ_x χ(f(x))| ≤ (m−1)√p for f with m distinct roots that is not a square.

## 1. Task 1: literature pin (CITED)

Sources: text extractions in the scratchpad `structured/` directory (arXiv numbers below) and `sources/sigma-chang-character-sums-survey-DubProc.txt`.

**Question.** Is it published that for every ε > 0 there is δ > 0 with |T_H(c)| ≤ p^{−δ}|H| for all subgroups H, |H| > p^ε, and all c ≠ 0?

**Answer: no.** The statement is posed as an open problem, and every source examined that gives an individual bound gives only the Weil-range bound √p.

- **Chang, "Character sums in finite fields", Contemp. Math. 518 (2010), Problem 5 (attributed to Bourgain):** "Obtain nontrivial bound on Σ_{x∈H} χ(a+x) for H < F_p^*, |H| ∼ √p, and a ∈ F_p^*." Problem 4 (Sarnak) is the case A = B = H, |H| ∼ √p, of Problem 3 (arbitrary A, B of size ∼ √p); Problem 4.3 is the two-set conjecture itself. The remark after Theorem 4.4 (Karacuba: |A| > p^{1/2+δ}, |B| > p^δ ⇒ |S(A,B)| ≪ p^{−0.05δ²}|A||B|) reads: "It is unknown if there is non-trivial bound on the character sum Σ_{x∈A,y∈B} χ(x+y) for |A| = |B| ∼ p^{1/2}, not even for the special case when A = B = H < F_p^*." So even |H| ∼ √p is open; p^ε is far beyond.
- **Gong, arXiv:1401.4618 (2014), Problem 1 (Bourgain) and Theorem 2:** "For any H < F_p^*, max_{a∈F_p^*} |Σ_{x∈H} χ(x+a)| < p^{1/2}. Thus, for any ε > 0 and H < F_p^* with |H| > p^{1/2+ε}, max_a |Σ_{x∈H}χ(x+a)| < p^{−ε}|H|." Nothing below √p.
- **Chang–Shparlinski, "Double character sums over subgroups and intervals" (arXiv, dated Feb. 23, 2014), §1.1:** the sums Σ_{λ∈G} χ(a+λ) "have been resisting all attempts to improve the classical bound (2) |Σ_{λ∈G} χ(a+λ)| ≤ √p", which "is instant from the Weil bound … (but can also be obtained via elementary arguments)". Their Theorem 1 (H > p^ε, T > p^{1/2−δ}) and Theorem 2 (T > p^ε, H > p^{1/4−δ}) concern *double* sums over an interval and a subgroup; their Theorem 4 bounds T_χ(a,G) = Σ_{λ,μ∈G} χ(a+λ+μ) nontrivially for T = |G| ≥ p^{13/33+ε} (using Shkredov's additive-energy bounds for subgroups). Theorem 4 is a two-set result with A = a+G, B = G, both structured, below √p; it is an average of T_G over a+G and gives no individual bound.
- **Volostnov–Shkredov, arXiv:1606.00358 (2016), §1:** "An affirmative answer was obtained just in the situation |A| > p^{1/2+δ}, |B| > p^δ … Even in the case |A| ∼ |B| ∼ p^{1/2} inequality (2) is unknown." Their Theorem 1 (Chang, Duke 2008, item (2) of its abstract: |A|, |B| > p^{4/9+ε}, |B+B| < K|B|) and their Theorem 2 (|A|,|B| > p^{12/31+δ}, |A+A| < K|A|, |A+B| < L|B|) are for **additive** doubling only; Volostnov arXiv:1712.09355 likewise ("sets having small additive doubling"). Neither paper treats multiplicative doubling.
- **Bourgain–Garaev–Konyagin–Shparlinski, arXiv:1110.0812 (hidden shifted power problem):** the only character-sum inputs over subgroups are Lemma 16, Σ_{x=1}^p χ(x^f + h) = O(f p^{1/2}) (Weil), and interval mixing (Lemmas 20–21, Corollary 33); no individual shifted-subgroup bound below √p.
- **Kowalski, arXiv:2401.04756 (2024)** is an exposition of Bourgain–Glibichuk–Konyagin, Theorem 1.1: for |H| ≥ p^γ, |Σ_{x∈H} e(ax/p)| ≪ |H|p^{−ν}. This is the **additive**-character theorem; no multiplicative analogue is mentioned. Bourgain's JAMS 2005 Mordell paper (local copy is a 5.6 KB HTML error page, so UNVERIFIED as a direct quote) is used in Chang–Shparlinski's Lemma 8 for additive characters only.
- **Kim–Yip–Yoo, arXiv:2309.09124 (Canad. J. Math., 2025) and arXiv:2602.20919 (2026):** the character-sum input they state is Vinogradov's bilinear bound |Σ_{a∈A,b∈B} ((ab+n)/p)| ≤ (p|A||B|)^{1/2}, their (1.1); their results are Stepanov-method structural theorems (AB+λ ⊂ S_d ∪ {0} ⇒ |A||B| ≤ |S_d| + …), not character-sum bounds.
- **Shparlinski, "Open problems on exponential and character sums" (2016), Problem 1:** "Obtain analogues of the results of J. Bourgain, A. A. Glibichuk and S. V. Konyagin for multiplicative character sums Σ_{x_1,…,x_k∈X} χ(x_1⋯x_k + a) and Σ_{x=1}^N χ(g^x + a) with very small values of N relative to p." (The consecutive-powers sum is the geometric-progression case of Task 3, see Prop. 4.3.)
- Not examined (not in the scratchpad, not fetched): Bourgain–Chang's multilinear Burgess note, Konyagin's 2002 survey, Shkredov's 2014 "medium size" paper (which Chang–Shparlinski cite for *additive-energy* bounds). A web search on 2026-09-05 returned only Gong 2014 and Kim–Yip–Yoo 2025 for the individual sum.

Numerical context (§6, scan of all 6753 subgroups of all 525 primes 101 ≤ p ≤ 4000): for every prime there are subgroups with M_H = |H| (a shift c+H monochromatic), the largest such having |H| up to 1.86·log p (p = 3121, |H| = 15, c = 1); max M_H/|H| over |H| > 2 log p is 0.944, over |H| > 6 log p is 0.52. So any uniform saving must come with a growth condition on |H| at least of order log p, consistent with the p^ε hypothesis.

## 2. Task 2 and its converse: exact equivalence for subgroups (PROVED)

**Proposition 2.1 (coset covariance).** For H ≤ F_p^*, c ∈ F_p, h ∈ H: T_H(ch) = χ(h) T_H(c). Hence |T_H| is constant on cosets of H, and for every θ the exceptional set E_H(θ) := {x ∈ F_p^* : |T_H(x)| > θ} is a union of cosets of H.

*Proof.* T_H(ch) = Σ_{h'∈H} χ(h(c + h'/h)) = χ(h) Σ_{h''∈H} χ(c + h''), since h' ↦ h'/h permutes H. ∎ (Verified for all subgroups of order ≤ 400 at p ∈ {101,197,229,257,1009,1093}.)

**Theorem 2.2 (equivalence).** Let H ≤ F_p^*.
(a) For every A ⊆ F_p and every union B = c_1H ∪ … ∪ c_kH of k distinct cosets,
    |S(A,B)| ≤ k( M_H |A| + |H| ).
In particular M_H ≤ p^{−δ}|H| implies |S(A,B)| ≤ (p^{−δ} + |A|^{−1}) |A||B| ≤ 2p^{−min(δ,ε)}|A||B| whenever |A| > p^ε.
(b) Conversely, if c ≠ 0 attains |T_H(c)| = M_H, then A_c := cH⁺ (so |A_c| ∈ {|H|, |H|/2}) satisfies S(A_c, H) = |H⁺| T_H(c), i.e. |S(A_c,H)| = M_H |A_c|.
Consequently, for all ε' > ε: [ |S(A,H)| ≤ p^{−δ}|A||H| for all |A| > p^ε and all subgroups |H| > p^ε ] implies [ M_H ≤ p^{−δ}|H| for all subgroups |H| > p^{ε'} ] (p large), and the converse holds with δ replaced by min(δ,ε) − log2/log p.

*Proof.* (a) For a ≠ 0, T_{cH}(a) = Σ_h χ(a + ch) = χ(c) Σ_h χ(a/c + h) = χ(c) T_H(a/c), so |T_B(a)| ≤ k M_H; for a = 0, |T_B(0)| ≤ k|H|. Sum over a ∈ A. (b) S(A_c,H) = Σ_{h∈H⁺} T_H(ch) = Σ_{h∈H⁺} χ(h) T_H(c) = |H⁺| T_H(c) by Prop. 2.1. ∎

Both directions were checked exactly (§6): e.g. p = 101, |H| = 20, c = 2, M_H = 6, |A_c| = 10, S(A_c,H) = −60. Thus **Task 2 is trivial, and — more importantly — the subgroup case of the two-set conjecture is not weaker than the shifted-subgroup problem: they are the same statement.** Note that T_H(c) = (Σ_y χ(y^d + c) − χ(c))/d with d = (p−1)/|H|, so the open problem is the multiplicative-character analogue of Bourgain's Mordell-type sparse-polynomial bound, for the binomial y^d + c with d = p^{1−ε}.

**Lemma 2.3 (top-m reformulation, any B).** For B ⊆ F_p and m ≥ 1 let τ_m(B) be the mean of the m largest values of |T_B(r)|, r ∈ F_p. Then (i) |S(A,B)| ≤ |A| τ_{|A|}(B) ≤ |A| τ_m(B) for all |A| ≥ m; (ii) if τ_m(B) > p^{−δ}|B| there is A with m p^{−δ}/2 ≤ |A| ≤ m and |S(A,B)| > p^{−δ}|A||B|/2. For a subgroup, τ_m(H) = M_H for all m ≤ |H|.
*Proof.* (i) is the rearrangement inequality. (ii) Let R be the set of the m largest positions and split R by the sign of T_B; the side A with the larger Σ|T_B| has S(A,B) = ±Σ_{A}|T_B| ≥ mτ_m/2 > m p^{−δ}|B|/2 ≥ |A|p^{−δ}|B|/2, and |A| ≥ mτ_m/(2|B|) ≥ m p^{−δ}/2 because |T_B| ≤ |B|. The last claim is Prop. 2.1. ∎ (Brute-force check of max_{|A|=m} S(A,B) = Σ of m largest T_B at p ∈ {13,17,19}.)

## 3. What elementary methods give, and exactly where they stop (PROVED)

**Proposition 3.1 (Weil moments).** For B ⊆ F_p, n = |B|, k ≥ 1: M_{2k}(B) ≤ (2k−1)!! n^k p + (2k−1) √p n^{2k}.
*Proof.* M_{2k} = Σ_{b_1..b_{2k}∈B} Σ_x χ(Π_i (x+b_i)). Tuples in which every value has even multiplicity number at most (2k−1)!! n^k (choose a perfect matching of the 2k positions and a value per pair) and contribute ≤ p each; all other tuples give a non-square polynomial with ≤ 2k distinct roots, hence ≤ (2k−1)√p by Weil. ∎ (Exact checks for k = 1,2,3 on subgroups and random sets at p ≤ 12289; the ratio LHS/RHS reached 0.98 at k = 1, so the k = 1 case is essentially the identity Σ T_B² = pn − n².)

**Corollary 3.2 (Chebyshev; Karatsuba's range).** #{x : |T_B(x)| > θ n} ≤ (2k−1)!! p θ^{−2k} n^{−k} + (2k−1) √p θ^{−2k}. For n ≥ p^ε, θ = p^{−δ} and kε ≥ 1/2 + 2kδ this is ≤ 2k p^{1/2+2kδ}. For a subgroup the exceptional set is a union of cosets (Prop. 2.1), so at most 2k p^{1/2+2kδ}/|H| cosets are exceptional; hence |S(A,H)| ≤ p^{−δ}|A||H| + |H|·|A ∩ E_H| is nontrivial as soon as |A| ≫ p^{1/2+2kδ+δ}. This is exactly Karacuba's range (Chang's survey, Theorem 4.4), and it is the same square-root barrier the brief describes. At p = 12289 the exceptional-coset counts at threshold |H|^{3/4} are 10/768, 8/384, 3/256 for |H| = 16, 32, 48 and 0 for |H| ≥ 64, versus Chebyshev (k = 2) bounds 129, 36, 17, 8.8.

**Proposition 3.3 (symmetric-set bound).** For B ⊆ F_p^*, n = |B|, integer t ≥ 0, let Sym_t(B) := {s ∈ F_p^* : |B ∩ sB| ≥ n − t}. Then for every r ≠ 0,
    (|T_B(r)| − 2t)_+² · |Sym_t(B)| ≤ p n − n²,   i.e.   M(B) ≤ 2t + ((pn − n²)/|Sym_t(B)|)^{1/2}.
*Proof.* T_B(rs) = Σ_b χ(s(r + b/s)) = χ(s) T_{s^{−1}B}(r) and |T_{s^{−1}B}(r) − T_B(r)| ≤ |s^{−1}B △ B| = 2(n − |B ∩ sB|) ≤ 2t for s ∈ Sym_t(B). Since s ↦ rs is injective, Σ_x T_B(x)² ≥ |Sym_t(B)|·(|T_B(r)| − 2t)_+², and the left side equals pn − n². ∎
For B = H: Sym_0(H) = H and M_H ≤ (p − |H|)^{1/2}, the elementary form of Gong's Theorem 2. Since Σ_s |B ∩ sB| = n², |Sym_t(B)| ≤ n²/(n−t), so the bound is nontrivial only when n ≳ p^{1/2}: it cannot leave the Weil range. (Exact checks at p ≤ 12289 for t ∈ {0, n/8, n/4}; and Sym_0(H) = H exactly.)

**Proposition 3.4 (multiplicative amplification and its obstruction).** For A ⊆ F_p, B, S ⊆ F_p^*, ν(y) := #{(a,s) ∈ A×S : a = ys}, E^×(A,S) := Σ_y ν(y)², and any k ≥ 1:
    |S|·|S(A,B)| ≤ |A| Σ_{s∈S} |B △ sB| + Σ_y ν(y)|T_B(y)|,
    (Σ_y ν(y)|T_B(y)|)^{2k} ≤ (|A||S|)^{2k−2} · E^×(A,S) · M_{2k}(B).
*Proof.* χ(s)T_B(a/s) = T_{sB}(a), so |T_B(a) − χ(s)T_B(a/s)| ≤ |B △ sB|; sum over (a,s) ∈ A×S and regroup by y = a/s. For the second line apply Hölder with exponents 2k/(2k−1) and 2k, then Σν^q ≤ (Σν)^{2−q}(Σν²)^{q−1} for q = 2k/(2k−1). ∎ (Exact integer checks, k = 1,2,3.)
If B is nearly S-invariant and E^×(A,S) ≤ L|A||S|, this gives |S(A,B)| ≤ |A|·avg_s|B△sB| + |A|(L/(|A||S|))^{1/2k} M_{2k}(B)^{1/2k}, nontrivial once |A||S| > L p^{1/2+ε} — a Burgess-type gain |A| → |A||S| — **but only under a multiplicative-energy hypothesis on A**. The hypothesis fails exactly on the extremal sets of Theorem 2.2: if A is a union of cosets of H and S ⊆ H then ν ≡ |S| on A and E^×(A,S) = |A||S|² (verified exactly), so the inequality collapses to Karatsuba's. Since the conjecture for (A,H) must in particular handle A = cH⁺ (Theorem 2.2(b)), no amplification of this kind can prove the subgroup case; this is the precise form of the brief's remark that the Hölder step must sum over all of F_p.

## 4. Task 3: A arbitrary, |BB| ≤ K|B| (verdict: OPEN; reductions PROVED)

Let 𝒞_K(ε) := {B ⊆ F_p^* : |B| > p^ε, |BB| ≤ K|B|}.

**Lemma 4.1.** |BB| = |B| iff B is a coset of a subgroup. (Standard: with b_0 ∈ B, |Bb_0| = |B| = |BB| forces BB = b_0B, hence bB = b_0B for all b ∈ B, so H = b_0^{−1}B satisfies hH = H for all h ∈ H and is a subgroup. Brute-force check over all subsets of size 2–4 of F_13^* and F_17^*: exactly the 13 resp. 12 cosets have |BB| = |B|.)

**Proposition 4.2 (the class contains two open problems).** (i) 𝒞_1(ε) is the set of cosets, and by Theorem 2.2 the two-set conjecture on {(A, B) : B ∈ 𝒞_1(ε)} is equivalent to the uniform shifted-subgroup bound (OPEN, §1). Hence any complete proof for 𝒞_K, K ≥ 1, would settle Bourgain's Problem 5, and none is offered here. (ii) The geometric progression P = {g, g², …, g^N} has |PP| = 2N − 1 ∈ 𝒞_2(ε) for N = p^ε; by Prop. 4.3 the instance (aP, P) is a Fejér-weighted consecutive-powers sum, which is Shparlinski's Problem 1 (§1).

**Proposition 4.3 (geometric progressions).** If g has multiplicative order > 2N and χ(g) = 1, then for every a ∈ F_p^*:
    S(aP, P) = Σ_{|k|<N} (N − |k|) χ(1 + a g^k).
*Proof.* χ(ag^i + g^j) = χ(g^j)χ(1 + ag^{i−j}) = χ(1 + ag^{i−j}); the number of (i,j) ∈ [1,N]² with i − j = k is N − |k|. ∎ (Exact checks at p ∈ {1009, 12289}, N ∈ {8,30,100}.) So the two-set conjecture for B = P and A = aP with |A| = |B| = p^ε implies |Σ_{|k|<N}(N−|k|)χ(1+ag^k)| ≤ p^{−δ}N², a bound for very short sums over consecutive powers that is not known.

**Proposition 4.4 (single-shift hypothesis suffices, trivially).** If M(B) ≤ p^{−δ}|B| then |S(A,B)| ≤ p^{−δ}|A||B| + |B| for every A. Conversely, for K close to 1 the propagation inequality |T_B(rs)| ≥ |T_B(r)| − |B △ sB| (Prop. 3.3's proof; verified for all r, s) shows that a single bad shift r spreads to r·Sym_t(B), so that Lemma 2.3(ii) produces a violating A of size ≥ |Sym_t(B)|·p^{−δ}/2; for general K ≥ 2 no such converse holds, since Sym_t(B) can be trivial (random B: |Sym_{|B|/4}(B)| = 1 at p = 12289).

**Proposition 4.5 (containers do not pass to dense subsets; REFUTES the Freiman route).** Multiplicative Freiman–Ruzsa in the cyclic group F_p^* places B ∈ 𝒞_K inside a coset progression of size ≤ f(K)|B|, but a bound for the container does not control its dense subsets, even when the subset is itself a subgroup: at p = 433 the subgroup H of order 24 has M_H = 8 = |H|/3, while its index-3 subgroup H' (order 8, density 1/3 in H) has M_{H'} = 8 = |H'|, attained at c' = 6, and A = 6·H'⁺ (|A| = 8) gives |S(A,H')| = 64 = |A||H'|. Further witnesses: (p, |H|, M_H, index, |H'|, M_{H'}, c') = (151,15,5,3,5,5,10), (211,15,5,3,5,5,17), (457,24,8,3,8,8,14), (463,21,7,3,7,7,38), (1777,24,8,2,12,12,108). Hence a reduction "𝒞_K → 𝒞_1 by containers" would need the two-set statement for all dense subsets of cosets, which is not implied by the subgroup hypothesis; the route is not closed.

**What is known for multiplicative doubling in the literature.** Nothing beyond the Weil range for K = 1 (§1). Chang's Duke 2008 Theorem (2) and Volostnov–Shkredov's Theorem 2 use additive doubling and sizes > p^{4/9}, > p^{12/31}; Chang–Shparlinski's Theorem 4 (A = a+G, B = G, |G| ≥ p^{13/33+ε}) is the strongest both-structured result below √p that was found.

**Numerical profile at p = 12289 = 3·2^12 + 1** (from the verifier; K = |BB|/|B|, Sym at t = |B|/4):

| B | |B| | K | M(B)/|B| | top-|B| mean/|B| | #{r: |T_B(r)| > 2√|B|} | |Sym_t(B)| |
|---|---|---|---|---|---|---|
| coset of |H| = 64 | 64 | 1.000 | 0.3125 | 0.3232 | 449 | 64 |
| subgroup |H| = 128 | 128 | 1.000 | 0.2188 | 0.2249 | 385 | 128 |
| subgroup |H| = 256 | 256 | 1.000 | 0.1406 | 0.1440 | 257 | 256 |
| union of 2 cosets, |H| = 64 | 128 | 1.500 | 0.2344 | 0.2404 | 257 | 64 |
| geometric progression, N = 128, χ(g) = 1 | 128 | 1.992 | 0.2969 | 0.2505 | 563 | 65 |
| geometric progression, N = 128, g primitive | 128 | 1.992 | 0.2812 | 0.2526 | 484 | 65 |
| subgroup 128 ∪ 8 random | 136 | 7.794 | 0.2647 | 0.2198 | 447 | 128 |
| random, 128 | 128 | 47.05 | 0.3906 | 0.2500 | 493 | 1 |
| random, 256 | 256 | 44.82 | 0.2344 | 0.1656 | 469 | 1 |

HEURISTIC reading: at this size the structured sets are not distinguishable from random ones by M(B)/|B| or by the top-|B| mean; the multiplicative structure shows up only in K and in Sym_t (65 = 2·32 + 1 for the progression, as predicted). "Subgroup ∪ few random points" is in 𝒞_K only for K ≈ number of added points, and then T_B = T_H + O(#points) trivially inherits the subgroup bound.

## 5. Task 4: bootstrapping by coset decomposition (PROVED obstruction, with witnesses)

**Identity.** Let H ≤ F_p^*, |H| = m, and for each coset C = c_C H put B_C := c_C^{−1}(B ∩ C) ⊆ H. Then for every A ⊆ F_p:
    S(A,B) = Σ_C χ(c_C) · S(c_C^{−1}A, B_C).
(Verified exactly: p = 12289, |H| = 64, |B| = 150, |A| = 200, both sides = 35, 106 cosets occupied.) The bootstrap would bound each S(c^{−1}A, B_C) by p^{−δ}|A||B_C| using the "structured" hypothesis for subsets of the subgroup H, then sum.

**Obstruction 5.1 (generic B is not reduced).** For random B ⊆ F_p^* of size N, the expected number of pairs in a common coset is C(N,2)(m−1)/(p−2) ≈ N²m/(2p). With N = p^ε, m = p^{ε/2} this is p^{5ε/2−1}/2 → 0 for ε < 2/5, so with probability → 1 every occupied coset is a singleton and the identity reads S(A,B) = Σ_{b∈B} χ(b) S(b^{−1}A, {1}), i.e. nothing was reduced. Witness: p = 1000003, m = 6, N = 251: 251 cosets met, maximum occupancy 1 (expected collisions 0.19); p = 100003, m = 21, N = 100: 99 cosets, max occupancy 2. For singletons no uniform hypothesis can hold: S(A', {1}) = |A'| for A' = {a : χ(a+1) = 1}, |A'| = (p−1)/2.

**Obstruction 5.2 (the input needed is the full conjecture on subsets of H).** When B is concentrated on few cosets the reduced sets B_C are arbitrary subsets of H of the occupation sizes. The required input "|S(A', B')| ≤ p^{−δ}|A'||B'| for every A' and every B' ⊆ H with |B'| ≥ |H|/K" is not implied by the shifted-subgroup bound for H (witnesses of Prop. 4.5, e.g. p = 433, |H| = 24, M_H = |H|/3, B' = H' of density 1/3 with |S(A,H')| = |A||H'|), and for |B'| ranging down to 1 it is the original conjecture with |B'| in place of |B|. So the decomposition changes the problem from (|A|, |B|) = (p^ε, p^ε) to (p^ε, |B ∩ C|) with |B ∩ C| ≤ p^{ε/2}, i.e. it only lowers the size of the second set. The bootstrap fails for a precise reason: **the two-set conjecture is not closed under passing from a set to its intersections with cosets, because the hypothesis is size-monotone in the wrong direction.**

**Obstruction 5.3 (smoothing is a different sum).** The Burgess-type alternative replaces S(A,B) by the smoothed sum Σ_{a∈A,b∈B,h∈H} χ(a+bh) = Σ_{b∈B} χ(b) S(b^{−1}A, H), which *is* ≤ p^{−δ}|A||B||H| + |A||B| under the shifted-subgroup hypothesis (Theorem 2.2(a)), but it is unrelated to S(A,B) unless B ≈ BH. Witness: p = 12289, A = {3}, B = {b : χ(3+b) = 1} (|B| = 6144), |H| = 64: S(A,B) = 6144 = |A||B| while the smoothed sum equals 6115 ≪ |A||B||H| = 393216.

## 6. Verifier record

`experiments/sigma_structured_2026_09_05.py` (numpy only, exact integer arithmetic; 16.3 s) performs 220,917 checks with 0 failures:
- §1: covariance T_H(ch) = χ(h)T_H(c), H-invariance of |T_H|, Σ T_H² = p|H| − |H|², M_H² ≤ p − |H|, the equivalence witnesses of Theorem 2.2(b), and the forward bound on unions of 1–3 cosets against random A (p ∈ {101,197,229,257,1009,1093}, all subgroups of order ≤ 400).
- §2: brute-force top-m/bottom-m identities (p = 13,17,19) and |BB| = |B| ⇔ coset (3281 subsets).
- §3: Weil moment bound and Chebyshev counts (k = 1,2,3), propagation lemma for all r, s (≈ 2·10^5 checks), symmetric-set bound, amplification inequalities (both lines, k = 1,2,3), E^×(A,H) = |A||H|² for H-invariant A.
- §4: geometric-progression identity (24 instances).
- §5: all 6753 subgroups of all 525 primes in [101, 4000]: M_H < √p, M_H² ≤ p − |H|, monochromatic-shift maxima, band maxima of M_H/|H|, index-2/3 dense-subset witnesses.
- §6: coset occupancy witnesses, decomposition identity, dense-subset witnesses re-verified, smoothing witness.
- §7: the profile table and the exceptional-coset table at p = 12289.
All witnesses are stored in `results/sigma_structured_2026_09_05.json`.

## 7. Remaining obligations

1. (OPEN, the direction's core) Prove max_{c≠0}|Σ_{h∈H} χ(h+c)| ≤ p^{−δ}|H| for |H| > p^ε — equivalently (Theorem 2.2) the two-set conjecture on (A, H). By §3 no combination of Weil moments, second-moment symmetry, or ratio-amplification can do it, because the extremal sets A = cH⁺ are H-invariant; a genuinely new input is needed (an analogue for χ(y^d + c), d = p^{1−ε}, of Bourgain's sum-product treatment of e(f(y)/p)).
2. (OPEN) The consecutive-powers sum Σ_{k}(N−|k|)χ(1 + ag^k) for N = p^ε (Prop. 4.3), needed already for 𝒞_2.
3. (OPEN) Even conditionally on 1, the passage from cosets to dense subsets of cosets (Prop. 4.5) and hence to 𝒞_K for K ≥ 2.
4. (Not examined) Bourgain–Chang's multilinear Burgess note and Konyagin's 2002 survey; Bourgain's JAMS 2005 text could not be fetched locally.

## Sources

- Chang, M.-C., *Character sums in finite fields*, Contemp. Math. 518 (2010) — `sources/sigma-chang-character-sums-survey-DubProc.txt`, Problems 3–5, Theorem 4.4 and its remark.
- Chang, M.-C., *On a question of Davenport and Lewis and new character sum bounds in finite fields*, Duke Math. J. 145 (2008) — scratchpad `chang-DL4.txt`, abstract item (2).
- Chang, M.-C., Shparlinski, I. E., *Double character sums over subgroups and intervals* (2014) — scratchpad `chang-151-intergroup.txt`, §1.1, Theorems 1–4.
- Gong, K., *An elementary approach to character sums over multiplicative subgroups*, arXiv:1401.4618 — scratchpad, Problem 1, Theorems 2 and 4 ([arXiv](https://arxiv.org/abs/1401.4618)).
- Volostnov, A. S., Shkredov, I. D., *Sums of multiplicative characters with additive convolutions*, arXiv:1606.00358 — scratchpad, §1, Theorems 1–2.
- Bourgain, J., Garaev, M. Z., Konyagin, S. V., Shparlinski, I. E., *On the hidden shifted power problem*, arXiv:1110.0812 — scratchpad, Lemmas 3, 16, 20, 21.
- Kowalski, E., *Exponential sums over small subgroups, revisited*, arXiv:2401.04756 — scratchpad, Theorem 1.1, Remark 1.2.
- Kim, S., Yip, C. H., Yoo, S., arXiv:2309.09124 (Canad. J. Math. 2025) and arXiv:2602.20919 (2026) — scratchpad, (1.1)–(1.2), Theorem 1.1 ([arXiv](https://arxiv.org/abs/2309.09124)).
- Shparlinski, I. E., *Open problems on exponential and character sums* (2016) — scratchpad, Problem 1.
- Shkredov, I. D., Vyugin, I. V., arXiv:1102.1172; Macourt–Shkredov–Shparlinski, arXiv:1701.06192 — scratchpad; examined, contain no shifted-subgroup character-sum bound.
