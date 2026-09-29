# Complete bicliques of the Paley relation: literature, exact data, and the sum-product route

**Status: PROVED — (a) for every complete biclique `A + B ⊆ Q ∪ {0}` the Hanson–Petridis bound `|A||B| ≤ (p−1)/2 + |B ∩ (−A)|` (CITED, Thm 1.2 of arXiv:1905.09134) refines to `|A||B| ≤ (p−1)/2 + min(|A|, |B|, ω(G_p))`, because `A ∩ (−B)` is a clique of the Paley graph; hence the balanced biclique number `b(p)` and the sum-clique number `s(p)` obey the same bound `(1+√(2p−1))/2` as `ω(G_p)`; (b) `M_2(p) = max_{|A|=2}|B(A)| = (p+3)/4` exactly, so `P(p) := max |A||B| ≥ (p+3)/2`, and this is a case of equality in the Hanson–Petridis inequality for every prime; (c) a Weil count bound `|B(A)| ≤ k + 2^{−k}[p − C(k,2) + c_k√p]` for `|A| = k`, which alone excludes `|A||B| > (p+3)/2` when `min(|A|,|B|) = 3` for `p ≥ 116` and `= 4` for `p ≥ 213`; (d) a reduction (via (a)) that makes `P(p)` finitely computable with candidate sets of size `≈ p/8`, and with it the exact values `P(p) = (p+3)/2` for all primes `43 ≤ p ≤ 1000`, `p ≡ 1 (mod 4)`, while `P(13) = 9`, `P(37) = 21`, `P(41) = 25` are clique-type; (e) `b(p) ≥ max(ω(G_p), s(p))`; (f) two no-go lemmas for the sum-product route on the complete case: the ratio-set constraint `(A+B)/(A+B) ⊆ Q` contributes only `T ≥ 2n⁴/(p−1)` to the weighted multiplicative energy `T` of the sumset, which is implied by the trivial diagonal count `T ≥ E⁺(A,B)² ≥ n²` whenever `2n² ≤ p−1`, and no upper bound on `T` valid for all pairs can push the chain below `|A||B| ≤ (2/9+o(1))p`. REFUTED (exact witnesses) — `b(p) = ω(G_p)` (`p = 89`), `b(p) ≤ ω(G_p)+1` (`p = 233`: `ω = 7`, `b = s = 9`), `P(p) = (p+3)/2` for all `p` (`p = 13, 37, 41`), and `gp_k(p) = M_k(p)` for `k ≥ 4` (`p = 997`, `k = 4`). OPEN — every bound `|A||B| ≤ p^{1−c}` (or `ω(G_p) ≤ p^{1/2−c}`), the identity `P(p) = (p+3)/2` for all `p ≥ 43`, exact `b(p)` beyond the computed range, and whether `b(p) − ω(G_p)` is bounded. The literature pin (§1) finds no improvement of the Hanson–Petridis clique bound for primes through August 2026; Sárközy's conjecture itself was settled by Kalmynin (2025).**

Worker slug `biclique`. Verifier: `experiments/sigma_biclique_2026_09_05.py` (standard library + sympy), results `results/sigma_biclique_2026_09_05.json`; search table produced by the C++ helper `experiments/sigma_biclique_2026_09_05.cpp` and stored in `results/sigma_biclique_2026_09_05_search.json`. Every number quoted below is in one of those two files.

---

## 0. Definitions, symmetries, and the quantities computed

Throughout `p ≡ 1 (mod 4)` is prime, `Q ⊆ F_p^*` is the set of nonzero squares, `χ` the Legendre symbol, `G_p` the Paley graph (`x ~ y` iff `χ(x−y) = 1`; symmetric since `χ(−1) = 1`).

- A **complete biclique** is a pair `A, B ⊆ F_p` with `A + B ⊆ Q ∪ {0}`. For `A ⊆ F_p` put
  `N'(a) := (Q ∪ {0}) − a = {b : a + b ∈ Q ∪ {0}}` and `B(A) := ∩_{a∈A} N'(a)`; `(A, B(A))` is the largest complete biclique with left side `A`.
- Writing `B' = −B`, the condition reads `A − B' ⊆ Q ∪ {0}`: every `a ∈ A` is adjacent or equal to every `b' ∈ B'` in `G_p`. A clique `C` gives the biclique `(C, −C)`; a set with `A + A ⊆ Q ∪ {0}` (a **sum-clique**) gives `(A, A)`.
- Symmetries preserving the relation: `(A,B) ↦ (uA + t, uB − t)` for `u ∈ Q`, `t ∈ F_p`, and `(A,B) ↦ (B,A)`. Hence any `A` with `|A| ≥ 2` may be normalised to contain `0` and either `1` or the least non-residue `n_0`; a clique may be normalised to contain `0, 1`.
- Quantities: `ω(p) = ω(G_p)`; `M_k(p) := max_{|A|=k} |B(A)|`; `b(p) := max min(|A|, |B(A)|)` (balanced biclique number); `s(p) := max{|A| : A + A ⊆ Q ∪ {0}}` (sum-clique number); `P(p) := max{|A||B| : A + B ⊆ Q ∪ {0}, |A|, |B| ≥ 2}`; `gp_k(p) := max |B(c·{r^i : 0 ≤ i < k})|` over `r ≠ 0, 1` with `r^0, …, r^{k−1}` distinct and `c ∈ {1, n_0}` (geometric progressions, the construction of `research/sigma-crux-2026-09-05.md` Thm 6.2).

Trivially `s(p) ≤ b(p)`, `ω(p) ≤ b(p)`, `b(p)² ≤ P(p)`, and `k·M_k(p) ≤ P(p)` whenever `M_k(p) ≥ 2`.

---

## 1. Literature pin (CITED; sources fetched as indicated)

**1.1 Hanson–Petridis.** B. Hanson, G. Petridis, *Refined estimates concerning sumsets contained in the roots of unity*, Proc. London Math. Soc. (3) 122 (2021) 353–358, DOI 10.1112/plms.12322; arXiv:1905.09134v3 (local copy `sources/sigma-hanson-petridis-1905.09134.txt`). Verbatim statements:

> **Theorem 1.2.** Let `p` be a prime and suppose `A, B ⊆ F_p` satisfy `A + B ⊆ Z_d ∪ {0}` for some `d` properly dividing `p − 1`. Then `|A||B| ≤ d + |B ∩ (−A)|`.

> **Corollary 1.3.** Let `p` be a prime and suppose `A, B ⊆ F_p` satisfy `A + B = Z_d` for some `d` properly dividing `p − 1`. Then `|A||B| = d` and all sums `a + b` are distinct. In particular, if `d` is prime and neither `A` nor `B` is a singleton, no such decomposition is possible.

> **Corollary 1.4 (Shakan).** As `x` tends to infinity, the number of primes `p ≤ x` such that there exist non-singleton sets `A, B ⊆ F_p` with `A + B = Z_{(p−1)/2}` is `o(π(x))`.

> **Corollary 1.5.** Let `p` be a prime, `d` properly dividing `p − 1` and suppose `A ⊆ F_p` is such that `A − A ⊆ Z_d ∪ {0}`. Then `|A|(|A| − 1) ≤ d`. In particular, for `p ≡ 1 (mod 4)`, we have `ω(G_p) ≤ (√(2p−1) + 1)/2`.

Here `Z_d` is the set of `d`-th roots of unity; with `d = (p−1)/2`, `Z_d = Q`. The paper also states (p. 4 of the arXiv version) that "The method developed to prove Theorem 1.2 only works in prime fields", and quotes [MP] (Maistrelli–Penman, Discrete Math. 306 (2006)) for `ω(G_p) ≤ √p − 4` "for certain primes". The constant in the clique bound is therefore exactly `ω(G_p) ≤ (1 + √(2p−1))/2 = √(p/2) + 1/2 + O(p^{−1/2})`, and for bicliques `|A||B| ≤ (p−1)/2 + |B ∩ (−A)|`.

**1.2 Sárközy's conjecture and its history.** A. Sárközy, *On additive decompositions of the set of quadratic residues modulo p*, Acta Arith. 155 (2012) 41–51, DOI 10.4064/aa155-1-4 (primary text not fetched; statements below are taken from the introductions of Chen–Xi and Kalmynin, both fetched): conjecture that `R_p := Q ≠ A + B` for `|A|, |B| ≥ 2` and large `p`; proof of the ternary version `Q ≠ A + B + C` for large `p`; and, if `Q = A + B`, `√p/(3 log p) ≤ |A|, |B| ≤ √p log p` (Chen–Xi, eq. (1.1)). I. Shparlinski, *Additive decompositions of subgroups of finite fields*, SIAM J. Discrete Math. 27 (2013) 1870–1879, and I. D. Shkredov, *Sumsets in quadratic residues*, Acta Arith. 164 (2014) 221–243 (arXiv:1305.4093, abstract fetched) removed the logarithms; Shkredov proved `A + A ≠ Q` for every `A` with `|A| ≥ 2` and obtained constants `1/6 − o(1)` and `3 + o(1)` (as quoted by Chen–Xi and by Hanson–Petridis [Shk1, Cor. 2.6]). I. D. Shkredov, *Any small multiplicative subgroup is not a sumset*, Finite Fields Appl. 63 (2020), arXiv:1702.01197 (abstract fetched): for every `ε > 0` and any multiplicative subgroup `Γ ⊆ F_p` with `1 ≪ |Γ| ≤ p^{2/3−ε}` there are no `B, C` with `|B|, |C| > 1` and `Γ = B + C`. Y.-G. Chen, X.-H. Yan, *A conjecture of Sárközy on quadratic residues*, J. Number Theory 229 (2021) 100–124 (as quoted in Chen–Xi §1, eq. (1.2)): no additive 3-decomposition of `Q` for any prime, and `((7−√17)/16)√p + 1 ≤ |A|, |B| ≤ ((7+√17)/4)√p − 6.63` for a 2-decomposition. Y.-G. Chen, P. Xi, *A conjecture of Sárközy on quadratic residues, II*, arXiv:2202.02780 (local copy `sources/sigma-chen-xi-2202.02780.txt`), Theorem 1.1: if `A + B = Q` with `|A|, |B| ≥ 2` then `√p/4 + 1/8 ≤ |A|, |B| < 2√p − 1`, and at least `(log 2)^{−1}√p − 1.6` elements of `A + B` have a unique representation.

**The conjecture is now a theorem.** A. Kalmynin, *On additive irreducibility of multiplicative subgroups*, arXiv:2504.10202v2 (28 May 2025; local copy `sources/sigma-kalmynin-2504.10202v2.txt`):

> **Theorem 3 (Sárközy's conjecture).** For all primes `p`, the set `R_p = μ_{(p−1)/2}` of quadratic residues modulo `p` is additively irreducible, i.e. `R_p ≠ A + B` with `|A|, |B| > 1`.

> **Theorem 2 (α = β theorem).** If `μ_d = A + B` with `|A|, |B| > 1`, then `|A| = |B| = √d`.

> **Theorem 1 (Lev–Sonn conjecture).** Suppose that `p` is a prime number and `d | p − 1`, `1 < d < p − 1`. If `μ_d` is the set of all `d`-th roots of unity in `F_p` and `A − A = μ_d ∪ {0}` for some `A ⊂ F_p` then `d = 2` or `6`.

together with Corollary 1 (no proper subgroup is `A + B + C` with all `|·| > 1`). M. Rudnev, F. Tyrrell, *Multiplicative subgroups of prime fields are not sumsets*, arXiv:2607.24270 (v2, 12 Aug 2026; abstract fetched): if `H ≤ F_p^*` is a proper multiplicative subgroup and `H = A + B`, then one summand is a singleton or `|A| = |B| = 2` and `|H| = 4`. C. H. Yip, S. Yoo, *Additive decompositions of multiplicative subgroups in prime fields: a self-contained approach*, arXiv:2608.02568 (3 Aug 2026; HTML fetched), Theorem 1.1: same classification, with the exceptional example `H = {1, −1, i, −i} = {0, −1−i} + {1, i}`.

**1.3 Clique bounds.** C. Bachoc, M. Matolcsi, I. Z. Ruzsa, *Squares and difference sets in finite fields*, Integers 13 (2013) A77 (arXiv:1305.0577, abstract fetched): `|B| ≤ √p − 1` for `B − B ⊆ Q ∪ {0}` "for approximately three quarters of the primes `p = 4k+1`". D. Di Benedetto, J. Solymosi, E. P. White, *On the directions determined by a Cartesian product in an affine Galois plane*, Combinatorica 41 (2021) 755–763 (arXiv:2001.06994, abstract fetched): `A × B ⊂ AG(2,p)` with `|A|, |B| ≥ 2`, `|A||B| < p` determines at least `|A||B| − min{|A|,|B|} + 2` directions, which yields the same clique bound as Hanson–Petridis. C. H. Yip, *On the clique number of Paley graphs of prime power order*, Finite Fields Appl. 77 (2022) 101930 (arXiv:2004.01175, abstract fetched): for `q = p^{2s+1}`, `p ≡ 1 (mod 4)`, `ω(G_q) ≤ min(p^s⌈√(p/2)⌉, √(q/2) + (p^s+1)/4 + (√(2p)/32)p^{s−1})`. C. H. Yip, *Exact values and improved bounds on the clique number of cyclotomic graphs*, Des. Codes Cryptogr. 93 (2025) 5131–5142 (arXiv:2304.13213, abstract fetched): `ω(Cay(F_q^+, S)) ≤ √|S/S| + √(q/p)`, which for `q = p`, `S = Q` is `√((p−1)/2) + 1`, i.e. Hanson–Petridis strength. D. Kunisky, X. Yu, *A degree 4 sum-of-squares lower bound for the clique number of the Paley graph*, CCC 2023 (arXiv:2211.02713, abstract fetched): the degree-4 SOS relaxation is `≥ Ω(p^{1/3})`. Randomstrasse101, *On the clique number of the Paley graph (problems 25–29)* (web page fetched) and *Open Problems of 2025* (arXiv:2603.29571, March 2026, HTML fetched): "Recently, Hanson and Petridis [HP (21)] showed that `ω(G_p) ≤ (1+o(1))√(p/2)`" is recorded as the best known bound; Conjecture 25 there is `ω(G_p) = O(polylog p)`. Two further web searches performed for this note (September 2026) surfaced improvements only for generalized Paley graphs (a search summary mentions an August 2026 preprint with `0.769√q` for cubic Paley graphs; not fetched, UNVERIFIED CITATION), none for the prime Paley graph.

**Best known bounds, as pinned.** For primes `p ≡ 1 (mod 4)`: `ω(G_p) ≤ (1+√(2p−1))/2` (HP Cor. 1.5; also DSW 2021, Yip 2025), `ω(G_p) ≤ √p − 1` on a set of primes of density about 3/4 (BMR 2013, weaker), and `ω(G_p) ≫ log p · log log log p` infinitely often (Graham–Ringrose, as in the brief). For complete bicliques: `|A||B| ≤ (p−1)/2 + |B ∩ (−A)|` (HP Thm 1.2), sharpened in Proposition 2.1 below to `(p−1)/2 + min(|A|,|B|,ω(G_p))`; the lower bound `(p+3)/2` of Proposition 2.3 shows the constant `1/2` cannot be improved. Nothing below `√(p/2)` for cliques or below `p/2` for products is known, and this note proves nothing of that kind either.

**1.4 Published clique numbers.** A. E. Brouwer's page *Paley graphs* (https://aeb.win.tue.nl/graphs/Paley.html, fetched) lists `ω(G_q)` for `q ≤ 197` (Shearer; `q = 125, 173, 197` due to Exoo) and states that Shearer computed the independence numbers of Paley graphs of prime order below 7000 and Exoo extended the table beyond 16000. The 21 prime entries there (`5:2, 13:3, 17:3, 29:4, 37:4, 41:5, 53:5, 61:5, 73:5, 89:5, 97:6, 101:5, 109:6, 113:7, 137:7, 149:7, 157:7, 173:8, 181:7, 193:7, 197:8`) agree with the values computed here (checked by the verifier). Brouwer also notes that the Hanson–Petridis clique bound is attained for `q = 5, 13, 41` (indeed `ω(ω−1) = (p−1)/2` there); the primes at which the product `P(p)` exceeds `(p+3)/2` are `13, 37, 41` (Theorem 2.7) — at `p = 37` through a 3-clique with a 7-element `B`-side, not through a maximum clique.

---

## 2. Elementary structure of complete bicliques (PROVED)

**Proposition 2.1 (Hanson–Petridis with the clique refinement).** Let `A + B ⊆ Q ∪ {0}`, `r := |B ∩ (−A)|` and `C := A ∩ (−B)`. Then `|C| = r`, `C` is a clique of `G_p`, and

    |A||B| ≤ (p−1)/2 + r ≤ (p−1)/2 + min(|A|, |B|, ω(G_p)).

*Proof.* `b ↦ −b` is a bijection `B ∩ (−A) → A ∩ (−B)`, so `|C| = r`. For distinct `c, c' ∈ C` we have `c ∈ A` and `−c' ∈ B`, hence `c − c' ∈ Q ∪ {0}`, and `c − c' ≠ 0`; so `χ(c − c') = 1` and `C` is a clique, `r ≤ ω(G_p)`. The first inequality is Theorem 1.2 of Hanson–Petridis with `d = (p−1)/2`, which properly divides `p − 1` for `p ≥ 3`, and `Z_d = Q`. ∎

**Corollary 2.2.** `b(p) ≤ (1+√(2p−1))/2`, `s(p) ≤ (1+√(2p−1))/2`, and `P(p) ≤ (p−1)/2 + ω(G_p) ≤ (p−1)/2 + (1+√(2p−1))/2`. (For a `k × k` biclique, `k² ≤ (p−1)/2 + k`.)

**Proposition 2.3 (the two-element side).** For `a ≠ 0`,
`#{b ∈ F_p ∖ {0, −a} : χ(b) = χ(b+a) = 1} = (p − 3 − 2χ(a))/4`. Consequently `|B({0,a})| = (p+3)/4` if `a ∈ Q` and `(p−1)/4` if `a ∉ Q`; so `M_2(p) = (p+3)/4`, `P(p) ≥ (p+3)/2`, and the pair `({0,1}, B({0,1}))` attains equality in Theorem 1.2 (`r = 2`).

*Proof.* `Σ_b (1+χ(b))(1+χ(b+a)) = p + 0 + 0 + Σ_b χ(b)χ(b+a) = p − 1`. The terms `b = 0` and `b = −a` contribute `1 + χ(a)` each (`χ(−a) = χ(a)`), the remaining terms are `4` on the set counted and `0` otherwise. For `a ∈ Q` the set `B({0,a})` is the counted set together with `b = 0` (`0 + 0 = 0`, `0 + a ∈ Q`) and `b = −a` (`−a ∈ Q`, `−a + a = 0`), giving `(p−5)/4 + 2`. For `a ∉ Q` neither `0` nor `−a` qualifies and the count is `(p−1)/4`. Equality in HP: `|A||B| = (p+3)/2 = (p−1)/2 + 2` and `B ∩ (−A) = {0, −1}`. ∎

**Proposition 2.4 (Weil count for `|B(A)|`).** For `|A| = k ≥ 2` put `c_k := Σ_{j=3}^{k} C(k,j)(j−1) = (k−2)2^{k−1} + 1 − C(k,2)`. Then

    #{b ∉ −A : χ(b+a) = 1 for all a ∈ A} ≤ 2^{−k} [ p − C(k,2) + c_k √p ],   |B(A)| ≤ k + 2^{−k}[ p − C(k,2) + c_k √p ].

*Proof.* Expand `Σ_b Π_{a∈A}(1 + χ(b+a)) = p + Σ_{∅≠I⊆A} Σ_b χ(Π_{a∈I}(b+a))`. Terms with `|I| = 1` vanish; for `|I| = 2`, `Σ_b χ((b+a)(b+a')) = −1`; for `|I| ≥ 3` the polynomial `Π_{a∈I}(x+a)` has `|I|` distinct roots, is not a constant times a square, and Weil gives `|Σ_b| ≤ (|I|−1)√p`. Every summand on the left is `≥ 0`, and each `b ∉ −A` with all `χ(b+a) = 1` contributes `2^k`. Finally `B(A) ⊆ {b ∉ −A : all χ(b+a) = 1} ∪ (−A)`. ∎

**Corollary 2.5 (small min-side).** If `min(|A|,|B|) = m ≥ 3` then `|A||B| ≤ m(m + 2^{−m}[p − C(m,2) + c_m√p])`, which is `< (p+3)/2` for `m = 3` and `p ≥ 116` (`p − 6√p − 51 > 0`), for `m = 4` and `p ≥ 213` (`p − 11√p − 52 > 0`), and in general for `m ≤ m*(p)` with `m*(p) → ∞` (the verifier records, for each `p ≤ 3000`, the set of `m ≤ 40` excluded this way; along the primes `p ≡ 1 (mod 4)` it is an initial segment `3 ≤ m ≤ m*(p)` with `m*(p) = 3, 4, 5, 6, 7, 8` from `p = 137, 229, 433, 853, 1597, 2713` on, respectively). So a biclique with `|A||B| > (p+3)/2` needs `min(|A|,|B|) ≥ 4` once `p ≥ 116`, `≥ 5` once `p ≥ 213`, etc.

**Proposition 2.6 (the reduction).** Let `A + B ⊆ Q ∪ {0}` with `|A||B| ≥ (p+5)/2` and, without loss of generality, `|A| ≤ |B|`. Then (i) `r = |B ∩ (−A)| ≥ 3`; (ii) after applying a symmetry `(A,B) ↦ (uA+t, uB−t)`, `u ∈ Q`, one has `{0, 1, c} ⊆ A` and `{0, −1, −c} ⊆ B` for some `c` with `χ(c) = χ(c−1) = 1`; (iii) `A ⊆ N'(0) ∩ N'(−1) ∩ N'(−c)`, `B ⊆ N'(0) ∩ N'(1) ∩ N'(c)`, and `|B| ≥ ⌈√((p+5)/2)⌉`, `|A| ≤ |B|`.

*Proof.* (i) is Proposition 2.1. For (ii) pick three elements `c_1, c_2, c_3` of the clique `C = A ∩ (−B)`; the map `x ↦ (x − c_1)/(c_2 − c_1)` on `A` and `x ↦ (x + c_1)/(c_2 − c_1)` on `B` preserves the relation because `c_2 − c_1 ∈ Q`; it sends `c_1, c_2, c_3` to `0, 1, c := (c_3−c_1)/(c_2−c_1)` and `−c_1, −c_2, −c_3` to `0, −1, −c`; `χ(c) = χ(c_3−c_1)χ(c_2−c_1) = 1` and `χ(c−1) = χ(c_3−c_2)χ(c_2−c_1) = 1`. (iii) `A ⊆ B(B) ⊆ B({0,−1,−c})` and symmetrically; `|B| ≥ √(|A||B|)`. ∎

Since the two candidate sets in (iii) have about `p/8` elements each, `P(p)` becomes a small finite computation: enumerate `c` (there are `(p−5)/4` of them, `c ~ 1−c`), then search `A` inside `N'(0) ∩ N'(−1) ∩ N'(−c)` with the pruning `|B(A)| ≥ ⌈√((p+5)/2)⌉`, `|A| ≤ |B(A)|`, and the product bound `max_r (|A|+r)·v_r` (`v_r` the `r`-th largest `|B(A) ∩ N'(a)|` over remaining candidates `a`). The search is exact: every biclique with `|A||B| ≥ (p+5)/2` is visited up to symmetry.

**Theorem 2.7 (computed, PROVED for the stated range).** For every prime `p ≡ 1 (mod 4)` with `43 ≤ p ≤ 1000`, `P(p) = (p+3)/2`, attained by `({0,1}, B({0,1}))`; and `P(5) = 4 = (5+3)/2`, `P(13) = 9` (`A = {0,1,4}`, `B = −A`, a 3-clique, `r = 3`), `P(37) = 21` (`A = {0,1,11}`, `|B(A)| = 7`, `r = 3`), `P(41) = 25` (`A = {0,1,2,10,33}`, `B = −A`, a 5-clique, `r = 5`). In each exceptional case `P(p) = (p−1)/2 + r`, i.e. Theorem 1.2 holds with equality. (C++ search for all `p ≤ 1000`; Python re-run of the same reduced search for `p ≤ 500`; brute-force search without the reduction for `p ≤ 100`, all agreeing.)

**Proposition 2.8 (bicliques from cliques and sum-cliques).** `b(p) ≥ max(ω(G_p), s(p))`. (If `C` is a clique, `C + (−C) ⊆ Q ∪ {0}`; if `A + A ⊆ Q ∪ {0}`, take `B = A`.) The extension `B(C)` of a maximum clique satisfies `|B(C)| ≥ ω`; in the table `|B(C)| = ω(G_p)` for every maximum clique found, i.e. the diagonal biclique `(C, −C)` was never extendable on the `B` side.

---

## 3. Exact data (`p ≡ 1 (mod 4)`, `p ≤ 1000`; `ω`, `s` to 3000)

All values are exact unless marked `*` (time-limited search; the value is then a lower bound with an explicit witness). Columns: `ω = ω(G_p)`; `s = s(p)`; `M_k = M_k(p)` for `k = 2..8`; `b = b(p)`; `b_gp = max{k : gp_k(p) ≥ k}` (balanced geometric-progression bicliques); `√(p/2)`; `L = log p · log log log p` (Graham–Ringrose scale); `log₂p`.

<!-- TABLE -->

**What the table shows (HEURISTIC, exact numbers).**

1. *Products.* `P(p) = (p+3)/2` from `p = 43` on: the product maximum is always attained by the most unbalanced biclique, and `3M_3, 4M_4, …` are all below `(p+3)/2` (`3M_3(997) = 399 < 500`, `4M_4(997) = 312`). The polynomial-method constant `1/2` is therefore optimal for the product problem *because of* the `|A| = 2` family; the interesting quantity is the balanced number.
2. *Profile.* `M_k(p)/(p/2^k)` grows with `k`: at `p ≈ 1000` the ratios are `≈ 1.07, 1.25, 1.65, 2.4` for `k = 3, 4, 5, 6` and `≈ 4–6, 7–10` for `k = 7, 8` at `p ≈ 400–650`. This is the expected "maximum of many fluctuations" effect: `|B(A)|` has mean `≈ p/2^k` and the maximum over `≈ p^{k−2}` normalised `A` adds `≈ √(2(p/2^k)(k−2)log p)`.
3. *Balanced bicliques versus cliques.* `b(p) − ω(G_p) ∈ {0, 1, 2}` for every `p` with exact `b(p)` (all `p ≤ 509`: `25` primes with difference `0`, `15` with `1`, `5` with `2`); the value `2` occurs at `p = 233, 257, 281` (`ω = 7`, `b = 9`), `449` and `509` (`ω = 9`, `b = 11`). The greedy lower bounds for `p ≤ 1000` also stay within `2` of `ω`. So the *balanced* product `b(p)²` exceeds `ω(G_p)²` by `2ω+1` or `4ω+4` at those primes — more than a constant, but by a quantity of order `ω = O(log p)`, and `b/ω ≤ 9/7` in the data.
4. *Sum-cliques.* `s(p)` (sets with `A + A ⊆ Q ∪ {0}`) is frequently larger than `ω(G_p)`: over the `211` primes `p ≤ 3000`, `s − ω` takes the values `−3` (once, `p = 1741`), `−2` (once, `p = 461`), `−1` (`23` primes), `0` (`90`), `1` (`53`), `2` (`28`; e.g. `p = 233, 281, 449, 1009`) and `3` (`5` primes: `1609, 2341, 2377, 2593, 2801`); `s > ω` at `25` of the `80` primes `p ≤ 1000`. Where `b(p)` is exact and `b > ω` (`20` primes `≤ 509`), `s(p) = b(p)` in `12` cases, i.e. the extremal balanced biclique is then of the form `(A, A)`; in the other `8` cases `b = s + 1` (e.g. `p = 257`: `ω = 7`, `s = 8`, `b = 9`). Sum-cliques satisfy the same Hanson–Petridis-type bound (Corollary 2.2), and `s(p) ≤ b(p)` always.
5. *Growth.* `ω(G_p)/log₂p` lies in `[0.86, 1.34]` for `200 ≤ p ≤ 1000` and `b(p)/log₂p` (greedy lower bounds where `b` is not exact) in `[0.99, 1.34]`; both are far below `√(p/2)` (`22.3` at `p = 997`) and above `L = log p · log log log p` (`≈ 4.5` at `p = 997`). The random-graph benchmark `2log₂p − 2log₂log₂p + O(1)` is about `15` at `p = 1000`; Paley cliques and bicliques are smaller. The data are consistent with `b(p) = ω(G_p) + O(1)` and with `ω(G_p) ≍ log p` up to a slowly varying factor; they say nothing about the asymptotic exponent.
6. *Geometric progressions.* `gp_3(p) = M_3(p)` for most `p` (explained: up to the symmetries, `{c, cr, cr²}` is `{0, 1, r+1}` scaled, so every 3-set except the arithmetic progression `{0,1,2}` is a geometric progression), while `gp_k < M_k` for `k ≥ 4` at most `p` (e.g. `p = 997`: `76 < 78`, `47 < 51`, `31 < 37`). The balanced GP number `b_gp(p)` is `0`–`3` below the exact `b(p)`; GPs are good but not extremal.

---

## 4. Refuted hypotheses (REFUTED, with witnesses in the results file)

- **H1: `b(p) = ω(G_p)`.** Witness `p = 89`: `A = {0,1,3,10,41,70}`, `B(A) = {1,8,39,68,87,88}` (`|B| = 6`, `ω(G_89) = 5`); exact `b(89) = 6`.
- **H2: `b(p) ≤ ω(G_p) + 1`.** Witness `p = 233`: `A = B = {0,1,7,8,25,30,170,195,207}` (`A + A ⊆ Q ∪ {0}`, `ω(G_233) = 7`); exact `b(233) = s(233) = 9`. Also `p = 257, 281`.
- **H3: `P(p) = (p+3)/2` for every prime.** Witnesses `p = 13, 37, 41` (Theorem 2.7). H3 holds for `43 ≤ p ≤ 1000`.
- **H4: the maximum clique's biclique extension is larger, `|B(C)| > |C|`.** For every maximum clique found (`p ≤ 3000`), `|B(C)| = |C|`.
- **H5: `gp_k(p) = M_k(p)`.** Fails for `k ≥ 4` at most primes (`p = 997`, `k = 4`: `76 < 78`).
- **H6: `s(p) ≤ ω(G_p)`.** Fails at `p = 89` (`A = {0,1,20,44,67,79}`, `A + A ⊆ Q ∪ {0}`, `s = 6 > 5 = ω`) and at `24` further primes `≤ 1000` (`86` of the `211` primes `≤ 3000`).

---

## 5. The sum-product route on the complete case (Task 3): OPEN, with two no-go lemmas (PROVED)

Fix a complete biclique and remove the zero sums: `B_1 := B ∖ (−A)`, so that `A + B_1 ⊆ Q`; this costs at most `|A|` elements of `B` (in the balanced case at most `ω`). Write `n := |A||B_1|`, `r(s) := #{(a,b) ∈ A×B_1 : a+b = s}`, `S := supp r ⊆ Q`, `E := E⁺(A,B_1) = Σ_s r(s)²`.

**5.1 The Weil-amplified moment (PROVED; the `√p` floor).** With `F(x) := Σ_{b∈B_1} χ(x+b)` one has `F(a) = |B_1|` for `a ∈ A`, hence for every `k ≥ 1`

    |A||B_1|^{2k} ≤ Σ_x F(x)^{2k} ≤ p·N_k + (|B_1|^{2k} − N_k)(2k−1)√p,   N_k ≤ (2k−1)!!·|B_1|^k,

where `N_k = (2k)!·[x^{2k}] cosh(x)^{|B_1|}` is the number of `2k`-tuples from `B_1` in which every element occurs with even multiplicity (for those tuples `Π(x+b_i)` is a square and the complete sum is `≤ p`; for the others Weil gives `(2k−1)√p`). Therefore `|A| ≤ (2k−1)√p + (2k−1)!!·p/|B_1|^k`; for `|B_1| ≥ p^{β}` with `kβ ≥ 1/2` this is `|A| ≤ (2k−1+o(1))√p`, and no choice of `k` gives an exponent below `1/2` because the second term of the Weil bound is present for every `k`. The verifier checks both inequalities and `N_k ≤ (2k−1)!!|B_1|^k` exactly for `k = 1,2,3` on every extremal biclique of the table with `|A||B| ≤ 400`.

**5.2 The ratio-set / multiplicative-energy chain.** For `λ ∈ F_p^*` let `ρ(λ) := Σ_t r(λt) r(t) = #{(a,b,a',b') : a'+b' = λ(a+b)}`, and let

    T := #{(a,b,a',b',a'',b'',a''',b''') ∈ (A×B_1)^4 : (a+b)(a'+b') = (a''+b'')(a'''+b''')}

be the multiplicative energy of the sumset counted with multiplicity.

**Lemma 5.1 (PROVED).** (i) `Σ_λ ρ(λ) = n²`, `ρ(λ) = 0` unless `λ ∈ S/S ⊆ Q`, and `ρ(1) = E ≥ n`.
(ii) `T = Σ_λ ρ(λ)²`.
(iii) `T ≥ n⁴/|S/S| ≥ 2n⁴/(p−1)`.
(iv) `T ≥ T_diag := 2E² − Σ_s r(s)⁴ ≥ E² ≥ n²`.
(v) `T ≤ n²·max_λ ρ(λ)`.
Consequently: if `2n² ≤ p − 1` then `2n⁴/(p−1) ≤ n² ≤ E² ≤ T_diag ≤ T`, so (iii) — the only consequence of `S/S ⊆ Q` for `T` — is implied by the trivial diagonal count; and (iii)+(v) give `max_λ ρ(λ) ≥ 2n²/(p−1)`, which is implied by `ρ(1) ≥ n` whenever `n ≤ (p−1)/2`.

*Proof.* (i) `Σ_λ ρ(λ) = Σ_{s,t} r(s)r(t) = n²`; `ρ(λ) ≠ 0` forces `λ = s/t` with `s, t ∈ S ⊆ Q`; `ρ(1) = Σ_t r(t)² = E ≥ Σ_t r(t) = n`. (ii) For `s, t, s', t' ∈ S ⊆ F_p^*`, `st = s't'` iff `s/s' = t'/t =: λ`; grouping by `λ` gives `Σ_λ (Σ_{s'} r(λs')r(s'))(Σ_t r(λt)r(t)) = Σ_λ ρ(λ)²`. (iii) Cauchy–Schwarz over the support of `ρ`, which has at most `|S/S| ≤ |Q| = (p−1)/2` elements. (iv) The solutions with `(s,t) = (s',t')` contribute `Σ_{s,t} r(s)²r(t)² = E²`, those with `(s,t) = (t',s')` contribute `E²`, the overlap `s = t = s' = t'` contributes `Σ r(s)⁴ ≤ E²`, and all other solutions contribute non-negatively. (v) `Σ_λ ρ(λ)² ≤ max ρ · Σ ρ`. ∎

**Lemma 5.2 (the arithmetic-progression obstruction, PROVED).** Let `A = B = {1, …, N} ⊆ F_p` with `2N < p`. Then `E⁺(A,A) = (2N³+N)/3` and `T ≥ E² ≥ (2N³/3)² = (4/9)n³` with `n = N²`. Hence any inequality `T ≤ U(|A|,|B|)` valid for all pairs of subsets of `F_p` — in particular every bound obtained from point–plane or point–line incidence theorems, which do not see whether `S ⊆ Q` — satisfies `U(N,N) ≥ (4/9)N⁶`, and the chain "`2n⁴/(p−1) ≤ T ≤ U`" can yield a contradiction only when `2n⁴/(p−1) > (4/9)n³`, i.e. `n > (2/9)(p−1)`. The chain therefore proves at best `|A||B| ≤ (2/9)(p−1)`, and nothing of the form `|A||B| ≤ p^{1−c}`.

*Proof.* `r(s) = s − 1` for `2 ≤ s ≤ N+1` and `2N+1−s` for `N+1 ≤ s ≤ 2N`, so `E = 2Σ_{j=1}^{N} j² − N² = (2N³+N)/3`; the rest is Lemma 5.1(iv). ∎

**5.3 Verdict on Task 3 (OPEN; the obstruction stated exactly).** The route "`A + B ⊆ Q` ⇒ `(A+B)/(A+B) ⊆ Q` (or `(A+B)(A+B) ⊆ Q`) ⇒ lower bound on the multiplicative energy of the sumset ⇒ contradiction with an incidence-theoretic upper bound" cannot give `|A||B| ≤ p^{1−c}` for any `c > 0`, for two independent reasons, each checked exactly on the extremal bicliques of the table:

1. *The lower bound is either dominated by the trivial solutions or already attained.* Membership of `S` in the index-2 subgroup enters only through `|S/S| ≤ (p−1)/2` (or `|SS| ≤ (p−1)/2`), giving `T ≥ 2n⁴/(p−1)`. (a) For `2n² ≤ p−1` the trivial solutions already give more: `T ≥ T_diag ≥ E² ≥ n² ≥ 2n⁴/(p−1)`. (b) For `2n² > p−1` — the case of every extremal biclique in the table, from `n = b(p)² ≈ 10²` up to `n ≈ p/2` for the `|A| = 2` extremal — the bound is attained by these very sets up to a few percent: for them `S/S = Q` exactly and `ρ` is nearly flat on `Q`, so `T·(p−1)/(2n⁴)` lies between `1.0015` and `1.17` on all extremal bicliques checked (values below), with the non-diagonal solutions alone accounting for `0.84`–`0.99` of the bound. Hence no upper bound `U ≥ T` that is true for these sets can contradict (iii). In the `λ`-form the chain's entire output is `max_λ ρ(λ) ≥ 2n²/(p−1)`, which `ρ(1) = E ≥ n` satisfies for every `n ≤ (p−1)/2` (and `max_λ ρ(λ) = ρ(1)` on every extremal set checked).
2. *No admissible upper bound is strong enough.* Any upper bound on `T` that does not use `S ⊆ Q` is at least `(4/9)n³` on arithmetic progressions (Lemma 5.2), so the chain cannot beat `|A||B| ≤ (2/9)(p−1)`, weaker than Hanson–Petridis' `(p−1)/2 + r`.

The quantity that is provably too large is thus the diagonal part `T_diag ≥ E⁺(A,B_1)²` of the energy (and, for structured pairs, `E⁺` itself, `≈ (2/3)n^{3/2}` for intervals). An input that could work must distinguish `Q` from an arbitrary subset of `F_p^*` of size `(p−1)/2` and must fail over `F_{p²}` (where `A = B = F_p` is a complete biclique with `|A||B| = q`); the only known input with both properties is the Stepanov/Rédei polynomial method, which stops at `p/2` (§1). The dilation-invariance used for Sárközy's problem (`Q·(A+B) = A+B` when `A + B = Q`) is unavailable here: for the extremal bicliques `|A + B| ≤ |A||B| ≪ |Q|`.

**Recorded values (exact integers; ratios rounded).** With `B_1 = B(A) ∖ (−A)`, `n = |A||B_1|`, `Qb = 2n⁴/(p−1)`:

| p | A | n | |S/S| | E = ρ(1) = max ρ | T_diag | T | T/Qb | (T − T_diag)/Qb |
|---|---|---|---|---|---|---|---|---|
| 997 | `{0,1}` (M₂) | 496 | 498 | 728 | 1 057 848 | 121 711 072 | 1.0015 | 0.993 |
| 997 | M₃ extremal `{0,1,93}` | 390 | 498 | 586 | 684 742 | 46 762 180 | 1.0066 | 0.992 |
| 997 | M₄ extremal `{0,1,5,456}` | 308 | 498 | 496 | 489 448 | 18 307 426 | 1.013 | 0.986 |
| 997 | M₆ extremal `{0,1,10,325,345,909}` | 216 | 498 | 274 | 149 458 | 4 451 342 | 1.018 | 0.984 |
| 997 | greedy 12×12 biclique | 144 | 498 | 164 | 53 508 | 896 576 | 1.038 | 0.976 |
| 509 | exact 11×11 biclique | 110 | 254 | 144 | 41 052 | 601 228 | 1.043 | 0.972 |
| 233 | `A = B` sum-clique, 9×9 | 72 | 116 | 136 | 36 400 | 252 832 | 1.091 | 0.934 |
| 101 | exact 6×6 biclique | 30 | 50 | 42 | 3 378 | 18 998 | 1.173 | 0.964 |
| 89 | exact 6×6 biclique | 30 | 44 | 54 | 5 598 | 21 044 | 1.143 | 0.839 |

(For the diagonal pair `(C, −C)` of a maximum clique every row contains a zero sum, so `B_1 = ∅`; that pair is not an input to this chain.) The verifier recomputes `T`, `T_diag`, `E`, `Qb`, `|S/S|`, `max ρ` for every profile extremal, balanced witness and maximum clique in the table and checks (i)–(v) and both consequences exactly.

---

## 6. Open obligations

1. **Below `√p`.** Nothing in this note bounds `ω(G_p)`, `b(p)`, `s(p)` below `(1+√(2p−1))/2`, or `P(p)` below `(p−1)/2 + ω(G_p)`. The balanced problem is at least as hard as the clique problem (`b ≥ ω`); a proof of `b(p) ≤ p^{1/2−c}` would be new.
2. **`P(p) = (p+3)/2` for all `p ≥ 43`.** By Propositions 2.1 and 2.6 this is equivalent to: no complete biclique with `min(|A|,|B|) ≥ 3` (indeed `≥ m*(p)+1`) has `|A||B| ≥ (p+5)/2`; it would imply `ω(G_p)² < (p+5)/2`, i.e. only Hanson–Petridis strength for cliques, but it also excludes all "near-extremal" configurations of Theorem 1.2 with `3 ≤ r ≤ ω`. Verified for `p ≤ 1000`.
3. **Exact `b(p)` and `s(p)` further out.** `b(p)` is exact only where the branch-and-bound finished (see the `exact` flags); the descending-`p` runs give lower bounds with witnesses. `s(p)` is exact to `3000`.
4. **Is `b(p) − ω(G_p)` bounded?** The random model predicts `b = ω + O(1)`: a `k×k` biclique with `A ∩ (−B) = ∅` needs `k²` independent conditions, a clique only `C(k,2)`, a sum-clique `C(k,2) + k`, and intermediate overlaps interpolate; the expected counts are `≥ 1` up to `k ≈ 2log₂p − 2log₂log₂p + O(1)` in all three cases, with the constant differing. The data show differences `0, 1, 2` only. A proof either way is open.
5. **Sum-cliques.** `s(p)` deserves a direct study: `A + A ⊆ Q ∪ {0}` is the "sum graph" analogue of the clique problem, has the same polynomial-method bound, and in the table exceeds `ω(G_p)` at `86` of the `211` primes `p ≤ 3000`, equals it at `90`, and falls below it at `25`.

---

## Verification

`experiments/sigma_biclique_2026_09_05.py` (runs in about [[TIME]] s on this machine; `[[N_CHECKS]]` checks, `[[N_FAIL]]` failures) recomputes `ω(G_p)` and `s(p)` exactly in Python for all `p ≤ 1300` (bitset branch and bound with greedy colouring, cliques normalised to contain `0, 1`, sum-cliques normalised to contain `1` or `n_0`), checks the 21 published prime entries of Brouwer's table, checks every clique, sum-clique and biclique witness of the C++ table, recomputes `M_k(p)` exactly for `k = 2` (`p ≤ 3000`), `k = 3, 4` (`p ≤ 1000`), `k = 5` (`p ≤ 400`), `k = 6` (`p ≤ 200`) and compares with the C++ values, verifies Proposition 2.3 for all `p ≤ 3000`, Proposition 2.4 on every profile extremal, Corollary 2.5's exclusion sets, re-runs the reduced product search of Proposition 2.6 for `p ≤ 500` and the unreduced brute force for `p ≤ 100` (Theorem 2.7), computes exact `b(p)` in Python for `p ≤ 120`, and checks all inequalities of §5 (moment inequality for `k ≤ 3`, Lemma 5.1 (i)–(v) and the two consequences, Lemma 5.2 for four `(p,N)` pairs) exactly with integers and `Fraction`s. The C++ helper `experiments/sigma_biclique_2026_09_05.cpp` (compile with `clang++ -O3 -std=c++17`) produced the table: modes `omega`, `sumclique`, `profile k_min k_max t`, `threshold k T t`, `balanced t`, `product t`, `product3`, `gp k_max`, `greedy iters`; its per-prime outputs are merged into `results/sigma_biclique_2026_09_05_search.json` together with the `exact` flags and run times.
