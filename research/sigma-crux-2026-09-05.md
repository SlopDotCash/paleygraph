# The local-to-average crux: dilation moments, level sets, and near-extremal rectangles

**Status: PROVED — (a) the two-set conjecture implies Vinogradov's conjecture `n_p ≤ 2⌊p^ε⌋` and, more strongly, that no run of `2⌊p^ε⌋−1` consecutive residues or non-residues exists in `[1,p−1]` for large `p`, via an exact triangle-weight identity; (b) an exact identity for the dilation second moment `Σ_t g(t)²` in terms of the signed ratio energy; (c) the exact square/non-square decomposition of `Σ_t g(t)^{2k}`, the generating function for the even-multiplicity count `N_k`, the bound `N_k ≤ (2k−1)!!·R^k`, and the exact value `N_2 = 3n²−2n`; (d) the level-set bound, the fact that it can never fall below `(2k−1)√p·η^{−2k}` (and never below `√p` in any `k` when `|A||B| ≤ √p`), the equivalence `C(ε,δ) ⟺ T(ε,δ)`, an unconditional density statement `Λ_{1/2+2kδ}`, and the invariance of every termwise-Weil dilation-moment bound under `(A,B) ↦ (uA,B)`, which is the precise form of the local-to-average obstruction; (e) a geometric-progression construction showing that complete rectangles with `≥ ¼log₂p − O(log log p)` biased dilates exist for every large prime. REFUTED (with witnesses) — H1 (small additive doubling), H2 (small multiplicative doubling), H3 (biased dilate set is a union of subgroup cosets), H4 in its fixed-`k` form (the maximum bias at `|A|=|B|=k` is exactly 1 for `k=4,5,6` at all 214 primes `101 ≤ p ≤ 1499`), and the tuple bound `N_k ≤ k!(|A||B|)^k`. SUPPORTED, not proved — sumset compression of near-extremal rectangles and Conjecture SI (spike isolation: the biased dilate set has size `O(log p)`, tested with zero violations on 1243 stored rectangles). OPEN — the conjecture itself, Conjecture SI, and any transfer from level-set density to the single dilate `t = 1`.**

Throughout `p` is an odd prime, `χ = χ_p` the Legendre symbol with `χ(0) = 0`, `A, B ⊆ F_p`, `m = |A|`, `n = |B|`, and

    S(A,B) = Σ_{a∈A} Σ_{b∈B} χ(a+b),        g(t) = g_{A,B}(t) = Σ_{a∈A} Σ_{b∈B} χ(ta+b)   (t ∈ F_p),

so that `S(A,B) = g(1)` and `g(t) = S(tA,B)`. The two-set conjecture is Conjecture 7 of Satake
(arXiv:2011.02907, `sources/satake-2011.02907.txt`, lines 204–209, stated there in the `χ(s−t)` form: "for
`0<α≤1` and `β>0`, property `P(α,β)` holds if for every pair `S,T ⊂ F_p` with `|S|,|T| > p^α`,
`|Σ_{s∈S,t∈T} χ(s−t)| ≤ p^{−β}|S||T|`"). We write `C(ε,δ)` for: there is `p₀` such that for all `p > p₀`
and all `A,B` with `|A|,|B| > p^ε`, `|S(A,B)| ≤ p^{−δ}|A||B|`. The `a−b` and `a+b` forms are equivalent
(`B ↦ −B`).

Verifier: `experiments/sigma_crux_2026_09_05.py` (standard library + numpy, ≈ 30 seconds, 59 738 exact checks), output
`results/sigma_crux_2026_09_05.json`. Stored search witnesses (exhaustive maxima, annealing lower bounds,
`k*` lower bounds): `results/sigma_crux_2026_09_05_search.json`. The verifier re-verifies every stored witness
and every identity below; it does not redo any search. Every numerical claim in this note is one of its checks
unless explicitly marked HEURISTIC.

---

## 1. Interval rectangles: Vinogradov and runs (PROVED)

For `N ≥ 1` put `w_N(s) = min(s−1, 2N+1−s)` for `2 ≤ s ≤ 2N` (the triangle weight); `Σ_{s=2}^{2N} w_N(s) = N²`.

**Proposition 1.1.** Let `2N < p` and `A = B = {1,…,N}`. Then

    S(A,A) = Σ_{s=2}^{2N} w_N(s) χ(s),

no `χ(0)` term occurs, and `S(A,A) = N²` if and only if `χ(s) = 1` for all `2 ≤ s ≤ 2N`, if and only if
`n_p > 2N`, where `n_p` is the least positive quadratic non-residue mod `p`.

*Proof.* For `a,b ∈ [1,N]` the integer `a+b` lies in `[2,2N] ⊂ [1,p−1]`, so `a+b ≢ 0` and `χ(a+b)` is the
Legendre symbol of a nonzero residue; `w_N(s)` is the number of representations `s = a+b`. Since `w_N > 0` on
`[2,2N]`, `|Σ w_N χ| ≤ Σ w_N = N²` with equality exactly when all `χ(s)` (`2 ≤ s ≤ 2N`) are equal to `+1`
(the value `−1` throughout would give `−N²`). As `χ(1) = 1`, "all of `2,…,2N` are residues" is "`n_p > 2N`". ∎

**Corollary 1.2 (the conjecture implies Vinogradov and a run bound).** Fix `0 < ε' < ε < 1` and assume
`C(ε',δ)`. Then for all sufficiently large `p`, with `N = ⌊p^ε⌋` (so `N > p^{ε'}` and `2N < p`):

1. `n_p ≤ 2N`;
2. `|Σ_{s=2}^{2N} w_N(s) χ(s)| ≤ p^{−δ} N²` (cancellation `o(N²)` in the triangle-weighted short sum);
3. the longest run of consecutive quadratic residues in `[1,p−1]`, and likewise of consecutive non-residues,
   has length at most `2N − 2`.

*Proof.* (1)–(2): apply `C(ε',δ)` to `A = B = {1,…,N}` and use 1.1; if `n_p > 2N` then `S = N² > p^{−δ}N²`.
(3): if `u+2, …, u+2N` are `2N−1` consecutive residues in `[1,p−1]`, take `A = {1,…,N}`, `B = {u+1,…,u+N}`;
every `a+b` lies in this run, `0 ∉ A+B`, and `S(A,B) = N²`; for non-residues `S(A,B) = −N²`. Both contradict
`C(ε',δ)`. ∎

The verifier checks the identity of 1.1 and the three-way equivalence at every prime `5 ≤ p < 3000` for
`ε ∈ {0.3, 0.45}` (856 cases), and constructs the run rectangles with `S = ±N²` from the actual longest runs
(the check `run_rectangle_cases`). Vinogradov's conjecture is open (brief, background facts); hence any proof
of `C` must in particular beat Burgess-range short character sums, and `1.2(3)` says it must do so for runs
starting anywhere, not only at `1`.

---

## 2. The dilation second moment (PROVED)

Let `Π = (A×B) ∖ {(0,0)}`. For `(a,b) ∈ Π` with `a ≠ 0` its *class* is `r = b/a ∈ F_p`; for `a = 0` (then
`b ≠ 0`) its class is `∞`. Let `μ_r`, `μ_∞` be the class multiplicities, and

    R(A,B) = Σ_r μ_r² + μ_∞²          (ratio energy: #{(π,π') ∈ Π² : ab' = a'b}, with (0,b),(0,b') counted),
    σ_r = Σ_{(a,b) in class r} χ(a),   R_χ(A,B) = Σ_{r∈F_p} σ_r²,
    s_{A*} = Σ_{a∈A, a≠0} χ(a),        s_B = Σ_{b∈B} χ(b).

**Theorem 2.1.** `Σ_{t∈F_p} g(t)² = p·R_χ(A,B) − |B|²·s_{A*}² + p·1_{0∈A}·s_B².`

*Proof.* `Σ_t g(t)² = Σ_{(a,b),(a',b')∈A×B} Σ_t χ(ta+b)χ(ta'+b')`. Three cases. If `a = a' = 0` the inner
sum is `p·χ(b)χ(b')`; summed over `b,b' ∈ B` this gives `p·s_B²`, present only if `0 ∈ A`. If exactly one of
`a,a'` is `0`, say `a' = 0`, the inner sum is `χ(b')·Σ_t χ(ta+b) = χ(b')·Σ_x χ(x) = 0`. If `a,a' ≠ 0`, then
`Σ_t χ(ta+b)χ(ta'+b') = χ(aa')·Σ_t χ((t+b/a)(t+b'/a'))`, which is `χ(aa')(p−1)` if `b/a = b'/a'` and
`−χ(aa')` otherwise (background fact `Σ_x χ(x)χ(x+c) = −1`, `c ≠ 0`). So the third case contributes
`p·Σ_{ab'=a'b} χ(a)χ(a') − Σ_{a,a'≠0}Σ_{b,b'} χ(a)χ(a') = p·R_χ − |B|²s_{A*}²`. ∎

**Corollary 2.2.** (i) `Σ_t g(t)² ≤ p·R(A,B)`, since `|σ_r| ≤ μ_r` and, when `0 ∈ A`, `|s_B| ≤ |B∖{0}| = μ_∞`.
(ii) If `0 ∉ A` then `μ_r ≤ min(m,n)` for every class, so `R ≤ mn·min(m,n)`; in general `R ≤ mn·n + n²`.
(iii) (Chebyshev at `k = 1`) `#{t : |g(t)| ≥ η mn} ≤ pR/(η mn)²`; since `R ≥ |Π| ≥ mn − 1`, this is `< 1` only
if `mn > p/η² − O(1)`: the square-root barrier of `research/sigma-barriers-2026-09-05.md` §B1 in dilation form.

The identity is checked exactly (integers) on 360 cases: 320 random pairs at `p ≤ 29` including all
combinations of `0 ∈ A`, `0 ∈ B`, and 40 stored witnesses at `101 ≤ p ≤ 1499`.

---

## 3. The dilation `2k`-th moment (PROVED)

Expand `g(t)^{2k} = Σ_{π ∈ Π^{2k}} Π_{i=1}^{2k} χ(a_i t + b_i)` (tuples containing `(0,0)` contribute `0` and
are excluded). For a tuple `π` put `P_π(t) = Π_i (a_i t + b_i)`, `d_π = #{i : a_i ≠ 0} = deg P_π`, and
`c_π = Π_{a_i≠0} a_i · Π_{a_i=0} b_i ≠ 0`, so `P_π = c_π · Π_{r} (t + r)^{μ_r(π)}` with `μ_r(π)` the number of
positions of `π` in the finite class `r`. Call `π` *square* if every `μ_r(π)` is even (the class `∞` is
unconstrained: it only contributes to `c_π`).

**Theorem 3.1.** For every `k ≥ 1`,

    Σ_{t∈F_p} g(t)^{2k} = Q_k + W_k,     Q_k = Σ_{π square} Σ_t χ(P_π(t)),     W_k = Σ_{π not square} Σ_t χ(P_π(t)),

and:

1. `Q_k = Σ_t Q_k(t)` with the exact pointwise formula, for `0 ∉ A`,
   `Q_k(t) = (2k)!·[x^{2k}] Π_{r ≠ −t} cosh(σ_r x)`, hence `0 ≤ Q_k(t) ≤ N_k^χ ≤ N_k`, where
   `N_k^χ = (2k)!·[x^{2k}] Π_r cosh(σ_r x)` and

       N_k = #{square tuples} = (2k)!·[x^{2k}] e^{μ_∞ x} Π_{r∈F_p} cosh(μ_r x).

   (For `0 ∈ A` the factor `e^{s_{B*}x}` with `s_{B*} = Σ_{b∈B,b≠0}χ(b)` multiplies the product and `Q_k(t)`
   may be negative; the upper bound `Q_k(t) ≤ N_k` persists.) In all cases `Q_k ≤ p·N_k`.
2. `|W_k| ≤ (2k−1)√p · T_k`, where `T_k = |Π|^{2k} − N_k` is the number of non-square tuples.

Consequently `Σ_t g(t)^{2k} ≤ p·N_k + (2k−1)√p·(|Π|^{2k} − N_k)`.

*Proof.* If `π` is square then `P_π = c_π Q(t)²` with `Q = Π_r (t+r)^{μ_r(π)/2}`, so
`Σ_t χ(P_π(t)) = χ(c_π)·#{t : Q(t) ≠ 0} = χ(c_π)(p − ρ_π)`, `ρ_π` the number of distinct roots; in absolute
value `≤ p`. The pointwise formula: for fixed `t`, `Π_i χ(a_i t + b_i) = Π_i χ(a_i)·Π_r χ(t+r)^{μ_r(π)}`, which
for a square tuple equals `Π_i χ(a_i)` if no class `r = −t` is present and `0` otherwise; summing
`Π_i χ(a_i)` over square tuples with prescribed even class counts `j_r` gives the multinomial
`(2k)!/Π j_r! · Π_r σ_r^{j_r}`, i.e. the stated coefficient of `Π_{r≠−t} cosh(σ_r x)`; all coefficients of
`cosh` are nonnegative and increase with `|σ_r| ≤ μ_r`, giving `0 ≤ Q_k(t) ≤ N_k^χ ≤ N_k`. Setting all signs to
`+1` and admitting arbitrary counts in class `∞` gives the generating function for `N_k`. If `π` is not square,
some finite class has odd multiplicity, so `P_π` is not a constant times a square, `1 ≤ d_π ≤ 2k`, and the Weil
bound (brief, background facts) gives `|Σ_t χ(P_π(t))| ≤ (d_π − 1)√p ≤ (2k−1)√p`. ∎

**Proposition 3.2 (counting square tuples).** If `0 ∉ A` then `N_k ≤ (2k−1)!!·R(A,B)^k`. For `k = 2` and
`n` classes each of multiplicity `1` (so `R = n`) the exact value is `N_2 = 3n² − 2n`; the bound
`(2k−1)!!R^k = 3n²` is sharp up to the lower-order term, while `k!·(|A||B|)^k = 2n²` is **false** for
`n ≥ 3` (`n = 3`: `N_2 = 21 > 18`).

*Proof.* Every square tuple arises from a perfect matching of the `2k` positions (`(2k−1)!!` choices) together
with, for each matched pair, an ordered pair of elements of `Π` in the same finite class (`Σ_r μ_r² = R`
choices): pair the positions of each class arbitrarily. This map is onto the square tuples, so
`N_k ≤ (2k−1)!!R^k`. For `k = 2` with `n` singleton classes, a square 4-tuple is either constant (`n` tuples)
or uses two classes twice each (`C(n,2)·C(4,2) = 3n(n−1)` tuples), total `3n²−2n`. ∎

By 2.2(ii), for `0 ∉ A`: `N_k/(mn)^{2k} ≤ (2k−1)!!·(R/(mn)²)^k ≤ (2k−1)!!/max(m,n)^k`.

Verification: `N_k` from the generating function agrees with brute-force enumeration of all `2k`-tuples for
`k ≤ 2` on 80 tiny cases (`p ≤ 17`, `0` allowed in `A` or `B`); the pointwise decomposition
`Q_k(t) + W_k(t) = g(t)^{2k}` and `Q_k(t) ≤ N_k` (and `≥ 0` when `0 ∉ A`) hold on all of them; the bound
`N_k ≤ (2k−1)!!R^k` on the 40 of them with `0 ∉ A`; the termwise Weil inequality is checked (by squaring, exact
integers) for 47 820 non-square tuples; the inequality
of 3.1 is checked for `k = 1,2,3` on 300 stored witnesses at `101 ≤ p ≤ 1499` with `mn` up to `400`; the
exact `N_2 = 3n²−2n` for `n ≤ 8`.

---

## 4. Level sets, the `√p` floor, and the transfer statement

For `η ∈ (0,1]` let `L(η) = #{t ∈ F_p : |g(t)| ≥ η mn}`.

**Corollary 4.1 (PROVED).** For `0 ∉ A` and every `k ≥ 1`,

    L(η) ≤ η^{−2k} [ p·N_k/(mn)^{2k} + (2k−1)√p ] ≤ η^{−2k} [ (2k−1)!!·p/max(m,n)^k + (2k−1)√p ],

and for `k = 1` the exact form `L(η) ≤ pR/(η mn)²`. (Chebyshev: `L(η)(η mn)^{2k} ≤ Σ_t g(t)^{2k}`, then 3.1,
3.2, and `T_k ≤ (mn)^{2k}`.)

**Proposition 4.2 (the floor; PROVED).** (i) For every `k ≥ 2` and `η ≤ 1` the right side of 4.1 is at least
`(2k−1)√p·η^{−2k} ≥ 3√p`. (ii) For `k = 1`, `pR/(η mn)² ≥ p/(η² mn) ≥ √p/η²` whenever `mn ≤ √p`. Hence, when
`mn ≤ √p`, no choice of `k` gives a level-set bound below `√p`; for arbitrary sizes no `k ≥ 2` does; and the
`k = 1` bound drops below `1` only past the square-root barrier `mn > p/η²`. In particular the termwise-Weil
moment method never certifies `L(η) = 0`, i.e. never controls an individual dilate.

The verifier evaluates all three bounds (`k = 1,2,3`) on 99 random pairs with `mn ≤ √p` at
`p ∈ {101, 401, 1009, 1499}` and confirms `min_k bound ≥ √p` in every case, and checks the Chebyshev step
`L(η)(η mn)^{2k} ≤ Σ_t g(t)^{2k}` for `η ∈ {½, ¾, 1}`, `k ≤ 3`, on 300 stored witnesses (2700 cases).

**Proposition 4.3 (unconditional density; PROVED).** Let `0 < ε < 1`, `k ≥ 1/(2ε)`, `0 ∉ A`,
`|A|,|B| > p^ε`, and `0 ≤ δ < 1/(4k)`. Then

    #{t ∈ F_p : |g(t)| ≥ p^{−δ}|A||B|} ≤ ((2k−1)!! + 2k − 1)·p^{1/2 + 2kδ} = o(p).

*Proof.* In 4.1 take `η = p^{−δ}`; `max(m,n)^k > p^{kε} ≥ p^{1/2}` so the first term is at most
`(2k−1)!!√p`. ∎ (If `0 ∈ A`, write `g = s_B + g_{A∖{0},B}` with `|s_B| ≤ |B|`; for `|A| ≥ 2p^{δ}+1` the count
is bounded by that of `(A∖{0},B)` at level `½p^{−δ}|A∖{0}||B|`, i.e. the same bound with `2^{2k}` in front.)

**Definition 4.4.** For `0 < ε < 1`, `δ > 0`, `κ ∈ ℝ`:

* `T(ε,δ)` (transfer): there is `p₀` such that for all `p > p₀`, all `A,B` with `|A|,|B| > p^ε`, and all
  `t ∈ F_p^*`: `|g_{A,B}(t)| ≤ p^{−δ}|A||B|`.
* `Λ_κ(ε,δ)` (density): there is `p₀` such that for all `p > p₀` and all such `A,B`:
  `#{t ∈ F_p^* : |g_{A,B}(t)| > p^{−δ}|A||B|} ≤ p^κ`.

**Proposition 4.5 (PROVED).**

1. `C(ε,δ) ⟺ T(ε,δ)`: `g_{A,B}(t) = S(tA,B)` and `|tA| = |A|`, so `T` is `C` applied to the pairs `(tA,B)`;
   conversely `t = 1`.
2. `C(ε,δ) ⟺ Λ_κ(ε,δ)` for every `κ < 0` (the count must be `0`), and `C(ε,δ) ⟹ Λ_κ(ε,δ)` for every `κ`.
3. `Λ_{1/2+2kδ}(ε,δ)` holds unconditionally (up to the constant of 4.3) for `k ≥ 1/(2ε)`, `δ < 1/(4k)`.
4. For `κ ≥ 0`, the hypothesis `Λ_κ(ε,δ)` is invariant under `(A,B) ↦ (uA,B)`, `u ∈ F_p^*`, and is satisfied
   by any pair whose dilation profile has at most `p^κ` biased dilates, in particular by a pair with
   `|S(A,B)| = |A||B|` whose other dilates are unbiased. Hence `Λ_κ` for `κ ≥ 0` carries no information
   about the dilate `t = 1` beyond membership in a set of size `≤ p^κ`; the passage from `Λ_κ` (`κ ≥ 0`) to
   `Λ_κ` (`κ < 0`), i.e. to `C`, is the whole content of the conjecture.

**Proposition 4.6 (invariance of the moment bound; PROVED).** The quantities `N_k`, `N_k^χ`, the profile
`t ↦ Q_k(t)` up to the relabelling `t ↦ t/u`, and the bound of Theorem 3.1 depend on `(A,B)` only through the
multiset of signed class sums `{σ_r}` (and `μ_∞`, `s_{B*}`). Under `(A,B) ↦ (uA,B)` the classes are
relabelled `r ↦ r/u` and `σ_r ↦ χ(u)σ_{r/u}`, so this multiset (up to sign) and hence the whole bound is
unchanged, and so is the level-set bound of 4.1. Therefore every estimate of this family assigns the *same*
bound to all `p−1` members of the dilation orbit and cannot single out `t = 1`.

*Proof.* Immediate from the formulas of 3.1 and `Σ_{a∈uA} χ(a) = χ(u)Σ_{a∈A}χ(a)` classwise. ∎

**What this isolates.** A proof of `C` must use information about `(A,B)` that is not a function of the ratio
multiplicities `{μ_r}` and signed sums `{σ_r}`. The dilation family sees `B/A = {b/a}`; the sum `S(A,B)` is
`Σ_s r_{A+B}(s)χ(s)`, a function of the additive representation function `r_{A+B}`. The obstruction is that
the moment method controls `t`-averages of a family indexed by *multiplicative* structure, while the target is a
single value determined by *additive* structure; nothing in §§2–4 relates `r_{A+B}` to `{σ_r}`.

---

## 5. Empirics from the stored search data

### 5.1 Data and validation (PROVED consistency)

`results/sigma_crux_2026_09_05_search.json` holds: (E) exact maxima `max_{|A|=m,|B|=n} |S(A,B)|` for
`n = 4` (`m ∈ {4,5,6,8,12,20,40}`) at all 214 primes `101 ≤ p ≤ 1499`, `n = 5` (`m ∈ {5,6,8,10,20}`) at the 37
primes `≤ 293`, `n = 6` (`m ∈ {6,8,12}`) at the 10 primes `≤ 149` — 1713 `(p,m,n)` triples, each obtained by
enumerating all `B ⊇ {0,1}` (affine normalisation, WLOG for `|B| ≥ 2` since `(A,B) ↦ (A−c, B+c)` and
`(sA,sB)` preserve `|S|`) and taking `A` to be the `m` largest or `m` smallest values of `x ↦ Σ_b χ(x+b)`
(exact for fixed `B`); (A) simulated-annealing lower bounds for `k×k`, `6 ≤ k ≤ 20`, and seven lopsided shapes
at 26 primes `101 ≤ p ≤ 1499` (404 triples, up to 3 distinct local optima each); (K) for each of the 26
primes the largest `k` at which a complete `k×k` rectangle (`|S| = k²`) was found.

Validation performed by the verifier: (i) the exhaustive routine reproduces all seven independently stored
extrema of `results/exact_checks.json` (`p = 5,7,13,17,29,37`; maxima `2,4,7,11,14,18,18`; the stored
witnesses there are in the `a−b` form and reproduce exactly under `B ↦ −B`); (ii) all 4601 stored witnesses
(`3362` exhaustive, `1039` annealing, `200` `k*`) have their `S` recomputed exactly and agree with the
stored values and maxima; (iii) 167 exhaustive maxima (`n = 4`, `p ≤ 200`; `n = 5`, `p ≤ 110`) are recomputed
from scratch and agree.

### 5.2 H4 (max bias for `|A| = |B| = k` decays like `p^{−c(k)}`): REFUTED for fixed `k`

For `(m,n) ∈ {(4,4),(5,4),(6,4),(8,4),(5,5),(6,5),(8,5),(6,6)}` the exact maximum equals `mn` (bias `1`) at
*every* prime of the grid; `(12,4)` has bias `1` at every grid prime `p ≥ 127`, `(10,5)` for `p ≥ 137`, `(8,6)` for
`p ≥ 149`, `(20,4)` for `p ≥ 223`, `(20,5)` for `p ≥ 293`, `(40,4)` for `p ≥ 479` (the verifier's
`complete_for_all_grid_primes_from`); only `(12,6)`, enumerated at `p ≤ 149` alone, never reaches bias `1` there
(maximum `65/72`). So `c(4) = c(5) = c(6) = 0` on `101 ≤ p ≤ 1499`, and the annealing lower
bounds at fixed `k ∈ {12,14,16,18,20}` *increase* with `p` (least-squares slopes of `log(bias)` against `log p`:
`+0.16, +0.19, +0.23, +0.26, +0.30`; recorded in the results file as HEURISTIC). This is
not an artefact of the range:

**Proposition 5.1 (PROVED; elementary).** Let `d = (p−1)/2`. If `k ≥ 1` and `p·C(d,k)/C(p,k) ≥ k` then there
exist `A,B` with `|A| = |B| = k` and `S(A,B) = ±k²`. (From `Σ_{|U|=k} |N(U)| = p·C(d,k)`,
`N(U) = {b : χ(a−b) = 1 ∀a∈U}`, proved in `research/moments-and-obstructions.md` §5; choose `U` with
`|N(U)| ≥ k` and `B = −(k elements of N(U))`.) Since `C(d,k)/C(p,k) ≥ 2^{−k}(1 − k²/(p−k+1))` (loc. cit.),
this holds for `p ≥ k·2^{k+1}` and `k ≥ 4`, say. Hence for every fixed `k` the maximum bias is exactly `1` for
all large `p`, and no decay `p^{−c(k)}` with `c(k) > 0` is possible.

The meaningful question is `k` growing with `p`. HEURISTIC: the `k*` lower bounds satisfy
`k* ≈ 1.1·log₂p` (ratio range `0.86`–`1.27`), below the random-model threshold `k₀(p) = max{k : C(p,k)²2^{−k²} ≥ 1}`
(`k₀ = 9,…,15` versus `k* = 6,…,13` across the 26 primes; both recorded by the verifier), consistent with annealing being incomplete rather than
with a genuine deficit; the annealing biases at `x = k/log₂p ∈ [1.2, 2.7]` fit `log(bias) ≈ −0.44(x − 1.2)`
(143 points; recorded as HEURISTIC), but these are lower bounds from an incomplete search at `p ≤ 1499` and no exponent `c` is
extracted from them. H4 in the form "bias `≤ exp(−c·k/log p)`" is OPEN.

### 5.3 H1 (small additive doubling) and H2 (small multiplicative doubling): REFUTED

For every stored witness the verifier computes `|A+A|/|A|`, `|B+B|/|B|`, `|A·A|/|A|`, `|B·B|/|B|`, minimal
AP/GP lengths, ratio energy `R/mn`, `|A+B|/mn`, `E_+(A,B)/mn` (`E_+ = Σ_s r_{A+B}(s)²`), and compares each with
the median over 3 random pairs of the same sizes in the same field. Paired ratios (witness / baseline median):

| population | N | add. doubling `A` | add. doubling `B` | mult. doubling `A` | mult. doubling `B` | `R/mn` | `|A+B|/mn` | `E_+/mn` |
|---|---|---|---|---|---|---|---|---|
| annealing, `k×k`, bias `< 1` | 434 | 1.005 | 1.01 | 1.00 | 1.00 | 1.00 | **0.90** | **1.14** |
| annealing, `m ≥ 3n`, bias `< 1` | 195 | 1.00 | 1.00 | 1.00 | 1.00 | 0.98 | **0.89** | **1.16** |
| annealing, bias `= 1`, `mn ≥ 64` | 266 | 0.99 | 1.00 | 1.00 | 1.00 | 1.00 | 0.94 | 1.12 |

(Medians of the paired ratios, from the verifier's fixed-seed baselines; the 10–90 % ranges of the four doubling
ratios lie within `0.91`–`1.07`.) Near-extremal rectangles have the additive and
multiplicative doubling of random sets: H1 and H2 are refuted. The explicit witness recorded in the results file
is the `20×20` annealing rectangle at `p = 1009` with `S = 333` (bias `0.83`), whose doubling constants
`(|A+A|/|A|, |B+B|/|B|, |A·A|/|A|, |B·B|/|B|) = (9.15, 9.25, 9.65, 10.0)` compare with random medians
`(9.65, 9.45, 9.65, 9.65)`. Caveat on the exhaustive witnesses: they are the lexicographically first optimal `B`
among typically hundreds (the count is capped at 400 per triple), so their `B ≈ {0,1,2,…}` and the forced
`A ⊂ ±Q` (from `0 ∈ B` and bias `1`) are enumeration artefacts, not structure; they are excluded from the
H1/H2 verdict.

What *is* systematically different is the sumset: `|A+B|/mn` is below the random baseline in 99–100 % of the
near-extremal annealing rectangles and `E_+(A,B)/mn` above it in 1158 of 1170 witnesses (99.0 %; median ratio
`1.20`). The multiplicity histogram of `A+B` has a longer tail
(median maximal multiplicity `4`–`5` versus `3`–`4` at `20×20`). SUPPORTED, with the caveat that within a single
annealing task the top local optima show no monotone relation between `|S|` and `E_+` (mean Spearman
`−0.11`), so the effect is a property of the near-extremal *class*, not a monotone predictor of bias. HEURISTIC
explanation: `S = Σ_s r_{A+B}(s)χ(s)` is a signed sum over `|A+B|` "coins" with total weight `mn`; fewer coins
make full alignment likelier. It is not a necessary feature: the construction of Theorem 6.2 gives bias `1` with
`|A+B| = mn(1 − o(1))`.

### 5.4 H3 (biased dilate set is a union of subgroup cosets): REFUTED

Let `D = D_{1/2}(A,B) = {t ∈ F_p^* : |g(t)| ≥ ½|S(A,B)|}` and `stab(D) = {h ∈ F_p^* : hD = D}`; `D` is a union of
cosets of a nontrivial subgroup iff `|stab(D)| > 1`. Over the 1500 stored rectangles (one exhaustive witness per triple) with `mn ≥ 64` and bias `≥ ½`,
`|stab(D)| > 1` in exactly 6 (stabiliser orders `2, 2, 2, 2, 4, 5`, all with `|D| ≤ 15`), and the stabiliser is
trivial for every rectangle with `|D| ≥ 16`. Explicit witness (in the results file): the `k*` rectangle at `p = 1009`,
`11×11`, `S = 121`, with `|D| = 44` and trivial stabiliser. (A trivial PROVED converse: if `HA = A` for a
subgroup `H ⊆ Q`, then `g(th) = g(t)` for `h ∈ H` and `D` is `H`-invariant; the data contain no such
rectangle with `|H| > 1`.)

### 5.5 Spike isolation and the exceptional class

Removing `0` from `A` and `B` (which otherwise adds the constant `χ(t)s_A`, resp. `s_B`, to every dilate and
inflates `D` for the `n = 4` exhaustive witnesses), 1243 stored rectangles have `mn ≥ 64` and bias `≥ ½`. For
them `|D_{1/2}|` has median `1`, is exactly `{1}` in 95.4 % (96.2 % of the annealing, 92.0 % of the exhaustive,
96.7 % of the `k*` witnesses), is `≤ log₂p` in 99.52 %, and its maximum is `53 = 5.31·log₂p` (the `p = 1009` rectangle above with `0` removed:
`10×11`, both `A` and `B` exact geometric progressions inside a coset of the order-36 subgroup). Exactly six
rectangles have `|D| > log₂p`: in three of them (`p = 1009`, `271`, `241`) both `A` and `B` are exact
geometric progressions, in a fourth (`p = 127`, `9×9`, `|D| = 9`) `A` is, the fifth (`p = 409`, `19×19`,
`|D| = 16 = 1.84·log₂p`) has `A` inside a GP of length `2.05|A|`, and the sixth is the `p = 151`, `9×9`
rectangle with `|D| = 10 = 1.38·log₂p` (three of the six came from the annealer's GP seeds). Conversely
GP-compactness does not force a large `D`: among the 39 rectangles with `min(gp(A)/|A|, gp(B)/|B|) ≤ 2`, 4
have `|D| > log₂p`.

Comparison with §4: the moment method allows `|D| ≲ 3√p·2^4 ≈ 1.9·10³` at `p = 1499`, `η = ½`, `k = 2`; the
observed `|D|` is `1` in the bulk and `≤ 53` in the worst case. The spike at `t = 1` is, empirically, an
isolated point of the dilation profile, and the only mechanism producing several biased dilates is
multiplicative self-similarity (`rU ≈ U`), exactly the mechanism of Theorem 6.2.

---

## 6. One structural conjecture, its test, and its proved part

**Conjecture SI (spike isolation).** Let `p` be a prime and `A, B ⊆ F_p^*` with `|A||B| ≥ 64` and
`|S(A,B)| ≥ ½|A||B|`. Put `D_{1/2}(A,B) = {t ∈ F_p^* : |g_{A,B}(t)| ≥ ½|S(A,B)|}`. Then

1. (universal) `|D_{1/2}(A,B)| ≤ 8·log₂p`;
2. (generic) if neither `A` nor `B` is contained in a geometric progression `{c r^i : 0 ≤ i < L}` with
   `L ≤ 2|A|`, resp. `L ≤ 2|B|`, then `|D_{1/2}(A,B)| ≤ 2·log₂p`.

**Test.** Both clauses hold for all 1243 stored rectangles described in §5.5 (zero violations; maximal ratios
`5.31` for clause 1 and `1.84` — the `p = 409`, `19×19` annealing rectangle with `gp(A) = 2.05|A|` — for clause 2).
The verifier records the eight largest cases with their `A`, `B`, GP lengths and stabilisers.

**What is proved about it.**

**Theorem 6.1 (trivial upper bound, PROVED).** `|D_{1/2}(A,B)| ≤ 4pR(A,B)/S(A,B)²`, by 2.2(i); for
`|S| ≥ ½mn` and `R ≤ mn·min(m,n)` this is `≤ 16p·min(m,n)/(mn)`, which is `O(√p)` only when `max(m,n) ≥ √p`.
Nothing in §§2–4 improves `O(√p)` in general (Proposition 4.2).

**Theorem 6.2 (the logarithm is necessary, PROVED).** Let `k ≥ 4`, let `r` be a primitive root mod `p`,
`U = {r^i : 0 ≤ i < k}`, and `N(U) = {b ∈ F_p : χ(u+b) = 1 for all u ∈ U}`. Then

    |N(U)| ≥ 2^{−k}[ p − √p·((k−2)2^{k−1} + 1) ] − k/2,

and if `N(U) ≠ ∅` the rectangle `(U, N(U))` satisfies `S(U,N(U)) = k|N(U)|` (bias `1`) and
`D_{1/2}(U,N(U)) ⊇ {r^j : |j| ≤ ⌊k/4⌋}`, so `|D_{1/2}| ≥ 2⌊k/4⌋ + 1`. Taking
`k = ⌊½log₂p − log₂log₂p⌋` gives, for every sufficiently large prime, a complete rectangle with
`|A||B| ≥ 64` and `|D_{1/2}| ≥ ¼log₂p − O(log log p)`. In particular the `log p` in Conjecture SI cannot be
replaced by a constant.

*Proof.* `Σ_b Π_{u∈U}(1+χ(u+b)) = Σ_{I⊆U} Σ_b Π_{u∈I}χ(u+b)`; the term `I = ∅` is `p`, terms with `|I| = 1`
vanish, and for `|I| ≥ 2` the polynomial `Π_{u∈I}(x+u)` has `|I|` distinct roots, is not a constant times a
square, and Weil gives `|Σ_b| ≤ (|I|−1)√p`; `Σ_{j≥2} C(k,j)(j−1) = (k−2)2^{k−1} + 1`. On the other hand the
same sum equals `2^k|N(U)|` plus the contribution of the `k` values `b = −u₀`, each at most `2^{k−1}`
(one factor is `1+χ(0) = 1`); note `N(U) ∩ (−U) = ∅`. This gives the lower bound. Every `u+b` with `b ∈ N(U)`
is a nonzero residue, so `S = k|N(U)|`. For `0 < j ≤ k/4`, `r^jU ∩ U = {r^j,…,r^{k−1}}` has `k−j` elements and
`r^jU ∖ U` has `j` (the powers `r^0,…,r^{k+j−1}` are distinct since `k+j < p−1`), so
`g(r^j) = S(r^jU, N(U)) ≥ (k−j)|N(U)| − j|N(U)| = (k−2j)|N(U)| ≥ ½k|N(U)| = ½S`; symmetrically for
`−k/4 ≤ j < 0`. For the choice of `k`: `k²4^k ≤ p` gives `√p ≥ k2^k`, hence
`p − √p((k−2)2^{k−1}+1) ≥ √p·k2^{k−1} ≥ k²2^{2k−1}`, so `|N(U)| ≥ k²2^{k−1} − k/2 ≥ 64/k`. ∎

The verifier checks the exact identity `Σ_b Π(1+χ(u+b)) = 2^k|N(U)| + Σ_{u₀}Π_{u≠u₀}(1+χ(u−u₀))`, the lower
bound, `S = k|N(U)|`, and the inclusion `{r^j : |j| ≤ ⌊k/4⌋} ⊆ D_{1/2}` for `p ∈ {101,…,4001}`, `k ∈ {4,…,7}`
(28 cases).

**Theorem 6.3 (what SI would give and would not give; PROVED implications).** (i) Conjecture SI implies that
for every subgroup `H ⊆ Q` (the residues) with `|H| > 8log₂p` and every `B ⊆ F_p^*` with `|B| ≥ 64/|H|`,
`|S(H,B)| < ½|H||B|`: indeed `g_{H,B}` is constant on cosets of `H`, so `D_{1/2}` is `H`-invariant and
nonempty as soon as `|S(H,B)| ≥ ½|H||B|`, forcing `|D_{1/2}| ≥ |H|`. Thus SI contains a weak form of the
subgroup (one-set-a-coset) case of the conjecture and is at least as hard. (ii) SI does not imply `C`: it is
consistent with `|S(A,B)| = |A||B|` for a single dilate. (iii) HEURISTIC: if the values `Σ_{h∈H}χ(h+b)` on
cosets behaved like sums of independent signs, then `|D_{1/2}| ≥ |H|` with `|H|` up to
`log₂p/(1−h(¼)) ≈ 5.3·log₂p` would occur, which is why clause 1 carries the constant `8` rather than `2`; the
observed maximum `5.31` is of this size.

---

## 7. Remaining obligations (OPEN)

1. Prove or refute Conjecture SI. Even its clause 1 implies the open subgroup estimate of 6.3(i), so a proof
   requires new input; a refutation would come from a rectangle with bias `≥ ½` and `> 8log₂p` biased
   dilates — the natural candidates are cosets of subgroups of order `≈ 6log₂p` (6.3(iii)), which the present
   grids (`p ≤ 1499`, `|H| ≤ 36`) do not reach.
2. Identify an estimate that breaks the invariance of Proposition 4.6: any bound on `S(A,B)` valid below the
   square-root barrier must depend on `r_{A+B}`, not only on `{σ_r}`. The sumset compression of §5.3 is the
   only empirical signal in that direction, and it is neither necessary (Theorem 6.2) nor monotone.
3. The `k ~ log p` regime of H4: exact maxima for `k ≥ 7` at `p ≥ 101` are out of reach of enumeration; the
   annealing gap `k* < k₀(p)` should be closed or shown to be real.
4. Section 4.3 with `0 ∈ A` is stated with a crude factor `2^{2k}`; a clean statement would carry the class
   `∞` through Proposition 3.2 (the bound `N_k ≤ Σ_j C(2k,j)μ_∞^j(2k−j−1)!!R_fin^{(2k−j)/2}` over even `2k−j`,
   provable by the same matching argument, is not needed above and is not checked).

No claim in this note bears on the truth of the conjecture; the proved statements delimit what dilation-moment
arguments can and cannot do, and the data locate the near-extremal rectangles as random-like sets with an
isolated spike at `t = 1`.
