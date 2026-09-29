# Referee report `referee3` on `stepanov-sharpen` and on Proposition 2.9 of `stepanov-robust` (2026-09-29)

**Status:** PROVED (independent re-derivation in §§1–6 of this note, taking the Hankel inequality
(★) of [R, Thm 2.1] as given and, where stated, Weil's bound): Proposition 2.9 of [R] ((★)_k) and
all five claims of [Sh] = `research/stepanov-sharpen-2026-09-27.md` are **mathematically correct**.
Verdicts: Prop. 2.9 of [R] CORRECT; claim 1 (Theorem 1.5, with Lemma 1.1, Cor. 1.2, Thm 1.3,
Thm 1.4) CORRECT; claim 2 (Prop. 2.1, LP optimality) CORRECT WITH FIXES (scope qualifier, and a
continuity step in (iii)); claim 3 (Cor. 3.3) CORRECT; claim 4 (Thm 4.3, Prop. 4.1) CORRECT WITH
FIXES (the restricted threshold needs a precise definition; two data statements are inaccurate);
claim 5 (Thm 5.1, Cor. 5.2) CORRECT WITH FIXES (paper §9.5 uses `θ₁` without defining it). No fix
changes a theorem. REFUTED: nothing in the claims. There is no counterexample among 5,555,096
checks, all but about 5,200 of them (float LPs and scans) in exact arithmetic. These cover every pair `(A,B)` up to affine equivalence for `p ≤ 23` (about 1.8·10¹³ pairs),
and every `A ∋ 0` at `p ≤ 19` for (★)_k (§9). Remark (not an error): for `κ > (1+√3)/2 ≈ 1.366`
the elementary Chung bound saves more than `η★(κ)` (§2.6). OPEN: the same questions as in [Sh]
(linear saving; the balanced window).

Verifier (independent; written from scratch, the worker's verifier was not imported):
`experiments/stepanov_referee3_2026_09_29.py` → `results/stepanov_referee3_2026_09_29.json`
(`/opt/miniconda3/bin/python3`; numpy, scipy, sympy; 38 s wall; 5,555,096 checks, 0 failures).
Everything is exact integer or rational arithmetic, except the LP values of §3 (scipy/HiGHS) and
the scans explicitly marked "float".

References: [R] = `research/stepanov-robust-2026-09-26.md`; [Sh] = the note under review;
[P] = `paper/robust-hanson-petridis.md`, Section 9; [Ref2] = `research/stepanov-referee2-2026-09-27.md`.

---

## 0. Notation

`p` odd prime, `χ` Legendre (`χ(0)=0`), `d=(p−1)/2`, `m=|A|`, `n=|B|`, `r=|B∩(−A)|`,
`δ_b=[b∈−A]`, `e_b=#{a∈A: χ(a+b)=−1}`, `N_−=Σ_{b∈B}e_b`, `t̄=N_−/(mn)`, `F_A(b)=Σ_{a∈A}χ(a+b)`.
Then `S(A,B)=mn−r−2N_−`. For `0<w≤1`: `u(w)=(√(12w−3w²)−w)/2`, `η★(κ)=1−u(1/(1+2κ))`.
Affine normalisation (used by the verifier): for `c≠0` and `t`, `S(c(A−t), c(B+t))=χ(c)S(A,B)`,
with the same `m,n,r`. So every `|S|` statement, and every statement about `e_b` for both `A` and
its non-residue dilate, reduces to `A ⊇ {0,1}` (`m ≥ 2`). This reduction was validated against brute
force over all `(A,B)` at `p = 7, 11` (check `B0_normalisation_bruteforce`).

## 1. Proposition 2.9 of [R] (the order-`k` Hankel inequality): CORRECT

**Statement.** Let `k | p−1`, `d_k=(p−1)/k`, `D=d_k+m−1 ≤ p−1`, and
`e_b=#{a∈A: a+b≠0, (a+b)^{d_k}≠1}`. Then for every integer `e ≥ 0` with `2e ≤ min(m−1,d_k)`,
`Σ_{e_b≤e}(e+1−e_b)(m−(3e+e_b)/2−δ_b) ≤ (e+1)(d_k−e)`.

**Re-derivation (PROVED).** Put `F=−1+Σ_l c_l(x+a_l)^D`, with `c_l=Π_{l'≠l}(a_l−a_{l'})^{−1}`. Fix
`b`, and let `y_l=b+a_l`, `ω_l=y_l^{d_k}` (for `y_l≠0`) and `E=E(b)`. For `1 ≤ j ≤ m−1` we have
`D−j ≥ d_k`, so `y_l^{D−j}=ω_l y_l^{m−1−j}` when `y_l≠0`, and `y_l^{D−j}=0` when `y_l=0`. Put
`R=Σ_{l∈E}λ_l(x+a_l)^D` with `λ_l=c_l(1−ω_l^{−1})`. Then

`F^{(j)}(b)−R^{(j)}(b) = (D)_j Σ_{y_l≠0} c_l y_l^{m−1−j}[ω_l−(ω_l−1)[l∈E]] = (D)_j Σ_{y_l≠0} c_l y_l^{m−1−j}`.

The bracket is `1` in both cases: `ω_l=1` for `l∉E`. By the partial-fraction identity
`Σ_l c_l y_l^{i}=[i=m−1]` (`0≤i≤m−1`), the right side is `0` for `j≥1` when `b∉−A`. When
`b=−a_{l₀}` it is `0` for `j≤m−2`, since the missing term is `c_{l₀}·0^{m−1−j}`. The case `j=0` is
the same, and the constant `−1` cancels the `1`. So `ord_b(F−R) ≥ m−δ_b`, and `R` is a combination
of `e_b` pure `D`-th powers. These are the only local inputs of Steps 1–3 of [R, Thm 2.1] (rank
`≤ e_b` of `[ρ_{i+j}]` and the row expansion). Step 4: `deg F=d_k` with leading coefficient
`C(D,m−1)`, so `u_s=F^{(s)}/(D)_s` has degree `d_k−s` and leading coefficient `C(D−s,m−1)`; this
needs `s ≤ 2e ≤ d_k`. The top coefficient of `H_{e+1}` is `Λ=det[C(D−i−j,m−1)]_{0≤i,j≤e}`. In
Lemma 2.2 of [R] one has `N'−c=d_k−e ≥ 0` and `c=m−1−e ≥ e`, and every factorial argument is
`≤ D−e ≤ p−1`, so `p∤Λ`. All weights are `≥ 0` because `m−2e−1 ≥ 0`. ∎ The hypothesis `2e ≤ d_k`
is stronger than the degree argument needs (`e ≤ d_k` would do), but it is a valid sufficient condition.

**Checks.** (★)_k holds for **every** `A ∋ 0` (translation invariance) at
`(p,k) ∈ {(7,3),(7,6),(11,5),(11,10),(13,3),(13,4),(13,6),(13,12),(17,4),(17,8),(17,16),(19,3),(19,6),(19,9),(19,18)}`.
This is 2,719,703 pairs `(A,e)`, all `m ≤ p−d_k`, all admissible `e`; `k=2` controls add 9,946.
It also holds for 42,946 random and structured instances (intervals, APs, geometric progressions,
subgroups, cosets) at `p ∈ {29,…,97}`. For `k ≥ 3`, equality at `e ≥ 1` occurs only when `d_k=2`
(at `e=1`). At `e=0` it is HP equality. **Direct algebra (A3):** for 413 triples `(p,k,A,e)`
with `p ≤ 41` and `k ∈ {2,…,6}`, the verifier computes `H_{e+1}` over `F_p`. In every case
`deg H_{e+1}=(e+1)(d_k−e)` exactly, the leading coefficient is `Λ mod p ≠ 0`, and
`ord_b H_{e+1} ≥ T_b` at all 666 points with `e_b ≤ e`. The order exceeds `T_b` at 210 of them, so
the Step-3 bound is not always sharp.

## 2. Claim 1: the sharpened bias bound (Theorem 1.5 of [Sh]): CORRECT

### 2.1 Lemma 1.1 (mean form)
With `e=E−1` (so `1≤E≤(m+1)/2`), the (★)-term of `b` is
`(E−e_b)(m+3/2−(3E+e_b)/2) − δ_b(E−e_b) = W_E(e_b) − δ_b(E−e_b)`, and it is `≥ 0`. Restricting to
`b∈B` and moving the `δ`-terms right gives `Σ_{B}W_E(e_b) ≤ E(d−E+1)+Er`. **Convexity:** in
`y=E−x`, `W_E=y²/2+y(m+3/2−2E)` on `x≤E`. This is a convex quadratic in `x`. Its left slope at `E`
is `−(m+3/2−2E) ≤ 0` (as `2E ≤ m+1`), and `W_E=0` on `[E,∞)`. So the slopes are nondecreasing,
`W_E` is convex, and Jensen gives `nW_E(x̄) ≤ Σ_B W_E(e_b)`. PROVED.

### 2.2 Corollary 1.2
Put `s=E/m`, `β=3/(2m)`, `a=1−2s+β ≥ 1/(2m)` and `y=s−t̄`. Then `W_E(x̄)=m²y(a+y/2)`, and the
right side over `n` is `≤ m²sw'`. For `y>0`, `y ≤ −a+√(a²+2sw')`, hence
`1−2t̄ ≤ −a−β+2√(a²−aw'+(1+β)w')`. The identity `2s=1+β−a` gives `a²+2sw'=a²−aw'+(1+β)w'`.
For `y ≤ 0` the bound `1−2t̄ ≤ a−β` is weaker. HP line: (★) at `e=0` restricted to `B` gives
`mn₀ ≤ d+r`. PROVED.

### 2.3 Theorem 1.3 (one-variable inequality)
I re-derived (F1)–(F6) (sympy: `3v²−3vw+w²−w=0`, `P(v)=(2v−w)²`, `h₀(v)=u`, the identity for
`v−w`, `u''=−2√3/(w(4−w))^{3/2}`). I also re-derived the Jensen step:

* `√(P+βw) ≤ √P+βw/(2√P)` gives `h_m ≤ h₀(a_E)−β(1−w/√P(a_E))`.
* Taylor at the minimum with `h₀'' ≤ 2/√c` gives `h₀(a_E) ≤ u+1/(√c m²)` when `|a_E−v| ≤ 1/m`.
  Such an `E` exists because the `a_E` have spacing `2/m`, run from `1−1/(2m)` down to `≤ 3/(2m)`,
  and `v ∈ [0.4045,1]`.
* `1−w/Q ≥ Q−w` for `w ≤ Q ≤ 1` is equivalent to `(Q−w)(1−Q) ≥ 0`.
* Hence `h_m−u ≤ (3/(2m))[(1+2/(3√c))/m−(q−w)] ≤ 0` once `m ≥ M_J`.

(F6) reduces to `ε²+2ε−2 ≤ 0`. The case split is correct:

* On `[11/20,1)`, `εM_J` is increasing: `c` decreases in `ε`, and `(3−2ε)/(1−ε)` increases. It is
  `≤ 3.7570 < 4` at `ε = 9/20`.
* On `[1/4,11/20]`, `M_J ≤ 10.0845`: the first factor is maximal at `ε=3/4`. The second factor's
  only critical point in `(0,1)` is `3/2−√3/2`, a minimum, so its maximum is at an endpoint.

Exact checks (C1): `εM_J(9/20) ≤ 3.757` and `M_J ≤ 10.085` hold with rational enclosures of `√c`.
**Independent interval check (C2):** for `1 ≤ m ≤ 10` I split `[1/4,11/20]` into 40 equal
intervals (a different subdivision from [Sh]'s 64). On each, some option at `w₂` has a rational
upper bound `≤` the rational lower bound of `u(w₁)`. All 400 boxes pass at the first level, with
minimum margin 0.0447. Monotonicity in `w` of every option, and of `u`, justifies the endpoint
evaluation. Cross-checks: the same interval method for `11 ≤ m ≤ 60` on `[1/4,0.9]` (3,000 boxes);
400 exact point checks with `m` up to 20,000 and `w` up to 0.999999; a float scan (`m ≤ 3000`,
`w ≤ 0.999`) with `max(min option − u) = −3.33·10⁻⁷`. PROVED (computer-assisted in the finite
range, as [Sh] says).

### 2.4 Theorems 1.4 and 1.5
Theorem 1.4 follows from 2.2 and 2.3 and the monotonicity of the options in `w'`. The `−S` case
follows by the non-residue dilation. In Theorem 1.5 I checked every step:

* `m₀ ≤ m₁ ≤ √(2p)+1 ≤ (p+1)/2` for `p ≥ 11` (`√22+1 = 5.69 ≤ 6`).
* `r(A',B) ≤ m₀`, so `w'(A',B) ≤ w₁`, and `w₁ ≥ w ≥ 1/4`.
* The subset average is exact: `S(A,B) = (m/m₀)·avg_{A'} S(A',B)`.
* Second inequality: `u' ≤ u'(1/4)=1.06525`. Also
  `w₁−w=(2m₁−1)/((1+2κ)p) ≤ 2(√(2p)+1)/p`. The case `w₁>1` uses
  `u(w)+u'(1/4)(1−w) ≥ u(1)=1`.

PROVED.

### 2.5 Exhaustive checks
For `p ∈ {7,11,13,17,19,23}`, every `A ⊇ {0,1}` (and `A={0}`) and every triple `(m,n,r)`, the
verifier computes the exact maximum and minimum of `S` over all `B` with `|B|=n` and
`|B∩−A|=r`: take the top or bottom `r` values of `F_A` on `−A` and the top or bottom `n−r` off it.
Since every bound depends on `(m,n,r)` only, this checks **all pairs** `(A,B)`: 1.77·10¹³ pairs up to
affine maps. Results, all exact and with no failures:

* Lemma 1.1: 979 `(m,E,r)` maxima, for `A` and its non-residue dilate.
* Corollary 1.2, every `E` plus the HP line: 6,984.
* Theorem 1.4: 6,446. The largest `S/(u(max(w',1/4))mn)` is 0.52, 0.70, 0.81, 0.81, 0.85, 0.84 for
  `p = 7,…,23`.
* Theorem 1.5: 4,959, at every critical `κ` (the right ends of the intervals of constant `m₁`,
  where the bound is tightest); 4,429 of them non-trivial (`w₁<1`). Worst ratio 0.81 (`p=23`, 4×4).
* Structured and adversarial sets at `p ∈ {29,…,197}`, with the exact worst `B` for every `(n,r)`:
  intervals, APs, QR-subsets, greedy sum-cliques and random sets; 1.79·10⁶ checks. The largest
  Theorem-1.4 ratio is 0.980 (`p=197`).

The exhaustive data quoted in [Sh] §3 are reproduced exactly: at `p=23` the maximal bias is
`0.917, 0.8125, 0.75, 0.64, 0.528`, and `0.833` at `p=19`, `κ=0.05`.

### 2.6 Remark: comparison with the second moment (not an error)
Chung's bound (brief, background) gives `|S| ≤ √(mn(p−n)) < √(p/(mn))·mn`. So for
`mn ≥ (1/2+κ)p` the saving is at least `1−√(2w)`. Now `u(w)=√(2w)` exactly at `w=2−√3`: there
`12w−3w²=3`, so `u=√3−1=√(4−2√3)`. Hence for `κ > (1+√3)/2 ≈ 1.366` the elementary bound is
better than Theorem 1.5. At `κ=3/2` it gives `1−1/√2 = 0.2929 > η★ = 0.2865` (exact enclosure,
check H3). The table rows `κ=1.5` in [Sh] §1 and [P] §9.1 are correct as values of `η★`, but they
are not the best known saving there. The same applies to the sentence "for `κ>3/2` the case `κ=3/2`
applies" in [R, Thm 2.5] / [P, Thm D], which the second moment supersedes. That is outside this
referee's scope.

## 3. Claim 2: optimality within (★) (Proposition 2.1): CORRECT WITH FIXES

**(i)** The LP relaxation keeps only `b∈B`. Dropping the nonnegative terms of `b∉B` is without loss,
because (★) alone is consistent with `e_b = m` off `B`. Lemma 1.1 and Corollary 1.2 use only the
LP constraints, with `Σn'_j ≤ min(m,n)` standing for `r`. So `V/(mn) ≤ u(max(w'',1/4))`. PROVED.
Float LPs at `p ∈ {101,1009}` (182 cases) give `max V/(mn)/u = 0.9915`.

**(ii)** I re-derived every step:

* `∂_tK=−(1−s−t)`.
* `τ(s)` is the root of `K(s,·)=sw`.
* `1−2τ(s)=h₀(1−2s)`, so `max τ = t* = (1−u)/2` by (F1).
* `K(s,t*) ≤ sw` in the three cases `τ(s) ≥ 0`, `τ(s) < 0`, and `s ∈ (1/2, 1/2+1/(2m)]`; the last
  uses `w ≥ 1/4`.
* `K(s,t*)−K(s,t) ≥ (t−t*)(1−s−t)`.
* `E(d−E+1) ≥ smd−sm²/2`. For `s=1/2+1/(2m)` this uses `s²m² ≤ sm²/2+sm/2`.

PROVED.

**(iii)** The constants (`t* ≤ 0.1432`, `δ=5(3/(2m)+m/(2d))`, `0.2685 ≥ 0.2`) are right. The
minimum of the concave product over `t ∈ [t*+δ, t*+δ+1/m]` is at an endpoint.

**Checks.**

* Exact (D1): for 8 triples `(w,m,d)` I check the all-`e_b=J` profile against every (★)_e in
  integers. The note's `J=⌈(t*+δ)m⌉` is feasible in all cases. The smallest feasible `J` is much
  closer to `t*m`: at `w=9/10`, `m=4000` it gives bias 0.9965 against `u=0.99655`; at `w=1/4`,
  `m=400`, 0.710 against 0.7135.
* Float LP with `d=∞` (D2): `m(u−V) ≈ 0.52, 0.20, 0.10, 0.05` at `w = 0.5, 0.8, 0.9, 0.95`
  (`m` up to 400), as [Sh] states. So `V = u−Θ(1/m)`.

**Fixes.**

1. *(Continuity.)* The profile of (ii) has `n=⌊d/(mw)⌋`, hence `mn ≤ (1/2+κ)(p−1) < (1/2+κ)p`. It
   sits just **below** the threshold of Theorem 1.5. To state (iii) for pairs with
   `mn ≥ (1/2+κ)p`, apply (ii) with some `w' < w`, `w'→w` (for example `w'=w(1−p^{−1/2})`), and use
   the continuity of `u`. This is a one-line fix and the conclusion is unchanged.
2. *(Scope.)* The optimality statement holds for (★) applied to the given `A`. By symmetry it also
   holds for `B`, via the doubly regular count profile, which (ii) handles on each side. It does
   **not** cover (★) applied to subsets `A' ⊆ A` or `B' ⊆ B`; [Sh] labels that part HEURISTIC.
   [Sh]'s status line ("exactly the best saving that the family (★) can certify") and [P] §9.2
   ("no non-negative combination of the inequalities (★) …") should say "(★) for `A` (or `B`)".
   The refutation of "combining (★) over several `e` gives a linear saving" is correct as worded
   (fixed `A`, all `e`).

## 4. Claim 3: upper bounds on the true saving (Cor. 3.3): CORRECT

Lemma 3.1 re-derived: the expansion of `[e_b=j]` over `I ⊆ A` gives
`Σ_{I≠∅}|T(I)| ≤ Σ_{k≥1}C(m,k)[(k−1)√p+m−k] ≤ m2^{m−1}(√p+1)`. Proposition 3.2 then follows by
choosing `B` greedily; the same Weil counts hold for every `A` of size `m`. The `β`-formulas were
checked exactly:

* `β₂=w` on `[1/3,1]`, `β₅=3/5+w/8` on `[8/15,1]` and `β₁=2w−1` on `[1/2,1]` (303 checks).
* The crossover `3/5+w/8=w ⇔ w=24/35 ⇔ κ=11/48`.
* `β_m(1/(2mw)) ≤ w` for all `3 ≤ m < 200` at all 518 breakpoints in `[24/35,1]`. `β_m − w` is
  piecewise linear in `w`: on each piece `β = 2mw·(α+γq)` with `q = 1/(2mw)`. The maximum of
  `β_m−w` is `0`, attained only at `m=5`, `w=24/35`.
* Hoeffding for `m ≥ 200`: `β ≤ 1/2+2me^{−m/8}`.

Table: sup over `m` of `β_m`, computed exactly for `m ≤ 400`. The tail `m > 400` uses Hoeffding
with threshold `0.3`, because [Sh]'s threshold `1/2` does not bound the tail at `κ = 3/2`, where
the sup is 0.458. The resulting savings are
`κ=1/4: 19/60 (m=5)`, `κ=1/2: 29/80 (m=5)`, `κ=1: 107/224 = 0.4777 (m=7)`,
`κ=3/2: 13/24 = 0.5417 (m=6)`, and `2κ/(1+2κ)` (`m=2`) for `κ ∈ {0.01,0.05,0.1}`. All match [Sh].
Realisation at `p=100003` (E4): random `A` with `m ∈ {2,5,6,7}` and greedy `B` give `S/(mn)`
within `2·10⁻³` of `β_m`. PROVED given Weil.

## 5. Claim 4: balanced sets (Prop. 4.1, Lemma 4.2, Thm 4.3): CORRECT WITH FIXES

Re-derived.

* **Lemma 4.2:** tuples in `A⁴` split by `deg R ∈ {0,2,4}`. There are `≤ 3m²` even tuples, each
  contributing `≤ p`. If `deg R = 2` the sum is `−1−χ(R(−x)) ∈ [−2,0]`. If `deg R = 4`, Weil gives
  `≤ 3√p`.
* **Theorem 4.3(a):** `S ≤ (m−2)n+2N₀`, with `m/2^m ≤ k/2^k` for `m ≥ k ≥ 2`. Also
  `(√p+5)/(2p) ≤ 1/√p` for `p ≥ 25`, and `M²/(τ√p) ≤ κ/(2τ)`. The saving is `≥ κ/(τm) ≥ κ/(τM)`.
* **Theorem 4.3(b):** `(3/(τm)+3√p/n)^{1/4} ≤ (1/4+1/4)^{1/4} = 0.8409 < 0.841`.

The two regions and the window `W={M<|A|≤|B|<12√p}` partition all pairs with `k ≤ |A| ≤ |B|`.
`|A| ≤ |B|` is without loss of generality. `M(k,κ)=⌈12/τ⌉` with `τ=k/2^k+κ` is stated identically
in [Sh] and [P]. Sharpness: the witnesses of Prop. 4.1 have `|A|=k ≤ M`, so they lie outside `W`.

**Checks (exact):**

* Lemma 3.1: 390 checks at `p ∈ {10007, 100003}`; max deviation/bound = 0.135.
* Lemma 4.2: 48 checks, with intervals, squares, APs and random sets up to `m=80`; max ratio 0.71.
* Theorem 4.3 at primes that satisfy its hypotheses, namely `p = 16411, 26249, 40009, 614657`,
  for `(k,κ) = (2,1), (3,1), (4,1), (3,1/2)`: 805 checks of (a) with the worst `B`, max ratio 0.70;
  138 checks of (b), max bias 0.52.
* Prop. 4.1 witnesses at `p=100003`: `k=2,…,6` give `mn/p ≈ k/2^k` and `S/(mn) ≥ 0.99997`.
* The LP-consistency of complete bicliques (`n₀=d/m` satisfies every (★)_e for `m ≤ 3d/2`): 897
  exact checks.

**Fixes.**

3. The "exact threshold outside the window" needs a definition. `W` depends on `κ` through `M`, so
   `τ_k` restricted to `W^c` should be defined as `τ_k^{out}` = inf of `τ` such that for every
   `κ>0` there are `c,p₀` with `|S| ≤ (1−c)|A||B|` for all `p ≥ p₀` and all pairs with
   `min ≥ k`, `|A||B| ≥ (τ+κ)p` and `(A,B) ∉ W(k,κ)`. With this definition Thm 4.3 and Prop. 4.1
   give `τ_k^{out}=k/2^k`, as claimed. The saving is `min(κ/(τM), 0.159)`, under `k ≥ 2`,
   `0<κ≤1` and `p ≥ max(25,(2M²/κ)²)`. [P] §9.4 states the savings but omits these hypotheses.
4. *(Data statements, [Sh] §§3–4.)*
   * "the extremisers are 3×3, 3×4 and 4×4 rectangles" is imprecise. At `p=23` the listed maxima
     are attained at 4×4, 4×5, 5×5 and 6×6. At `p=19` they are attained at 3×4, 3×5, 4×4, 3×6, 4×5
     and 5×6. Ties were not enumerated.
   * "largest complete biclique with `min ≥ 3` has `mn/p = 0.52, 0.47, 0.53, 0.69`" counts
     zero entries (`b ∈ −A`), i.e. `B ⊆ {e_b=0}` with `S = mn−r`. Without zeros (a genuine
     `K_{m,n}` in the graph `χ(a+b)=1`) the values are `0.39` (3×3, `p=23`) and `0.47` (`p=19`),
     and there is none at `p = 13, 17`. I reproduced `M_m(23) = [12,6,4,3,2,2,1,1]` exactly.

## 6. Claim 5: characters of order `k` (Thm 5.1, Cor. 5.2): CORRECT WITH FIXES

Re-derived. The dilation by `g^{−1}` with `ψ(g)=ω` maps ω-bad pairs to the `e_{b'}` of
Proposition 2.9. `r` is unchanged. `m ≤ d_k+1` gives `2e ≤ m−1 ≤ d_k` and `m+d_k−1 ≤ 2d_k ≤ p−1`.
Lemma 1.1, Corollary 1.2 and Theorem 1.3 go through verbatim with `d_k`, so
`#bad_ω ≥ θmn` and `N_ω = M−#bad_ω ≤ (1−θ)M`.

**Vertex lemma.** On `{N ≥ 0, ΣN=M, N_ω ≤ (1−θ)M}` with `θ ≤ 1/2`, a vertex has at most one
fractional coordinate, and two coordinates at `(1−θ)M` would exceed `M` unless `θ = 1/2`. So the
vertices are `(1−θ)M` and `θM` on two roots, and the maximum of `|ΣN_ωω|` is `M|1−θ+θζ|`. The
identity `|1−θ+θζ|²=1−2θ(1−θ)(1−cos 2π/k)` holds.

**Cor. 5.2.** `p ≥ (1+λ)k+1` gives `d_k ≥ 1+λ`, hence `m₁ ≤ d_k+1`. Also `λ ≤ 3` gives
`w₁ ≥ 1/4`. The bound is decreasing in `θ ∈ [0,1/2]`. Series: `θ(λ)=λ²/6+O(λ³)` (sympy).
Sharpness: `A={0}`, `B=H_k` gives `S_ψ=|A||B|=d_k`, and adding `λd_k` points of `gH_k` gives bias
`|1+λζ|/(1+λ)` (exact at `p ∈ {1009,10009}`, `k ∈ {3,4,6}`; this needs `λ ≤ 1`).

**Checks.**

* Vertex lemma: 140 exact cases.
* Theorem 5.1 exhaustively over **all** `(A,B)` at `p=7` (`k=3,6`) and `p=13` (`k=3,4,6`): exact
  `|S_ψ|²` in `Z` and exact comparison with the irrational bound. This is about 1.0·10⁸ `(A,B)` pairs,
  giving 482 `(m,n,r)` maxima for the modulus bound and 1,341 worst-`B` bad-pair counts.
  Every `A ∋ 0` with `m ≤ d_k+1` is included at `p = 19` (`k = 3, 6, 9`).
* Cor. 5.2 at every critical `λ`: 844 exhaustive checks. With random `A` and direction-greedy `B`
  at `p ∈ {31,37,61,73}`: 130,512 Cor. 5.2 checks and 130,194 Theorem 5.1 modulus checks.
* 699k exact worst-`B` bad-pair counts.
* Largest `|S_ψ|/bound`: 0.91 (`p=13`, `k=6`).

**Fixes.**

5. [P] §9.5 uses `θ₁` without defining it. Insert `m₁=⌈√((1+λ)d_k)⌉`, `w₁=(d_k+m₁)/((1+λ)d_k)`,
   `θ₁=(1−U(min(1,w₁)))/2` as in [Sh, Cor. 5.2].
6. [Sh] Cor. 5.2 says "For `k=2` this is Theorem 1.5 with `λ=2κ`". This holds only up to
   normalisation: the thresholds are `(1/2+κ)(p−1)` against `(1/2+κ)p`, and `w₁` differs
   (`d` against `(1/2+κ)p` in the denominator). Say "the same statement up to `O(1/p)` in the
   threshold".

## 7. Transcription check of [P] Section 9

* §9.1: identical to [Sh, Thm 1.5] (`m₁`, `w₁`, `p ≥ 11`, `0<κ≤3/2`, the constant 2.14, the series
  `4/3κ²−40/9κ³`). The quoted numbers (`9.838e−3` vs `7.591e−3`, `0.10436` vs `0.08579`,
  `0.28647` vs `0.25`) are correct, and the ratio → 4/3 (sympy). See §2.6 for `κ = 3/2`.
* §9.2: correct up to fix 2 (scope qualifier, "for `A`").
* §9.3: correct (`29/80`, `11/48`).
* §9.4: correct up to fix 3 (hypotheses of Thm 4.3).
* §9.5: correct up to fix 5 (`θ₁` undefined). The (★)_k hypotheses are transcribed correctly.
* §9.6 is from [R] §1.2 and §4 and was not in my brief; I did not referee it.

## 8. Consolidated fix list (none changes a theorem)

1. [Sh] Prop. 2.1(iii): apply (ii) with `w'<w`, `w'→w`, to land in `mn ≥ (1/2+κ)p`.
2. [Sh] status (2), §2, and [P] §9.2: restrict the optimality claim to (★) for `A` (and `B`).
   Subsets `A' ⊆ A` are only HEURISTIC.
3. [Sh] §4 and [P] §9.4: define the restricted threshold `τ_k^{out}` and state the hypotheses
   `k ≥ 2`, `0<κ≤1`, `p ≥ max(25,(2M²/κ)²)`.
4. [Sh] §3 and §4 data sentences: extremiser shapes at `p=23`, and "complete biclique" includes
   zero entries.
5. [P] §9.5: define `m₁`, `w₁`, `θ₁`.
6. [Sh] Cor. 5.2: the `k=2` comparison is "up to `O(1/p)`".
7. (Recommended, §2.6.) Record that for `κ > (1+√3)/2` the Chung bound beats `η★`, and that the
   best known saving is `max(η★(κ), 1−(1/2+κ)^{−1/2})`.

## 9. Checks performed (all exact unless marked)

| section | content | checks |
|---|---|---|
| A1 | (★)_k, every `A ∋ 0`, 15 `(p,k)` pairs plus `k=2` controls | 2,729,649 |
| A2 | (★)_k, random and structured `A`, `p ≤ 97` | 42,946 |
| A3 | `H_{e+1}` over `F_p`: degree, `Λ`, orders | 1,079 |
| B0–B5 | Lemma 1.1, Cor. 1.2, Thm 1.4, Thm 1.5 on all `(A,B)`, `p ≤ 23`; normalisation vs brute force | 24,377 |
| B6 | the same, structured and adversarial `A`, worst `B`, `p ≤ 197` | 1,790,242 |
| C | Thm 1.3: constants, interval boxes (400 + 3,000), 400 exact points, float scan | 672 |
| D | LP profiles (exact), LP values (float) | 218 |
| E | `β_m` (exact), table, realisation | 837 |
| F | Lemma 3.1, Lemma 4.2, Thm 4.3, Prop. 4.1, biclique LP | 2,283 |
| G | vertex lemma, Thm 5.1, Cor. 5.2, sharpness | 962,775 |
| H | quoted numbers, Chung crossing | 18 |
| **total** | | **5,555,096, 0 failures** |

## 10. Remaining obligations

* The whole of [Sh] remains conditional on (★) ([R, Thm 2.1], refereed in [Ref2] and Lean-checked
  per [P] §10.2). Theorem 1.3's finite range is computer-assisted. It has now been verified by two
  independent rational-interval implementations.
* OPEN, unchanged: whether the true saving is linear in `κ`; any threshold `< 1/2` inside the
  balanced window; order-`k` Weil-range thresholds.
