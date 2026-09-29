# Stepanov wave, worker `padic`: lifting Hanson–Petridis to ℤ/p² and ℤ_p (2026-09-29)

**Status:** Verdict: a precise obstruction, not a new constraint. The lift to `ℤ/p²` or `ℤ_p`
adds no usable input for complete bicliques `A + B ⊆ Q ∪ {0}`.
PROVED (complete elementary proofs below, all checked exactly):
(i) the refined Euler criterion `x^{d+j(p−1)} ≡ χ(x)(1 + p(j+½)q(x)) (mod p²)`, the
Fermat-quotient rules, and an exact integer cocycle for `q` of a sum. It gives two
decompositions of `q(ã+b̃)`. One is a Witt/Teichmüller part `£₁(a/(a+b))` plus a separable term
`(u(a)+u(b))/(a+b)`. The other is a residue part `φ(a+b)` minus a carry term `ε(a,b)/(a+b)`
(Cor. 1.3).
(ii) Exact `mod p²` Taylor data of the lifted Hanson–Petridis polynomial `F̃_A` at every point
(Prop. 2.2). At complete points they are Fermat-quotient combinations. For `j ≤ m−2` they do not depend
on the integer lift of `b`, and their dependence on the lift of `A` is a fixed polynomial of
degree `≤ m−2` (Lemma 2.3, Prop. 2.4).
(iii) An exact criterion for extra `mod p²` vanishing at complete points: all sums `b̃ + ã_k`
must lie in `ker q`. For `A = {0,1}`, `b = 1`, this holds iff `p` is a Wieferich prime
(Prop. 2.5). A second family of extra vanishing comes from the Cauchy–Mirimanoff factor
`(x²+x+1)²`.
(iv) **Root counting over `ℤ_p` is exactly root counting over `F_p`.** For every `P̃ ∈ ℤ_p[x]`,
the number of `p`-adic roots in the residue disc of `b` equals `ord_b(P̃ mod p)` (Prop. 2.6).
So no `ℤ_p` or `ℤ/p²` multiplicity count can beat `mn ≤ d + r` with the same polynomial.
Complete points are generically Eisenstein clusters of `m` conjugate roots at distance
`p^{−1/m}` (Prop. 2.7).
(v) Completeness imposes no condition on the second `p`-adic digits (Prop. 2.8).
(vi) Order-`m` vanishing `mod p²` at all complete points is available, but only for the
Frobenius-exponent polynomial `F̃_*`, whose reduction has degree `pd` (Prop. 2.9).
(vii) `Σ_{a,b}(ã+b̃)^d ≡ S + (p/2)T (mod p²)` with the carries handled exactly. `S` is
recovered from its residue and parity when `mn − r < p`. `T` splits into Witt, separable,
residue and carry parts. Its translation average is fixed by `mod p` data (Props. 3.1–3.4).
(viii) The Witt part has a verbatim analogue in `W₂(F_{p²})` up to a Frobenius twist
(Prop. 1.5). The only lift ingredient with no `F_{p²}` analogue is the canonical section
`[0,p) ⊂ ℤ`, and it enters only through a separable coboundary (Cor. 1.3(3)).
REFUTED:
- "complete points are roots of `F̃_A` modulo `p²`". Witness: `p = 17`, `A = {0,1}`: the
  second digits at the three complete points `1, 8, 15` are `13, 2, 13 ≢ 0`. Only 8 of 1,139
  tested complete points vanish `mod p²`.
- "the `mod p²` data are functions of the configuration's profile". Witness: `p = 101`,
  `A = {0,1,4}`, `B = B(A)` has no complete point with vanishing second digit. Its image
  `(4A + 85, 4B − 85)` has the same profile and one such point, `b = 21`. `T` takes between
  5 and `p` values on every one of 29 tested symmetry orbits.
HEURISTIC: on complete bicliques the Fermat quotients of the sums behave like those of random
rectangles (ranks, zero counts, exponential sums). The second-digit function has full degree
`p − 1` in 25/25 cases.
OPEN: nothing here narrows the balanced Hanson–Petridis question. The LP-tight profiles of
`research/stepanov-balanced-2026-09-29.md` are untouched, because the lift yields no
count-level (profile) constraint (§4).

Verifier: `experiments/stepanov_padic_2026_09_29.py` → `results/stepanov_padic_2026_09_29.json`
(run with `/opt/miniconda3/bin/python3`; standard library + numpy + sympy; about 8 s;
**1,314,056 exact checks, 0 failures**).

---

## 0. Notation

`p` is an odd prime (`p ≥ 5` where stated) and `d = (p−1)/2`. `χ` is the Legendre symbol and
`Q` the nonzero squares. For `x ∈ F_p`, `x̃ ∈ [0,p)` is its **canonical** integer
representative; sets `A, B ⊆ F_p` are lifted canonically unless said otherwise. For an integer
`x` with `p ∤ x`, the Fermat quotient is `q(x) := (x^{p−1} − 1)/p mod p`. We write
`f(x) := (x^p − x)/p ∈ ℤ` for every integer `x`, so `f(x) = x·q(x)`. On `F_p` put

- `u(x) := f(x̃) mod p` (so `u(0) = 0`),
- `φ(z) := q(z̃)` for `z ≠ 0` (so `u(z) = zφ(z)`),
- `ε(a,b) := [ã + b̃ ≥ p] ∈ {0,1}` (the **carry**),
- `£₁(s) := Σ_{k=1}^{p−1} s^k/k ∈ F_p[s]` (truncated logarithm),
- `Γ(X,Y) := Σ_{k=1}^{p−1} (−1)^{k−1} X^k Y^{p−k}/k ∈ F_p[X,Y]`.

`ω` is the Teichmüller map. Modulo `p²`, `ω(x) ≡ x̃^p`. `P^{[j]} = P^{(j)}/j!` is the Hasse
derivative. For `A = {a_1,…,a_m}` with `2 ≤ m ≤ (p+1)/2`: `D = d + m − 1 ≤ p − 1` and
`c_k = Π_{l≠k}(a_k − a_l)^{−1}`. As in the robust note, `F_A(x) = −1 + Σ_k c_k(x + a_k)^D`, with
`deg F_A = d`. For `b ∈ F_p`: `y_k = b + a_k`. `B₁ = B ∖ (−A)`, `B₀ = B ∩ (−A)`, `r = |B₀|`.
`B(A)` is the full partner set `{b : b + A ⊆ Q ∪ {0}}`.
`S = S(A,B) = Σ_{a,b} χ(a+b)`.

## 1. Exact congruences

**Lemma 1.1 (Fermat quotients; PROVED).** Let `x, y, k` be integers with `p ∤ xy`.

(a) `q(xy) ≡ q(x) + q(y)` and `q(x + pk) ≡ q(x) − k/x (mod p)`. In particular `q` depends
only on `x mod p²`.

(b) `q : (ℤ/p²)^* → F_p` is a surjective homomorphism. Its kernel is
`μ_{p−1} = {x^p mod p²}`, and `ω(x) ≡ x^p ≡ x(1 + p q(x)) (mod p²)`.

(c) **Refined Euler criterion.** For every `j ≥ 0`,
`x^{d + j(p−1)} ≡ χ(x)·(1 + p(j + ½)q(x)) (mod p²)`. In particular
`x^d ≡ χ(x)(1 + p q(x)/2)` and `x^{pd} ≡ χ(x) (mod p²)`.

*Proof.* (a) `(xy)^{p−1} = (1+pq(x))(1+pq(y)) ≡ 1 + p(q(x)+q(y))`. Also
`(x+pk)^{p−1} ≡ x^{p−1} + (p−1)x^{p−2}pk ≡ x^{p−1} − pk x^{p−2} (mod p²)`, and
`x^{p−2} ≡ 1/x (mod p)`.

(b) The homomorphism property is (a). Surjectivity: `q(1 + pk) ≡ −k`. Next,
`q(x^p) = p q(x) ≡ 0`, and the `p − 1` classes `x^p` (for `x = 1,…,p−1`) are distinct mod `p`.
So they fill the kernel, which has order `p(p−1)/p`. Finally `x^p = x·x^{p−1} = x(1+pq(x))`.

(c) Write `x^d ≡ χ(x) + pt`. Squaring gives `1 + pq ≡ x^{p−1} ≡ 1 + 2χ(x)pt`, so
`t ≡ χ(x)q/2`. Then multiply by `x^{j(p−1)} = (1+pq)^j ≡ 1 + jpq`. For `j = d`,
`j + ½ = p/2 ≡ 0`. ∎

**Lemma 1.2 (exact cocycle; PROVED).** For all integers `x, y`,

    f(x + y) = f(x) + f(y) + C_p(x, y),     C_p(X,Y) := Σ_{k=1}^{p−1} (1/p)·C(p,k)·X^k Y^{p−k} ∈ ℤ[X,Y],

and `C_p ≡ Γ (mod p)` coefficientwise. For `x, y ∈ F_p` with `x + y ≠ 0`:
`Γ(x,y) = (x+y)·£₁(x/(x+y))`.

*Proof.* `f(x+y) − f(x) − f(y) = ((x+y)^p − x^p − y^p)/p`. Also
`(1/p)C(p,k) = (p−1)⋯(p−k+1)/k! ≡ (−1)^{k−1}/k`.

For `s ∈ F_p` take any integer lift `s̃`. Expanding `(1−s̃)^p`,
`C_p(s̃, 1−s̃) = (1 − s̃^p − (1−s̃)^p)/p = −Σ_{k=1}^{p−1}(C(p,k)/p)(−1)^k s̃^k ≡ £₁(s)`.
`Γ` is homogeneous of degree `p`, and `λ^p = λ` in `F_p`. With `λ = x+y` and `s = x/λ`:
`Γ(x,y) = λ^pΓ(s,1−s) = λ£₁(s)`. ∎

**Corollary 1.3 (the Fermat quotient of a sum; PROVED).** Let `a, b ∈ F_p` with `a + b ≠ 0`.
Let `ã + b̃ ∈ [1, 2p−2]` be the integer sum of the canonical representatives. Then modulo `p`:

1. *(Witt form)* `q(ã+b̃) ≡ £₁(a/(a+b)) + (u(a) + u(b))/(a+b)`;
2. *(carry form)* `q(ã+b̃) ≡ φ(a+b) − ε(a,b)/(a+b)`;
3. *(the carry is a coboundary plus Γ)* `ε(a,b) ≡ u(a+b) − u(a) − u(b) − Γ(a,b)`;
4. *(Teichmüller)* `q(ω(a) + ω(b)) ≡ £₁(a/(a+b))`.

*Proof.* (1): Lemma 1.2 with `x = ã`, `y = b̃`, divided by `a+b`. (2): Lemma 1.1(a) with
`ã + b̃ = (a+b)~ + εp`. (3): multiply (1) and (2) by `a+b` and use `u(z) = zφ(z)`.
(4): Lemma 1.2 with integer representatives of `ω(a)` and `ω(b)`, which have `q = 0`. ∎

**Lemma 1.4 (facts about `£₁` and `u`; PROVED).**
- `£₁(s) = £₁(1−s)`.
- `£₁(1/s) = −£₁(s)/s`.
- `Σ_{s∈F_p} £₁(s) = 1`.
- `£₁(s) ≡ −(s̃^p + (1−s̃)^p − 1)/p`.
- `Σ_{x∈F_p} u(x) = ½`.
- `u(x+1) − u(x) = −£₁(−x)` for `0 ≤ x ≤ p−2`.
- As functions `F_p → F_p`, both `u` and `£₁` have **degree exactly `p − 1`**.

*Proof.*
- The first is the symmetry of `C_p`.
- The second: `x^{−k} = x^{p−1−k}` for `x ≠ 0`; reindex.
- The third: `Σ_s s^k = −[p−1 | k]`, and the `k = p−1` coefficient is `1/(p−1) = −1`.
- The fourth is the computation in Lemma 1.2.
- `Σ_{x=0}^{p−1} x^p ≡ 0 (mod p²)`: pair `x` with `p − x`, using `(p−x)^p ≡ −x^p (mod p²)`.
  Hence `Σ u = (0 − p(p−1)/2)/p ≡ ½`.
- The difference formula is `C_p(x,1) ≡ Γ(x,1) = −£₁(−x)`.
- The coefficient of `x^{p−1}` in the interpolation polynomial of any `g : F_p → F_p` is
  `−Σ_x g(x)`, which is `−½ ≠ 0` for `u` and `−1 ≠ 0` for `£₁`. ∎

**Proposition 1.5 (the Witt part is field-agnostic; PROVED).** Let `q' = p^f` and let
`R = W₂(F_{q'})` be the Galois ring of order `q'²`, with `ω(x) ≡ x̂^{q'} (mod p²)` for any lift
`x̂`. For `z ∈ R^*` put `q_{q'}(z) := (z^{q'−1} − 1)/p mod p ∈ F_{q'}`. Then
`z^{(q'−1)/2} ≡ χ_{q'}(z̄)(1 + p q_{q'}(z)/2)` and, for `a + b ≠ 0`,
`q_{q'}(ω(a) + ω(b)) ≡ Γ(a,b)^{1/p}/(a+b)`.

*Proof.* The first statement is Lemma 1.1(c) verbatim in `R`. For the second, write
`c(a,b) = ω(a) + ω(b) = ω(a+b)(1 + p w(a,b))`, with `w(a,b) ∈ F_{q'}`, modulo `p²`.

- Raise to the `p`-th power. On one side, `c^p ≡ ω((a+b)^p) = ω(a^p + b^p)`, because
  `(1+pw)^p ≡ 1` modulo `p²`.
- On the other side, `c^p = ω(a^p) + ω(b^p) + p·C_p(ω(a),ω(b))`, and the first two terms
  equal `ω(a^p + b^p)(1 + p w(a^p,b^p))`.
- Comparing the two gives `w(a^p, b^p)(a+b)^p ≡ −Γ(a,b)`, so `w(a,b) = −Γ(a,b)^{1/p}/(a+b)`.
- Finally `(1+pw)^{q'−1} ≡ 1 − pw`. ∎

For `f = 1` this is Cor. 1.3(4). The verifier checks it in `W₂(F_{p²})` for `p ≤ 11` (§H). So
the Teichmüller half of every computation below has a verbatim `F_{p²}` analogue, and it
cannot be the prime-specific input asked for (barrier B4 of `research/sigma-barriers-2026-09-05.md`).
The canonical section `x ↦ x̃ ∈ [0,p)` has no analogue. By Cor. 1.3(3) it enters only through
the separable term `u(a) + u(b)` (equivalently, the carry bits).

## 2. The Hanson–Petridis polynomial over ℤ_p and ℤ/p²

**Setup.** Take any integer lifts `ã_k` (canonical unless stated) and
`c̃_k := Π_{l≠k}(ã_k − ã_l)^{−1} ∈ ℤ_(p)^*`. Put `F̃(x) := −1 + Σ_k c̃_k(x + ã_k)^D ∈ ℤ_(p)[x]`.

**Lemma 2.1 (PROVED).** Over `ℚ`, `Σ_k c̃_k (z + ã_k)^l = [l = m−1]` for `0 ≤ l ≤ m−1` and every
`z`. So `F̃` has degree exactly `d`, with leading coefficient `C(D, m−1)`, a `p`-adic unit.
Its reduction is `F_A`.

*Proof.* `Σ_k g(ã_k)c̃_k` is the `x^{m−1}`-coefficient of the Lagrange interpolant of `g` on
`{ã_k}`. For `g(s) = (z+s)^l` with `l ≤ m−1`, that coefficient is `[l = m−1]`. The
coefficient of `x^i` in `F̃ + 1` is `C(D,i)Σ_k c̃_k ã_k^{D−i}`. It vanishes for
`D − i ≤ m−2` and equals `C(D,m−1)` for `i = d`. Finally `D < p`. ∎

**Proposition 2.2 (Taylor data mod p²; PROVED).** Let `z` be an integer with
`p ∤ z + ã_k` for all `k`, and `0 ≤ j ≤ m−1`. With `y_k = z + ã_k`,

    F̃^{[j]}(z) ≡ C(D,j) Σ_k c̃_k χ(y_k)(1 + p q(y_k)/2) y_k^{m−1−j} − [j = 0]      (mod p²).

If `b = z mod p` is complete (`b + A ⊆ Q`), then

    F̃^{[j]}(z) ≡ (p/2)·C(D,j)·Σ_k c_k q(y_k) y_k^{m−1−j}      (mod p²).

If `b = −a_{k₀}` with `b + A ⊆ Q ∪ {0}`, `z = −ã_{k₀}` and `0 ≤ j ≤ m−2`, the same formula holds
with the sum over `k ≠ k₀`.

*Proof.* `F̃^{[j]}(z) = C(D,j)Σ_k c̃_k y_k^{D−j} − [j=0]`, and
`y^{D−j} = y^d·y^{m−1−j}`. Apply Lemma 1.1(c).

- *Complete case.* Here `Σ_k c̃_k y_k^{m−1−j} = [j=0]` by Lemma 2.1, and `C(D,0) = 1` cancels
  the `−1`.
- *`B₀` case.* `y_{k₀} = 0` kills that term, since `D − j ≥ 1`. For `j ≤ m−2`,
  `Σ_{k≠k₀}c̃_k y_k^{m−1−j} = Σ_k c̃_k y_k^{m−1−j} = [j=0]`. ∎

**Lemma 2.3 (which digits are intrinsic; PROVED).** Suppose `P̃ ∈ ℤ_(p)[x]` and
`ord_b(P̃ mod p) ≥ μ`. Then `P̃^{[j]}(z + pt) ≡ P̃^{[j]}(z) (mod p²)` for all `t` and all
`j ≤ μ − 2`. Hence the second digits of `F̃` at complete points `b ∈ B₁` are independent of the
lift of `b` for `j ≤ m−2`, and at `B₀` for `j ≤ m−3`. The last digit is not: at `b ∈ B₁`,
`F̃^{[m−1]}(z+pt) ≡ F̃^{[m−1]}(z) + pt·m·F̃^{[m]}(z)`, and `F̃^{[m]}(z)` is a unit
(balanced note, Theorem 1.1(1)).

*Proof.* `P̃^{[j]}(z + h) = Σ_i C(j+i,i) P̃^{[j+i]}(z) h^i`, and Hasse derivatives of `p`-integral
polynomials are `p`-integral. The `i = 1` term is `≡ 0 (mod p²)` when `j + 1 < μ`, and the
terms with `i ≥ 2` are divisible by `p²`. ∎

**Proposition 2.4 (what the canonical lift contributes; PROVED).** Let `b ∈ B₁` be complete,
`0 ≤ j ≤ m−2`, and read `y_k` mod `p`. Then

    F̃^{[j]}(b̃)/p ≡ (C(D,j)/2)·[ Σ_k c_k £₁(a_k/y_k) y_k^{m−1−j}  +  Σ_k c_k u(a_k) y_k^{m−2−j} ]
                 ≡ (C(D,j)/2)·[ Σ_k c_k u(y_k) y_k^{m−2−j}  −  Σ_k c_k ε(a_k,b) y_k^{m−2−j} ]   (mod p).

Now build `F̃_T` from the Teichmüller nodes `ω(a_k)`. It has the same reduction, and
`F̃_T^{[j]}(ω(b))/p ≡ (C(D,j)/2)Σ_k c_k £₁(a_k/y_k)y_k^{m−1−j}` for all `j ≤ m−1`.

*Proof.* Substitute Cor. 1.3(1), (2) and (4) into Prop. 2.2. In the Witt form, the `u(b)`
term is `u(b)Σ_k c_k y_k^{m−2−j}`, which vanishes by Lemma 2.1 because `m−2−j ∈ [0, m−2]`. ∎

*Reading.* The second digits at complete points split into two parts:

- a Teichmüller part, which is the same formula over `W₂(F_{q'})` for every `q'`
  (Prop. 1.5);
- the Taylor data at `b` of the **fixed** polynomial `h_A(x) = Σ_k c_k u(a_k)(x+a_k)^{m−2}`
  (degree `≤ m−2`), up to the unit factors `C(m−2, j)`.

`h_A` does not depend on `B` or on completeness. It is exactly the effect of lifting `A`
canonically instead of by Teichmüller (`F̃ − F̃_T ∈ pℤ_(p)[x]`). **So the prime-specific part
of the lift carries no `B`-dependent information about the Hanson–Petridis data.** In the
carry form the same data read as "`u` evaluated at the residues `a_k + b`" minus a threshold
term: `ε(a_k,b) = [b̃ ≥ p − ã_k]` is piecewise constant in `b̃`. Both are fixed functions,
and `u` has full degree `p−1` (Lemma 1.4).

**Proposition 2.5 (extra vanishing mod p²; PROVED).** Let `b ∈ B₁` be complete.

1. For `1 ≤ M ≤ m−1`: `F̃^{[j]}(b̃) ≡ 0 (mod p²)` for all `j < M` iff
   `Σ_k c_k q(y_k) y_k^t ≡ 0` for `m−M ≤ t ≤ m−1`.
2. Some lift `z ≡ b` has `F̃^{[j]}(z) ≡ 0 (mod p²)` for all `j ≤ m−1` iff some integer
   `z ≡ b (mod p)` has `z + ã_k ∈ ker q` (that is, `(z+ã_k)^{p−1} ≡ 1 (mod p²)`) for every `k`.

*Proof.* (1) is Prop. 2.2, since `C(D,j)` is a unit.

(2) The conditions for `j ≤ m−2` form an `(m−1) × m` Vandermonde-type system in
`w_k = c_k q(y_k)`. Its kernel is spanned by `w_k = c_k/y_k`, because
`Σ_k c_k y_k^{t−1} = 0` for `1 ≤ t ≤ m−1`. So `q(y_k) = κ/y_k` for some `κ`.

At the lift `z = b̃ + pt`, `q(z + ã_k) = (κ − t)/y_k` by Lemma 1.1(a). Using
`Σ_k c_k/y_k = (−1)^{m−1}/Π_k y_k` (partial fractions),
`F̃^{[m−1]}(z) ≡ (p/2)C(D,m−1)(κ−t)(−1)^{m−1}/Π_k y_k`. This vanishes iff `t = κ`, iff every
`q(z+ã_k) = 0`. The converse follows from Prop. 2.2. ∎

*Example 1 (Wieferich).* For `A = {0,1}`, `F̃(x) = −1 − x^D + (x+1)^D` and
`F̃(1) = 2^{d+1} − 2`. For `p ≡ ±1 (mod 8)` the point `b = 1` is complete and
`F̃(1) ≡ p·q(2) (mod p²)`. So `F̃_{0,1}` vanishes `mod p²` at `b = 1` **iff `p` is a
Wieferich prime**. Checked for all 244 primes `p ≡ ±1 (mod 8)` with `7 ≤ p ≤ 3600`: the only
one is `p = 3511`; `1093 ≡ 5 (mod 8)` is excluded.

*Example 2 (Cauchy–Mirimanoff).* If `D ≡ 1 (mod 6)` then `(x²+x+1)²` divides
`(x+1)^D − x^D − 1` in `ℤ[x]`. Proof: at a primitive cube root of unity `ρ`,
`(ρ+1)^D − ρ^D − 1 = −ρ^{2D} − ρ^D − 1 = 0`, and `(ρ+1)^{D−1} − ρ^{D−1} = ρ^{2(D−1)} − ρ^{D−1} = 0`,
using `D` odd and `3 | D−1`.

For `p ≡ 1 (mod 12)`, `D = d+1 ≡ 1 (mod 6)`. The two points `b = −ζ` (`ζ` a primitive sixth
root of unity) are complete: `−ζ = ζ⁴` and `1 − ζ = ζ^{−1}` are squares. At them `F̃_{0,1}`
has an **exact double root in `ℤ_p`**, at `−ω(ζ)`. This explains all six `mod p²` zeros of
the `A = {0,1}` rows at `p = 13, 37, 61` (verifier C12: 23 primes `≤ 600`, 46 points).

*Data* (C8, C9, C11, C12). Among 1,139 complete points `b ∉ −A` in 56 configurations
(tight examples, the balanced worker's extremal witnesses at `p = 101, 197, 509, 997`, greedy
cliques and balanced pairs up to `p = 1009`, random `A` with `B = B(A)`), only 8 have
`F̃(b̃) ≡ 0 (mod p²)`. At rate `1/p` one would expect 4.9. Six of the eight are the
Cauchy–Mirimanoff points; `C_summary` in the results file lists all of them. The rows with
`m = |A| ≥ 3` show no extra vanishing of order `≥ 2`.

**Definition and Proposition 2.6 (multiplicity over ℤ_p; ℤ_p-counting = F_p-counting;
PROVED).** Let `P̃ ∈ ℤ_(p)[x]` with `P̄ := P̃ mod p ≠ 0`. Let `K/ℚ_p` be a splitting field with
valuation `v(p) = 1`, ring of integers `O_K` and residue field `k`. For `b ∈ F_p` define
`μ_b(P̃) := #{roots ρ of P̃ in K, with multiplicity : v(ρ − b̃) > 0}`. Then
`μ_b(P̃) = ord_b P̄` for every `b`, and `Σ_b μ_b(P̃) ≤ deg P̄`.

The "order `mod p²`", `ν_b(P̃;z) := min{j : P̃^{[j]}(z) ≢ 0 (mod p²)}`, satisfies
`ν_b ≤ ord_b P̄`.

*Proof.* Write `P̃ = c·Π_i(x − ρ_i)` and set:

- `U = Π_{v(ρ_i)≥0}(x − ρ_i) ∈ O_K[x]`, monic;
- `V = Π_{v(ρ_i)<0}(1 − x/ρ_i) ∈ O_K[x]`, with constant term 1 and all other coefficients
  in `m_K`;
- `c' = c·Π_{v(ρ_i)<0}(−ρ_i)`.

Then `P̃ = c'UV`. By Gauss's lemma over the discrete valuation ring `O_K` (reductions of
polynomials with a unit coefficient are nonzero in `k[x]`, a domain), `UV` has a unit
coefficient. So does `P̃`, because `P̄ ≠ 0`. Hence `c' ∈ O_K^*` and `P̄ = c̄'·Π_{v(ρ_i)≥0}(x − ρ̄_i)`.

Now `ord_b P̄` counts the integral `ρ_i` with `ρ̄_i = b`, i.e. `v(ρ_i − b̃) > 0`. Non-integral
roots have `v(ρ_i − b̃) < 0`. The bound on `ν_b` is trivial: vanishing `mod p²` implies
vanishing `mod p`. ∎

**Consequences.** Any Stepanov-type count over `ℤ_p` (roots near `F_p`-points, with
multiplicity) is literally the `F_p` count for the normalized reduction. Counting over `ℤ/p²`
itself is not available at all: `ℤ/p²` is not a domain, `x²` has the `p` roots `pℤ/p²`, and
by Lemma 2.3 `F̃` has either 0 or `p` roots above each `b` with `ord_b F_A ≥ 2`. So a lift
can help only by **certifying `mod p` vanishing of a new `F_p`-polynomial of small degree**,
never by counting.

**Proposition 2.7 (complete points are Eisenstein clusters; PROVED, using the standard
Eisenstein criterion over ℤ_p).** Let `b ∈ B₁` be complete with `F̃(b̃) ≢ 0 (mod p²)`. Then the
`m` roots of `F̃` in the disc `v(x − b̃) > 0` are conjugate over `ℚ_p`, generate a totally
ramified extension of degree `m`, and all satisfy `v(ρ − b̃) = 1/m`.

*Proof.* Set `U_b := Π_{v(ρ_i−b̃)>0}(x − ρ_i)`. The index set is Galois-stable, so
`U_b ∈ ℤ_p[x]`, and `deg U_b = m` by Prop. 2.6 and Theorem 1.1(1) of the balanced note. The
non-leading coefficients of `U_b(x + b̃)` have positive valuation, so they lie in `pℤ_p`.

In `F̃ = c'·U_b·U_rest·V`, the values `U_rest(b̃)` and `V(b̃)` are units. So
`v(U_b(b̃)) = v(F̃(b̃)) = 1`, and `U_b(x + b̃)` is Eisenstein. ∎

Verifier C9: every complete point has first unit Taylor coefficient at index exactly `m`, and
1,131 of 1,139 are Eisenstein. So "multiplicity `m` over `F_p`" becomes, `p`-adically, a
ramified cluster of `m` distinct roots, not a higher-order root. The exceptions are the
Wieferich/Cauchy–Mirimanoff type points of Prop. 2.5.

**Proposition 2.8 (the second digit is an F_p-polynomial with no forced zeros; PROVED).**
Let `A + B ⊆ Q ∪ {0}`, so `F_A = λP₁^mP₀^{m−1}G_A` (balanced note, Theorem 1.1(3)). Put
`F₀ := C(D,m−1)·Π_{B₁}(x − b̃)^m Π_{B₀}(x − b̃)^{m−1}·G̃`, with any monic lift `G̃` of `G_A`.
Then `F₁ := (F̃ − F₀)/p ∈ ℤ_(p)[x]` has degree `≤ d−1`. For `b ∈ B₁`, `j < m` (and `b ∈ B₀`,
`j < m−1`), `F̃^{[j]}(b̃) = p·F₁^{[j]}(b̃)`.

So "`F̃` vanishes to order `M` `mod p²` at `b`" means exactly `ord_b F̄₁ ≥ M`, and the values
`F̄₁^{[j]}(b)` are the Fermat-quotient combinations of Prop. 2.2. The Taylor data of `F̄₁` at
`B` do not depend on `G̃`. Conversely, any prescribed Taylor data at `B` (`mn − r ≤ d`
numbers) are realized by exactly one polynomial of degree `< mn − r` (CRT). So the reduction
`F_A` imposes no condition on the second digits: they are free parameters, and they are in
bijection with the Fermat-quotient matrix `(q(ã+b̃))_{a∈A, b∈B}` (Vandermonde, Prop. 2.2).

*Proof.* `F̃` and `F₀` have degree `d`, the same leading coefficient (Lemma 2.1) and the same
reduction. `F₀` vanishes exactly to the stated orders at `b̃`. Changing `G̃` changes `F₁` by a
multiple of `P̃₁^mP̃₀^{m−1}`. ∎

A Stepanov count with `F̄₁` would give `Σ_{b∈B} ord_b F̄₁ ≤ d − 1`. That needs `F̄₁(b) = 0` on
`B`, which by Prop. 2.5 is a Wieferich-type condition that fails at almost all complete
points (data above). **Task 2's answer: no combination built from `F̃_A` vanishes `mod p²` to
higher order on complete points except under such arithmetic accidents, and even extra
vanishing would not raise the root count (Prop. 2.6).**

**Proposition 2.9 (free extra vanishing costs a factor p in degree; PROVED).**
For `j ≥ 0` let `F̃_j := −1 + Σ_k c̃_k(x+ã_k)^{D+j(p−1)}`. At complete `b ∉ −A`, for
`i ≤ m−1`,

    F̃_j^{[i]}(z) ≡ (2j+1)·(p/2)·C(D+j(p−1), i)·Σ_k c_k q(y_k) y_k^{m−1−i}   (mod p²).

So `2j+1 ≡ 0`, i.e. `j ≡ d (mod p)`, makes the second digits vanish at every complete point
of every `A`. For `j ≢ d` they do not vanish in general (the `p = 17` witness of §5). The
first case is `F̃_* := F̃_d = −1 + Σ c̃_k(x+ã_k)^{pd+m−1}`. It vanishes `mod p²` to order `m` at
every lift of every complete point, but `F̃_* mod p` has degree exactly `pd`: its leading
coefficient is `C(pd+m−1, m−1) ≡ 1` by Lucas. It agrees with `F_A` as a function on `F_p`.
Root counting with it gives `mn ≤ pd + r`.

*Proof.* Lemma 1.1(c) with exponent `d + j(p−1)`, then the argument of Prop. 2.2. ∎

(Verifier C10: 3,151 Taylor coefficients.)

**Remark 2.10 (height; HEURISTIC, with exact instances).** A product-formula argument would
play the `p`-adic clustering (Props. 2.6–2.7) against the archimedean size of `Disc(F̃)`. The
clusters force only `v_p(Disc) ≳ Σ_b(μ_b − 1) ≤ d`. Here `Disc` is the discriminant of the
integral polynomial `L·F̃`, with `L` the common denominator of the `c̃_k`, which is prime
to `p`. The archimedean size is enormously larger: for `(p, A) = (17,{0,1}), (29,{0,1}), (29,{0,1,5}), (37,{0,1,11}), (13,{0,1,4})`
we find `v_p(Disc) = 3, 6, 7, 12, 4` against `log_p|Disc| = 12.8, 29.4, 126.6, 271.5, 32.0`
(verifier §G; `Disc ≠ 0` in all cases). No Hensel or height argument of this shape comes
close.

## 3. S(A,B) modulo p²

**Proposition 3.1 (PROVED; `p ≥ 5`).** For `A, B ⊆ F_p`, lifted canonically,

    Σ_{a∈A,b∈B} (ã + b̃)^d ≡ S(A,B) + (p/2)·T   (mod p²),     T := Σ_{a+b≠0} χ(a+b)·q(ã + b̃),

where `q` is taken of the **integer** `ã + b̃ ∈ [1, 2p−2]`. When `ã + b̃ > p`,
`q(ã+b̃) = φ(a+b) − 1/(a+b)`.

*Proof.* Lemma 1.1(c) termwise. The terms with `ã + b̃ ∈ {0, p}` give `0^d = 0` and
`p^d ≡ 0 (mod p²)` (`d ≥ 2`). The carry rule is Lemma 1.1(a). ∎

**Lemma 3.2 (PROVED).** `S ≡ mn − r (mod 2)` and `|S| ≤ mn − r`. So if `mn − r < p`, `S` is
the unique integer in `[−(mn−r), mn−r]` with the right residue mod `p` and the right parity.
(If `mn < p/2`, the residue alone suffices.) (Verifier D2: 95 cases.)

**Proposition 3.3 (decompositions of T; PROVED).**

    T = T_W + T_u = T_φ − T_ε,

where
- `T_W = Σ χ(a+b)£₁(a/(a+b))`,
- `T_u = Σ χ(a+b)(u(a)+u(b))/(a+b)`,
- `T_φ = Σ_{a,b} ψ(a+b)` with `ψ := χφ`, a function of the residue `a+b` alone,
- `T_ε = Σ χ(a+b)ε(a,b)/(a+b)`.

This is Cor. 1.3 summed. For a complete biclique, all `χ = 1` off `a + b = 0`.

**Proposition 3.4 (behaviour under the symmetries of the problem; PROVED).** For `s ∈ Q` and
`t ∈ F_p`, the map `(A,B) ↦ (sA + t, sB − t)` preserves `A + B ⊆ Q ∪ {0}`, the sizes, `r`,
every bad-partner count, and hence every LP profile variable. Under it:

- `T_φ` depends only on `s`;
- `T_W` depends only on `t/s`;
- the second digits `F̃^{[j]}(b̃)/p` at complete points depend only on `s`, the point, and the
  carry vector `(ε(a'_k, b'))_k`;
- `Σ_{t∈F_p} T(A + t, B − t) ≡ S(A,B) + Σ_{a+b≠0} χ(a+b)/(a+b) (mod p)`.

*Proof.*
- `T_φ` depends only on the multiset `{s(a+b)}`.
- For `T_W`: `£₁((sa+t)/(s(a+b))) = £₁((a + t/s)/(a+b))` and `χ(s) = 1`.
- For the digits: by Cor. 1.3(2), `q(ã'+b̃') = φ(s(a+b)) − ε'/(s(a+b))`, and
  `c'_k = s^{−(m−1)}c_k`.
- For the average: by Cor. 1.3(1), `Σ_t £₁((a+t)/z) = Σ_{s'} £₁(s') = 1` and
  `Σ_t u(a+t) = Σ_t u(b−t) = ½` (Lemma 1.4). So `Σ_t q = 1 + 1/z` for each pair with sum
  `z = a+b`. ∎

**Data (exact, verifier §E; 29 orbits).** `T` is **not** a function of the profile.

- On full orbits (`p ≤ 197`), `T` takes between 5 values (`p = 13`, tight 3-clique) and all
  `p` values (e.g. every extremal witness at `p = 101` with `|A| ≠ 4`, and at `p = 197` with
  `|A| ∈ {3,5,6,7}`).
- On the partial orbits at `p = 509, 997, 1009` it takes 133–806 values.
- The vanishing pattern of the second digits also changes within an orbit. Example:
  `p = 101`, `A = {0,1,4}`, `B = B(A)` (15 elements) has no complete point with vanishing
  second digit. Its image under `s = 4`, `t = 85` has one, at `b = 21`.
- Recorded witnesses exist for 20 of the 29 orbits (`E_orbits` in the results file).

The translation average, the only natural orbit-invariant built from `T`, is determined by
`mod p` data (Prop. 3.4).

**HEURISTIC evidence that Fermat quotients over complete bicliques carry no structure**
(verifier §F). Because `q(ã+b̃)` depends only on the integer sum, statistics are taken over
distinct integer sums, and exponential sums are normalized by the `ℓ²` norm of the
representation function `r(z)`. Over 46 complete bicliques with `|B₁| ≥ 2`:

- the matrix `(q(ã+b̃))_{a∈A,b∈B₁}` has full rank `min(m, n₁)` in every case;
- 14 zero values occur among the distinct sums, against 11.6 expected at rate `1/p` (random
  rectangles of the same shapes: 17);
- `max_h |Σ_{pairs} e(h q/p)| / ‖r‖₂` has mean 2.38, against 2.31 for random rectangles;
- the Teichmüller second-digit function
  `σ_A(x) = Σ_k c_kχ(x+a_k)£₁(a_k/(x+a_k))(x+a_k)^{m−1}` has degree exactly `p − 1` in all 25
  tested cases, against `deg F_A = d`;
- `max_h|Σ_{r<p} e(hφ(r)/p)|/√p` and `max_h|Σ_{r<p} χ(r)e(hφ(r)/p)|/√p` stay in
  `[1.5, 6.1]` for `29 ≤ p ≤ 383`.

The last bullet concerns a character of order `2p` modulo `p²` over an interval of length `p`,
the natural "Paley problem modulo `p²`". It is no easier than the original.

## 4. Verdict

**No new constraint. The obstruction is precise:**

1. *Counting.* By Prop. 2.6, multiplicities over `ℤ_p` are exactly the `F_p` multiplicities
   of the normalized reduction. Over `ℤ/p²` there is no root-counting bound at all. A lift can
   only help by producing a new `F_p`-polynomial that vanishes on `B` with degree-to-order
   ratio below `d/m`.
2. *The candidate polynomials.*
   - `F̃_A` itself acquires extra `mod p²` vanishing only under Wieferich-type conditions
     (Prop. 2.5), which completeness does not imply. Its second-digit polynomial `F̄₁`
     (degree `≤ d−1`) has, at complete points, the values of Fermat-quotient combinations
     (Prop. 2.8).
   - The polynomial that vanishes for free, `F̃_*`, has reduction of degree `pd`
     (Prop. 2.9).
   - All second-digit functions have full degree `p − 1` (Lemma 1.4; data §F).
3. *Prime-specificity.*
   - The Teichmüller/Witt half of the lift exists verbatim over `W₂(F_{p²})`, up to a
     Frobenius twist (Prop. 1.5), so it cannot be the prime-specific input.
   - The genuinely prime-specific half (the section `[0,p)`) enters as the separable term
     `u(a) + u(b)`, equivalently the carries (Cor. 1.3(3)). In the Hanson–Petridis data this
     term is `B`-independent (Prop. 2.4).
4. *Freedom.* Completeness imposes no condition on the second digits (Prop. 2.8). They are in
   bijection with the Fermat-quotient matrix, which on complete bicliques behaves like that
   of random rectangles (§3, HEURISTIC).

**The LP-tight profiles.** The profiles of the balanced worker (Theorem 3.5 there) are vectors
of bad-partner counts `n_j, n'_j`, not sets. A Fermat-quotient constraint needs the integers
`ã, b̃`, so it can only be evaluated on explicit sets. It would feed the LP only through a
consequence invariant under the symmetry group `(sA + t, sB − t)`, which fixes every profile.
The `mod p²` data are not invariant (Prop. 3.4 and data). The natural invariant, the average
of `T` over translations, is a function of `mod p` data. We found no orbit-invariant
consequence and none is suggested by the structure above (HEURISTIC). So nothing can be added
to the LP, and the tight profiles are not excluded.

**The `F_{p²}` test.** Applied to this programme, the test shows that the lift contributes only
two things. One is the Witt part, which exists over `F_{p²}`. The other is a coboundary. The
only prime-specific step in anything above remains `D < p`, exactly as in Hanson–Petridis.

## 5. Refuted, open, obligations

- REFUTED: "complete points are automatically roots of the lifted polynomial `mod p²`". Witness:
  `p = 17`, `A = {0,1}`, second digits `13, 2, 13` at `b = 1, 8, 15`. More generally only 8 of
  1,139 complete points vanish, and 6 of those are explained by Cauchy–Mirimanoff.
- REFUTED: "the `mod p²` data (`T`, second-digit zero patterns) are determined by the
  configuration's profile". Witnesses in §3 and `E_orbits`.
- OPEN (unchanged): a prime-specific input coupling `A` and `B` for balanced bicliques. The
  one route left open here: a **new** `F_p`-polynomial whose vanishing on `B` is certified by
  `mod p²` information and whose degree-to-order ratio on `B` is below `d/m`. We know of no
  candidate. Every
  second-digit object we could build has degree `p − 1` or `pd`.
- Not done: a literature search. We make no novelty claim. The statements are elementary,
  Lemma 1.1 and Example 2 are classical, and Prop. 2.6 is a standard Newton-polygon fact
  with a self-contained proof here.

## 6. Verification (`results/stepanov_padic_2026_09_29.json`)

| section | what is checked | checks |
|---|---|---|
| A | refined Euler, exponent family, `x^{pd} ≡ χ`, shift rule, homomorphism, Teichmüller kernel, all `x < p²`, `p ≤ 61` | 119,399 |
| B | exact integer cocycle; `Γ = (a+b)£₁`; Witt and carry forms of `q(ã+b̃)` and the carry identity for all pairs, `p ≤ 151`; `£₁` identities; `Σu = ½`; `deg u = deg £₁ = p−1` | 1,093,354 |
| C | degree/leading coefficient; Hasse vs coefficients; Taylor formula `mod p²` at every point (57,990); complete and `B₀` formulas; lift independence; decomposition; Teichmüller polynomial; extra-vanishing criterion; Newton polygons; `F̃_*`; Wieferich (`p < 3600`); Cauchy–Mirimanoff | 82,158 |
| D | `S mod p²` identity, recovery of `S`, `T` decompositions, translation average (120 cases incl. random rectangles) | 575 |
| E | 29 symmetry orbits: `T = T_φ − T_ε`, `T_φ = T_φ(s)`, `T_W = T_W(t/s)`, second digits = f(s, carries), non-invariance | 145 |
| F | ranks, degree of `σ_A` | 71 |
| G | `deg_ℚ F̃ = d`, `Disc ≠ 0` | 10 |
| H | refined Euler and the Witt second-digit formula in `W₂(F_{p²})`, `p ≤ 11` | 18,344 |

Total **1,314,056 checks, 0 failures**; runtime about 8 s on one process.
