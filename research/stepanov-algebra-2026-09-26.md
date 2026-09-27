# Stepanov wave, direction `algebra`: exact algebraic reformulations of HP and RHP (2026-09-26)

**Status:** PROVED: the lifting lemma (below `p/2`, and with parity below
`p`), the central-binomial kernel `binom(d,j) ≡ (−1/4)^j binom(2j,j)`
with `(1+u)^d ≡ (1+u)^{−1/2} mod u^p`, the gcd/point-count forms of HP and
RHP, the identity `F(x) = h_d(A+x) − 1` (the HP polynomial is a complete
homogeneous symmetric polynomial, derivative = shift, jets = RS syndromes
of the bad sets), a sign-pattern HP (every fibre of `b ↦ (χ(a+b))_a` has
size `≤ (d+m−1)/m`), a Hankel auxiliary polynomial giving
`N_1 ≤ (2d−2)/(m−2)` (and `(e+1)(d−e)/(m−2e)` for general `e` when the
Hankel determinant is nonzero), and precise reasons (§3) why no
list-decoding bound (Johnson, Guruswami–Sudan, capacity-achieving,
linear-algebraic folded) beats the second moment's factor 2. REFUTED (with witnesses): the
least-absolute-residue lift above `p/2`; the HP bound `d/m` for mixed sign
patterns; multiplicity for the unweighted power-sum polynomial; HP for
multisets; large minimum distance of the translate code; any
field-uniform Gram/Johnson bound with constant below `2p/(p+1)`. OPEN:
whether any use of the syndrome structure reaches `f(η) → 1` (RHP).
Nothing here proves RHP or the conjecture.

Verifier: `experiments/stepanov_algebra_2026_09_26.py`, writing
`results/stepanov_algebra_2026_09_26.json`.

Notation as in `research/stepanov-brief-2026-09-26.md`: `p` odd prime,
`d=(p−1)/2`, `χ` Legendre with `χ(0)=0`, `Q` the nonzero squares,
`A={a_1,…,a_M}` (`m=M=|A|`), `n=|B|`, `S=S(A,B)=Σ_{a,b}χ(a+b)`,
`r=|B∩(−A)|`, `e_b=#{a: χ(a+b)=−1}`, `c_k=1/Π_{l≠k}(a_k−a_l)`,
`D=d+M−1`, `F(x)=−1+Σ_k c_k(x+a_k)^D`. `P_j(X)=Σ_{x∈X}x^j` with `0^0=1`.
`h_i(y_1,…,y_M)` is the complete homogeneous symmetric polynomial
(`h_0=1`, `h_i=0` for `i<0`), `e_j` the elementary one.

## 1. Three exact reformulations (PROVED)

### 1a. Lifting from one element of `F_p`

**Lemma 1.1.** Let `R = Σ_{a∈A,b∈B}(a+b)^d ∈ F_p`.
(i) `R ≡ S (mod p)`.
(ii) If `mn < p/2` (equivalently `mn ≤ d`), `S` is the least absolute
residue of `R`.
(iii) If `mn < p`, `S` is the unique integer in `[−mn, mn]` with
`S ≡ R (mod p)` and `S ≡ mn − r (mod 2)`.

*Proof.* (i) Euler: `x^d ≡ χ(x)` for `x ∈ F_p^×`, and `0^d = 0 = χ(0)`
since `d ≥ 1`. (ii) `|S| ≤ mn < p/2`, and the integers in `(−p/2,p/2)`
are a complete residue system. (iii) Each nonzero sum contributes `±1`,
so `S = (mn − r) − 2·#{(a,b): χ(a+b) = −1}`, hence
`S ≡ mn − r (mod 2)`; with (i) and `p` odd, `S` is fixed modulo `2p`, and
`[−mn,mn]` has length `2mn < 2p`. ∎

So up to the square-root scale (`mn < p/2`) the conjecture is a statement
about one element of `F_p`, and up to the Chung scale (`mn < p`) about one
element of `F_p` together with the parity of `r`.

*Warning (checked).* For `mn > p/2` the least absolute residue of `R` is
in general **not** `S`: witness `p=101`, `A={0,1,2}`, `B` = the 33 values
of `b` with the largest `Σ_aχ(a+b)`, `mn=99`, `S=58`, least absolute
residue `−43`; the rule (iii) recovers `58` (also at `p=103,197,401`).
Consequently the hypothesis "the least absolute residue `v` of `R`
satisfies `v ≥ (1−2η)mn`" forces `mn ≤ d/(1−2η)` for the trivial reason
`v ≤ d`; it says nothing about pairs with `mn > p/2`, where the target
consequence lives (§2.6).

For complete bicliques (`A+B ⊆ Q∪{0}`), `S = mn − r ≤ d` by HP, so there
the least absolute residue always equals `S`.

### 1b. The central-binomial kernel

**Lemma 1.2.** For `0 ≤ j ≤ d`,
`binom(d,j) ≡ binom(−1/2, j) = (−1/4)^j·binom(2j,j) (mod p)`; for
`d < j < p`, `binom(2j,j) ≡ 0 (mod p)`; and
`β_j := (−1/4)^j binom(2j,j)` satisfies `β_j ≡ β_{d−j}`.

*Proof.* `2d = p−1 ≡ −1`, so `d ≡ −1/2` in `F_p`. `binom(d,j) = (d)_j/j!`
with `j! ` invertible for `j ≤ d < p`, and `(d)_j` is an integer polynomial
in `d`, so `binom(d,j) ≡ (−1/2)_j/j! = binom(−1/2,j)`. The identity
`binom(−1/2,j) = (−1)^j (1·3···(2j−1))/(2^j j!) = (−1/4)^j (2j)!/j!^2`
is standard. For `d<j<p`, `p ≤ 2j < 2p` divides `(2j)!` once and does
not divide `j!^2`. Symmetry: `binom(d,j) = binom(d,d−j)`. ∎

**Corollary 1.3 (kernel formula).** For all `A, B ⊆ F_p`,

    S(A,B) ≡ Σ_{j=0}^{d} β_j · P_j(A) · P_{d−j}(B)   (mod p).

*Proof.* Expand `(a+b)^d = Σ_j binom(d,j) a^j b^{d−j}` (with `0^0=1`,
which is consistent when `a=0` or `b=0`), sum, apply Lemmas 1.1(i), 1.2. ∎

**Lemma 1.4 (power series).** In `F_p[[u]]` put
`(1+u)^{−1/2} := Σ_{j≥0} β_j u^j` (the coefficients `(−1/4)^j binom(2j,j)`
lie in `Z[1/2]`, so this is the reduction of the binomial series). Then
`((1+u)^{−1/2})^2 (1+u) = 1` and

    (1+u)^d ≡ (1+u)^{−1/2}   (mod u^p).

*Proof.* The first identity is `Σ_{i+j=n} binom(2i,i)binom(2j,j) = 4^n`.
For the second, `(1+u)^{2d+1} = (1+u)^p = 1+u^p`, so
`((1+u)^d)^2 ≡ (1+u)^{−1} (mod u^p)`. Two series with constant term 1
and equal squares modulo `u^p` are congruent modulo `u^p`:
`(f−g)(f+g) ≡ 0` and `f+g` is a unit because `p` is odd. ∎

*Meaning for `χ(1+u)`.* Let `T_d(u) = Σ_{j=0}^{d} β_j u^j`. Since `β_j ≡ 0`
for `d<j<p`, Lemma 1.4 says `(1+u)^d = T_d(u)` as polynomials, hence

    χ(1+u) = T_d(u)   for every u ∈ F_p,
    χ(x_0+u) = χ(x_0)·T_d(u/x_0)   for x_0 ≠ 0.

The Legendre symbol is the degree-`d` truncation of the Taylor series of
`1/√(1+u)`, and (derivative form) `(x^d)^{(j)}(x_0) = (−1/2)_j χ(x_0) x_0^{−j}`
with `(−1/2)_j = (−1/2)(−3/2)···(−1/2−j+1)`, nonzero for `j ≤ d`. This is
the formal sense in which `χ(x)` "is" `x^{−1/2}` locally at every nonzero
point, with the sign `χ(x_0)` as the only global datum. It is not an
identity of functions for the full series: the series has degree `≥ p`
terms that vanish only modulo `u^p`.

*The kernel is field-agnostic.* The same proof gives
`(1+u)^{(q−1)/2} ≡ (1+u)^{−1/2} (mod u^q)` over `F_q`, `q = p^k`
(checked by integer arithmetic for `q ∈ {9,27,25,125,49,121,169}`). So
Corollary 1.3 and Lemma 1.4 cannot by themselves distinguish `F_p` from
`F_{p²}`. The prime-field input of HP is elsewhere: `D = d+M−1 ≤ p−1`,
so every `binom(D,l)` and every falling factorial `(D)_j` is nonzero
mod `p`; over `F_{p^k}`, `k ≥ 2`, the exponent `d_q+M−1` exceeds `p`
and Lucas' theorem kills binomials and derivatives.

### 1c. The gcd form of HP and the point count `N_e`

For `a ∈ F_p` put `f_a(x) = (x+a)^d − 1` and `h_a(x) = (x+a)^{d+1} − (x+a)`.
`f_a` has the `d` simple roots `Q − a`; `h_a` has the `d+1` simple roots
`(Q∪{0}) − a`.

**Proposition 1.5 (gcd form).** For every nonempty `T ⊆ F_p`,

    deg gcd_{a∈T} f_a ≤ d/|T|,       |T|·deg gcd_{a∈T} h_a ≤ d + |B'_T ∩ (−T)|,

where `B'_T` is the common root set of the `h_a`. Conversely the first
inequality for all `T` is HP restricted to `A+B ⊆ Q`, and the second is
HP in full.

*Proof.* Each `f_a` is squarefree and splits, so the gcd is
`Π_{b∈B_T}(x−b)` with `B_T = {b : T+b ⊆ Q}`, of degree `|B_T|`. HP with
`A=T`, `B=B_T`, `r=0` gives `|T||B_T| ≤ d`. Likewise for `h_a` with
`B'_T = {b: T+b ⊆ Q∪{0}}`. Conversely any `B` with `A+B ⊆ Q` lies in
`B_A`. ∎

**Definition.** `N_e(A) = #{b ∈ F_p : e_b ≤ e}` and
`N_e^×(A) = #{b ∉ −A : e_b ≤ e}`. Since each `b ∉ −A` has `χ(a+b)=±1`,
and `h_a(b) = 0 ⟺ χ(a+b) ∈ {0,1}`,

    N_e(A) = #{b ∈ F_p : b is a root of at least m − e of h_{a_1},…,h_{a_m}},
    N_e^×(A) = #{b ∉ −A : b is a root of at least m − e of f_{a_1},…,f_{a_m}}.

**RHP restated.** Because `B ↦ m|B| − |B∩(−A)|` is monotone, RHP(f,g)
is equivalent to: for all `p`, all `A` with `m ≤ √p`, all `0 ≤ e < m/2`,

    m·N_e(A) − #{b ∈ −A : e_b ≤ e}  ≤  f(e/m)·d + g(m),

i.e. *the number of points of `F_p` lying on at least `m−e` of the `m`
translated curves `h_a = 0` is at most `(f(e/m)+o(1))·d/m`.* HP is the
case `e=0`, `f=1`, `g=0`. (Checked: 510 random `T` for each gcd identity
and bound, 606 point-count identities.)

## 2. HP as a statement about power-sum sequences (PROVED)

### 2.1. The auxiliary polynomial is a complete homogeneous polynomial

**Lemma 2.1 (divided differences).** For distinct `y_1,…,y_M` and
`N ≥ 0`, `Σ_k c_k y_k^N = h_{N−M+1}(y)` where `c_k = 1/Π_{l≠k}(y_k−y_l)`.
(Partial fractions of `Π_k(1−y_k t)^{−1} = Σ_k c_k y_k^{M−1}/(1−y_k t)`.)
The `c_k` depend only on differences, so they are the same for `A` and
`A+x`.

**Proposition 2.2.** `F(x) = h_d(a_1+x, …, a_M+x) − 1` and, for `j ≥ 1`,
`F^{(j)}(x) = (D)_j · h_{d−j}(A+x)`, i.e. `d/dx` acts on the sequence
`i ↦ h_i(A+x)` as the shift `i ↦ i−1` times `(i+M−1)`. Leading coefficient
`binom(D, M−1) ≢ 0` because `D ≤ p−1`.

*Proof.* Lemma 2.1 with `N = D−j`, and `d/dx (x+a)^N = N(x+a)^{N−1}`. ∎

### 2.2. What HP uses about the sequences

The power sums `P_s(A+x) = Σ_k (x+a_k)^s` and the complete homogeneous
sums `h_i(A+x) = Σ_k c_k (x+a_k)^{i+M−1}` are two solutions of the **same**
order-`M` linear recurrence, with characteristic polynomial
`Π_k(T − a_k − x) = Σ_j (−1)^j e_j(A+x) T^{M−j}` (Newton; `i h_i = Σ_{k≤i} P_k h_{i−k}`,
invertible for `i<p`). They differ only in initial conditions:
`(h_{−M+1},…,h_{−1},h_0) = (0,…,0,1)` (the impulse response), while the
power sums have `P_{−j}(A+x) = Σ_a (a+x)^{−j} ≠ 0` in general.

**Proposition 2.3 (HP in sequence language).** For `b ∉ −A`:
(i) `A+b ⊆ Q` iff the impulse response is `d`-periodic,
`h_{i+d}(A+b) = h_i(A+b)` for all `i ≥ −M+1`;
(ii) in general the periodicity defect is
`h_{i+d}(A+b) − h_i(A+b) = −2 Σ_{k∈E(b)} c_k (a_k+b)^{i+M−1}`,
where `E(b) = {k: χ(a_k+b) = −1}`; it is a solution of an order-`e_b`
recurrence (roots `a_k+b`, `k∈E(b)`);
(iii) the `M`-jet of `F` at `b` is the window of the defect:
`F^{(j)}(b)/(D)_j = −2 Σ_{k∈E(b)} c_k (a_k+b)^{M−1−j}` for `0 ≤ j ≤ M−1`
(the brief's formula, verified).

*Proof.* `(a_k+b)^{d+s} = χ(a_k+b)(a_k+b)^s`; write `χ = 1 − 2·1_E`, apply
Lemma 2.1 and Proposition 2.2. For (i) "if" is (ii) with `E=∅`; "only if":
the defect sequence at `i = −M+1,…,0` is `−2V c_E` with `V` a
Vandermonde in the nonzero distinct `a_k+b`, so it vanishes only if
`E = ∅`. ∎

So the exact algebraic content of HP is:

1. **(A-side)** the characteristic polynomial of the recurrence satisfied
   by `(P_j(A))_j` is squarefree (`disc = det(P_{i+j}(A))_{0≤i,j<M} = Π_{k<l}(a_k−a_l)^2 ≠ 0`);
   this is exactly what allows the impulse-response solution `h_i(A+x)`
   with `M−1` vanishing initial terms (weights `c_k`);
2. **(Euler)** a complete biclique point `b` is a point where the impulse
   response is `d`-periodic, so the jet of `h_d(A+x)−1` at `b`, which by
   Proposition 2.2 is the window `(h_d−1, h_{d−1}, …, h_{d−M+1})(A+b)`,
   equals `(h_0−1, h_{−1}, …, h_{−M+1}) = 0`;
3. **(prime field)** `D = d+M−1 ≤ p−1`, so the shift factors `(i+M−1)` and
   the leading coefficient `binom(D,M−1)` are nonzero;
4. **(B-side)** only that `B` is a set of `n` distinct points: a nonzero
   polynomial of degree `d` has at most `d` roots with multiplicity. No
   property of `(P_j(B))_j` beyond squarefreeness of its characteristic
   polynomial is used.

The power sums themselves carry no multiplicity: for
`Φ_A(x) = Σ_a (x+a)^d = P_d(A+x)` one has
`Φ_A^{(j)}(b) = (d)_j P_{−j}(A+b)` at a complete point, and
`Φ_A'(b) = d Σ_a (a+b)^{−1}` is generically nonzero. Observed: at 389
complete points, `Φ_A − m` vanished to order exactly 1 at 382, order 2
at 6 and order 3 at 1. The unweighted route gives only the trivial `n ≤ d`.

Squarefreeness is essential on both sides: for the multiset
"`A = {0}` with multiplicity 10", `B = Q` at `p=101`, all sums are squares
and `mn = 500 > p` (witness in the results file). HP is false for
multisets, so any proof of RHP must use distinctness, as HP does through
the Vandermonde.

### 2.3. Consequence for bad points

**Proposition 2.4.** If `b ∉ −A` and `e_b ≥ 1` then `ord_b F ≤ e_b − 1`.

*Proof.* By Proposition 2.3(iii), if `ord_b F ≥ e_b` then
`Σ_{k∈E} c_k y_k^{s} = 0` for the `e_b` consecutive exponents
`s = M−1, …, M−e_b` (`y_k = a_k+b ≠ 0`); this is a nonsingular
`e_b×e_b` Vandermonde system (after factoring `y_k^{M−e_b}`) in the
nonzero unknowns `c_k`, contradiction. ∎

So one bad partner kills the root, and there is no lower bound on the
order at a bad point (observed orders at 5,309 bad points: 0 at all but
43, which have order 1).

### 2.4. Sign-pattern HP: every fibre of `b ↦ (χ(a+b))_a` is small

**Proposition 2.5.** Let `m ≤ (p+1)/2` and `u ∈ {±1}^A`. Then

    μ(u) := #{b ∉ −A : χ(a+b) = u_a for all a ∈ A}  ≤  (d + m − 1)/m,

and `μ(u) ≤ d/m` when `u` is constant.

*Proof.* Put `F_u(x) = −1 + Σ_k c_k u_k (x+a_k)^D`. At `b` in the fibre,
`(b+a_k)^{D−j} = u_k (b+a_k)^{M−1−j}`, so
`F_u^{(j)}(b) = (D)_j Σ_k c_k u_k^2 (b+a_k)^{M−1−j} − [j=0] = 0` for
`0 ≤ j ≤ M−1` (Lemma 2.1). `F_u ≢ 0`: otherwise
`Σ_k c_k u_k a_k^l = 0` for `0 ≤ l ≤ M−1` (all `binom(D,l) ≠ 0` as
`D < p`), and the Vandermonde forces `c_k u_k = 0`. `deg F_u ≤ D`. Hence
`m μ(u) ≤ d + m − 1`. For constant `u = ±1`, `F_u = ±(F+1) − 1` has
degree `d` (HP). ∎

Exhaustive check: all `A ∋ 0` with `|A| ∈ {3,4}` for `11 ≤ p ≤ 47`
(59,006 sets) and 120 random sets up to `p = 1009`. The `+m−1` is
necessary: at `p = 11`, `A = {0,1,3,4}` a mixed pattern has `μ = 2`,
`mμ = 8 = d+m−1 > d = 5`. (Proposition 2.5 is the natural extension of
HP's argument; I have not found it stated in the literature, but it is
routine and I do not claim novelty.)

Consequence: `N_e^×(A) = Σ_{wt(u) ≤ e} μ(u) ≤ V(m,e)·(d+m−1)/m` with
`V(m,e) = Σ_{i≤e} binom(m,i)`. This is useless for RHP (`V(m,1) = m+1`),
but it isolates the content of RHP exactly: **RHP says the Hamming ball of
radius `e` around the all-ones word carries at most `f(η)+o(1)` fibres'
worth of mass, although it contains `V(m,e)` patterns each of which may
carry up to one fibre.**

### 2.5. Linear complexity of the defect and the Hankel auxiliary polynomial

By Proposition 2.3(iii) the jet `z_j(b) := F^{(j)}(b)/(D)_j`,
`0 ≤ j ≤ M−1`, is an RS/BCH *syndrome*: error locators `a_k+b`,
`k ∈ E(b)`, error values `−2c_k`. A syndrome of weight `e` has linear
complexity `e` (Berlekamp–Massey), so all `(e+1)×(e+1)` Hankel minors
`det[z_{i+l+s}(b)]_{0≤i,l≤e}` vanish when `e_b ≤ e` and `2e+s ≤ M−1`.
This gives the natural "robust" auxiliary polynomials
`W_e(x) = det[z_{i+l}(x)]_{0≤i,l≤e}` with `z_j(x) = F^{(j)}(x)/(D)_j`.

**Proposition 2.6 (e = 1).** `W_1 = z_0 z_2 − z_1^2` has degree exactly
`2d−2`, with leading coefficient
`binom(D,M−1)binom(D−2,M−1) − binom(D−1,M−1)^2 = −binom(D−1,M−1)^2 (M−1)/(d(D−1)) ≢ 0`.
Its order is `≥ 2M−2` at points with `e_b = 0` and `≥ M−2` at points with
`e_b = 1` (`b ∉ −A`). Hence

    (M−2)·#{b∉−A: e_b = 1} + (2M−2)·#{b∉−A: e_b = 0} ≤ 2d − 2,
    so N_1^×(A) ≤ (2d−2)/(M−2).

*Proof.* `z_j' = (D−j) z_{j+1}`, so `W_1^{(s)}` is a fixed linear
combination (depending on `D, s` only) of products `z_i z_j` with
`i+j = s+2`. At `e_b = 1` with bad index `k`, `z_i(b) z_j(b) = 4c_k^2 y_k^{2M−2−i−j}`
depends only on `i+j` as long as `i,j ≤ M−1`, i.e. for `s ≤ M−3`; the same
combination evaluated for the pure translate `(x+a_k)^D` (whose `W_1`
vanishes identically, being a rank-one Hankel form) is zero, so
`W_1^{(s)}(b) = 0` for `s ≤ M−3`. At `e_b = 0` every product with
`min(i,j) ≤ M−1` vanishes, which covers `s ≤ 2M−3`. The leading
coefficient: `z_j` has leading coefficient `binom(D−j,M−1)` and degree `d−j`. ∎

Observed minimal orders (75 sets, `29 ≤ p ≤ 139`) are exactly `2M−2` and
`M−2`, so Proposition 2.6 is sharp for this polynomial. **It recovers
exactly the second-moment constant**: `m·N_1^× ≤ 2d·m/(m−2) ≈ 2d`, the
same as `m·N_1^× ≤ m^2(p−m)/(m−2)^2 ≈ 2d`.

**Proposition 2.7 (all `e`, conditional on `W_e ≢ 0`).** `deg W_e ≤
(e+1)(d−e)`, and `ord_b W_e ≥ M−2e` at every `b ∉ −A` with `e_b ≤ e`
(`≥ (e+1)(M−e)` if `e_b = 0`). Hence, if `W_e ≢ 0`,

    (M − 2e)·N_e^×(A) ≤ (e+1)(d − e).

*Proof.* Degree: each term of the determinant is a product of `e+1`
entries `z_{i+σ(i)}` of degrees `d−i−σ(i)`, total `(e+1)d − e(e+1)`.
Order: `W_e^{(s)}` is a universal combination (depending on `D, e, s`) of
products `z_{i_0}···z_{i_e}` with index sum `s+e(e+1)` and largest index
`≤ 2e+s`. At `b` with bad set `E`, `|E| ≤ e`, put
`G = Σ_{k∈E} α_k (x+a_k)^D` with `α_k = −2c_k (b+a_k)^{M−1−D}`; then
`z_i[G](b) = z_i(b)` for all `i ≤ M−1` (Proposition 2.3(iii)), the jets of
`G` satisfy the same rule `z_i' = (D−i)z_{i+1}`, and `W_e[G] ≡ 0` because
the Hankel matrix of a sum of `|E| ≤ e` geometric sequences has rank
`≤ e`. So `W_e^{(s)}(b) = (W_e[G])^{(s)}(b) = 0` whenever `2e+s ≤ M−1`.
If `e_b = 0` all `z_i(b)`, `i ≤ M−1`, vanish and a product is nonzero only
if every index is `≥ M`, which needs `s + e(e+1) ≥ (e+1)M`. ∎

`W_1 ≢ 0` is Proposition 2.6. For `e = 2` the verifier finds `deg W_2 =
3d−6` exactly in all 40 sets tested and minimal orders exactly `M−4`,
`2M−6`, `3M−6` at `e_b = 2, 1, 0` (773 points), so Proposition 2.7 is
sharp there and gives `N_2^× ≤ 3(d−2)/(M−4)`, worse than the second
moment. The determinantal use of linear complexity loses a factor `e+1`
against the target `d/M`; at `e = 1` this is exactly the factor 2 of the
second moment.

### 2.6. The robust question in this language

*Question.* If `S ≡ v (mod p)` with lift `v ≥ (1−2η)mn`, do the recurrence
structure and the kernel force `mn ≤ (1+o_η(1))p/2`?

*Answer (PROVED, in three parts).*
(a) As literally posed with the least absolute residue it is trivially
yes and empty: `v ≤ d` gives `mn ≤ d/(1−2η)`; but for `mn > p/2` the
least absolute residue is not `S` (witness in §1a), so the statement says
nothing about the pairs the target concerns.
(b) With the correct lift (Lemma 1.1(iii), `mn < p`), the statement is
exactly the brief's target consequence in the window `p/2 < mn < p`, and
by Markov it is implied by RHP; nothing is gained or lost by passing to
`F_p`, since each column sum `Σ_a χ(a+b) ≡ Φ_A(b) = P_d(A+b)` is itself
exactly liftable (`|Σ_aχ(a+b)| ≤ m < p/2`).
(c) What the structure gives: the scalar `Σ_b P_d(A+b)` uses the
power-sum solution of the recurrence, which carries no multiplicity
(§2.2), while the only known multiplicity comes from the impulse-response
solution, whose jets are syndromes of the bad sets (§2.3). The robust
content is a statement about how many points `b` have a jet of linear
complexity `≤ e`. The determinantal use of linear complexity gives the
constant `2` at `e = 1` and `e+1` in general (Propositions 2.6, 2.7), not `1+o(1)`. Whether a
non-determinantal use of the syndrome structure reaches `1+o(1)` is OPEN.

*Symbolic vanishing and the Hamming ball (PROVED).* HP imposes vanishing
at the unknown points *symbolically*: the conditions hold as identities in
`x`, so they hold at every `b` with the right pattern. For
`Q(x) = A_0(x) + Σ_k A_k(x) Y_k(x)`, `Y_k = (x+a_k)^d`, symbolic vanishing
(order 1) at all points with pattern `u` means `A_0 + Σ_k u_k A_k ≡ 0`.
If this holds for `u = 1` and for each `u = 1 − 2δ_k`, subtracting gives
`2A_k ≡ 0`, so `Q ≡ 0`. Hence no nonzero auxiliary polynomial linear in
the `Y_k` vanishes symbolically on the patterns of the Hamming ball of
radius 1; linear symbolic vanishing on a pattern set `U` is possible only
when the affine span of `U` is proper, and for the patterns with bad set
inside a fixed `E_0` it is HP on `A ∖ E_0`. Higher degree in `Y` pays `d`
per factor: `Π_{k∈S}(1−Y_k)` vanishes to order `≥ |S|−e` at good points
but has degree `|S|d`.

## 3. The list-decoding view

### 3.1. The dictionary (PROVED)

Fix `A`. For `b ∉ −A` the word `w_b = (χ(a+b))_{a∈A} ∈ {±1}^A` is the
evaluation on `A` of the translate `(x+b)^d`. `N_e^×(A)` is the number of
`b` whose word lies in the Hamming ball of radius `e` around the all-ones
word, counted **with multiplicity**: the map `b ↦ w_b` is far from
injective (fibres `μ(u)`, Proposition 2.5). HP says the fibre over `1` has
size `≤ d/m`; RHP says the ball carries `≤ (f(η)+o(1))d/m`.

**Lemma 3.1 (the family is canonical).** A polynomial `g ∈ F_p[x]` of
degree `≤ d` with values in `{±1}` on `F_p ∖ {−a}` and `g(−a) = 0` is
`±(x+a)^d`; with values in `{±1}` everywhere it is `±1`.

*Proof.* `g^2 − (x+a)^{p−1}` (resp. `g^2 − 1`) has degree `≤ p−1` and
vanishes on `F_p`, so it is zero; factor. ∎ (Exhaustively confirmed at
`p = 7, 11`: 1,773,962 polynomials, exactly `2p` resp. `2` solutions.)
So "list decoding of `{0,±1}`-valued low-degree codewords" is literally
the translate family; a bound for general `±1`-valued codewords is no
weaker than the problem itself.

### 3.2. Why the Reed–Solomon list-decoding theorems do not apply (PROVED)

The two natural ways to read the translates as RS codewords both fall
outside every list-decoding regime.

(i) *Evaluation set `A`.* The translates have degree `d ≥ m` (as
`m ≤ √p`), and `RS_A[d+1] = F_p^A`: every word is a codeword, so no bound
(unique or list, Johnson, Guruswami–Sudan, capacity) says anything. The
same holds for the bivariate reading `(x+y)^d` on the grid `A×B` when
`m+n ≤ d+2`. For GS-type interpolation `Q(x,y)` with multiplicity `r` at
the points `(a,1)`, `a ∈ A`, and `(1,d)`-weighted degree `< D'`: if
`D' ≤ d` then `Q = Q(x)` needs `D' > mr`, and the root-counting step needs
`t·r > D'` with agreement `t ≤ m`, impossible; if `D' > d` then
`r > d/m ≥ √p/2` and the monomial count `≈ D'^2/(2d) < m^2r^2/(2d)` is far
below the `m·r(r+1)/2` conditions. GS needs roughly `m > d`.

(ii) *Evaluation set `F_p` (or any `E` with `|E| > d`).* Now
`RS_E[d+1]` has dimension `k = d+1` and rate `(d+1)/|E| ≥ 1/2`, but a
translate agrees with the all-ones word in at most `d < k` positions
(exactly `d` on `F_p`): `(x+a)^d − 1` has at most `d` roots. Every
list-decoding theorem requires agreement `> √(kn) ≥ k` (Johnson,
Guruswami–Sudan: "e < n − √(kn)", ECCC TR98-043) or `≥ (R+γ)n > k`
(capacity). At agreement `d` the RS list around the all-ones word
contains all `1 + λΠ_{y∈Y}(x−y)`, `|Y| = d`, `λ ≠ 0`: `(p−1)binom(p,d)`
codewords. This applies equally to the capacity-achieving algorithms for
RS codes over prime fields on arbitrary evaluation sets claimed in
September 2026 (Jeronimo, arXiv:2609.05870; Brakensiek–Chen–Putterman–
Zhang–Zheng, arXiv:2609.08005; abstracts checked, results not refereed):
although they are prime-field statements and so not excluded by the
`F_{p²}` test, the parameter `agreement ≤ d < k` puts the translates
strictly below capacity. Moreover the question is transposed: RHP counts
*positions* at which many of the `m` given codewords agree with `1`, not
codewords agreeing with one word, and for general codewords the common
agreement of `m` codewords can be `d` (take `1 + λ_iΠ_{y∈Y}(x−y)`).

### 3.3. Johnson-type bounds equal the second moment (PROVED)

Let `L` be a set of `b ∉ −A` with `e_b ≤ e`, `τ = 1 − 2e/m`. Cauchy–Schwarz
on `Σ_{b∈L} w_b` against `1_A` (the Johnson argument) gives
`|L|^2 τ^2 m ≤ Σ_{b,b'∈L} ⟨w_b,w_{b'}⟩`. The right side can be bounded
only through the Gram matrix of the words, and there are three options.

(J1) *List-variable second moment:*
`Σ_{b,b'∈L}⟨w_b,w_{b'}⟩ = Σ_{a∈A} F_L(a)^2 ≤ Σ_{x∈F_p} F_L(x)^2 = |L|(p−|L|)`,
so `|L| ≤ pm/((m−2e)^2 + m)`.

(J2) *Spectral:* `≤ λ_max(W^TW)|L|` with `W` the `(p−m)×m` word matrix.
Since `W^TW = pI − J − Σ_{b∈−A} v_b v_b^T` with `v_b = (χ(a+b))_a`, one has
`W^TW ⪯ pI`, so `|L| ≤ pm/(m−2e)^2` (observed `λ_max ∈ [p−1, p)`). This is
the brief's second-moment bound `N_e^× ≤ m(p−m)/(m−2e)^2` up to lower order.

(J3) *Distinct-word Johnson:* if all pairs of distinct words have
correlation `≤ λ` and `τ^2 > λ`, the ball contains `≤ (1−λ)/(τ^2−λ)`
distinct words, each with fibre `≤ (d+m−1)/m`. This needs the translate
code on `A` to have minimum distance `> 2e`-ish. In all 71 sets tested
(`p ∈ {101,197,401,809}`, `m` up to 28, random, intervals and greedy
near-bicliques) the minimum distance between distinct words was 1, 2 or 3
(heuristically `O(log p/log m)` at `m ≤ √p`, since there are `p^2` pairs
among `2^m` words), and (J3) is vacuous for every `e ≥ 1`.

At `e = 0` the best of (J1), (J2) and the second moment gives
`m N_0^×/d ∈ [1.60, 1.94]` over all 71 sets (`≥ 1.90` once `m ≥ 20`),
against HP's `1` and a largest observed value `0.73`. For `e ≥ 1`, (J1) is
`≈ 2m^2/((m−2e)^2+m)` and the others `≈ 2/(1−2η)^2`; all tend to `2` as
`η → 0` and none is below `2 − O(1/m)` at `e = 0` (table `S3_list_decoding` in the
results file; at `p = 809`, `m = 28`, `e = 1`: (J1) 2.23, (J2) 2.32, SM
2.24, Hankel (Prop. 2.6) 2.15, truth for the greedy sets 0.35).

**Proposition 3.2 (the factor 2 is forced for Gram-only arguments).** Let
`q = p^2`, `A = B = F_p ⊂ F_q`, `χ_q` the quadratic character of `F_q`.
Then `A + B ⊆ Q_q ∪ {0}`, `m = n = √q`, `r = √q`, and

    mn − r = q − √q = (2p/(p+1))·(q−1)/2,

while `Σ_{x∈F_q} F_A(x)^2 = m(q−m)` with the `p` points of `F_p` carrying
the fraction `(p−1)/p` of it. Hence the second-moment count
`n(m−1)^2 ≤ m(q−m)` is attained within one point (`n ≤ p+1`, actual `p`),
and any inequality of the RHP shape that is valid over every finite field
of odd order, with `g = o(q)` at `m = √q`, has `f(0) ≥ 2p/(p+1) → 2`;
(J1)–(J3), Chung and Vinogradov are deduced from Gram identities that hold
in every such field, so none of them can do better.

*Remark on the brief's `g`.* At `m = √q` (and at `m = √p` over `F_p`) the
brief's example `g(m) = O(m^2)` is of order `p`, not `o(p)`, and would
absorb this example; the `F_{p²}` test and the target consequence
(which allows `|A|, |B|` near `√p`) both need `g(m) ≤ κ'p` with `κ'` small
at the relevant `m`, e.g. `g(m) = O(m^2)` only for `m ≤ p^{1/2−δ}`.

*Proof.* Every element of `F_p` is a square in `F_{p²}` (`x^{(p²−1)/2} =
(x^{p−1})^{(p+1)/2} = 1`). The Gram identity `Σ_x χ_q(x+a)χ_q(x+a') = −1`
holds in every `F_q`. ∎ (Checked for `p = 7,…,23`.)

So Johnson-type list decoding is the second moment and is sharp for the
subfield; improving the constant 2 needs a prime-field input, which in
HP is the degree bound `D ≤ p−1` (§1b). The Hankel polynomial of
Proposition 2.6 does use that input and still returns 2 at `e = 1`.

### 3.4. What the analogy does say (HEURISTIC)

HP's `F = Σ_k c_k(x+a_k)^{M−1}·Y_k − 1` is exactly an interpolation
polynomial of the linear-algebraic form `A_0(x) + Σ_k A_k(x)Y_k` used by
Guruswami for folded RS codes (arXiv:1106.0436, "a linear polynomial" in
the `y` variables), with the "folding" given by translation by the
arbitrary set `A` and the multiplicities supplied by the differential
equation `(x+a)Y_a' = dY_a`. In folded-RS decoding, errors are absorbed
because interpolation constraints are imposed at the known received
symbols; here the positions are unknown and the constraints must be
symbolic, and the end of §2.6 shows symbolic linear constraints cannot tolerate a
single error. In list-recovery language, RHP asks that a multiplicity code
(the `M`-jets of the degree-`D` polynomials `F_u`) have its coordinates in
the `V(m,e)`-element input lists `{jet of F_u : wt(u) ≤ e}` at only
`(1+o(1))d/m` positions, whereas the union bound gives `V(m,e)(d+m−1)/m`.

## 4. Witnesses (all in `results/stepanov_algebra_2026_09_26.json`)

| claim refuted | witness |
|---|---|
| least absolute residue of `R` equals `S` for `mn < p` | `p=101`, `A={0,1,2}`, `n=33`: `S=58`, residue `−43` (also `p=103,197,401`) |
| every sign-pattern fibre is `≤ d/m` | `p=11`, `A={0,1,3,4}`: a mixed fibre of size 2, `mμ = 8 = d+m−1 > d = 5` |
| `Φ_A − m = Σ_a(x+a)^d − m` vanishes to high order at complete points | order exactly 1 at 382 of 389 complete points, `≤ 3` at all |
| HP for multisets | `A = {0}×10`, `B = Q`, `p = 101`: all sums square, `mn = 500` |
| translate code on `A` has large minimum distance | minimum distance `≤ 3` in all 71 tested sets, e.g. `p=809`, `m=28` |
| field-uniform Gram/Johnson constant `< 2p/(p+1)` | `F_{p²}`, `A=B=F_p`, `p = 7,…,23` |

## 5. Open obligations

1. Whether the syndrome structure (Proposition 2.3(iii): at a point with
   `e_b ≤ e` the `M`-jet of `F` is a weight-`e_b` syndrome with *known*
   error values `−2c_k`) can be exploited non-determinantally to give
   `N_e^× ≤ (1+o_η(1))d/m`. The determinantal route gives
   `(e+1)d/(M−2e)`; the union over sign patterns gives `V(m,e)(d+m−1)/m`;
   neither is RHP. (OPEN)
2. Whether the fibres `μ(u)` for `u` in a Hamming ball can be
   simultaneously near `d/m` for many `u` (this is where a counterexample
   to RHP would have to live; `adversary`'s domain). (OPEN)
3. Whether the September 2026 prime-field capacity results (§3.2) have a
   formulation in which the translates have agreement above rate, e.g.
   after a change of code (not of roles). None of the encodings in §3.2
   do. (OPEN; the claims themselves are unrefereed.)
4. Proposition 2.5 is a routine extension of HP's argument; I found no
   statement of it in a brief search (the abstract of Kim–Yip–Yoo,
   arXiv:2309.09124, concerns product sets in shifted subgroups) and do not
   claim novelty.

## 6. Verifier

`experiments/stepanov_algebra_2026_09_26.py` (standard library + numpy)
checks every identity and inequality above by exact modular or integer
arithmetic (floating point only for the observed eigenvalue in (J2)).
Final run: 2,118,597 checks (including 1,773,962 enumerated polynomials
for Lemma 3.1 and 59,006 exhaustive sign-pattern sets), 0 failures,
about 2 minutes on a heavily loaded machine. Counts per identity and all
witnesses are in `results/stepanov_algebra_2026_09_26.json`.
