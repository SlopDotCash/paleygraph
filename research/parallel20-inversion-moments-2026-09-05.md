# Inversion removes short additive relations without discarding the input

Status: a proved transfer between restricted and unrestricted moment
bounds, with explicit finite hypotheses. The uniform moment bound itself
remains unproved. This strengthens the pass-19 reduction but proves no
new Paley cancellation exponent or official prize result. No novelty
claim is made for the elementary ingredients.

Let p be an odd prime, χ its quadratic character with χ(0)=0, and
C⊂F_p have size k≥1. For signs w_c∈{−1,1}, write

\[
F_{C,w}(x)=\sum_{c\in C}w_c\chi(x-c),\qquad
M_{2r}(C,w)=\sum_x|F_{C,w}(x)|^{2r}.
\]

Omit w when all signs are positive. A B_h set has no nontrivial
equal h-term sums, allowing repeated summands. Throughout h≥2 and p>h.

## 1. An explicit bound for bad inversion poles

For z∉C put D_z={1/(c−z):c∈C}. Let

\[
N_h(k)=\binom{k+h-1}{h},\qquad
Q_h(k)=k+(2h-2)\binom{N_h(k)}2.                         \tag{1}
\]

Then at most Q_h(k) elements z∈F_p either lie in C or make D_z
fail the B_h property. In particular, if p>Q_h(k), a good pole exists.

To prove the assertion, fix two distinct h-element multisets from C.
After cancellation their reciprocal-sum difference is

\[
R(z)=\sum_{c\in T}\frac{\nu_c}{c-z},\qquad
\nu_c\in\{-h,\ldots,h\}\setminus\{0\},\quad
\sum_c\nu_c=0,\quad |T|\le2h.                         \tag{2}
\]

Because p>h, every displayed coefficient remains nonzero in F_p.
Multiplying by ∏_{c∈T}(c−z) gives

\[
P(z)=\sum_{c\in T}\nu_c\prod_{d\in T\setminus\{c\}}(d-z).
                                                               \tag{3}
\]

This polynomial is nonzero: at any c∈T its value is
ν_c∏_{d≠c}(d−c)≠0. Its leading coefficient vanishes because
Σν_c=0, so deg P≤|T|−2≤2h−2. A nonzero polynomial has at most
its degree many roots. There are \(\binom{N_h(k)}2\) unordered
pairs of distinct multisets. A union bound and exclusion of C prove
(1), including pairs that share summands. No character-sum estimate
or probabilistic independence assumption enters this count.

The bound is sufficient, not a sharp count. Its hypothesis matters:
for C={0,1,2} in F_5, neither available pole produces a B_2 set.
This finite failure without a size condition does not concern Paley.

## 2. The exact character-moment change, including the missing row

Give d=1/(c−z) the sign w_d=χ(c−z). For x≠z put t=1/(x−z).
The elementary identity

\[
t-d=\frac{c-x}{(x-z)(c-z)}
\]

implies

\[
F_{D_z,w}(t)=\chi(-1)\chi(x-z)F_C(x).                  \tag{4}
\]

The remaining row t=0 has value χ(−1)k. Thus the exact finite-field
moment identity is

\[
\boxed{M_{2r}(D_z,w)=M_{2r}(C)-|F_C(z)|^{2r}+k^{2r}.}  \tag{5}
\]

Since |F_C(z)|≤k, (5) gives M_(2r)(C)≤M_(2r)(D_z,w).
Dropping the missing row before making the change of variables would
lose this exact correction. The signs in (4) work for both congruence
classes of odd primes; no χ(−1)=1 assumption is made.

## 3. Removing the signs using equal-size B_h completions

Use the [pass-19 extension bound](parallel19-relation-free-reduction-2026-09-05.md)

\[
f_h(s)=s+\sum_{j=1}^h s^{2h-j}.
\]

If p>f_h(2k−1), every B_h set D of size k can be enlarged to a
B_h set U of size 2k, by applying the same forbidden-point rule in
the whole field. Write U=D⊔E, |E|=k.

Let P,N partition D according to the positive and negative signs,
a=|P| and s=|P|−|N|=2a−k. Choose uniformly a (k−a)-element
subset J of E and an a-element subset K of E. Then
C_+=P⊔J and C_−=N⊔K are B_h sets of size exactly k. Pointwise
as functions on F_p,

\[
\boxed{w1_D=\mathbb E1_{C_+}-\mathbb E1_{C_-}
                     +\frac{s}{k}1_E.}                \tag{6}
\]

On E the first two expectations have difference (k−2a)/k=−s/k;
on D they already give w. The two random choices need not be coupled.
This identity also covers the all-positive and all-negative cases.

Suppose every B_h set V of size exactly k satisfies M_(2r)(V)≤U₀.
The triangle inequality in ℓ^(2r), applied to (6) after convolution
with χ, proves

\[
\boxed{M_{2r}(D,w)\le\left(2+\frac{|s|}{k}\right)^{2r}U_0
                              \le3^{2r}U_0.}           \tag{7}
\]

Only moments at the original cardinality k are used; no estimate
on smaller subsets or on size 2k is assumed. Size 2k is needed only
for the relation-free completion. Only the subset choices in (6) are
probability averages; the signs w remain coefficients.

## 4. A finite transfer with no partition remainder

Combining (1), (5) and (7) proves the following theorem.

**Theorem.** Let r≥1, h≥2, p>h be prime, and k≥1. If

\[
p>Q_h(k),\qquad p>f_h(2k-1),                            \tag{8}
\]

then a bound M_(2r)(V)≤U₀ uniform over B_h sets of size k implies

\[
\boxed{M_{2r}(C)\le3^{2r}U_0
                    \quad\text{for every }|C|=k.}      \tag{9}
\]

Restriction gives the converse bound with constant one. Thus, under
(8), restricted and unrestricted uniform moment bounds are equivalent
up to the displayed factor. This is a theorem about transferring a
bound, not a bound on either supremum without an input U₀.

An entirely explicit sufficient size threshold is available. Define

\[
K_h=\max\{2,\ 4h(h-1),\ (h+1)2^{2h-1}+1\}.             \tag{10}
\]

If k≥K_h and p≥k^(2h), both strict inequalities in (8) hold.
For the second,
f_h(2k−1)≤(h+1)(2k)^(2h−1)<k^(2h).
For the first,

\[
Q_h(k)\le k+\frac{h-1}{(h!)^2}(k+h-1)^{2h}.
\]

The coefficient (h−1)/(h!)² is at most 1/4 for h≥2. Also
\((1+(h-1)/k)^{2h}\le2\) when k≥4h(h−1): expand by the
binomial theorem, bound each binomial coefficient by (2h)^j,
and compare with the geometric series whose ratio is at most 1/2.
Consequently Q_h(k)≤k+k^(2h)/2<k^(2h), since k≥2 and h≥2.

## 5. Stronger scope for the restricted SS criterion

Fix r≥3 and 2≤h≤⌊(r+1)/2⌋. At k=⌊p^(1/(r+1))⌋, for every
prime p≥K_h^(r+1) we have k≥K_h and p≥k^(r+1)≥k^(2h).
The finite transfer therefore applies on the precise SS slice.

It follows that the sufficient restricted hypothesis can now use an
unbounded sequence r_j with

\[
2\le h_j\le\left\lfloor\frac{r_j+1}{2}\right\rfloor,
\quad\beta_j\ge0,\quad\beta_j\longrightarrow0,
\]

and demand only

\[
M_{2r_j}(V)\le C_j p^{1+\beta_j}k_j(p)^{r_j}
\quad\text{for }B_{h_j}\text{ sets }|V|=k_j(p).
\tag{SS-B*}
\]

Equation (9) supplies the unrestricted SS estimate with constant
3^(2r_j)C_j. The [original single-size sampling transfer](parallel5-classical-2026-09-04.md)
then proves full classical Paley. For completeness, given |B|≥k,
uniform sampling of k-subsets B′⊂B and convexity imply
M_(2r)(B)≤(|B|/k)^(2r) E M_(2r)(B′). Hölder consequently gives

\[
\frac{|S(A,B)|}{|A||B|}
\le 3\left(\frac{C_jp^{1+\beta_j}}{|A|k_j^{r_j}}\right)^{1/(2r_j)}.
                                                               \tag{11}
\]

For |A|,|B|>p^ε choose one j with θ_j+β_j<ε/2, where
θ_j=1/(r_j+1). Then k_j≤|B| and k_j≥p^θ_j/2 eventually.
The right side of (11) is at most
3(C_j2^r_j)^(1/(2r_j))p^{−ε/(4r_j)}; absorbing the fixed
constant into the threshold permits δ=ε/(8r_j)>0.

Every moment order, relation order and constant is fixed before p
grows. The possibly enormous explicit K_h threshold is not an
estimate at r=r(p). There is no remaining h_j/r_j→0 requirement.

The parameter boundary 2h=r+1 is included. As k grows there,
Q_h(k)/k^(2h) tends to (h−1)/(h!)²≤1/4, while the completion
bound has degree 2h−1. The factor from counting unordered multisets
matters at that boundary; replacing it by a constant greater than
one times k^(2h) would not justify the same conclusion.

For example, an unproved sixth-moment SS bound only for Sidon sets
at k=⌊p^(1/4)⌋ would now give the unrestricted sixth-moment bound
and hence Paley for ε>1/4. A fourth-moment Sidon bound at
k=⌊p^(1/3)⌋ lies outside (8)'s proved asymptotic range and is not
claimed to transfer. No such new character-moment estimate has been
proved in this pass.

## 6. Validation and the remaining estimate

The [exact verifier](../experiments/parallel20_inversion_moments_2026_09_05.py)
enumerates multiset relations and their nonzero numerator polynomials,
checks every pole in bounded cases, verifies (4)–(5) at every field
point, constructs B_h completions, and tests (6)–(7) with exact integer
moments. It also checks rational polynomial threshold certificates.
These finite checks complement the proofs and do not establish SS-B*.

All new arguments are elementary and supplied here; no new external
theorem is imported. The earlier projective graph calculations
likewise retained their missing-point terms, but the setwise signed
moment and completion transfers above are separate arguments.
Independent mathematical review remains outstanding.

The next required input is an actual uniform upper bound in SS-B*,
not merely the absence of short additive relations. Classical Paley,
the uniform subgroup square-root target, the spectral edge and the
quantitative bridge to the official Reed–Solomon prize remain open.
