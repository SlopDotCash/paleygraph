# Distinct word sheaves and a quadratic aggregate over all label sequences

**Status: source-dependent proof, checked locally; independent review remains
outstanding.** This extends the individual-word theorem to simultaneous
quadratic bounds for every coefficient vector on a complete word family.
Its useful depth is still only a bounded multiple of log p. It does not
prove the extreme-eigenvalue aggregate or the full Paley conjecture.

Fix p≡1 mod 4, distinct anchors A⊂F_p of size a≥1, and y∉A. Use
the principal sheaves F_w, raw perverse objects K_w, and raw paths R_w
from the [ninth-pass proof](parallel9-all-degrees-2026-09-05.md).
A word w=(T₁,…,T_m) has nonempty T_i⊂A. Set U=F_p\(A∪{y}),
d_w=rank F_w, and r_w=generic rank of H^(−1)(K_w).

The empty word has F_empty=L_χ(x−y), d_empty=r_empty=1.
For every word, the principal sheaf is geometrically irreducible,
tame, and pure of usual weight |w|. Its finite singularities lie
in A∪{y}. All statements below are uniform in the anchor positions.

## 1. Recovering the word, and arithmetic self-duality

At each v∈A let t_v(F) be the number of trivial Jordan blocks minus
the number of quadratic Jordan blocks. The preceding rank argument
proves t_v(F_w)>0 for every word, including the empty word, where
it is one. Apply geometric MC_χ to F_w for w nonempty. Involution
recovers G=ME(F_prefix⊗L_dT), where T is the last label. Consequently

\[
 T=\{v\in A:t_v(G)<0\}.                                      \tag{1}
\]

No difference is zero. Twisting G by the recovered label gives
F_prefix. Repetition recovers the entire word, stopping at rank one;
all nonempty words have rank at least two. Hence

\[
 F_u\cong F_v\text{ geometrically}\quad\Longleftrightarrow\quad u=v.
                                                               \tag{2}
\]

This uses middle-convolution inversion, not an assumption that local
monodromy by itself classifies arbitrary sheaves.

The arithmetic normalization also satisfies

\[
 F_w^\vee\cong F_w(|w|).                                      \tag{3}
\]

Indeed, Verdier duality exchanges compact and ordinary convolution,
so it commutes with the middle image. It commutes with intermediate
extension and sends L_χ[1] to L_χ[1](1). Starting with the
self-dual rank-one quadratic sheaf, induction gives
D(F_w[1])≅F_w[1](|w|+1), which is (3). This is an arithmetic
statement with the displayed Tate twist. On U, purity and (3)
make p^(−|w|/2)Tr(Frob_x|F_w) real: the eigenvalue multisets are
stable under α↦p^|w|/α, which equals complex conjugation there.

## 2. A simultaneous Gram bound for all principal words

For any finite set W of distinct words, including words of different
lengths, put

\[
 h_w(x)=p^{-|w|/2}\operatorname{Tr}(\operatorname{Frob}_x|F_w),
 \quad D_W=\sum_{w\in W}d_w^2,
 \quad \epsilon_W=\frac{aD_W+1}{\sqrt p}.
\]

For every complex coefficient vector z,

\[
 \left|\frac1p\sum_{x\in U}\left|\sum_w z_wh_w(x)\right|^2
       -\|z\|_2^2\right|\le\epsilon_W\|z\|_2^2.          \tag{4}
\]

For u≠v, (2) and (3) exclude both invariants and coinvariants of
F_u⊗F_v on U. The tame tensor has rank d_ud_v and a+2 punctures,
so H_c¹ has dimension a d_ud_v. Its trace sum, after normalization,
has magnitude at most a d_ud_v√p. On the diagonal, geometric
irreducibility leaves exactly one invariant line. By (3), its
H_c² eigenvalue is p^(|w|+1); H_c¹ has dimension a d_w²+1.
Thus the normalized Gram entries differ from the identity by at
most (a d_ud_v+1_(u=v))/√p. Cauchy-Schwarz against the rank vector
(d_w) proves (4). No probabilistic independence is assumed.

If B⊂F_p\(A∪{y}) is nonempty, s=|B|, then similarly

\[
 \left\|\left(\frac1p\sum_{x\in U}
     \prod_{b\in B}\chi(x-b)h_u(x)h_v(x)\right)_{u,v}\right\|
 \le\frac{(a+s)D_W}{\sqrt p}.                            \tag{5}
\]

At any b∈B, the untwisted factors are lisse and the new twist has
scalar quadratic inertia, excluding H_c² even on the diagonal.
There are a+s+2 punctures. Extension by zero at B supplies exactly
the displayed zero character values. This proves (5).

## 3. Exact raw-rank accounting with all boundary constituents

For a tame perverse object K on this stratification, define

\[
 c_v(K)=r(K)-\dim H^{-1}(K_v)+\dim H^0(K_v),
 \qquad v\in A\cup\{y\}.
\]

These nonnegative integers and the generic rank are additive on
perverse exact sequences. For a simple middle extension, c_v is its
invariant codimension; a one-dimensional punctual object contributes
one at its support. This is a stalk Euler characteristic, not a
claim that the two ordinary stalk dimensions are separately additive.

Compact convolution with L_χ[1] preserves every finite c_v and
has generic rank Σ_v c_v. This follows for ordinary irreducible
middle extensions by Katz's finite local-monodromy and rank formulas:
the constant infinity kernel restores the subtracted infinity rank.
The exceptional cases are explicit. A constant object goes to zero;
a quadratic Kummer translate goes to a punctual middle object plus
a rank-one constant kernel; a punctual object goes to a Kummer
translate. They obey the same formulas. Additivity proves the claim
for all the raw objects used here.

Extension by zero with mask T sets c_v to r(K) at v∈T and leaves
the other c_v unchanged. Since initially c_y=1 and c_v=0 for v∈A,
the exact recurrence is particularly simple:

\[
 c_y=1,\quad r=1+\sum_{v\in A}c_v;\qquad
 c'_v=\begin{cases}r,&v\in T,\\c_v,&v\notin T,\end{cases}
 \quad r'=1+\sum_{v\in A}c'_v.                           \tag{6}
\]

In particular r_w≥d_w. These r_w include the generic boundary
contributions; they need not equal the principal ranks. For a=2,
the first three labels have principal ranks 2,2,2 but raw ranks
2,2,3. Discarding the last constant correction would give a wrong
rank recurrence.

For all length-m words write W_m=(2^a−1)^m, R_m=Σ_w r_w,
Q_m=Σ_w r_w², and Z_m=Σ_wΣ_v c_v(w)². For a≥2, set
N=2^a−1, h=2^(a−1), and A₂=2^(a−2)(a²+3a−1)−1. Then

\[
 \begin{aligned}
 W_{m+1}&=NW_m,\\
 R_{m+1}&=(2^{a-1}(a+1)-1)R_m+hW_m,\\
 Q_{m+1}&=A_2Q_m+h(a+2)R_m+(h/2)Z_m+(h/2)W_m,\\
 Z_{m+1}&=ahQ_m+(h-1)Z_m,
 \end{aligned}                                            \tag{7}
\]

with (W₀,R₀,Q₀,Z₀)=(1,1,1,0). To verify the square recurrence,
write r'=r+Σ_(v∈T)(r−c_v). The sum of r−c_v is (a−1)r+1,
and the sum of their squares is (a−2)r²+2r+Σc_v². Summing
over subsets, and then deleting the empty subset, gives (7).
For a=1 there is a single word and r_w=d_w=m+1, Q_m=(m+1)².

For fixed a≥2, Q_m=Θ_a(Λ_a^m), where

\[
 \Lambda_a=2^{a-3}\left(a^2+3a+1+
     \sqrt{(a^2+3a-3)^2+8a}\right)-1.                    \tag{8}
\]

This is the Perron eigenvalue of the positive 2×2 block acting on
(Q,Z) in (7). It exceeds the eigenvalues driving R and W; the
nonnegative forcing and the positive initial Q component give the
stated upper and lower comparisons. For a=2,
Λ₂=(9+√65)/2≈8.531, and Q₀,Q₁,Q₂,Q₃=1,17,199,2001.
The exact recurrence supplies computable constants at every depth.

## 4. The quadratic bound for the original character paths

Define the actual raw trace function on U by

\[
 q_w(x)=(-1)^{|w|}p^{-|w|/2}R_w(x,y).
\]

The sign follows from Tr(K_w)=(-1)^(|w|+1)R_w and
Tr(F_w[1])=−Tr(F_w). The ninth-pass weight-gap sequence gives

\[
 |q_w(x)-h_w(x)|\le(r_w-d_w)/\sqrt p.                   \tag{9}
\]

Let L_W=Σ_w(r_w−d_w)² and γ_W=√(L_W/p), so L_W≤Q_m
for the complete length-m family. For the norm
||f||_U²=p^(−1)Σ_U|f(x)|², (9) gives
||Σz_w(q_w−h_w)||_U≤γ_W||z||₂. Combining with (4) proves

\[
 \left|\frac1p\sum_U\left|\sum_w z_wq_w(x)\right|^2
       -\|z\|_2^2\right|
 \le\left(\epsilon_W+2\gamma_W\sqrt{1+\epsilon_W}
                         +\gamma_W^2\right)\|z\|_2^2.   \tag{10}
\]

This is a signed quadratic aggregate for the original matrix paths,
valid simultaneously for every complex coefficient vector. It
does not replace them by their principal terms.

There is also uniform restriction to any fresh adjacency cell. For
B as above, s≥1, choose arbitrary signs σ_b and set
C={x∈U\B:χ(x−b)=σ_b for all b∈B}, β=2^(−s). Expanding the
indicator and using (4)-(5), then removing the soft half-valued
boundary at B, gives the principal Gram error

\[
 E_{W,s}=\frac{(a+s/2)D_W+\beta}{\sqrt p}
                     +\frac{sD_W}{2p}.                 \tag{11}
\]

Indeed, each point of B contributes at most D_W||z||²/2 before
division by p, since |h_w|≤d_w there. For the raw functions the
corresponding estimate is

\[
 \left|\frac1p\sum_C|\sum_wz_wq_w|^2-\beta\|z\|^2\right|
 \le\left(E_{W,s}+2\gamma_W\sqrt{\beta+E_{W,s}}
                              +\gamma_W^2\right)\|z\|^2. \tag{12}
\]

## 5. Useful range and unresolved spectral step

For fixed a≥2, the complete length-m family in (10) is an
asymptotic isometry whenever Λ_a^m=o(√p). For example,
m≤(1/2−ε)log(p)/log(Λ_a) works for any fixed ε>0. The relative
cell error in (12) tends to zero if
2^s(a+s)Λ_a^m=o(√p); polynomial factors in s must be retained.
This quantifies control of all label patterns and their next-anchor
weights simultaneously, instead of a span indexed only by one rank.

This does not reach the [spectral criterion](parallel2-spectral-transfer-2026-09-04.md).
That criterion needs the full cyclic trace, including its J terms
and exceptional anchor fibers, at lengths with k/log p→∞. The
present coefficient-to-function map is controlled only in the
displayed shorter range. Nor is there a proof that extremal
eigenvectors lie in, or are approximated by, this word span.
Geometric inequivalence yields correlation bounds, not independent
random signs or unrestricted cancellation of Frobenius traces.
The full Paley, subgroup and official prize targets remain open.

There is a stronger, elementary obstruction to extending (10) without
restricting the coefficient vectors. The evaluation map has (2^a−1)^m
columns and only p−a−1 rows. If

\[
 (2^a-1)^m>p-a-1,                                      \tag{13}
\]

it has a nonzero kernel, for either the principal or raw functions.
For a vector in that kernel, the left side of (10) is exactly ||z||².
Thus an error less than one is impossible in this range, regardless
of any improvement in the sheaf estimates. In particular a≥2 fixed
and m/log p→∞ cannot be handled by unrestricted coefficient isometry.
This does not exclude cancellation for the specific coefficients in
the cyclic spectral trace. It shows why that structure, or a suitable
quotient of the word space, is a necessary next ingredient for this
approach. The verifier records an exact integer raw-path null vector
for p=13, A={0,1}, y=2, m=3, where 27 columns map to 10 points.

## Sources and finite verification

The source-dependent inputs are those audited in the ninth pass.
Katz, [Rigid Local Systems](../sources/katz-rigid-local-systems.pdf),
Section 2.5.1 (complete PDF page 49, newly rendered and inspected)
gives convolution duality; Section 2.6.2 identifies the middle image.
The earlier inspected inversion, local-monodromy and compact/middle
exact sequences supply (1), (3) and (6). The previous BBD weight
and Katz-GKM curve-cohomology inputs supply (4), (5) and (9).
No new external theorem or claim of literature priority is assumed.

The [companion verifier](../experiments/parallel10_word_aggregate_2026_09_05.py)
checks backwards recovery of finite words,
distinct local signatures, raw ranks against full numerical inventories
of simple-constituent local types, and (7) against exhaustive subset
updates. It also evaluates the original character paths with integer
circular convolution, retaining every zero mask, and verifies Gram
and cell bounds by exact positive-definiteness certificates. It does
not computationally prove geometric inequivalence, duality or weights.
Independent mathematical review remains outstanding.
The [results](../results/parallel10_word_aggregate_2026_09_05.json)
record 7,541 backwards recoveries and distinct signatures, 773 raw
constituent-inventory comparisons, 114 integer convolutions, 573
literal path entries, 172 zero-mask rows and 70 exact positive-definiteness
certificates for 35 Gram or cell bounds. Six of those bounds have
relative error below one, so the checks include nonvacuous cases.
Small-field tests with larger errors are not presented as evidence
of asymptotic cancellation. The separate null vector certifies (13)
for an actual raw-path family, without a cohomological assumption.
