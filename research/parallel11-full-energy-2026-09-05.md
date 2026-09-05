# The complete spectral trace through a positive matrix energy

**Status: a source-dependent bound for the complete aggregate in an explicit
logarithmic depth range, checked locally.** Independent mathematical review
remains outstanding. The longer-depth estimate needed for the spectral edge,
the full Paley conjecture, and the official prize goal remain unproved.

Let p≡1 mod 4, let A₀ be a clique of size a≥1, and let C be its
common neighborhood with the anchors removed. Put m=|C|, b=2^a,
N=b−1, S_xy=χ(x−y), J=11ᵀ, and E=diag(1_C). Define

\[
 R=S/\sqrt p-J/p=2P-I,\qquad D=bE-I,\qquad
 T=RD/\sqrt N,\qquad H_j=\|T^j\|_F^2.                 \tag{1}
\]

Here ||·||_F is the Frobenius norm. R is an orthogonal involution;
D²=NI+(N−1)D. These are the exact matrices of the
[aggregate criterion](parallel2-spectral-transfer-2026-09-04.md):
T=W/(p√N). All rank-one J and anchor contributions are included.

## 1. An exact scalar recurrence and the even trace

Set K_j=||ET^j||_F², the energy in rows belonging to C. Orthogonality
of R gives

\[
 H_{j+1}=H_j/N+(N-N^{-1})K_j
        =H_j+\frac{N-1}{N}(bK_j-H_j).                  \tag{2}
\]

In particular H₀=p and H₁=mN+(p−m)/N. Since T is invertible,
H_j>0. With ρ_j=bK_j/H_j−1,

\[
 \log(H_j/p)=\sum_{\ell=0}^{j-1}
      \log\left(1+\frac{N-1}{N}\rho_\ell\right).        \tag{3}
\]

The factors are positive, between 1/N and N. The biases need not
be nonnegative, and neither monotonicity nor independence is assumed.
For every j≥1, Cauchy-Schwarz on matrix entries gives

\[
 |\operatorname{tr}T^{2j}|\le H_j.                    \tag{4}
\]

There is also an exact curvature identity:

\[
 H_{j+1}-2H_j+H_{j-1}
       =\frac{(N-1)^2}{N}\operatorname{tr}T^{2j}.       \tag{5}
\]

To prove it, use the orthogonal two-projection decomposition from
the earlier criterion. If λ is an eigenvalue of P_(C,C), set
c=√N/b and z=(λ−1/2)/c. Its two-dimensional block of T has
determinant one, trace 2z and squared Frobenius norm N+N⁻¹.
Writing its jth power using the Chebyshev polynomial U_(j−1) gives
squared norm 2+γ U_(j−1)(z)², where γ=(N−1)²/N. The residual
p−2m directions have squared norm N^(−j). Therefore

\[
 H_j=2m+\gamma\sum_{i=1}^m U_{j-1}(z_i)^2
                         +(p-2m)N^{-j}.               \tag{6}
\]

This includes j=0 with U_(−1)=0 and the degenerate endpoints of
the projection decomposition. The identity
U_j²−2U_(j−1)²+U_(j−2)²=2T_(2j), together with the residual
term, proves (5). For a=1, N=1 and T is orthogonal, so H_j=p
at every depth. This is a calibration case, not a new edge result.

The positive formulation does not make the open trace condition weaker.
Writing τ_j=tr(T^(2j)), there is a converse bound

\[
 H_j\le p+\gamma j^2\bigl(|\tau_j|/2+p\bigr).          \tag{6a}
\]

For |z|≤1, |U_(j−1)(z)|≤j. For |z|>1, put |z|=cosh t;
sinh(jt)/sinh t≤j cosh(jt), so U_(j−1)(z)²≤j²T_(2j)(z).
The trace formula and T_(2j)≥−1 on [−1,1] give
Σ_out T_(2j)≤|τ_j|/2+m_in. Substitution in (6), and 2m≤p,
proves (6a). Together, (4) and (6a) compare the two criteria with
explicit factors. At the required scales a=O(log p), j/log p→∞,
log(H_j+2p)/j and log(|τ_j|+2p)/j have the same vanishing
upper-bound requirement. Equation (14) below is a positive reformulation of
that difficulty, not a claim to have weakened it away.

## 2. The actual coefficient vector and its raw-rank budget

Let D_soft=Σ_(empty≠Z⊂A₀)D_Z, with D_Z(x)=∏_(z∈Z)χ(x−z).
The clique condition gives exactly

\[
 D=D_{\rm soft}-(b/2)E_{A_0}.                          \tag{7}
\]

Set B=(S/√p)D_soft/√N. Both ||T|| and ||B|| are at most √N.
In B^j(S/√p), every length-j word has coefficient one. This is
the particular coefficient vector needed here, rather than an
unrestricted isometry on a space that eventually has too many columns.

For y∉A₀, use the raw paths R_w(x,y) and principal traces of the
[word-aggregate proof](parallel10-word-aggregate-2026-09-05.md).
Write r_w for the generic raw rank. Its first-moment recurrence is

\[
 L_{j+1}=\lambda_a L_j+2^{a-1}N^j,\quad L_0=1,
 \qquad\lambda_a=2^{a-1}(a+1)-1,
\]

where L_j=Σ_(|w|=j)r_w. Hence, for a≥2,

\[
 L_j=\frac{a\lambda_a^j-N^j}{a-1},\qquad
 F_j=\frac{L_j^2}{N^j},\qquad
 \Theta_a=\frac{\lambda_a^2}{N},
 \quad F_j\le\left(\frac a{a-1}\right)^2\Theta_a^j.     \tag{8}
\]

For a=1, L_j=j+1 and F_j=(j+1)². The coefficient-specific budget
F_j is at most the sum of squared raw ranks by Cauchy-Schwarz.
For a=2, Θ₂=25/3, compared with the preceding unrestricted
second-moment rate (9+√65)/2.

On U_y=F_p\(A₀∪{y}), define the normalized sum
f_y=N^(−j/2)Σ_w(−1)^j p^(−j/2)R_w(x,y). The principal
off-diagonal correlations and diagonal normalization from pass ten,
applied before bounding this particular coefficient vector, give

\[
 \frac1p\sum_{x\in U_y}|f_y(x)|^2\le V_j,
 \quad \epsilon_j=\frac{aF_j+1}{\sqrt p},\quad
 g_j=\sqrt{F_j/p},\quad
 V_j=1+\epsilon_j+2g_j\sqrt{1+\epsilon_j}+g_j^2.         \tag{9}
\]

Indeed, if d_w is the principal rank, its error for this vector is
at most (a(Σd_w)²/N^j+1)/√p≤ε_j. The lower-weight raw correction
has normalized L² norm at most Σ(r_w−d_w)/√(pN^j)≤g_j.
No cancellation of those corrections is assumed.

## 3. All finite rows and the missing constant direction

At a finite point v let c_v=r_w−dim H^(−1)(K_w,v)+dim H^0(K_w,v).
The preceding exact conductor recurrence gives c_y=1 and 0≤c_v≤r_w
at each anchor. Also dim H^(−1)(K_w,v)≤r_w, by additivity bounds
over simple perverse constituents. Thus dim H^0(K_w,v)≤c_v.
The weight-gap sequence bounds both ordinary stalk traces by their
dimensions times p^(j/2): the principal object has no H^0, and
the error has perverse weights at most j. Consequently

\[
 |R_w(v,y)|\le(r_w+c_v)p^{j/2}\le2r_wp^{j/2}.           \tag{10}
\]

For generic y, each of the a+1 omitted x values therefore contributes
at most 4F_j/p to the squared norm of B^j(S/√p). Summing (9) over
generic columns and (10) over omitted rows costs at most
pV_j+4(a+1)F_j. The a exceptional columns y∈A₀ cost at most aN^j
by the operator norm. Finally S²/p=I−J/p gives

\[
 \|B^j\|_F^2=\|B^jS/\sqrt p\|_F^2
                       +\|B^j\mathbf1/\sqrt p\|_2^2.
\]

The last term is at most N^j. Therefore

\[
 \|B^j\|_F^2\le U_j:=pV_j+4(a+1)F_j+(a+1)N^j.        \tag{11}
\]

This step explicitly restores the constant direction lost by using S.
Neither that direction nor the exceptional fibers are dropped.

## 4. Restoring J and the anchor correction in the powers

From (1) and (7),

\[
 T-B=-\frac{b}{2\sqrt N}\frac S{\sqrt p}E_{A_0}
                       -\frac1{\sqrt N}\frac JpD.
\]

Their Frobenius norms imply

\[
 \|T-B\|_F\le d_0:=\frac b{2\sqrt N}\sqrt{a(1-1/p)}
        +\sqrt{\frac{mN^2+p-m}{pN}}.
\]

Telescoping the difference of powers, using ||T||,||B||≤√N,
gives ||T^j−B^j||_F≤j d₀N^((j−1)/2). Combining with (11) proves
the complete, nonasymptotic estimate

\[
 \boxed{|\operatorname{tr}T^{2j}|\le H_j
       \le\left(\sqrt{U_j}+j d_0N^{(j-1)/2}\right)^2.}  \tag{12}
\]

All constants are explicit and computable from p,a,m,j. No claim is
made that (12) is sharper than the trivial bound at every parameter.

## 5. What this proves, and the remaining depth requirement

Fix a≥2 and 0<ε<1/2. Uniformly over all a-anchor cliques, for

\[
 1\le j\le\frac{(1/2-\epsilon)\log p}{\log\Theta_a},
\]

equation (8) gives F_j=O_a(p^(1/2−ε)). Thus V_j=1+o(1),
U_j=p(1+o(1)), and j d₀N^((j−1)/2)=o(√p), since N<Θ_a
and d₀=O_a(1). The new complete-aggregate consequence is

\[
 \boxed{H_j\le(1+o(1))p,\qquad
  \left|\frac{\operatorname{tr}W^{2j}}
                    {[p\sqrt{b-1}]^{2j}}\right|
                   \le(1+o(1))p.}                    \tag{13}
\]

The result concerns the full trace, including J and anchors, at a
specified logarithmic range of even depths. The cohomological inputs
remain the locally checked source-dependent results of passes nine
and ten. This is not a proof at j/log p→∞.

Equation (3) isolates a concrete sufficient next estimate. At a possibly
growing anchor count, if j/log p→∞ and the actual trajectory satisfies

\[
 \sum_{\ell<j}\log\left(1+\frac{N-1}{N}\rho_\ell\right)
                    \le o(j),                         \tag{14}
\]

uniformly over the required anchor cliques, then (4) supplies the
earlier aggregate criterion with log(Bound+2p)=o(j). For
a=floor(σ log₂p), this gives its stated O(p^((1−σ)/2)) clique
consequence. Establishing (14) is still open. It uses actual row-energy
fractions, not a guessed random distribution or the unavailable
unrestricted coefficient isometry. Nor would the clique consequence
alone prove the arbitrary-two-set conjecture, uniform subgroup target,
or the quantitative bridge to the official Reed-Solomon prize.

## Verification and sources

The exact energy identities use only finite-dimensional real linear
algebra and the already verified Paley projection. The full bound
imports the rank, weight and correlation statements in
[pass nine](parallel9-all-degrees-2026-09-05.md) and
[pass ten](parallel10-word-aggregate-2026-09-05.md), with their
audited Katz and BBD sources. No additional external theorem is used.
The [primary discussion of the spectral depth issue](https://arxiv.org/html/2303.16475v1#S7)
was rechecked. No independent mathematical review, formal verification
or literature-priority claim is made.

The [companion verifier](../experiments/parallel11_full_energy_2026_09_05.py)
uses integer matrices and exact pairs u+v√p
to check (2), (4), (5), the complete anchor expansion, the soft/full
comparison, the missing constant direction, and the specific all-word
matrix sum. Squared pointwise bounds retain all diagonal and anchor
values. Rational upper square roots certify (12). An independent
two-dimensional Jordan-block case checks the polynomial j² behavior
at the spectral edge. The finite computations supplement the proof;
they do not establish (14) or extrapolate to untested primes.
The [results](../results/parallel11_full_energy_2026_09_05.json)
record 28 clique cases over five fields, 168 row-energy recurrences,
140 curvature identities, 168 trace/energy comparisons in each
direction, 168 complete upper bounds, 168 constant-direction checks,
672 boundary checks, 30 literal all-word matrix sums and 117,120
pointwise path bounds. A separate exact first-step case at p=10009
has a certified upper bound below 2p. Both positive and negative
row-energy biases occur in the checked cases.
