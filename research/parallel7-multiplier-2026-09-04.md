# A signed Fourier multiplier improves the subgroup exponent to 8/9

**Status: a stronger maximum-period bound on the existing density-one
quartic prime class is proved. The uniform thin-subgroup target and
the full Paley conjecture remain open.** This argument removes the
outer Cauchy–Schwarz square from the pass-six amplification without
requiring parity or the absence of zero triples. It uses a signed
Fourier multiplier on one side and a probability measure on the other.
No novelty, best-known-bound, or formalization claim is made.

For N dyadic, the same sixth-energy class K_N in
[pass five](parallel5-subgroup-2026-09-04.md), among
P_N={p prime:N⁴/4≤p≤N⁴, p≡1 mod N}, satisfies

\[
 \boxed{M(H_N)\le(17+\log N)^{1/9}N^{8/9},\qquad
         |K_N|/|P_N|=1-O(1/\log N).}                    \tag{1}
\]

This improves the preceding exponent 17/18 to 8/9 on that class.
The target exponent is still 1/2+o(1) at every eligible prime.

## 1. The multiplier and its exact origin contribution

Let H≤F_p* have size n≥2; no assumption that −1∈H is needed.
Write e_p(t)=exp(2πit/p), η(t)=Σ_(h∈H)e_p(th), and θ(t)=η(t)/n.
For a real r≥1 define the spectral energy

\[
 \mathcal E_r=\frac1p\sum_{t\in F_p}|\eta(t)|^{2r},
 \qquad\widetilde{\mathcal E}_r=\mathcal E_r-n^{2r}/p.
\]

For integer r this is the usual r-term additive energy E_r. For real
r the displayed spectral formula is the definition, not a tuple count.
Let μ_s be the probability distribution of the sum of s independent
uniform elements of H, where s≥1 is an integer. Put u_s=μ_s(0).
Then 0≤u_s≤1/n and Σ_x μ_s(x)e_p(tx)=θ(t)^s.

Define the possibly signed function

\[
 \alpha_r(x)=\frac1p\sum_t|\theta(t)|^r e_p(-tx),
 \qquad a_r=\alpha_r(0)=\frac1p\sum_t|\theta(t)|^r.
 \tag{2}
\]

Since |θ(−t)|=|θ(t)|, α_r is real and even. Fourier inversion and
Parseval give exactly

\[
 \sum_x\alpha_r(x)=1,\quad
 \sum_x\alpha_r(x)e_p(tx)=|\theta(t)|^r,\quad
 \sum_x\alpha_r(x)^2=\mathcal E_r/n^{2r}.               \tag{3}
\]

Also 0≤a_r≤1. For r≥2, |θ|≤1 and orthogonality imply

\[
 0\le a_r\le\frac1p\sum_t|\theta(t)|^2=1/n.             \tag{4}
\]

These facts do not assert that α_r(x)≥0. This distinction is essential:
the new function is used in the bilinear estimate, where signed weights
are allowed, and never as the measure in Jensen's inequality.

On Ω=F_p*, set

\[
 A_r(x)=\alpha_r(x)-(1-a_r)/(p-1),\quad
 B_s(x)=\mu_s(x)-(1-u_s)/(p-1).
\]

Both functions have sum zero on Ω. Their squared norms are

\[
 \begin{split}
 V_r^\alpha&=\sum_\Omega A_r(x)^2
   =\frac{\widetilde{\mathcal E}_r}{n^{2r}}
      -\frac p{p-1}(a_r-1/p)^2,\\
 V_s^\mu&=\sum_\Omega B_s(x)^2
   =\frac{\widetilde E_s}{n^{2s}}
      -\frac p{p-1}(u_s-1/p)^2.
 \end{split}                                           \tag{5}
\]

In particular both are nonnegative and no principal frequency has
been hidden in the centered quantities.

## 2. A mixed-order theorem with no outer square

**Theorem.** For every a≠0, Δ=|η(a)|/n, every real r≥1 and integer
s≥1,

\[
 \boxed{\Delta^{rs}\le
 a_r+u_s-a_ru_s-\frac{(1-a_r)(1-u_s)}{p-1}
                    +\sqrt{pV_r^\alpha V_s^\mu}.}       \tag{6}
\]

For r≥2 this implies the simpler bound

\[
 \boxed{\Delta^{rs}\le\frac2n+
   \frac{\sqrt{p\,\widetilde{\mathcal E}_r\widetilde E_s}}
        {n^{r+s}}
 \le\frac2n+\frac{\sqrt{p\,\mathcal E_r E_s}}{n^{r+s}}.} \tag{7}
\]

**Proof.** Subgroup invariance gives

\[
 \sum_y\mu_s(y)\theta(ay)
 =\frac1n\sum_{h\in H}\theta(ah)^s=\theta(a)^s.
 \tag{8}
\]

The function z↦|z|^r on C is convex for r≥1. Jensen with respect to
the nonnegative probability measure μ_s therefore proves

\[
 \begin{split}
 \Delta^{rs}
 &\le\sum_y\mu_s(y)|\theta(ay)|^r\\
 &=\sum_{x,y}\alpha_r(x)\mu_s(y)e_p(axy)=:T.
 \end{split}                                           \tag{9}
\]

The last identity is Fourier inversion. T is real and nonnegative
because of its middle expression, irrespective of the signs of α_r.

The terms with x=0 or y=0 sum exactly to
a_r+u_s−a_ru_s, since both full functions have total mass one. On
Ω×Ω, expand the constants from (5). The mixed constant/zero-sum
terms vanish because Σ_(y≠0)e_p(axy)=−1 for x≠0. The constant term
is −(1−a_r)(1−u_s)/(p−1). Thus

\[
 T=a_r+u_s-a_ru_s-\frac{(1-a_r)(1-u_s)}{p-1}
                  +\sum_{x,y\in\Omega}A_r(x)B_s(y)e_p(axy).
 \tag{10}
\]

For arbitrary complex functions U,V, Cauchy–Schwarz and additive
orthogonality give

\[
 \left|\sum_{x,y}U(x)V(y)e_p(axy)\right|
                  \le\sqrt p\,\|U\|_2\|V\|_2.          \tag{11}
\]

Apply this to A_r,B_s extended by zero. Taking the real part in (10)
and using (9) proves (6). Now (4), u_s≤1/n, and (5) imply (7). ∎

The difference from the preceding probability-only argument is the
Fourier transform |θ|^r. One Jensen step yields Δ^(rs) directly,
without squaring an average of |θ(axy)|. Positivity is needed only
for μ_s; assuming it for α_r would be an invalid extra step.

## 3. Sixth energy now gives exponent 8/9

At r=s=3, formula (7) becomes

\[
 \Delta^9\le2/n+\sqrt p\,E_3/n^6.                       \tag{12}
\]

Consequently, for every multiplicative subgroup with p≤n⁴ and
E₃≤Bn³,

\[
 \boxed{M(H)\le(B+2)^{1/9}n^{8/9}.}                    \tag{13}
\]

This conditional implication has explicit constants and no eventual
size restriction. It can be weaker than M≤n for small parameters.
It does not require fourth-energy control or zero-triple exclusion.

The pass-five prime budget gives E₃(H_N)≤(15+log N)N³ on K_N.
Substitution proves (1). The asymptotic density assertion imports
exactly the previously checked
[Thorner–Zaman prime count](https://arxiv.org/html/2108.10878v2#S3.SS1).
No additional external analytic theorem is needed for (2)–(13).

For every integer j≥3, positivity of the nonprincipal Fourier moments
then gives, on this class with B=15+log N,

\[
 \widetilde E_j\le M^{2j-6}\widetilde E_3
 \le B(B+2)^{2(j-3)/9}N^{(16j-21)/9}.                   \tag{14}
\]

In particular the centered eighth-energy exponent improves from
44/9 to 43/9, with logarithmic factor (log N)^(11/9). It is still
above the N⁴ scale required by the next proposed arithmetic budget.
The [parallel subgroup analysis](parallel7-subgroup-2026-09-04.md)
examines the limits of reinserting this feedback and the failure of
treating α₃ as a probability measure.

## 4. The payoff of a genuinely stronger eighth moment

For r=4,s=1, subgroup invariance makes the Jensen expression exact:
Σ_(h∈H)|θ(ah)|⁴/n=Δ⁴. Here u₁=0, a₄=E₂/n⁴, and

\[
 a_4-\frac{1-a_4}{p-1}
                 =\frac p{p-1}\frac{\widetilde E_2}{n^4}.
\]

Suppose p≤n⁴ and \(\widetilde E_4\le Dn^4\), D≥1. Cauchy–Schwarz
over nonzero frequencies gives \(\widetilde E_2\le\sqrt D\,n^2\).
Using V₁^μ≤1/n and V₄^α≤\(\widetilde E_4/n^8\) in (6),

\[
 \Delta^4\le\frac p{p-1}\frac{\sqrt D}{n^2}
                     +\sqrt D\,n^{-1/2}
                 \le2\sqrt D\,n^{-1/2}.                \tag{15}
\]

The last inequality holds for p≥3,n≥2, since p/(p−1)≤3/2 and
(3/2)n^(−3/2)≤1. Hence the explicit conditional consequence is

\[
 \boxed{M(H)\le2^{1/4}D^{1/8}n^{7/8}.}                 \tag{16}
\]

Thus the unproved averaged centered eighth budget would have a
stronger spectral payoff than recorded in pass six. Equations
(14)–(16) do not supply that budget themselves.

## Verification and remaining scope

The [verifier](../experiments/parallel7_multiplier_2026_09_04.py)
uses exact rational intervals for period magnitudes and Fourier
coefficients, exact integer additive energies, and rational squared
norms. It checks the centered gate and its coarser consequence for
fixed even and odd subgroup orders, and checks the explicit exponent
and constant calculations. Recorded input hashes accompany the
[results](../results/parallel7_multiplier_2026_09_04.json).

These are finite checks supporting an ordinary proof. They do not
prove prime-density asymptotics or uniform square-root cancellation.
The exceptional primes, logarithmic-depth centered moments, full
two-set Paley conjecture and exact prize reduction remain open.
