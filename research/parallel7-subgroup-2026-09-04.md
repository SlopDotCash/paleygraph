# Signed multipliers: a positivity obstruction and continuous-order audit

**Status: the main lane's multiplier argument is independently verified.
Its 8/9 exponent supersedes this lane's intermediate parity bound.
The averaged centered eighth-energy estimate remains unproved.**

The substantive safeguard here is an exact counterexample to treating
the absolute-Fourier multiplier as a probability measure, even on the
multiplicative quotient and even in the quartic prime window.
The second result is an optimization over every real multiplier order
r≥1 and every integer walk order s≥1: the present energy estimates,
inserted in the centered bilinear norm estimate, have best saving 1/9.
This is a limitation of that specified substitution, not a lower
bound for subgroup periods or a limit on methods retaining more information.

## Independent audit of the multiplier argument

Let H≤F_p* have size n≥2, write
\[
 \eta(t)=\sum_{h\in H}e_p(th),\quad
 \theta(t)=\eta(t)/n,\quad
 \delta=\max_{a\ne0}|\theta(a)|,\quad M=n\delta,
\]
and let μ_s be the distribution of a sum of s independent uniform H
elements. For a real r≥1 define
\[
 \alpha_r(x)=\frac1p\sum_{t\in F_p}|\theta(t)|^r e_p(-tx).
\]
This function is real because |\theta(-t)|=|\theta(t)|. It can have
negative entries. Nevertheless,
\[
 \sum_x\alpha_r(x)=1,\qquad
 a_r:=\alpha_r(0)\in[1/p,1],\qquad
 \|\alpha_r\|_2^2=\frac1p\sum_t|\theta(t)|^{2r}.
\]
For r≥2, |\theta|≤1 and orthogonality give a_r≤1/n.

Fourier inversion and subgroup invariance imply
\[
 \begin{aligned}
 \sum_{x,y}\alpha_r(x)\mu_s(y)e_p(axy)
 &=\sum_y\mu_s(y)|\theta(ay)|^r\\
 &\ge\left|\sum_y\mu_s(y)\theta(ay)\right|^r
 =|\theta(a)|^{rs}.                                      \tag{1}
 \end{aligned}
\]
Convexity of |z|^r holds on C for every real r≥1. The only measure
used in this Jensen step is the nonnegative μ_s; no positivity of
α_r is assumed. No condition −1∈H is needed.

Put u_s=μ_s(0), and center both functions on F_p*. Let
\[
 W_r=\sum_{x\ne0}\left(\alpha_r(x)-\frac{1-a_r}{p-1}\right)^2,
 \quad
 V_s=\sum_{x\ne0}\left(\mu_s(x)-\frac{1-u_s}{p-1}\right)^2.
\]
If
\[
 \widetilde{\mathcal E}_r
       =\frac1p\sum_{t\ne0}|\eta(t)|^{2r},
\]
then exact expansion gives
\[
 W_r=\frac{\widetilde{\mathcal E}_r}{n^{2r}}
       -\frac p{p-1}(a_r-1/p)^2,\qquad
 V_s=\frac{\widetilde E_s}{n^{2s}}
       -\frac p{p-1}(u_s-1/p)^2.                          \tag{2}
\]
The weighted bilinear estimate from
[pass 6](parallel6-subgroup-2026-09-04.md), valid for signed weights,
now proves
\[
 |\theta(a)|^{rs}
 \le C_{r,s}+\sqrt{pW_rV_s},\qquad
 C_{r,s}=a_r+u_s-a_ru_s-\frac{(1-a_r)(1-u_s)}{p-1}.        \tag{3}
\]
The constant has the useful equivalent form
\[
 C_{r,s}
 =\frac{p(a_r-1/p)}{p-1}
    +\frac{p(1-a_r)u_s}{p-1}\ge0.                        \tag{4}
\]
Thus discarding the exact constant is not a hidden cancellation.
The proof of (3) needs no positivity in the bilinear Cauchy–Schwarz
step; the centered weights have zero sum, and orthogonality bounds
their bilinear sum by the product of their L2 norms times √p.

For r=s=3, p≤n⁴, and E₃(H)≤Bn³ with B≥1, the simple bounds
a₃,u₃≤1/n and (2) give
\[
 \delta^9\le(B+2)/n,\qquad
 M\le(B+2)^{1/9}n^{8/9}.                                \tag{5}
\]
This confirms the main-lane conclusion on the sixth-energy class,
without a triple-free or fourth-energy restriction.

## A multiplier is not another probability measure

Take H={±1} in F₁₃. Then θ(t)=cos(2πt/13), and
\[
 \alpha_3(6)
 =\frac1{13}\sum_{t=0}^{12}
       |\cos(2\pi t/13)|^3\cos(12\pi t/13)
 \in(-0.004014894,-0.004014892).                           \tag{6}
\]
The companion checker proves a much narrower rational interval.
Since α₃ is H-invariant, α₃(7)=α₃(6). For the nonnegative,
H-invariant function g=1_{6H}, let b=Σ_x α₃(x)g(x). Then b<0 and
\[
 \sum_x\alpha_3(x)g(x)^2=b<b^2
   =\left(\sum_x\alpha_3(x)g(x)\right)^2.                 \tag{7}
\]
This explicitly violates the convex-square Jensen inequality with
α₃ as the averaging weights. It cannot be used as a probability
measure in another Jensen stage.

A second fixed quartic example is the order-four subgroup
H={1,27,46,72} in F₇₃:
\[
 \alpha_3(5)\in(-0.001688700,-0.001688698).
\]
Its entire H-coset has negative mass, giving the same failure with
an H-invariant indicator. Both examples satisfy n⁴/4≤p≤n⁴.

These examples refute positivity and generic signed-weight Jensen.
They do **not** refute (1)–(5), which use positive μ_s, nor prove that
a specially structured double-multiplier period inequality is false.
Establishing such an inequality would require additional structure,
not another application of Jensen to α_s.

## Every real multiplier order with the present moment feedback

Assume p≤n⁴ and E₃≤Bn³, with B=n^{o(1)}. The following is a
precise audit of replacing the variances in (3) by their available
centered energy upper bounds and estimating the bilinear term by L2
norms. All orders are fixed as n grows.

For 1≤r≤3, Hölder on the nonprincipal frequencies gives
\[
 \frac{\widetilde{\mathcal E}_r}{n^{2r}}
       \le B^{r/3}n^{-r}.
\]
For every real r≥3, positivity of the spectral sum gives instead
\[
 \frac{\widetilde{\mathcal E}_r}{n^{2r}}
       \le Bn^{-3}\delta^{2r-6}.                          \tag{8}
\]
The same bounds apply to integer s. At s=1, one may use the exact
bound \(\widetilde E_1/n^2\le1/n\). At s=2, Hölder gives
\(\widetilde E_2/n^4\le B^{2/3}/n^2\).
Using centered rather than raw energies in (8) is essential at
large orders; the principal frequency is absent.

Write T for the resulting upper-bound expression for √(pW_rV_s),
and suppress factors n^{o(1)}. The following table lists T and the
best saving supplied by its balance with δ^{rs}. Since C_{r,s}≥0,
adding that term cannot improve these limits.

| Multiplier order r | Walk order s | Available expression T | Balance saving |
|---|---|---|---|
| 1≤r≤3 | 1 | n^{(3−r)/2} | 0 |
| 1≤r≤3 | 2 | n^{(2−r)/2} | (r−2)_+/(4r) |
| 1≤r≤3 | s≥3 | n^{(1−r)/2}δ^{s−3} | (r−1)/(2[s(r−1)+3]) |
| r≥3 | 1 | δ^{r−3} | 0 |
| r≥3 | 2 | n^{-1/2}δ^{r−3} | 1/(2r+6) |
| r≥3 | s≥3 | n^{-1}δ^{r+s−6} | 1/(rs−r−s+6) |

The entries agree at r=3. For clarity, this is a constraint on the
displayed **upper-bound expression**, not a claim that T is an
actual lower bound on the bilinear sum. To see its implication,
consider an expression K n^{-κ}δ^λ with rs>λ. Even the stronger
scalar inequality δ^{rs}≤K n^{-κ}δ^λ permits
δ=c n^{-κ/(rs−λ)} for an appropriate positive constant c.
Consequently it cannot by itself force a larger power saving.
Including the nonnegative origin constant cannot repair that issue.

The global optimization is elementary and includes noninteger r:

- On 1≤r≤3, s=2, (r−2)_+/(4r) increases to 1/12.
- On 1≤r≤3, s≥3, put x=r−1∈[0,2]. The function
  x/(2(sx+3)) increases with x and decreases with s. Its maximum
  is 1/9 at r=s=3.
- On r≥3, s=2, 1/(2r+6) decreases from 1/12.
- On r≥3, s≥3, the denominator is
  (r−1)(s−1)+5≥9, with equality only at r=s=3.
- The remaining ranges have no positive balance saving.

Therefore
\[
 \boxed{\text{the best saving in this specified continuous-order
 substitution is }1/9,\text{ attained only at }r=s=3.}    \tag{9}
\]
This does not exclude orders varying with n, different probability
measures, an improved variance estimate in (2), or cancellation in
the centered bilinear sum before its absolute value is bounded.
It shows why merely adjusting the real multiplier order is not the
next missing input.

## Subordinate parity variant and the eighth-energy gap

Before the multiplier extension, this lane obtained the following
useful but weaker special case. If −1∈H and r is even, θ is real and
\[
 \sum_{x,y}\mu_r(x)\mu_s(y)e_p(axy)
   =\sum_y\mu_s(y)\theta(ay)^r
   \ge|\theta(a)|^{rs}.
\]
Thus the same unsquared centered gate holds with α_r replaced by
the positive μ_r. If H has no zero triple, r=4,s=3 gives
\[
 \delta^{12}\le E_2/n^4+\sqrt{pV_4V_3}
               \le B/n^2+B\delta/n.
\]
Here raw Hölder gives E₂≤E₃^{2/3}, and
\(\widetilde E_4\le M^2\widetilde E_3\).
Parseval gives
\[
 \delta^2\ge\frac{p-n}{n(p-1)}\ge\frac1{n^2},
\]
so \(\delta^{11}\le2B/n\). Its exponent 10/11 is superseded by
8/9 in (5).

The triple-free condition can be imposed at negligible prime density.
Let Z₃(p,n) count ordered zero triples in H. For dyadic n, p>3,
the action of H×C₃ by scaling and cyclic rotation is free: a
nontrivial stabilizer would give a cube root of unity in a group of
power-of-two order, then force all three entries equal. Thus a
nonzero Z₃ is at least 3n. The existing
[norm estimate (N')](cyclotomic-prime-average.md), with T₃=0, gives
\[
 \sum_{p\equiv1\ (n)}Z_3(p,n)\log p\le\frac{n^3\log3}{2}.
\]
At most n²log3/(6log(n⁴/4)) quartic primes have a zero triple.
The prime count established in
[pass 4](parallel4-subgroup-2026-09-04.md) makes their relative
proportion at most \((\log3/9+o(1))/n\).
This exclusion is unnecessary for (5). The fixed quartic witness
n=64, p=6700417, H=〈2〉 has Z₃=192=3n and is explicitly
excluded from the parity variant, not from the multiplier theorem.

The primary requested averaged centered eighth budget,
\[
 \sum_{p\in\mathcal P_N}(E_4(H_N)-N^8/p)\log p
       \ll N^7\operatorname{polylog}N,
\]
is still unproved. This pass does not obtain cancellation in the
arithmetic norm budget. The main-lane spectral improvement instead
implies on the sixth-energy class
\[
 \widetilde E_4
 \le B(B+2)^{2/9}n^{43/9},
\]
which improves the pass-6 exponent 44/9 but is larger than the
required n⁴ scale. The concrete next input is stronger control of
the centered variance or signed bilinear correlations, or an
arithmetic estimate for the centered eighth moment; multiplying
signed averaging measures is not a justified substitute.

## Bounded exact certificates

The checker uses the rational cosine intervals already certified in
[the pass-6 checker](../experiments/parallel6_subgroup_2026_09_04.py).
It verifies the two negative multiplier cosets, the strict Jensen
counterexamples for their indicators, and the valid Jensen step with
positive μ₃. It checks the continuous-order ledger on a bounded
rational grid and records the analytic monotonicity certificates
that cover the full parameter range. The n=64 witness is checked by
exact modular arithmetic, without factoring a large integer.

These finite checks complement the proofs above. They do not
establish the asymptotic claims by sampling.

Run [the checker](../experiments/parallel7_subgroup_2026_09_04.py).
[Results](../results/parallel7_subgroup_2026_09_04.json) include the
terminal status and source hashes.
