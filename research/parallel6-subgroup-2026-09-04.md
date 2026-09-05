# A weighted bilinear bound: exponent 17/18 on a growing class

**Status: a stronger spectral consequence is proved on a density-one
quartic class. The uniform subgroup conjecture remains open.**
For every prime in the sixth-moment class from pass 5,
\[
 \boxed{M_N\le 2^{1/6}(15+\log N)^{1/18}N^{17/18}.}        \tag{1}
\]
The exceptional proportion is O(1/log N). This improves the earlier
23/24 exponent on that class. The new argument uses a weighted
bilinear inequality and subgroup invariance, with no trilinear theorem
or dyadic selections. It does not merely interpolate the previous
maximum estimate. No novelty, optimality, or formalization claim is made.

## A centered mixed-order inequality

Let H≤F_p* have size n≥2, fix a≠0, and write
\[
 \eta(t)=\sum_{h\in H}e_p(th),\qquad
 \Delta=|\eta(a)|/n,\qquad f(t)=|\eta(at)|/n.
\]
For j≥1 let r_j(x) be the number of ordered j-term sums from H, and put
\[
 \mu_j(x)=r_j(x)/n^j,\quad u_j=\mu_j(0),\quad m_j=1-u_j.
\]
The probability measure μ_j is real and nonnegative. Choosing the
first j−1 entries determines the last, so
\[
 0\le u_j\le1/n.                                         \tag{2}
\]
Let Ω=F_p*, d=p−1, and define on Ω
\[
 \nu_j(x)=\mu_j(x)-m_j/d,\qquad V_j=\sum_{x\in\Omega}\nu_j(x)^2.
\]
Writing \(\widetilde E_j=E_j(H)-n^{2j}/p\), expansion gives
\[
 \boxed{V_j=\frac{\widetilde E_j}{n^{2j}}
             -\frac p{p-1}(u_j-1/p)^2\ge0.}              \tag{3}
\]
Thus V_j retains the principal-frequency subtraction and also removes
the separate origin discrepancy.

**Theorem A.** For every pair of positive integers r,s,
\[
 \boxed{
 (\Delta^{rs}-u_r-u_s+u_ru_s)_+^2
 \le \frac{p-n}{n(p-1)}m_r^2m_s^2
       +(1-1/n)m_rm_s\sqrt{pV_rV_s}.}                    \tag{4}
\]
In particular,
\[
 \begin{split}
 (\Delta^{rs}-u_r-u_s)_+^2
 &\le\frac1n+
        \frac{\sqrt p\,\sqrt{\widetilde E_r\widetilde E_s}}
             {n^{r+s}}\\
 &\le\frac1n+
        \frac{\sqrt p\,\sqrt{E_rE_s}}{n^{r+s}}.           \tag{5}
 \end{split}
\]
Both statements are uniform for all subgroups and all orders r,s.
Neither assumes a generalized selection argument.

### Proof: the weighted amplification

Subgroup invariance gives, for every t and j≥1,
\[
 \sum_x\mu_j(x)\frac{\eta(atx)}n
       =\left(\frac{\eta(at)}n\right)^j.                  \tag{6}
\]
Indeed averaging the j sums after multiplication by each h in H
permutes every chosen entry. Taking absolute values of the right side
in (6) yields
\[
 \sum_y\mu_s(y)f(xy)\ge f(x)^s.
\]
Average over μ_r and use Jensen for the nonnegative function f:
\[
 \sum_{x,y}\mu_r(x)\mu_s(y)f(xy)
 \ge\left(\sum_x\mu_r(x)f(x)\right)^s
 \ge\Delta^{rs}.                                        \tag{7}
\]
Here f(0)=1. Removing all pairs with x=0 or y=0 therefore removes
exactly the mass u_r+u_s−u_ru_s. If
\[
 S=\sum_{x,y\in\Omega}\mu_r(x)\mu_s(y)f(xy),
\]
then \(S\ge(\Delta^{rs}-u_r-u_s+u_ru_s)_+\).
When −1 belongs to H the periods are real, but an odd power is not
convex on the whole real line. Real-valuedness and odd orders alone
therefore do not justify replacing the absolute-value Jensen step
by a signed one.

### Proof: the weighted bilinear estimate and centering

For any functions A,B on F_p and c≠0, Cauchy–Schwarz and additive
orthogonality give
\[
 \begin{split}
 \left|\sum_{x,y}A(x)B(y)e_p(cxy)\right|^2
 &\le\|A\|_2^2\sum_x\left|\sum_yB(y)e_p(cxy)\right|^2\\
 &=p\|A\|_2^2\|B\|_2^2.                                 \tag{8}
 \end{split}
\]
This supplies the needed estimate directly, including all weights.
Extend ν_r,ν_s by zero at the origin. Their sums vanish. Since
\(\sum_{y\in\Omega}e_p(cxy)=-1\) for x≠0, expanding their constants
proves, for c≠0,
\[
 \begin{split}
 B_{r,s}(c)&:=\sum_{x,y\in\Omega}\mu_r(x)\mu_s(y)e_p(cxy)\\
           &=-m_rm_s/d+\sum_{x,y}\nu_r(x)\nu_s(y)e_p(cxy),\\
 \operatorname{Re}B_{r,s}(c)
           &\le-m_rm_s/d+\sqrt{pV_rV_s}.                 \tag{9}
 \end{split}
\]
The negative constant term in (9) is retained.

Weighted Cauchy–Schwarz, using total weight m_rm_s, gives
\[
 S^2\le m_rm_s
 \sum_{x,y\in\Omega}\mu_r(x)\mu_s(y)\frac{|\eta(axy)|^2}{n^2}.
\]
Expand the squared period. The n choices h=h′ contribute m_rm_s/n.
The n(n−1) remaining ordered pairs have nonzero coefficient
c=a(h−h′), so (9) applies. Consequently
\[
 S^2\le m_rm_s\left\{\frac{m_rm_s}{n}
       +(1-1/n)\left[-\frac{m_rm_s}{p-1}+\sqrt{pV_rV_s}\right]\right\}.
\]
Combining the first two terms proves (4). Since m_r,m_s≤1 and (3)
gives \(V_j\le\widetilde E_j/n^{2j}\le E_j/n^{2j}\), (5) follows. ∎

## The sixth-moment consequence

At r=s=3, (2) and (5) imply
\[
 \Delta^9\le2/n+\sqrt{1/n+\sqrt p\,E_3(H)/n^6}.            \tag{10}
\]
Suppose p≤n⁴ and E₃≤Bn³ with B≥1. Then
\[
 \Delta^9\le2/n+\sqrt{(1+B)/n}
       \le2\sqrt{(1+B)/n}.
\]
The last inequality uses n≥2 and B≥1. Squaring gives
\[
 \boxed{\Delta^{18}\le8B/n,\qquad
             M\le2^{1/6}B^{1/18}n^{17/18}.}              \tag{11}
\]
The constants in (11) are explicit and no eventual-size qualifier is
needed for this conditional implication. The displayed bound may be
weaker than M≤n in small cases.

Now let N be dyadic and
\[
 \mathcal P_N=\{p:N^4/4\le p\le N^4,\quad p\text{ prime},
                                        \quad p\equiv1\pmod N\}.
\]
Use exactly the sixth-moment class \(\mathcal K_N\) of pass 5,
defined by
\[
 \sum_{2\le t\mid N}
 \frac{E_3(H_t)-(15t^3-45t^2+40t)}{t^3}\le\log N.        \tag{12}
\]
The sum is over dyadic divisors. Its explicit bad-prime count is at most
\[
 \frac{(8N^3-8)\log6}{14\log N\,\log(N^4/4)}.
\]
Together with \(|\mathcal P_N|\sim3N^3/(8\log N)\), this proves
\[
 |\mathcal K_N|/|\mathcal P_N|=1-O(1/\log N).
\]
The norm argument and the powerful-modulus prime-count theorem behind
this statement are proved and sourced in
[pass 5](parallel5-subgroup-2026-09-04.md) and
[pass 4](parallel4-subgroup-2026-09-04.md); the primary analytic source is
[Thorner–Zaman, equation (3.2)](https://arxiv.org/html/2108.10878v2#S3.SS1).
The current argument adds no external analytic theorem.

For p in \(\mathcal K_N\), E₃(H_N)≤(15+log N)N³. Substitution
in (11) proves (1). No intersection with an extra fourth-energy
condition is required. All smaller levels still satisfy (12), but
p≤t⁴ need not hold there; applying (11) indiscriminately to them
would be invalid. Formula (10), with its actual p, remains available
at every level.

## Feedback through higher moments has a precise limitation

As a consequence, not as the main argument, positivity of nonprincipal
Fourier terms gives for every j≥3
\[
 \widetilde E_j\le M^{2j-6}\widetilde E_3
 \le2^{(j-3)/3}B^{(j+6)/9}n^{(17j-24)/9}.                \tag{13}
\]
In particular the centered eighth energy is bounded on this class by
\[
 \widetilde E_4\le2^{1/3}B^{10/9}n^{44/9}.                \tag{14}
\]
This improves the previous class exponent 59/12 to 44/9, but is still
larger than n⁴. It does not establish an averaged bound over the
exceptional primes.

Could inserting (13) back into the general mixed-order gate (5) improve
the maximum again? The answer is no for the exponent supplied by that
substitution, even though (4) contains information discarded by it.

On the current class, \(\widetilde E_1\le n\), and Hölder gives
\[
 \widetilde E_2\le(1-1/p)^{1/3}\widetilde E_3^{2/3}
                    \le B^{2/3}n^2.
\]
For fixed j define a_j by the available centered energy exponent
\(\widetilde E_j\le n^{2j-a_j+o(1)}\). Thus
\[
 a_1=1,\quad a_2=2,\quad a_j=(j+24)/9\quad(j\ge3).
\]
From (5), p≤n⁴, and u_r,u_s≤1/n, the saving supplied is
\[
 \sigma(r,s)=
 \frac{\min\{1,\max[0,(a_r+a_s-4)/2]\}}{2rs}.             \tag{15}
\]
The cap at one comes from the n^(-1) diagonal term.

For r,s≥3, (15) is at most 1/(2rs)≤1/18, with equality at r=s=3.
By symmetry the remaining cases have r=1 or r=2. For r=2,s≥3,
\[
 \sigma(2,s)=\frac{\min(1,(s+6)/18)}{4s}\le1/24;
\]
the expression decreases through s=12 and is then 1/(4s).
For r=1,s≥4,
\[
 \sigma(1,s)=\frac{\min(1,(s-3)/18)}{2s}\le1/42;
\]
it increases through s=21 and decreases afterwards. The remaining
pairs have zero saving. This proves the global bound
\[
 \boxed{\sigma(r,s)\le1/18\quad\text{for every fixed }r,s.} \tag{16}
\]
This is an obstruction to improving the result by feeding the present
moment estimates into the same coarse gate. It is not a lower bound
on M or a limit on all weighted bilinear arguments. In particular,
the signed terms before (9) and the sharp variance information in
(3) have not been estimated beyond Cauchy–Schwarz.

## A better eighth moment would have a verified spectral payoff

The general theorem provides a rigorous next bridge, without an
assumption about generalized selections. Suppose instead that
\[
 \widetilde E_4(H)\le Dn^4,\qquad D\ge1,\qquad p\le n^4.
\]
Use (5) with r=1,s=4. Here u₁=0, u₄≤1/n, and
\(\widetilde E_1\le n\). Therefore
\[
 \Delta^4\le1/n+\sqrt{1/n+\sqrt D/\sqrt n}
               \le3D^{1/4}n^{-1/4},
\]
which proves the explicit conditional consequence
\[
 \boxed{M\le3^{1/4}D^{1/16}n^{15/16}.}                   \tag{17}
\]
Thus the still-unproved averaged centered eighth budget
\[
 \sum_{p\in\mathcal P_N}\widetilde E_4(H_N)\log p
                       \ll N^7(\log N)^A
\]
would, by Markov and the prime count, give a density-one class with
D=(log N)^(A+1) up to a constant. Formula (17) would improve 17/18
to 15/16 on that class. The current feedback bound (14) supplies
neither that moment scale nor that stronger exponent.

## A zero-atom correction cannot simply be dropped

Take H={±1} in F₅, a=2, r=1, s=2. Then
\[
 \Delta=(1+\sqrt5)/4,\quad u_1=0,\ u_2=1/2,\quad
 V_1=1/4,\ V_2=1/16.
\]
The right side of (4) is (3+√5)/32, whereas omitting the atom
correction on the left would give
\[
 \Delta^4=(7+3\sqrt5)/32>(3+\sqrt5)/32.
\]
The correct left side is (3−√5)/32 and satisfies the bound.
This example is outside the endpoint quartic family, but it lies in
the general theorem and illustrates why origin removal is exact.
No failure of (4), (5), or the subgroup target follows.

## Bounded certificates

The companion checker uses fixed small prime fields, exact rational
walk distributions, and rational intervals for the real periods.
The intervals use a Machin arctangent certificate for π and a Taylor
remainder for cosine. It checks (3) and the zero-atom terms exactly,
then verifies the pointwise centered gate (4) at the specified mixed
orders and frequencies; full-group equalities are handled algebraically.
The sharp F₅, n=2, a=2, r=s=1 equality is also checked in Q(√5).
The completed run has 8 fixed prime/subgroup pairs, 32 exact walk moments,
216 pointwise gate certificates covering 1,440 nonzero-frequency/order
pairs by subgroup invariance, and 64 polynomial coefficient certificates
for (9). Of the 216 gate checks, 213 have strict rational interval
certificates and three are algebraic equalities. The checker separately
certifies the F₅ counterexample above. The moment and
feedback exponents, the explicit constants in (11) and (17), and the
finite ranges sufficient for the global optimization (16) are audited.
No prime scan, large integer factorization, or numerical asymptotic fit is
used.

The first verifier run stopped at the sharp pentagon equality because
outward intervals cannot certify equality. The terminal run passed
after adding the exact equality certificate; no tolerance or weakened
inequality was introduced.

Run [the checker](../experiments/parallel6_subgroup_2026_09_04.py).
[Results](../results/parallel6_subgroup_2026_09_04.json) record the
terminal status and source hashes.
