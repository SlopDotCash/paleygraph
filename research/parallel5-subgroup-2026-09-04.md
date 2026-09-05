# Centered sixth moments and a stronger maximum bound on a growing class

**Status: an unconditional improvement on a density-one quartic class is
proved, not a worst-case Paley or prize bound.** The sixth additive
energy has size N^(3+o(1)) on this class. Combining it with the previous
fourth-energy result and a checked fixed-order amplification argument gives

\[
 \boxed{\max_{b\ne0}|\eta_{H_N}(b)|
       \le C N^{23/24}(\log N)^{7/72}.}                    \tag{1}
\]

This improves the recorded uniform exponent 2849/2880 only on the class
specified below. It does not approach the desired exponent 1/2 at
logarithmic moment depth. No novelty, optimality, or formalization claim
is made. The literature input is identified separately from the new
prime-average deduction.

## A centered sixth-moment theorem at every dyadic level

Let N=2^m≥4 and
\[
 \mathcal P_N=\{p\text{ prime}:N^4/4\le p\le N^4,\ p\equiv1\pmod N\}.
\]
At p in this set let H_s denote the order-s subgroup for every dyadic
divisor s≥2 of N. Define
\[
 U_s=E_3(H_s)-T_6(s)\ge0,\quad
 T_6(s)=15s^3-45s^2+40s,\quad
 \widetilde E_r(H_s)=E_r(H_s)-s^{2r}/p.
\]
The last quantity is the centered energy, exactly
\[
 \widetilde E_r(H_s)
 =\sum_{x\in\mathbb F_p}(r_r(x)-s^r/p)^2
 =\frac1p\sum_{b\ne0}|\eta_{H_s}(b)|^{2r}\ge0,             \tag{2}
\]
where r_r counts ordered r-term sums.

**Theorem A.** For every w>0, at most
\[
 \boxed{\frac{(8N^3-8)\log6}{14w\log(N^4/4)}}              \tag{3}
\]
primes of \(\mathcal P_N\) fail
\[
 \sum_{2\le s\mid N}\frac{U_s}{s^3}\le w.                  \tag{4}
\]
Thus at all levels of every remaining tower,
\[
 \begin{split}
 E_3(H_s)&\le T_6(s)+ws^3\le(15+w)s^3,\\
 \widetilde E_3(H_s)&\le T_6(s)+ws^3-s^6/p.                \tag{5}
 \end{split}
\]
In particular, the mean over nonzero frequencies satisfies the exactly
centered bound
\[
 \frac1{p-1}\sum_{b\ne0}|\eta_{H_s}(b)|^6
 \le\frac{p[T_6(s)+ws^3]-s^6}{p-1}.                       \tag{6}
\]
The principal term is retained explicitly; it is not replaced by an
uncentered surrogate in (2), (5), or (6).

**Proof.** The normalized-tuple norm argument from
[the prime-average note](cyclotomic-prime-average.md) applies to six
terms and gives
\[
 \sum_{p\equiv1\ (s)}U_s(p)\log p
 \le\frac{s^6-T_6(s)}2
       \log\frac{6s^6}{s^6-T_6(s)}
 \le\frac12s^6\log6.                                    \tag{7}
\]
For the last step, y log(6/y) increases on 0<y≤1. Intrinsic tuples
remain zero modulo every splitting prime, so U_s≥0. Divide (7) by
s³, sum over all levels, and use
\(\sum_{2\le s\mid N}s^3=(8N^3-8)/7\) and
\(\log p\ge\log(N^4/4)\). Markov's inequality proves (3).
Equation (2) proves the centered forms. ∎

The verified prime-counting input from
[pass 4](parallel4-subgroup-2026-09-04.md) is
\[
 |\mathcal P_N|\sim 3N^3/(8\log N).                       \tag{8}
\]
It follows from Thorner–Zaman's
[powerful-modulus prime number theorem, (3.2)](https://arxiv.org/html/2108.10878#S3.SS1),
at q=N, rad(q)=2, x=N⁴, h=3N⁴/4. Its fixed analytic parameter is
1/12; the w in (3) is unrelated and may depend on N. All constants
along these moduli are absolute.

The bad proportion in (3) is at most
\[
 \frac{8\log6}{21w}(1+o(1)).                              \tag{9}
\]
For w=log N this tends to zero. Intersect this class with the
pass-4 class having relative fourth-energy error N^(-1/2), whose bad
proportion is O(N^(-1/2)). Call the intersection \(\mathcal J_N\).
Then
\[
 \frac{|\mathcal J_N|}{|\mathcal P_N|}
   =1-O(1/\log N),\quad
 E_2(H_s)\le4s^2,\quad E_3(H_s)\le(15+\log N)s^3           \tag{10}
\]
simultaneously at every dyadic level. These are full energies where
needed for support estimates; (5)–(6) give their separate centered
consequences. The eventual prime-count threshold is not numerically
computed here.

## Why the sixth moment improves the maximum through amplification

A sixth-moment maximum bound by itself still loses the number of cosets:
\(M^6\le(pE_3-N^6)/N\). In the quartic window, E_3=O(N³) only
gives a bound of order N. The additional ingredient is multiplicative
subgroup invariance, used at each step below.

We use the weighted trilinear estimate stated as Lemma 4.1 of
[Di Benedetto et al., arXiv:2003.06165](https://arxiv.org/html/2003.06165#S4):
for X,Y,Z contained in F_p*, and weights of modulus at most one,
\[
 \left|\sum_{x,y,z}\alpha_x\beta_y\gamma_z e_p(axyz)\right|
 \ll p^{1/4}|X|^{3/4}|Y|^{3/4}|Z|^{7/8}.                 \tag{11}
\]
The paper attributes this estimate to Petridis–Shparlinski and uses it
in its §5 three-selection proof. Here the same fixed selection orders
are used, with the energy parameters retained. We do not assume that
the abstract arbitrary-order selections in the earlier ledger exist.

The primary source
[Petridis–Shparlinski, Theorem 1.1](https://arxiv.org/html/1604.08469v4#S1.SS3)
orders the sizes X≥Y≥Z. This does not require our selected sets to
arrive in that order: permuting the sets and their factorwise weights
gives the bound \(p^{1/4}(|X||Y||Z|)^{3/4}
\min(|X|,|Y|,|Z|)^{1/8}\), at most the right side of (11) with
our original labels. Multiplying one set by the nonzero frequency
absorbs a. The primary theorem has no further size restriction.

**Amplification lemma.** Let H have size n and let p lie in its quartic
window. Suppose
\[
 E_2(H)\le A n^2,\qquad E_3(H)\le B n^3,\qquad
 1\le A\le4,\quad 1\le B\le15+\log n.
\]
For all sufficiently large n, with an absolute constant,
\[
 M\le C p^{1/72}n^{65/72}
             A^{1/144}B^{1/36}(1+\log n)^{5/72}.          \tag{12}
\]
In the quartic window this implies (1). This lemma concerns the
endpoint n=N in (10); no quartic-window claim is made for its smaller
subgroups.

We give the selections and the zero deletions explicitly.

**Elementary selection fact.** If f:V→[0,1] has average at least a>0,
discard f<a/2 and split the remaining values into dyadic intervals.
There are at most
\(\ell=\lceil\log_2(2/a)\rceil+1\) intervals. Some interval with upper
endpoint t has
\[
 t\ge a/2,\quad f>t/2\text{ on }G,\quad
 |G|\ge a|V|/(2\ell t).                                  \tag{13}
\]
The discarded mass is at most a|V|/2; summing the mass in the remaining
intervals proves the cardinality assertion. In particular, the lower
bound for t loses no logarithmic factor.

Fix a nonzero frequency a and write \(|\eta_H(a)|=n\Delta\).
If Δ<n^(-1/24), the desired quartic conclusion already holds.
Otherwise all selections below have \(\ell\ll L=1+\log n\).

For u=(u₁,u₂,u₃) in H³, subgroup invariance gives
\[
 \sum_{u\in H^3}|\eta_H(a(u_1+u_2+u_3))|
       \ge n|\eta_H(a)|^3=n^4\Delta^3.
\]
Apply (13) to f(u)=|\eta_H(a\sum u)|/n. There is a set G₁ of triples
and a level Δ₁≥Δ³/2, with
\[
 |G_1|\gg\frac{n^3\Delta^3}{L\Delta_1},\qquad
 |\eta_H(a\sum u)|>n\Delta_1/2\quad(u\in G_1).
\]
Let X be the distinct sums of G₁ after deleting zero. Cauchy–Schwarz
and the full energy E₃ give
\[
 |X|\gg\frac{n^3\Delta^6}{B L^2\Delta_1^2}.               \tag{14}
\]
Before deleting zero, the same bound holds without loss. Its lower
bound is at least a constant times \(n^{11/4}/(BL^2)\), so subtracting
one is harmless for all sufficiently large n.

For v in H³ now put
\[
 f_2(v)=\frac1{n|X|}\sum_{x\in X}|\eta_H(ax(v_1+v_2+v_3))|.
\]
For every x in X, subgroup invariance and the triangle inequality give
\(\sum_v|\eta_H(ax\sum v)|\ge n|\eta_H(ax)|^3\).
Thus f₂ has average at least Δ₁³/8. Select G₂ using (13) at a
level Δ₂≥Δ₁³/16. Let Y be its distinct sums, again omitting zero.
Then
\[
 |Y|\gg\frac{n^3\Delta_1^6}{B L^2\Delta_2^2},\qquad
 \sum_{x\in X}|\eta_H(axy)|>n|X|\Delta_2/2\quad(y\in Y).
                                                                  \tag{15}
\]
The support bound before deletion is at least a constant times
\(n^{9/4}/(BL^2)\), so this second deletion is also valid.

Summing (15) over y and applying Cauchy–Schwarz in (x,y) gives
\[
 \sum_{z_1,z_2\in H}
 \left|\sum_{x\in X,y\in Y}e_p(axy(z_1-z_2))\right|
       \ge n^2|X||Y|\Delta_2^2/4.
\]
Apply (13) to the normalized summand as a function of (z₁,z₂).
There is a level Δ₃≥Δ₂²/8 and a set G₃ of pairs. Let Z be their
distinct differences after deleting zero. The full energy E₂ yields
\[
 |Z|\gg\frac{n^2\Delta_2^4}{A L^2\Delta_3^2},\qquad
 \left|\sum_{x\in X,y\in Y}e_p(axyz)\right|
       >|X||Y|\Delta_3/2\quad(z\in Z).                   \tag{16}
\]
Before deletion the lower bound is at least a constant times
\(n^{1/2}/(AL^2)\). Indeed Δ₂≥cΔ⁹, hence Δ₂⁴≥cΔ³⁶.
This proves the final deletion, the most restrictive one.

Choose a phase γ_z of modulus one to make each inner sum in (16)
nonnegative real. In (11) set α_x=β_y=1. Every member of X,Y,Z is
nonzero, and all weights meet its hypotheses. Comparing (11) with
(16), then taking fourth powers, gives
\[
 p\gg |X||Y||Z|^{1/2}\Delta_3^4
 \gg\frac{n^7\Delta^6\Delta_1^4\Delta_3^3}
           {A^{1/2}B^2L^5}.
\]
Finally Δ₁≥cΔ³, Δ₂≥cΔ₁³, and Δ₃≥cΔ₂² imply
\[
 \boxed{p\gg\frac{n^7\Delta^{72}}{A^{1/2}B^2L^5}.}        \tag{17}
\]
Taking the 72nd root proves (12) in the nontrivial branch. In the
discarded branch Δ<n^(-1/24), the right side of (12) dominates
n^(23/24) up to an absolute constant since p≥n⁴/4, A,B≥1, L≥1.
Thus both branches are covered. This proves the lemma and Theorem (1).

The exponent comparisons are exact:
\[
 1+\frac4{72}-\frac7{72}=\frac{23}{24}
       <\frac{2849}{2880},\qquad
 \frac5{72}+\frac1{36}=\frac7{72}.
\]
This is a stronger maximum-period bound on \(\mathcal J_N\), whose
density tends to one. It is not a new uniform bound for its complement.

## Higher centered moments that now follow, and the remaining gap

On \(\mathcal J_N\), for every integer r≥3, positivity of the
nonprincipal Fourier terms gives
\[
 \widetilde E_r\le M^{2r-6}\widetilde E_3.
\]
Consequently (with an absolute constant C, enlarged if necessary)
\[
 \boxed{\widetilde E_r(H_N)
  \le C^{2r}N^{(23r-33)/12}
             (1+\log N)^{\,1+7(r-3)/36}.}                \tag{18}
\]
The constant's dependence on r is displayed and the principal term
has already been removed. For example, at r=4 this gives
\(\widetilde E_4\ll N^{59/12}(\log N)^{43/36}\), improving the
N⁵ scale delivered by the raw prime-average norm budget. At
r comparable to log N, (18) is still far above the required
\((Cr)^rN^r\). It recovers the exponent 23/24, not 1/2.

For a precise audit of the raw budget's limitation, fix r≥4 and let
\(V_r(p)=E_r(H_N)-T_{2r}(N)\). Its norm bound is
\[
 \sum_{p\in\mathcal P_N}V_r(p)\log p
 \le A_{N,r}:=\frac{N^{2r}-T_{2r}(N)}2
       \log\frac{2rN^{2r}}{N^{2r}-T_{2r}(N)}
 \sim\tfrac12\log(2r)N^{2r}.                             \tag{19}
\]
The exact centered conversion is
\[
 \sum_p\widetilde E_r(p)\log p
 \le A_{N,r}
      +T_{2r}(N)\sum_p\log p
      -N^{2r}\sum_p\frac{\log p}{p}.                      \tag{20}
\]
The prime theorem used in pass 4 and the bounds of the quartic annulus give
\[
 \sum_p\log p\sim\tfrac32N^3,\qquad
 \sum_p\frac{\log p}{p}=\Theta(1/N).                      \tag{21}
\]
The second assertion follows directly by bounding 1/p between 1/N⁴
and 4/N⁴ in the sum. No additional prime-counting hypothesis is used.

At fixed r≥4, the intrinsic term in (20) is O(N^(r+3)), while the
subtracted principal term has size Θ(N^(2r−1)). Both are lower order
than the N^(2r) budget in (19). Thus this specific centered upper-bound
expression still has leading term \(\tfrac12\log(2r)N^{2r}\).
Its average scale is N^(2r−3). Subtraction alone cannot improve the
exponent. This is a statement about the displayed upper bound, not
an impossibility theorem for cyclotomic methods.

A strictly stronger actionable averaged estimate is already meaningful
at the next fixed depth:
\[
 \boxed{\sum_{p\in\mathcal P_N}
       \left(E_4(H_N)-N^8/p\right)\log p
          \ll N^7(\log N)^A\quad\text{for a fixed }A.}     \tag{22}
\]
It would give \(\widetilde E_4\le N^4(\log N)^{A+1}\) on a
density-one class by Markov and (8). This asks for an averaged
eighth-moment bound, much less than a uniform logarithmic-depth
theorem. It improves the raw norm budget by a factor N up to logs,
and also improves (18)'s N^(59/12) bound. Estimate (22) is unproved.

## Exact fourth energy does not force exact sixth energy

A fixed previously recorded quartic field supplies a counterexample
to that potential implication:
\[
 p=262657,\quad N=32,\quad N^4/4\le p\le N^4.
\]
Literal pair and triple sums give
\[
 E_2=2976=3N^2-3N,\qquad
 E_3=458240>446720=T_6(N).
\]
Thus the exact fourth-energy class can contain extra six-term
relations. This does not refute Theorem A, a Gaussian upper coefficient,
or the subgroup period target. For the generator g=247126 of order 32,
the six exponents (17,17,8,16,15,29) give values
(15531,15531,118090,262656,222204,153959) summing to zero modulo p,
with no opposite pair. The size and signs of all further relations
still require arithmetic control.

## Bounded verification

The companion script passed 12 fixed prime certificates and 47 moment
levels. It computes exact
pair/triple/fourfold convolution energies as specified in its certificate
list, verifies (2) with denominator-cleared integer arithmetic and
independent normalized zero-word counts, and records the order-32
counterexample above. It reconstructs bounded complex intrinsic counts
at orders six and eight, checks the all-level geometric sums, the dyadic
selection inequality on rational examples, and the complete exponent
and logarithm ledger in (14)–(18). It does not numerically test the
asymptotic theorem or claim small certificates verify its large-order
conclusion.

Run [the checker](../experiments/parallel5_subgroup_2026_09_04.py);
[the result](../results/parallel5_subgroup_2026_09_04.json) records its
terminal status and the hashes of all local proof dependencies.
