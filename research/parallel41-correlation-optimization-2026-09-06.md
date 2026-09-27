# Pass 41: mixed-level correlation optimization without a dyadic logarithm

Date: 2026-09-06. This is a mathematical optimization of the inequalities in
[the level-collision audit](/Users/shawwalters/Desktop/paleygraph/research/parallel41-level-collision-bound-2026-09-06.md).
It retains all triangle and level-intersection multiplicities. No new
incidence theorem, subgroup energy theorem, or full-target proof is asserted.

The improvement is

\[
\boxed{W\ll E_*^3/n^2+n^{14/3}X^{1/3}B^{1/5}.}       \tag{1}
\]

With the already scoped weak cubic tail, the exact mass Σ_Ca(C)=n−1,
and B≪n²(1+log n), X≤B, this gives

\[
\boxed{W\ll n^{86/15}(1+\log n)^{8/15}.}             \tag{2}
\]

The prior bound had logarithmic power 53/15. The power of n is unchanged
and still exceeds 17/3 by 1/15. The sharper asymmetric intersection bound
and summation by the largest level remove the dyadic loss entirely. The
root lane supplied the simpler global-energy summation used below after
the initial distribution-sensitive optimization; it has been independently
checked here.

## 1. Quantities retained before inserting published estimates

Use the notation of the level-collision audit: H=−H, n=|H|, n²<p,
G=F_p^*/H, a(C)=r_(H−H)(C), b=a², and

\[
T(u,v)=\sum_C b(C)b(uC)b(vC),\qquad
A_2=\sum_Ca(C)^2=E_*/n.
\]

Partition the positive coefficients into levels D_i with
T_i<a(C)≤2T_i, where T_i=2^i and i may start at −1 to include
a(C)=1. Put d_i=|D_i| and define the actual dyadic tail constant

\[
Q=\max_{i:d_i>0} T_i^3d_i.                           \tag{3}
\]

The previously checked subgroup tail gives Q≪n². For the intermediate
inequalities Q remains the actual value in (3).

It is also useful to retain the weighted quotient energy

\[
K=\sum_q\left(\sum_C a(C)a(qC)\right)^2.             \tag{4}
\]

Let R₀=(H−1)\{0}, B=E×(R₀), and X be its nontrivial excess as
defined in the level-collision audit. Grouping ratios of elements of R₀
into H-cosets gives the exact identity

\[
\sum_Ca(C)a(qC)=\sum_{\lambda\in qH}r_{R_0/R_0}(\lambda).
\]

Cauchy–Schwarz on each n-element ratio coset therefore gives

\[
K\le nB.                                            \tag{5}
\]

The main correlation estimate proved below is the stronger data-sensitive
form

\[
\boxed{\|T\|_{3/2}\ll K^{1/5}Q^{26/15}.}            \tag{6}
\]

All norms here are counting norms. No norm is over a selected injective
subfamily of triangles.

## 2. Two bounds for each full mixed-level multiplicity

Set

\[
t_i(q)=|D_i\cap qD_i|,\qquad
m_{ijk}(u,v)=|\{C\in D_k:uC\in D_i,\ vC\in D_j\}|.
\]

Write P_T=T_iT_jT_k and P_d=d_id_jd_k. The exact identities and bound

\[
\|m_{ijk}\|_1=P_d,\qquad
\|m_{ijk}\|_2^2=\sum_qt_i(q)t_j(q)t_k(q),\qquad
\|m_{ijk}\|_\infty\le\min(d_i,d_j,d_k)               \tag{7}
\]

include coincident cosets and arbitrary fiber multiplicities.

The last bound and L¹–L∞ interpolation improve the symmetric Young
estimate on unequal levels:

\[
\boxed{\|P_T^2m_{ijk}\|_{3/2}
\le P_T^2P_d^{2/3}\min(d_i,d_j,d_k)^{1/3}.}          \tag{8}
\]

Indeed Σm^(3/2)≤||m||∞^(1/2)Σm. Since
min(d_i,d_j,d_k)≤P_d^(1/3), (8) is at least as strong as
P_T²P_d^(7/9), with equality of these bounds when the three sizes agree.

For completeness, the energy arm can also use K rather than nB. On D_i,

\[
T_i^2t_i(q)
\le\sum_Ca(C)\mathbf1_{D_i}(C)
             a(qC)\mathbf1_{D_i}(qC)
\le\sum_Ca(C)a(qC).
\]

Thus Σ_q t_i(q)²≤K/T_i⁴. Combining t_i≤d_i, Hölder in
(7), and L¹–L² interpolation yields

\[
\begin{aligned}
\|m_{ijk}\|_2^2&\le K P_d^{1/3}P_T^{-4/3},\\
\boxed{\|P_T^2m_{ijk}\|_{3/2}
&\le K^{1/3}P_T^{14/9}P_d^{4/9}.}
\end{aligned}                                      \tag{9}
\]

Alternatively one can retain the individual level energies throughout.
The simpler global K bound suffices for the improvement below; no assumption
about distributing the whole energy evenly among levels is used.

## 3. Largest-scale summation removes all dyadic logarithms

First order the three level endpoints so T_i≤T_j≤T_k=z and write

\[
T_i=z2^{-r},\quad T_j=z2^{-s},\qquad r\ge s\ge0.
\]

This covers all ordered triples after at most six permutations. Repeated
indices cause only harmless overcounting in an upper bound.

Since d_l≤Q/T_l³,

\[
P_d\le Q^3/P_T^3,\qquad
\min(d_i,d_j,d_k)\le d_k\le Q/z^3.
\]

The two arms (8) and (9) consequently become

\[
U(z)=Q^{7/3}z^{-1},\qquad
V(z,r,s)=K^{1/3}Q^{4/3}z^{2/3}2^{-2(r+s)/9}.         \tag{10}
\]

The crucial point is that the first arm depends on the largest endpoint,
not merely on the product of all three endpoints.

For any A,D>0, a two-sided geometric sum gives

\[
\sum_{z\in2^{\mathbb Z}}\min(Az^{-1},Dz^{2/3})
\ll A^{2/5}D^{3/5}.                                 \tag{11}
\]

To prove it, split at z₀=(A/D)^(3/5). Below z₀, sum the increasing
geometric arm Dz^(2/3); above z₀, sum the decreasing arm A/z.
Rounding the split to a dyadic endpoint changes only an absolute constant.
Extending a finite collection of levels to all dyadic z only increases the
sum.

For fixed gaps r,s, (10) and (11) give

\[
\sum_z\min(U(z),V(z,r,s))
\ll K^{1/5}Q^{26/15}2^{-2(r+s)/15}.                  \tag{12}
\]

The double sum over r≥s≥0 converges, since it is bounded by the product
of two convergent geometric series. Finally, each actual mixed-level block
of T is pointwise at most 64P_T²m_ijk, because a≤2T_i on D_i.
Triangle inequality, the finite permutation factor, and (12) prove (6).

There is no factor for the number of levels in (6). Merely replacing the
minimum of the two arms by one interpolated expression before summing all
triples would lose the geometric decay and reintroduce an unnecessary log.

## 4. Weighted triangle bound with actual energy and excess

The exact normalization and shifted-energy identities from the input audit
give

\[
W=n\sum_{u,v}\rho(u,v)T(u,v),\qquad
\sum_{\rho\ge3}\rho^3\le\frac92X,\qquad
\sum_{u,v}T(u,v)=A_2^3.
\]

Splitting at ρ=2 and applying Hölder and (6) therefore yields

\[
\boxed{W\le\frac{2E_*^3}{n^2}
+C n X^{1/3}K^{1/5}Q^{26/15}.}                      \tag{13}
\]

This is the distribution-sensitive bound before Q≪n² and K≤nB.
The earlier fractional-moment argument also gives
||T||_(3/2)≪A₂Q^(4/3). For example, split the weak-cubic tail
Qλ^(−3) against the A₂λ^(−2) tail, apply the layer integral at
18/7, and then the triple-correlation Young inequality. Hence one can retain
both alternatives:

\[
\boxed{W\le\frac{2E_*^3}{n^2}
+C nX^{1/3}
\min\{A_2Q^{4/3},\ K^{1/5}Q^{26/15}\}.}             \tag{14}
\]

Dyadic Q gives the weak-cubic tail up to an absolute geometric-series
constant, so no extra logarithm is needed here either. In particular,

\[
W\ll E_*^3/n^2+
X^{1/3}\min\{n^{8/3}E_*,\ n^{14/3}B^{1/5}\}.       \tag{15}
\]

The switch occurs at E_*≈n²B^(1/5), approximately n^(12/5)
when B has quadratic size. Substituting the published inputs already checked
in the input audit proves (1).

For (2), no MRSS 49/20 energy bound is required. The older bound
E(H)≪n^(5/2) already suffices for the low-incidence term, giving
W_low≪n^(11/2). It also follows directly from the inputs retained here:
with A₁=Σ_Ca(C)=n−1, the distribution function obeys

\[
N(\lambda)\le\min\{A_1/\lambda,\ C Q/\lambda^3\}.
\]

Splitting 2∫₀^∞λN(λ)dλ at λ=(CQ/A₁)^(1/2) gives
A₂≤4√(CA₁Q)≪n^(3/2), hence E_*≪n^(5/2).
The high-incidence term in (1) has the asserted n^(86/15)L^(8/15)
bound, and 11/2<86/15. Thus the strongest W theorem uses only the
weak cubic tail, the shifted-energy input, and the exact identities. The
MRSS estimate appears below only to compare the extremal level with the
stronger energy cap available in the project. This note imports no stronger
bound on B, X, Q, or E_*.

## 5. The extremal scale of these particular inequalities

The power of n in this optimization is sharp at the level of the stated
functional inequalities. It does not assert the existence of an actual
prime-field subgroup with extremal character or collision behavior.

For an exact functional witness, let D be a subgroup of an arbitrary finite
abelian quotient group, |D|=d, and put a=A·1_D, A>0. Then

\[
T(u,v)=A^6d\,\mathbf1_{D\times D}(u,v),\qquad
\|T\|_{3/2}=A^6d^{7/3},\qquad K=A^4d^3.            \tag{16}
\]

For one dyadic level with lower endpoint T₀ and T₀<A≤2T₀,
Q=T₀³d. Consequently

\[
K^{1/5}Q^{26/15}
=(T_0/A)^{26/5}A^6d^{7/3},
\]

within an absolute factor of the left side of (6). Also
A₂Q^(4/3)=(T₀/A)⁴A⁶d^(7/3). Thus both alternatives in
(14) can have the same scale. The level collision identity, L¹–L∞
bound, and L¹–L² interpolation are all attained for this constant
subgroup-supported multiplicity.

Ignoring logarithms, let one potentially dominant level have amplitude
n^τ and size n^δ. With K≤n³ and Q≤n², its two bound exponents
are at most

\[
\min\{6\tau+7\delta/3,
       1+14\tau/3+4\delta/3\},\qquad
3\tau+\delta\le2.
\]

They are consequently bounded by

\[
\min\{14/3-\tau,\ 11/3+2\tau/3\}.
\]

The first arm decreases and the second increases. Their crossing is

\[
\boxed{\tau=3/5,\quad\delta=1/5,\quad
\|T\|_{3/2}\text{ exponent }=61/15.}                \tag{17}
\]

At this scale A₂ has exponent 7/5, below the available cap 29/20,
and the level mass A₁ has exponent 4/5, below the available total n−1.
Thus neither the global upper mass budget nor the imported additive-energy
upper bound eliminates this potentially dominant level. This observation
concerns a level's contribution, not a claim that the one-level witness
already satisfies every exact identity of H−1.

For a homogeneous bound using only K and Q, the witness makes the same
constraint explicit. A putative bound ||T||_(3/2)≪K^αQ^β must
have 4α+3β=6 under rescaling a. The subgroup-size dependence in
(16) then requires 3α+β≥7/3, hence α≥1/5. Under K≤n³
and Q≤n², the smallest corresponding power is 61/15. Additional
information beyond these two functional budgets can change this conclusion.

Finally, Hölder alone permits its allowed high-incidence cubic mass to lie
where T is large. At the scale in (17), the support D² has size
n^(2/5); cubic mass X≈n² corresponds to incidence height
n^(8/15) on that support. The resulting upper-bound scale is

\[
n\,X^{1/3}\|T\|_{3/2}
=n^{1+2/3+61/15}=n^{86/15}.                          \tag{18}
\]

This is an extremizer of the Hölder budget, not a construction of the
actual incidence function ρ associated with H. It identifies precisely
why this optimization retains a 1/15 power gap. New actual-H structure,
a smaller K or X, or a constraint on their joint concentration could still
improve that power.

## 6. The exact incidence row couples the tail to the excess

The root lane identified an additional consequence of the existing exact
identities. It is verified here independently. Because −1∈H, every
H-coset is stable under negation, so

\[
\rho(1,C)=|(H+1)\cap C|=|(H-1)\cap C|=a(C).
\]

The zero element never lies in C. This single row is a nonnegative part
of X=Σ_(u,v)(ρ(u,v))₃. Writing A_j=Σ_Ca(C)^j and
(a)₃=a(a−1)(a−2), one obtains

\[
\boxed{X\ge\sum_C(a(C))_3=A_3-3A_2+2A_1.}           \tag{19}
\]

In particular, a³≤4a for a=0,1,2, while
a³≤(9/2)(a)₃ for integral a≥3. Since A₁=n−1,

\[
Q\le A_3\le4(n-1)+\frac92X,
\qquad Q\ll\min\{n^2,n+X\}.                         \tag{20}
\]

The more precise A₃≤X+3A₂−2A₁ from (19) can also be retained
when the actual A₂ is known. Combining (13), (20), and the exact identity
B=2(n−1)²−(n−1)+X gives

\[
\boxed{W_{\rho\ge3}\ll
nX^{1/3}(nB)^{1/5}\min\{n^2,n+X\}^{26/15}.}         \tag{21}
\]

For n≤X≤n², one has B≍n² and therefore

\[
\boxed{W_{\rho\ge3}\ll n^{8/5}X^{31/15}.}           \tag{22}
\]

For X<n, the uniform form
n^(8/5)X^(1/3)(n+X)^(26/15) handles the remaining case
and is already below the required triangle power. The case X=0 has
W_(ρ≥3)=0 exactly.

Consequently the unproved hypothesis

\[
\boxed{X\ll n^{61/31}\quad\text{up to fixed logarithmic factors}}
                                                               \tag{23}
\]

would suffice for W_(ρ≥3)≪n^(17/3), up to logarithmic factors:
8/5+(31/15)(61/31)=17/3. More explicitly, if
X≪n^(61/31)L^c for fixed c, then B≪n² eventually and
W_(ρ≥3)≪n^(17/3)L^(31c/15). The low term remains
O(n^(11/2)). Equivalently, for 0<δ<1, a saving
X≪n^(2−δ) produces high-part power
86/15−31δ/15; δ≥1/31 reaches the desired power.
No such uniform saving for the actual X is proved in this note.

This joint relation clarifies the extremal discussion in section 5. The
level a≈n^(3/5), d≈n^(1/5) has A₃ of order n² and
Σ(a)₃ of the same order. It therefore cannot occur with subquadratic X
in an actual-H configuration satisfying (19). It remains compatible in
power with the currently available quadratic bound on X. Thus (19)
does not improve the unconditional power in (2), and it does not change
the sharpness of the standalone K,Q functional inequality. It does make
the sufficient new excess estimate substantially weaker than treating X
and Q as independent upper-bound budgets.

The audit is algebraic and read-only outside this note. No Lean builds,
subgroup scans, or central-file edits were performed.
