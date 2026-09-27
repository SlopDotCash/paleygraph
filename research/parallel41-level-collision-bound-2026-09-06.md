# A bound for the full weighted triangle sum and its exponent barrier

This records the intermediate audit. The [subsequent correlation optimization](parallel41-correlation-optimization-2026-09-06.md)
removes the extra three logarithms, giving the final bound
W<<n^(86/15)(1+log n)^(8/15). The bounds below remain valid.

The weighted collision problem admits a direct upper bound for the full
nondegenerate triangle sum. With the currently imported inputs, the bound is

\[
\boxed{W\ll n^{86/15}(1+\log n)^{53/15}.}                         \tag{1}
\]

This improves the power `6` supplied by the previously used generic Young
bound by `4/15`. It remains above the required power `17/3` by `1/15`.
The root lane proposed the triple-correlation route and then its dyadic
energy refinement during this audit; both are independently verified
below. The direct route first gives power `347/60`; the final section
uses additional level energy information to obtain (1).

This is an ordinary mathematical argument using scoped published inputs.
No uniform energy improvement or full-conjecture proof is established.

## Notation, zero terms, and the shifted-energy input

Let `p` be odd, `H=-H`, `n=|H|`, and `n^2<p`. Put

\[
 f(x)=r_{H-H}(x),\quad E_*=\sum_{x\ne0}f(x)^2,\quad
 F_3^*=\sum_{x\ne0}f(x)^3,\quad
 W=\sum_{\alpha,\beta\ne0,\ \alpha\ne\beta}
 f(\alpha)^2f(\beta)^2f(\alpha-\beta)^2.
\]

Thus `E(H)=n^2+E_*`. On `G=F_p^*/H`, write `a(C)=f(C)` and
`A_r=sum_C a(C)^r`; in particular `A2=E_*/n` and `A3=F3*/n`.
The exact identity from pass 39 gives

\[
 a(C)=|(H-1)\cap C|.
\]

Set `R0=(H-1)\{0}`, `B=E^times(R0)`, and

\[
 X=B-\bigl(2(n-1)^2-(n-1)\bigr).
\]

For

\[
 \rho(u,v)=\#\{(x,y)\in uH\times vH:y-x=1\},
\]

[the existing exact shifted-energy identity](mixed-periods-and-shifted-energy.md)
is

\[
 X=\sum_{u,v}(\rho(u,v))_3\ge0.                                  \tag{2}
\]

Hence, for `H3={(u,v):rho(u,v)>=3}`,

\[
 \sum_{H_3}\rho^3\le\frac92 X.                                  \tag{3}
\]

The constant follows from `c^3<=9c(c-1)(c-2)/2` for integers `c>=3`.
One must restrict to these cells: the sum of `rho^3` over all cells also
contains a term of order `p` and cannot be bounded by `B` alone.

Shkredov's [Theorem 6](https://arxiv.org/html/1504.04522v1#Thmsatz6),
checked against the archived primary HTML
`sources/mixed-periods-2026-09-04/shkredov-1504.04522.html`, gives

\[
 B\le E^\times(H-1)\ll n^2(1+\log n),\qquad X\le B.               \tag{4}
\]

The source allows arbitrary nonzero shifts and requires `|Gamma||Pi|<p`;
here both subgroups are `H` and both shifts are `-1`. Its energy convention
includes zero, while `B` omits zero, so the displayed inequality has the
correct direction. No subtraction from a coarse upper bound is used to
claim a power saving for `X`.

## Exact full normalization and the fractional correlation bound

Put `b=a^2` and define

\[
 T(u,v)=\sum_C b(C)b(uC)b(vC).
\]

Normalizing `y-x=z` by `z`, with `alpha=y` and `beta=x`, gives

\[
 W=n\sum_{u,v}\rho(u,v)T(u,v),\qquad
 \sum_{u,v}T(u,v)=A_2^3.                                         \tag{5}
\]

Every edge is nonzero throughout. No injectivity, coset-distinctness,
dyadic selection, or three-set incidence hypothesis is needed in (5).

For any nonnegative function `b` on a finite group and `t>=1`, ordinary
counting norms give

\[
 \sum_{u,v}T(u,v)^t
 \le \left(\sum_C b(C)^{3t/(t+2)}\right)^{t+2}.                    \tag{6}
\]

Here the group is abelian; writing correlations instead of convolutions
only inserts inverse maps, which preserve the norms. A complete proof
uses Young twice. Set

\[
 s=\frac{3t}{t+2},\quad q=\frac{3t}{2t+1},\quad
 r=\frac tq=\frac{2t+1}{3},\quad v=\frac sq=\frac{2r}{r+1}.
\]

All four exponents are at least one. For fixed `u`, let
`h_u(C)=b(C)b(uC)`. Since `1+1/t=1/q+1/s`,

\[
 \|T(u,\cdot)\|_t\le\|h_u\|_q\|b\|_s.
\]

The second Young inequality applies to the correlation of `b^q` with
itself, with target exponent `r` and two source exponents `v`. Therefore

\[
 \sum_u\|h_u\|_q^t
 \le\|b^q\|_v^{2r}
 =\left(\sum_Cb(C)^s\right)^{2t/s}.
\]

Multiplying by `||b||_s^t` gives exponent `3t/s=t+2`, proving (6).
In particular, at `t=3/2`,

\[
 \|T\|_{3/2}\le\left(\sum_Ca(C)^{18/7}\right)^{7/3}.              \tag{7}
\]

Split (5) at `rho=2` and use (3), Holder, and (7). This proves

\[
 W\le2nA_2^3+
 \left(\frac92\right)^{1/3}nX^{1/3}
 \left(\sum_Ca(C)^{18/7}\right)^{7/3}.                            \tag{8}
\]

Using only `sum a^(18/7)<=A2^(3/7)A3^(4/7)` would already give
`W<<n^(347/60)(1+log n)^(28/15)`. The next step removes the cubic-moment
logarithm, without changing the power of `n`.

## Using the actual subgroup difference tail

The published Shkredov paper
`sources/parallel38-energy-source/shkredov-2013-published.pdf`, PDF page 10,
printed page 197, Theorem 4, orders values on distinct nonzero subgroup
cosets. Its `d=2` case, under `n<sqrt p`, states

\[
 a_j\ll n^{2/3}j^{-1/3}.
\]

The page was visually checked in this audit. Its convolution is `H*H`,
which agrees with `f` because `H=-H`. Thus, for an absolute `C0`,

\[
 N(\lambda):=|\{C:a(C)>\lambda\}|
 \le \min\{A_2\lambda^{-2},\ C_0n^2\lambda^{-3}\}.                \tag{9}
\]

The first arm is Markov's inequality; the second is the valid published
tail, not the corrected general three-set incidence statement.
For any `2<z<3`, split the layer integral at `C0 n^2/A2` to obtain

\[
 \sum_C a(C)^z
 \le\frac{z}{(z-2)(3-z)}
 A_2^{3-z}(C_0n^2)^{z-2}.                                       \tag{10}
\]

At `z=18/7`, the prefactor is `21/2`. Substitution in (8) gives

\[
 \boxed{W\le\frac{2E_*^3}{n^2}
       +C n^{8/3}X^{1/3}E_*}                                    \tag{11}
\]

for an absolute constant `C`. This is the useful bound before replacing
the actual shifted-energy excess by its published upper bound.

With `L=1+log n`, (4) turns (11) into

\[
 W\ll E_*^3/n^2+n^{10/3}E_*L^{1/3}.                              \tag{12}
\]

The project's scoped MRSS bound, Corollary 12, equation (32), in
`sources/sigma-subgroup-2026-09-05/mrss-1712.00410.txt`, is
`E(H)<<n^(49/20)L^(1/5)` for `n<=sqrt p`. The primary text was checked.
Its insertion in (12) gives powers `107/20` and `347/60`, with logarithmic
powers `3/5` and `8/15`, respectively. The second power dominates for
`n>=2`, proving the preliminary bound
`W<<n^(347/60)L^(8/15)`.

In particular, all terms supported on `rho<=2` satisfy
`W_low<<n^(107/20)L^(3/5)`, below the target `n^(17/3)` by `n^(19/60)`.
The remaining difficulty is concentrated at `rho>=3`, exactly where
the shifted-energy excess occurs. The positive quartic collision in
pass 40 had `rho=3` and belongs to this remaining part.

## What remains after the estimate

The actual dyadic collision correction from pass 40 is also bounded by
(11), up to an absolute factor. If `tau_i` is the upper endpoint of a
dyadic level, then `f<=tau_i<2f` on that level. Its weight
`tau_i^2 tau_j^2 tau_k^2` is at most `64` times the actual triangle
weight. The collision correction is a nonnegative part of this full
dyadic envelope. This implication retains all multiplicities and does
not assume any cancellation in the counts `C_ijk`.

The power gap in the preliminary bound is `347/60-17/3=7/60`. Even if one feeds
the energy dependence in (12) back into the old retained-triangle lower
bound `E(H)^6/(n^9 L)`, the two branches give only

\[
 E(H)\ll n^{7/3}L^{1/3}
 \quad\hbox{or}\quad
 E(H)\ll n^{37/15}L^{4/15}.                                     \tag{13}
\]

This calculation uses `E_*<=E(H)` and assumes the already-audited
zero-edge/retention alternative has been handled. The exponent `37/15`
is worse than the imported `49/20` by `1/60`. Thus (13) is a limitation
of this substitution, not an improved energy theorem.

Within the specific Holder/triple-correlation/tail substitution family,
changing the dual Holder exponent does not improve the power. Write
`E_*=n^(e+o(1))`, ignoring logarithms. For `2<r<=3`, the resulting
power is

\[
 3e-2+\frac{16-6e}{r};
\]

for `r>=3`, it is

\[
 3e-\frac43+\frac{14-6e}{r}.
\]

These follow from (6), (10), and respectively
`sum_(rho>=3)rho^r<<X` or `||rho||_r<=||rho||_3` on that set.
For `7/3<e<8/3`, the first decreases to `r=3` and the second increases
from `r=3`; the minimum is `e+10/3`. At `1<r<=2`, the higher moments
from (9) give power at least `6` (with endpoint logarithms), so that
range also does not improve the result at `e=49/20`.
This is not a no-go theorem for other uses of the subgroup structure.

One concrete sufficient new input is a power saving for the actual
excess `X`, rather than for the full shifted energy including its
quadratic trivial term. If `X<<n^(2-delta)` up to logarithms, the same
self-consistent calculation gives exponent `(37-delta)/15` in (13).
A saving `delta=1/3`, namely `X<<n^(5/3)` up to logarithms, reaches
power `22/9`; no such estimate is proved here. Direct insertion of the
current energy bound into (11) instead requires `X<<n^(33/20)` up to
logarithms to reach `W<<n^(17/3)` at once.

## Dyadic level energy: the stronger final bound

Partition the positive values of `a` into dyadic levels `D_i` with
`T_i<a(C)<=2T_i`, and put `d_i=|D_i|`. There are `O(L)` nonempty levels,
since `a` is integer valued and `a<=n`. Let

\[
 t_i(q)=|D_i\cap qD_i|,\quad
 m_{ijk}(u,v)=\#\{C\in D_k:uC\in D_i,\ vC\in D_j\}.
\]

These multiplicities include coincident edge cosets, as required for a
bound on the full `W`. Write `P_T=T_iT_jT_k`, `P_d=d_id_jd_k`.
Let `R_i` be the elements of `R0` lying in the cosets of `D_i`. The
correlation of the coset counts of `R_i` is at least `T_i^2 t_i(q)`.
Aggregating field ratios into `H` cosets and applying Cauchy gives

\[
 \sum_q t_i(q)^2\le\frac{nE^\times(R_i)}{T_i^4}
 \le\frac{nB}{T_i^4}.                                           \tag{14}
\]

The factor `n` is essential: each quotient ratio contains `n` field
ratios. Also `t_i(q)<=d_i`, so Holder and the exact collision identity
give

\[
 \begin{aligned}
 \|m_{ijk}\|_2^2
 &=\sum_qt_i(q)t_j(q)t_k(q)\\
 &\le nB\,P_d^{1/3}P_T^{-4/3}.
 \end{aligned}                                                   \tag{15}
\]

Since `||m_ijk||_1=P_d`, interpolation between `l1` and `l2` yields

\[
 \|P_T^2m_{ijk}\|_{3/2}
 \le(nB)^{1/3}P_T^{14/9}P_d^{4/9}.                              \tag{16}
\]

The mixed-function version of (6), with the same two-Young proof,
also gives the ordinary block estimate

\[
 \|P_T^2m_{ijk}\|_{3/2}\le P_T^2P_d^{7/9}.                        \tag{17}
\]

Taking `(16)^(3/5)*(17)^(2/5)` bounds the minimum of the two right-hand
sides from above. Its result is

\[
 \|P_T^2m_{ijk}\|_{3/2}
 \le(nB)^{1/5}\prod_{l\in\{i,j,k\}}(T_l^3d_l)^{26/45}
 \ll(nB)^{1/5}n^{52/15}.                                        \tag{18}
\]

The last step uses (9). The actual block contribution to `T` in (5)
is at most `64 P_T^2m_ijk`. Triangle inequality over `O(L^3)` blocks
therefore proves

\[
 \|T\|_{3/2}\ll(nB)^{1/5}n^{52/15}L^3.
\]

Applying (3) and Holder to the high-incidence part of (5), then (4),
gives

\[
 W_{\rho\ge3}\ll
 nX^{1/3}(nB)^{1/5}n^{52/15}L^3
 \ll n^{86/15}L^{53/15}.                                        \tag{19}
\]

Combining this with the previously checked low-incidence estimate
proves (1). This improvement uses the energy of the actual level
subsets, while bounding each such energy by the full shifted energy.
It does not make the failed injectivity assumption.

If one instead closes the old energy inequality using (19), the
high-incidence branch gives

\[
 E(H)\ll n^{221/90}L^{34/45}.
\]

Its exponent exceeds the imported `49/20` by `1/180`; it is therefore
still not an improvement to the project's energy input. The remaining
power gap for `W` is now `86/15-17/3=1/15`.

A more precise hereditary input would keep `E^times(R_i)` in (14),
instead of replacing it by `B`. A bounded primary-source search did
not verify a density-sensitive bound that improves this replacement
for these particular subsets. The inspected Shkredov Theorem 6 applies
to whole shifted subgroups. Its Theorem 7 requires invariant sets;
`R_i subset H-1` is generally not invariant. Its Corollary 5 includes a
subset `X subset H`, but gives a sumset estimate using the whole shifted
energy, not the needed density saving. No hereditary energy estimate
is imported on that basis.

The audit used exact norm identities, primary-source checks, a bounded
source search, and rational exponent arithmetic. It launched no Lean
builds or subgroup scans and changed no central files. The final bounded
result is (1), with actual-excess bound (11) also retained, and a remaining
power gap of `1/15` to the required triangle scale.
