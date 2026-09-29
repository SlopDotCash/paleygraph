# Pass 39: independent review of the single-degree comparison

**Verdict:** the ordinary proof of (L), (A), and (R) in [pass 35](parallel35-single-degree-comparison-2026-09-05.md) is valid as stated, using the cited derivative-root theorem. The review found no missing range restriction, normalization factor, principal-term error, or constant error in those three inequalities.

The verdict applies to every odd prime `p`, nonempty symmetric set `A=-A` contained in `F_p^*`, and integer `s>=1`. It does not prove the additional positive upper hypothesis on `B_s`, does not improve the project's accepted period exponent, and does not settle the Paley conjecture or the prize problem. No Lean verification or full build was performed.

## External theorem and scope checked

The input is Ravichandran, *Principal submatrices, restricted invertibility and a quantitative Gauss-Lucas theorem*, [arXiv:1609.04187v2, Theorem 4.4](https://arxiv.org/html/1609.04187v2#S4.Thmtheorem4). The exact statement was checked in the current version-pinned primary HTML and visually on PDF page 12 of the archived `sources/parallel35-single-degree/ravichandran-1609.04187v2.pdf`.

For a degree-`d` real-rooted polynomial whose roots are in `[-1,1]` and sum to zero, the theorem encloses the roots of the `k`th derivative in

\[
[-2\sqrt{c(1-c)},2\sqrt{c(1-c)}],\qquad c=k/d\ge 1/2.
\]

Pass 35 uses an integer `k` with `k<d`, so the source's notation `p^(cd)` entails no rounding or constant-polynomial ambiguity. Repeated roots are permitted: the source does not assume simple roots. This review accepts the cited theorem as an external mathematical input and checks its application; it is not a new formal proof of the source theorem.

The radius can also be checked directly from the source's barrier estimate, Proposition 3.2. For roots `lambda_i in [-1,1]` of mean zero and `b>1`, the convex chord bound gives `Phi(b)=sum_i 1/(b-lambda_i) <= db/(b^2-1)`. Thus the barrier estimate bounds the largest derivative root by

\[
b-k/\Phi(b)\le(1-c)b+c/b.
\]

Its infimum over `b>1` is `2sqrt(c(1-c))` for `1/2<=c<1`, attained at `b=sqrt(c/(1-c))` when `c>1/2` and approached as `b` decreases to `1` at `c=1/2`. Applying the same argument to reflected roots supplies the lower endpoint. This checks the source constant and its boundary case without relying only on a transcription of Theorem 4.4; the barrier estimate remains the external input.

Both local source hashes were recomputed and agree with `results/parallel35_source_scope_2026_09_05.json`:

- HTML: `92e9d9909a9383246f095f2884859a5bb63814435cda01fc404ece580af76b64`.
- PDF: `01389932e3b45075eea3155d035b3d6caea254ed043aa1f23ea907fc6c450f27`.

## Pointwise factorization, including the scale factor

For arbitrary real `y_i in [-2,2]`, let

\[
N\ge1,\quad n=2N,\quad r=2s,\quad
F=\sum_{i=1}^{N}y_i,\quad H=r!e_r(y).
\]

First take `r<=N/2`. In particular `N>=4`, and the derivative order `k=N-r` lies in `[N/2,N)`. Put `mu=F/N`, `z_i=y_i-mu`, and `P(X)=prod_i(X+z_i)`. Then `sum z_i=0`, `|z_i|<=4`, and

\[
Q(X)=4^{-N}P(4X)=\prod_i(X+z_i/4)
\]

satisfies all the source hypotheses. Its derivative of order `N-r` is `4^(-r)P^(N-r)(4X)`. Consequently the roots `lambda_j` of `P^(N-r)` satisfy

\[
|\lambda_j|\le8\sqrt{(r/N)(1-r/N)}.
\]

The identity

\[
e_r(X+z_1,\ldots,X+z_N)=\frac{P^{(N-r)}(X)}{(N-r)!}
\]

has leading coefficient `binom(N,r)`, not one. Evaluating at `X=mu` and multiplying by `r!` therefore gives exactly

\[
H=a\prod_{j=1}^{r}(F-\theta_j),\qquad
a=\frac{(N)_r}{N^r},\qquad\theta_j=N\lambda_j.
\]

Every factor `1-j/N`, for `0<=j<r`, is greater than `1/2`. Thus `2^(-r)<=a<=1`. Also

\[
|\theta_j|\le8\sqrt{r(N-r)}\le R,
\qquad R^2=64rN=64sn.
\]

This verifies both changes of scale and the falling-factorial normalization in pass 35.

## Signs and numerical constants

The degree `r=2s` is even. Thus, when `|F|>=R`, every factor has the same weak sign and the product is nonnegative, including the negative-`F` branch. When `|F|<=R`, every factor has magnitude at most `2R`. It follows that

\[
H\ge-(2R)^r=-(256sn)^s.
\]

The absolute bound follows from

\[
|H|\le(|F|+R)^r
\le2^{r-1}(|F|^r+R^r)
\le4^s[F^{2s}+(64sn)^s].
\]

For the reverse bound, `|F|>=2R` gives

\[
H\ge a(|F|/2)^r\ge F^r/4^r,
\quad\text{hence}\quad F^r\le16^sH.
\]

For `|F|<2R`, put `L=(2R)^r=(256sn)^s`. The already established bounds `H>=-L` and `F^r<=L` imply

\[
F^r\le16^sH+(1+16^s)L
\le16^sH+2(4096sn)^s.
\]

The last numerical inequality is `1+16^s<=2*16^s` and `16*256=4096`. It remains valid when `H` is negative; no positivity of `H` is silently assumed in this branch.

## All remaining values of N and s

The complementary range `r>N/2` implies `n<8s`. If `r<=N`, the triangle inequality gives `|H|<=2^r(N)_r<=n^r`. If `r>N`, the convention `e_r=0` gives `H=0`. In both cases,

\[
|H|\le n^{2s}<(8sn)^s,\qquad
F^{2s}\le n^{2s}<(8sn)^s.
\]

These imply the lower and absolute bounds directly, and imply the reverse bound by the same `H>=-L`, `F^r<=L` argument. This covers `N=1,2,3`, the transition `r=N/2`, `r=N`, and arbitrarily large `s` with `r>N`. No theorem about derivatives is used in this complementary range. When `r>N`, both the word count and the principal falling-factorial term are zero, so `B_s=0` as well.

## Fourier normalization and principal subtraction

Because `p` is odd and `0` is excluded, negation partitions `A` into exactly `N=n/2` pairs. Choose representatives `h_i` and write

\[
y_i(a)=e_p(ah_i)+e_p(-ah_i)\in[-2,2].
\]

The product defining `e_r(y(a))` chooses `r` distinct opposite classes and one sign in each class. Every resulting set of elements has exactly `r!` orderings. Additive character orthogonality therefore gives

\[
\frac1p\sum_{a\in\mathbb F_p}H(a)=O_r.
\]

At the principal frequency `a=0`, one has `y_i(0)=2` and `H(0)=2^r(N)_r`. Hence

\[
B_s=\frac1p\sum_{a\ne0}H(a)
\]

exactly, including `r>N`. The row sum `F(a)=eta(a)` is real, so `F(a)^(2s)=|eta(a)|^(2s)`. Averaging the pointwise inequalities over the nonprincipal frequencies yields (L) and (R) by linearity and yields (A) by the triangle inequality. The constant contribution is multiplied by `(p-1)/p<=1`. Thus all three constants in pass 35 are preserved. This argument makes no real-rootedness claim about the averaged polynomial.

## Consequences checked and hypotheses that remain

For `K>=0`, the single-degree implications follow from the integer-power inequalities

\[
(16K)^s+2\,4096^s\le(16K+8192)^s,
\qquad K^s+64^s\le(K+64)^s.
\]

The subgroup corollary is valid for a multiplicative subgroup satisfying the standing symmetry hypothesis, equivalently `-1 in H`. Its nonprincipal Fourier values are constant on cosets of cardinality `n`, giving `M^(2s)<=(p/n)T_s`. The displayed `sqrt(e)` factor uses the natural logarithm in `s>=log(p/n)`. In a quartic window `p` comparable to `n^4`, choosing an integer `s` of order `log n` is consistent with all the proved inequalities.

The degree-ten deduction is also valid: `sum_(j=0)^9 2j=90`, the product estimate is nonnegative for `n>=90`, and using `p<=n^4` therefore preserves the stated inequality direction. Exact integer arithmetic confirms `90+1280^5=3435973836800090`.

The unresolved ingredient is a **uniform positive** estimate `B_s<=(Ksn)^s` with an appropriate absolute `K` at the required growing degrees. Neither (L), (A), nor (R) supplies that input. A fixed-degree instance or the unconditional negative bound does not supply it either. Extending a symmetric-subgroup criterion to the full arbitrary-set Paley problem is an additional issue, already excluded from the pass-35 claim.

This lane independently checked the universal algebra and source application. It inspected the existing verifier's relevant routines but did not rerun the pass-35 finite test corpus or treat its successful finite checks as a proof of universality. Apart from source-hash verification and the displayed integer constant, no new finite enumeration was necessary. No central project file was changed.
