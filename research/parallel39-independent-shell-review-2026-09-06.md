# Independent mathematical review of shell inversion, sections 1-4

Date: 2026-09-06. Scope: ordinary mathematical review by a separate agent of [parallel36-shell-inversion-2026-09-05.md](parallel36-shell-inversion-2026-09-05.md), sections 1-4. This is not Lean verification, external peer review, or a proof of a uniform thin-annulus bound.

Audited source SHA-256: `986655882b1790acfd1a34e27edc982675ac59d74b51dc3a58f9bae151297d9b`. The hash identifies the source version; the review below checks its mathematics rather than accepting its prior hash-based audit.

## Verdict and precise domain clarification

**Sections 1-3 check out as stated. Section 4 checks out when its dyadic-subgroup argument is read with `A=H` and `a != 0`.** I found no incorrect centering constant, inverse coefficient, norm comparison constant, or lattice-norm exponent in that domain.

One local clarification is needed before quoting section 4 as a standalone theorem. At the start of its dyadic-subgroup paragraph (source lines 137-139), explicitly say:

> Set A=H and fix a in F_p^*.

The preceding definitions allow `a=0` and introduce a general symmetric set `A`; the new paragraph introduces `H` without explicitly making either restriction. The calculations clearly concern the nonzero subgroup cosets discussed immediately before it. Taken literally at `a=0`, however, the assertions `aN != 0`, rank one, nonzero norm, and equation (7) are false: `F_0=0` and `V(0)=0 < p^(-2/N)`. This is a domain clarification, not a counterexample to the intended nonzero-frequency argument. No source or central file was edited by this review.

## 1. Exact centering and the zero extension

Let `q=(p-1)/2`. Directly,

\[
\sum_{a\ne0}r_j(a)^2=2\sum_{t=1}^qt^2
=\frac{p(p^2-1)}{12}.
\]

Dividing by `(p-1)p^2` and summing the `N=n/2` representatives gives

\[
\overline V=\frac{N(p+1)}{12p}=\frac{n(p+1)}{24p}.
\]

For each nonzero `h`, the complete additive-character sum over `a` is zero, hence its sum over `a != 0` is `-1`. Thus `average eta=-n/(p-1)`. Pairing `h,-h` proves that `eta` is real. These arguments do not require `A` to be a subgroup.

Writing `s=n/(p-1)`, the elementary norm inequality gives `|M-M_f|<=s<=1`. More generally, for every normalized `l^r` norm with `1<=r<=infinity`,

\[
\big|\|\eta\|_r-\|f\|_r\big|\le s.
\]

The artificial assignments `f(0)=d(0)=0` in the source are essential and consistent: they apply to the operator argument, not to the physical distance `V(0)` or the true principal period `eta(0)=n`.

## 2. Fourier identity and inverse

I independently reopened [NIST DLMF 24.8.1](https://dlmf.nist.gov/24.8.E1). At its index `n=1`, it gives the cosine series for `B_2`. For `-1/2<=y<=1/2`,

\[
B_2(y+1/2)=y^2-1/12,
\qquad
\cos(2\pi k(y+1/2))=(-1)^k\cos(2\pi ky).
\]

Periodization yields the source's squared-distance series. Since the coefficient magnitudes sum to `zeta(2)`, it converges absolutely and uniformly; endpoints cause no problem.

Pairing the two elements of each representative pair first gives the uncentered identity

\[
V(a)=\frac N{12}+\frac1{2\pi^2}\sum_{k\ge1}\frac{(-1)^k}{k^2}\eta(ka).
\]

For `p` not dividing `k`, subtracting the average over nonzero `a` replaces `eta(ka)` by `eta(ka)+s=f(ka)`. For `p` dividing `k`, its centered contribution is zero, exactly represented by `f(0)=0`. Uniform absolute convergence permits this averaging. This proves equation (1), including its factor `2 pi^2`.

For the Dirichlet inverse, put `u_k=1/k^2` and `v_k=mu(k)/k^2`. Then

\[
(u*v)_k=k^{-2}\sum_{d\mid k}\mu(d)=\mathbf1_{k=1}.
\]

The claimed factorization `c=-u*(delta_1-delta_2/2)` is exact: for even `k`, its coefficient is `-1/k^2+(1/2)/(k/2)^2=1/k^2`; for odd `k`, it is `-1/k^2`. The second factor has its absolutely convergent geometric inverse because its nonidentity part has `l^1` norm `1/2`.

Consequently, for `k=2^t m` with `m` odd,

\[
b_k=-\sum_{j=0}^{t}2^{-j}
\frac{\mu(2^{t-j}m)}{(2^{t-j}m)^2}.
\]

For `t>=1`, only `j=t,t-1` can contribute, giving exactly `-mu(m)/(2^(t+1)m^2)`. For `t=0`, it gives `-mu(m)/m^2`. This establishes (2) for all integers, including nonsquarefree odd parts.

The Euler product and the independently rechecked [DLMF zeta values](https://dlmf.nist.gov/25.6.E1) give

\[
\sum_{m\text{ odd}}\frac{\mu(m)^2}{m^2}
=\frac{\zeta(2)}{\zeta(4)(1+2^{-2})}
=\frac{12}{\pi^2}.
\]

The powers-of-two weight is `1+sum_(t>=1)2^(-t-1)=3/2`. Therefore `sum |b_k|=18/pi^2`. Since `p` is odd, excluding its multiples removes exactly the Euler factor `1+p^(-2)`, yielding `18p^2/[pi^2(p^2+1)]`. All three masses in (3) are correct.

For completeness, the inverse remains valid after reducing the arguments modulo `p`. On functions extended by zero at zero, the operators `U_k` obey `U_k U_l=U_(kl)`; `U_k=0` if `p|k`, and otherwise `U_k` is an isometry for every norm under review. The double series is bounded by `||b||_1||c||_1||f||_r`. It can therefore be regrouped according to `k*l`, producing `T_b T_c=T_(b*c)=I`. This proves (4) without assuming a faithful representation of the integer semigroup.

## 3. Constants and every norm exponent

Define

\[
\|v\|_r=\left((p-1)^{-1}\sum_{a\ne0}|v(a)|^r\right)^{1/r},
\qquad 1\le r<\infty,
\]

and use the maximum norm when `r=infinity`. Minkowski, followed by the permutation property of `U_k` for `p` not dividing `k`, gives

\[
2\pi^2\|d\|_r
\le\frac{\pi^2}{6}(1-p^{-2})\|f\|_r,
\]
\[
\|f\|_r
\le2\pi^2\frac{18p^2}{\pi^2(p^2+1)}\|d\|_r.
\]

Thus the complete comparison, for every `1<=r<=infinity`, is

\[
\boxed{\frac{12p^2}{p^2-1}\|d\|_r
\le\|f\|_r
\le\frac{36p^2}{p^2+1}\|d\|_r.}
\]

If `2A=A`, both functions are invariant under multiplication of the argument by 2. In the forward series, each odd-part block has coefficient

\[
-m^{-2}+\sum_{t\ge1}(2^tm)^{-2}=-\frac2{3m^2}.
\]

Its mass over odd `m` not divisible by `p` is `pi^2(1-p^(-2))/12`. The lower constant consequently doubles to `24p^2/(p^2-1)`. The inverse upper constant `36p^2/(p^2+1)` remains valid. This proves both (5) and (6) with no dependence on `r`.

There is no hidden restriction to integer or fixed `r`: real exponents `r>=1` and exponents growing with `p` are covered. When phrased as unrooted moments, the corresponding comparison constants are raised to the power `r`; the uniform assertion is about the norms. Moving between `f` and `eta` incurs the additive shift `s`, as above. In a fixed quartic window, that shift is harmless for the stated asymptotic scale, but these comparisons provide no bound on either function by themselves.

## 4. Lattice and cyclotomic norm

The homomorphism `Z^N -> F_p`, `z -> sum z_j h_j`, is onto because at least one `h_j` is nonzero. Hence its kernel `L` has index `p`. The displayed candidate for `L^*` pairs integrally with `L`. Its added vector has exact order `p` modulo `Z^N`, so its covolume is `1/p`, equal to that of the dual. This proves the dual-lattice identity.

In the coset represented by `a(h_1,...,h_N)/p`, the minimization is over independent integer coordinates. Choosing the centered residue in each coordinate therefore gives shortest squared norm exactly `V(a)`. This part is valid for all `a`, including zero.

Now set `A=H=<g>`, with `n=2N` a power of two, `N>=2`, and **fix `a != 0`**. Then `p>n`, since `n|p-1`. The following steps of the norm argument are all valid:

1. `(X+1)^N+1` has leading coefficient 1, constant coefficient 2, and all intermediate coefficients even. Eisenstein at 2 proves irreducibility of `X^N+1` over `Q`.
2. The roots modulo `p` are exactly `g^k` for odd `k` with `1<=k<n`, and they are distinct. At each such root, the geometric sum for `F_a(g^k)` vanishes unless `k=n-1`; at that root it equals `aN != 0`. Evaluation identifies `R/pR` with the product of these `N` copies of `F_p`, so multiplication has rank exactly one.
3. `F_a` is nonzero and has degree less than `N`. Irreducibility makes its multiplication determinant over `Q` nonzero. Smith normal form then shows that at least `N-1` integer invariant factors are divisible by `p`, proving `p^(N-1)|Norm(F_a)`.
4. The orthogonality sum over roots is `N^(-1) sum_(k odd) zeta_n^(k(j-l))=1_(j=l)` for `0<=j,l<N`. Consequently, `sum_(k odd)|F_a(zeta_n^k)|^2=N S`, where `S=sum r_j^2`. Applying AM-GM to the `N` squared absolute values gives `|Norm(F_a)|^2<=S^N`.

Combining the last two facts yields

\[
S\ge p^{2-2/N},\qquad V(a)=S/p^2\ge p^{-2/N}.
\]

There is no missing power of `p` or `N` in (7). Also `g^2` has order `N>=2`, so `S=a^2 sum_(j=0)^(N-1)g^(2j)=0 mod p`. Thus the `1/p` lattice of possible values of `V` is correct. Keeping the ordered powers and their actual signed centered coordinates is necessary for this polynomial calculation.

For `p` in a fixed quartic window, `p^(-2/N)=exp(-2 log(p)/N)` tends to 1, whereas `vbar~N/12`. The target width `O(sqrt(n log n))` is `o(N)`. A lower bound tending to 1 supplies neither the required lower edge near `N/12` nor a corresponding upper edge. The review therefore confirms the source's stated limitation; it does not establish the missing annulus estimate.

## Finite cross-checks and provenance

I read the original verifier's analytic and norm-check implementations, its recorded results, the source-scope record, the pass audit, the NIST transcription, and the target statement. The old record reports 41,574 exact assertions, 4,096 recursively recovered inverse coefficients, six symmetric-set cases, and 28 dyadic norm cases. These were treated as prior finite evidence, not as a substitute for the derivations above. This review did not rerun the original verifier's `main`, invoke its C++ executable, or overwrite its result files.

Separate bounded calculations performed in this review were:

- Direct rational divisor convolution `sum_(d|k)c_d b_(k/d)=1_(k=1)` for every `1<=k<=256`, using trial-factorization Mobius values rather than the original sieve/recurrence.
- Exact nonzero-frequency means for all 105 nonempty symmetric subsets in the fields `p=3,5,7,11,13`.
- The stated `p=1153,n=8,g=75,a=1` example, with the determinant computed by the 24-term Leibniz formula for its 4-by-4 multiplication matrix. Its centered coordinates are `[1,75,-140,-123]`, its root evaluations modulo `p` are `[0,0,0,4]`, and its determinant is `1,532,808,577=1153^3`. The square sum is `40,355=35p`, so `V=35/1153` and `vbar=1154/3459`. The ratio to the norm lower bound is exactly `35/sqrt(1153)`, approximately `1.0307501121`, confirming the stated 3.1% figure.

All these checks passed. Their finite character does not certify an unbounded family; the universal conclusions reviewed here follow from the explicit arguments above. Sections 5-7 and the large finite certificate were outside this lane's mathematical audit. No Lean build, submission, external message, or full-goal claim was made.
