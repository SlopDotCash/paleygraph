# Polynomial-size sets: a quantified directional bound from Burgess

This pass obtains an actual upper bound in a polynomial-size range.
For each fixed \(0<\beta<1/15\), arithmetic progressions of size
\(p^{13/30+\beta}\le n\le\sqrt p\), and sets sufficiently close to
them, satisfy

\[
\boxed{|\langle W_B,E\rangle|
 \ll_\beta p^{1/2-\beta}n^4.}
\tag{A}
\]

This saves a factor \(p^\beta\) over the termwise elliptic Weil scale.
For example, \(n\ge p^{9/20}\) allows \(\beta=1/60\) and the
upper bound \(O(p^{29/60}n^4)\).
It uses the established Burgess theorem, rather than proving a new
character-sum theorem. The transfer to the directional pairing and the
perturbation/range audit are proved below. A rank-two progression
corollary is also recorded. These bounds do not establish (LM), the full
Paley conjecture, or a size-restricted counterexample to either.

All fields here are prime fields. No existence statement about
arbitrarily large primes is used as polynomial control of their size.

## 1. A nonnegative remainder in the small-set range

Retain the definitions from the
[exact elliptic reduction](parallel2-classical-2026-09-04.md):
\(n=|B|\), \(F_B(x)=\sum_{b\in B}\chi(x-b)\),
\(A_B=\sum_{b\in B}F_B(b)^2\), and

\[
M_4(B)=R_{p,n,B}+D_B,\qquad D_B=\langle W_B,E\rangle,
\]

where

\[
R_{p,n,B}=p(3n^2-2n)-6n^3+14n^2-9n-6A_B
             +\frac{3(n)_4}{p-2}.
\]

The earlier note proves \(R_{p,n,B}\le3pn^2\). There is an additional
lower bound in the relevant range:

\[
\boxed{0\le R_{p,n,B}\le3pn^2
       \quad(p\ge5,\ n\le\sqrt p).}
\tag{1}
\]

**Proof of the lower bound.** Since \(A_B\le n^3\) and \((n)_4\ge0\),

\[
R_{p,n,B}\ge n\bigl(p(3n-2)-12n^2+14n-9\bigr).
\]

For \(n\ge4\), use \(p\ge n^2\). The expression in parentheses
is at least

\[
3n^3-14n^2+14n-9
=(n-4)(3n^2-2n+6)+15>0.
\]

For \(n=3\), the same bound is \(3(7p-75)\ge6\), since the least
eligible prime is 11. Directly, \(R=0\) when \(n=0\), \(R=p-1\)
when \(n=1\), and \(R=8p-22\) when \(n=2\); in the latter case
\(A_B=2\).

Let \(U_B=\max_x|F_B(x)|\). The exact second moment
\(M_2(B)=pn-n^2\) now yields

\[
M_4(B)\le U_B^2(pn-n^2),\qquad
\boxed{|D_B|\le\max\{3pn^2,\ U_B^2(pn-n^2)\}.}
\tag{2}
\]

The maximum in (2) uses both signs correctly: \(D_B\le M_4(B)\)
and \(D_B\ge-3pn^2\). In particular, a pointwise translated
character bound gives directional cancellation without bounding the
entire coefficient vector or its Fourier-L1 norm.

## 2. Imported character bounds and their exact scope

We use two results in the pinned primary source
[Alsetri–Shao, arXiv:2509.07765v1](https://arxiv.org/html/2509.07765v1#S1.SS1):
the classical Burgess estimate stated in its introduction, and its
Theorem 1.1. Burgess gives, for an arithmetic progression \(P\), every
integer \(r\ge2\), and every \(\eta>0\),

\[
\max_x\left|\sum_{b\in P}\chi(x-b)\right|
 \ll_{r,\eta}|P|^{1-1/r}p^{(r+1)/(4r^2)+\eta}.
\tag{B}
\]

Theorem 1.1 states that a proper rank-two generalized arithmetic
progression of size at least \(p^{1/4+\varepsilon}\) has a character
sum bounded by \(C_\varepsilon |P|p^{-\delta(\varepsilon)}\), for
some positive \(\delta(\varepsilon)\). Proper means that its box
parametrization is injective. Translations and nonzero dilations
preserve properness and size, so the bound applies uniformly to
\(x-P\). We use the stated existence of \(\delta\), without claiming
a numerical value extracted from that proof.

Both statements are imported inputs, not results of the present
investigation. For arithmetic progressions in (B), an affine change
of variables reduces to the source's interval formulation; the scalar
character factor has absolute value one.

## 3. Arithmetic progressions and growing perturbations

Fix \(0<\beta<1/15\). Let \(B\subseteq\mathbb F_p\) have size

\[
p^{13/30+\beta}\le n\le\sqrt p,
\]

and suppose there is an arithmetic progression \(P\) such that

\[
|B\mathbin\triangle P|\le n p^{-1/30}.
\tag{3}
\]

The progression has nonzero step and distinct elements. Condition (3)
allows a number of arbitrary edits growing polynomially with \(p\).
The case \(B=P\) is included.

Put \(k=|B\mathbin\triangle P|\). For every \(x\),

\[
|F_B(x)|\le|F_P(x)|+k.
\]

Also \(|P|\le n+k\le2n\). Apply (B) with \(r=3\) and
\(\eta=\beta/3\) to obtain

\[
\begin{aligned}
U_B
&\ll_\beta n^{2/3}p^{1/9+\beta/3}+np^{-1/30}\\
&\ll_\beta np^{-1/30}.
\end{aligned}
\tag{4}
\]

The exponent in the second line follows exactly from
\(n^{-1/3}\le p^{-13/90-\beta/3}\). This gives the genuine
uniform discrepancy bound

\[
\boxed{|S(A,B)|\ll_\beta |A|n p^{-1/30}
       \quad\text{for every }A\subseteq\mathbb F_p.}
\tag{5}
\]

There is no size restriction on \(A\). Using (2) and (4),

\[
M_4(B)\ll_\beta p^{14/15}n^3,
\qquad
|D_B|\ll_\beta pn^2+p^{14/15}n^3.
\tag{6}
\]

Divide the last expression by \(\sqrt p\,n^4\). Its two terms are
bounded respectively by

\[
\frac{\sqrt p}{n^2}\le p^{-11/30-2\beta},\qquad
\frac{p^{13/30}}n\le p^{-\beta}.
\]

This proves (A), with all constants independent of \(p,B,P\).
The condition \(\beta<1/15\) is precisely what makes the stated
size window nonempty at the exponent level.

### Why the threshold is 13/30 for this transfer

For an unperturbed progression, combining (B) at general integer
\(r\ge2\) with (2) makes the Burgess term smaller by a power than
\(\sqrt p\,n^4\) when

\[
\log_p n>
\theta_r:=\frac{r^2+r+1}{2r(r+2)}
\]

by a fixed margin. The arbitrarily small \(\eta\) is absorbed by
that margin. Among integers \(r\ge2\), the minimum is attained
uniquely at \(r=3\), because

\[
\theta_r-\frac{13}{30}
=\frac{(2r-5)(r-3)}{30r(r+2)}\ge0.
\]

Thus 13/30 is the optimum of this particular second-moment times
pointwise-Burgess argument. It is not claimed to be a mathematical
barrier for other approaches.

## 4. A rank-two polynomial-size window

Apply Theorem 1.1 with \(\varepsilon=1/8\), and let its exponent be
\(\delta_*>0\). Set
\(\sigma=\min(\delta_*,1/16)>0\). Every proper rank-two progression
of size

\[
p^{1/2-\sigma}\le n\le\sqrt p
\]

satisfies \(U_B\ll np^{-\sigma}\): its size exceeds the required
\(p^{3/8}\) threshold, since \(1/2-\sigma\ge7/16\).
Equation (2) then proves

\[
\boxed{|D_B|\ll p^{1/2-\sigma}n^4.}
\tag{7}
\]

Indeed, division by \(\sqrt p\,n^4\) bounds the two terms by
\(O(p^{-1/2+2\sigma})\) and \(O(p^{-\sigma})\). The former is
smaller because \(\sigma\le1/16\). The same source also directly
gives \(|S(A,B)|\ll |A|n p^{-\delta_*}\) for every \(A\).

This is a rigorously quantified existence of a polynomial-size window;
its numerical width is unspecified because the imported exponent is.
No claim for arbitrary rank or arbitrary subsets of such a progression
follows from that theorem.

## 5. What has and has not advanced

There is now a directional power saving for explicit growing classes
below the square-root size threshold. The estimate survives the growing
arbitrary perturbations in (3), and it avoids the all-size norm
hypotheses refuted in the earlier passes. The cancellation itself comes
from existing Burgess results; the Kloosterman transform does not
independently improve their exponents here.

This still does not prove the moment hypothesis. In the interior of
\(n\le\sqrt p\), (LM) would essentially require \(M_4\ll pn^2\)
up to the stated logarithmic term. The upper bound in (6) can exceed
\(pn^2\) by the factor \(n p^{-1/15}\), which grows polynomially
throughout our window. A saving over termwise Weil is therefore not a
Gaussian fourth-moment estimate. Nor does one fixed fourth-moment bound
settle all exponents in the full Paley conjecture.

No counterexample in a fixed range \(n\ge p^\varepsilon\) was
constructed. The earlier prime-family obstructions remain limited to
their stated all-size quantifiers. The unresolved bridge is extending
controlled cancellation from these structured classes to arbitrary
sets, or obtaining a sufficiently strong directional estimate directly.

The accompanying script checks (1)–(2), exact cross-ratio pairings, the
perturbation transfer with actual finite row maxima, preservation of
proper progressions under translation, and the rational exponent
calculations. It does not use finite experiments to estimate the
unspecified constants or establish the imported asymptotic theorems.
The completed run contains 3,869 exhaustive remainder/directional
cases, 20 progression cases including eight perturbations, and checks
of 99 integer Burgess parameters and four rational exponent choices.
