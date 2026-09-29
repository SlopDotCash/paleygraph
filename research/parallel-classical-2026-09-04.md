# Classical Paley: biased rows sharpen the necessary logarithmic allowance

This bounded pass does **not** prove or refute the logarithmic moment
hypothesis (LM), or prove the Paley conjecture. It gives a new refinement
within this project of the existing forced-row obstruction. No claim of
novelty in the literature is made. Everything below is proved directly;
no external theorem about character-sum cancellation is imported.

Write

\[
F_B(a)=\sum_{b\in B}\chi(a-b),\qquad
M_{2r}(B)=\sum_{a\in\mathbb F_p}F_B(a)^{2r}.
\]

The hypothesis under examination, for every fixed integer \(r\ge2\), is

\[
M_{2r}(B)\le C_r\bigl(p|B|^r+|B|^{2r}\log p\bigr).
\tag{LM}
\]

The original all-positive common-neighborhood construction proves that the
coefficient of the second term must grow logarithmically in some ranges.
Allowing a small fraction of negative entries gives a strictly larger
asymptotic coefficient. This strengthens a necessary condition on (LM),
without giving an unbounded ratio that would disprove it.

## 1. An exact finite construction

Let \(p\) be an odd prime, \(d=(p-1)/2\), \(1\le k\le d\), and
\(k/2<j\le k\). For a row set \(U\) of size \(k\), let

\[
N_j(U)=\left\{b\notin U:
   |\{a\in U:\chi(a-b)=1\}|\ge j\right\}.
\]

Define

\[
P_{p,k,j}=
\frac{\sum_{i=j}^k\binom di\binom d{k-i}}{\binom pk}.
\]

Then some \(U\) has \(|N_j(U)|\ge\lceil pP_{p,k,j}\rceil\). For every integer
\(1\le n\le\lceil pP_{p,k,j}\rceil\), some set \(B\) of size \(n\) satisfies

\[
\boxed{M_{2r}(B)\ge
       k\left(\frac{2j-k}{k}\right)^{2r}n^{2r}.}
\tag{1}
\]

**Proof.** Each column \(b\) has exactly \(d\) positive, \(d\) negative,
and one zero entry. The number of row sets \(U\) avoiding the zero and
having exactly \(i\) positive entries is
\(\binom di\binom d{k-i}\). Double counting therefore gives

\[
\sum_{|U|=k}|N_j(U)|
   =p\sum_{i=j}^k\binom di\binom d{k-i}.
\tag{2}
\]

Choose a maximizing \(U\), then any \(n\)-element \(B\subseteq N_j(U)\).
For each \(b\in B\),
\(\sum_{a\in U}\chi(a-b)\ge2j-k\). Hence

\[
\sum_{a\in U}F_B(a)\ge(2j-k)n.
\]

Convexity of \(x\mapsto x^{2r}\) now gives

\[
M_{2r}(B)\ge\sum_{a\in U}F_B(a)^{2r}
 \ge k\left(k^{-1}\sum_{a\in U}F_B(a)\right)^{2r},
\]

which proves (1). Individual rows need not all have positive bias; only
their average is used. The construction works for both prime residue
classes modulo four, and more generally for any square sign matrix whose
columns have these same counts.

Taking \(j=k\) recovers the all-positive construction already in
[the earlier moment note](moments-and-obstructions.md). The case \(j<k\)
is the extra flexibility here.

## 2. A binomial estimate with an explicit finite-field error

For \(0\le j\le k\le d\),

\[
\frac{\binom dj\binom d{k-j}}{\binom pk}
\ge \left(1-\frac{k^2}{p}\right)2^{-k}\binom kj.
\tag{3}
\]

If the right side is negative the assertion is automatic. For a positive
right side, and in fact directly for all \(k\le d\), expand falling
factorials. The ratio of the left side to \(2^{-k}\binom kj\) is

\[
\frac{
 \prod_{a=0}^{j-1}(1-(2a+1)/p)
 \prod_{a=0}^{k-j-1}(1-(2a+1)/p)}
 {\prod_{a=0}^{k-1}(1-a/p)}.
\]

The denominator is positive and at most one. Applying
\(\prod(1-u_a)\ge1-\sum u_a\) to the numerator bounds it below by
\(1-[j^2+(k-j)^2]/p\ge1-k^2/p\).

Put

\[
I(t)=\frac{(1+t)\log(1+t)+(1-t)\log(1-t)}2,
\quad 0\le t\le1,
\]

where \(0\log0=0\). The elementary type bound

\[
\binom kj\ge\frac1{k+1}
 \exp\bigl(kH(j/k)\bigr)
\]

follows by considering the modal atom of a binomial variable with success
probability \(j/k\); its probability is at least \(1/(k+1)\).
Here \(H(q)=-q\log q-(1-q)\log(1-q)\). Thus, with
\(t_k=2j/k-1\), (3) yields

\[
\boxed{P_{p,k,j}\ge
 \frac{1-k^2/p}{k+1}\exp(-kI(t_k)).}
\tag{4}
\]

This estimate uses only one atom of the upper tail, so it is a lower
bound, not a claimed asymptotic evaluation of the whole tail.

There is also a quantitative fixed-density version. Put

\[
q_{k,j}=2^{-k}\sum_{i=j}^k\binom ki.
\]

Summing (3) shows \(P_{p,k,j}\ge(1-k^2/p)q_{k,j}\). Thus, for any
fixed \(0<c<q_{k,j}\), every odd prime satisfying

\[
p\ge\max\left(2k+1,\frac{k^2q_{k,j}}{q_{k,j}-c},\frac1c\right)
\]

admits \(B\) of exactly \(n=\lfloor cp\rfloor\ge1\) elements with

\[
M_{2r}(B)\ge k((2j-k)/k)^{2r}n^{2r}.
\tag{FD}
\]

This states the density/bias tradeoff without a power-scale assumption.

## 3. The sharper asymptotic obstruction

Fix \(r\ge2\) and \(t\in(0,1)\). Consider any sequence of odd primes
and integers \(1\le n\le p\) such that \(p/n\to\infty\). Put
\(L=\log(p/n)\to\infty\) and choose

\[
k=\left\lfloor
\frac{L-2\log L}{I(t)}
\right\rfloor,
\qquad j=\left\lceil\frac{1+t}{2}k\right\rceil.
\]

Eventually these are admissible: \(k\to\infty\),
\(k=O_t(\log p)=o(\sqrt p)\), and \(k/2<j\le k\). Since \(t\) is
fixed, \(t_k=t+O_t(k^{-1})\) and
\(kI(t_k)=kI(t)+O_t(1)\). Equation (4) implies

\[
pP_{p,k,j}\ge
 \frac{p(1-o(1))}{k+1}e^{-kI(t)-O_t(1)}
 \gg_t nL\ge n.
\]

Use (1), and note that
\(k=L/I(t)+O_t(\log L)\). This proves the existence of sets of exactly
\(n\) elements satisfying

\[
M_{2r}(B)\ge
 \left(\frac{t^{2r}}{I(t)}+o(1)\right)
 n^{2r}\log(p/n).
\tag{5}
\]

Define

\[
D_r=\max_{0<t\le1}\frac{t^{2r}}{I(t)}.
\]

The quotient tends to zero at \(t=0\) for \(r\ge2\), so a maximum
exists. Equation (5) holds with \(D_r\) in place of the quotient by
choosing a maximizing \(t\), which lies strictly inside \((0,1)\), as
shown below. In particular, the general tradeoff is

\[
\boxed{M_{2r}(B)\ge(D_r+o(1))n^{2r}\log(p/n)
\quad\text{for some }|B|=n,\quad p/n\to\infty.}
\tag{BT}
\]

Specializing to \(n=\lfloor p^\alpha\rfloor\), fixed \(0<\alpha<1\),
gives the coefficient \((1-\alpha)D_r\) in front of \(n^{2r}\log p\).
If \(\alpha r>1\), then
\(pn^r=o(n^{2r}\log p)\). Consequently any constant in (LM) must satisfy

\[
\boxed{C_r\ge(1-1/r)D_r.}
\tag{6}
\]

Indeed (5) gives \(C_r\ge(1-\alpha)D_r\) for every \(\alpha>1/r\);
let \(\alpha\) decrease to \(1/r\). This is a lower bound on the
possible constant, not a contradiction to its existence.

### Strict improvement and numerical size

At \(t=1\), the quotient is \(1/\log2\), the earlier all-positive
coefficient. Write \(t=1-u\). Direct expansion gives

\[
I(1-u)=\log2+\frac u2\log(u/2)-\frac u2+O(u^2),
\]

and hence

\[
(1-u)^{2r}\log2-I(1-u)
=u\left(-2r\log2-\frac12\log(u/2)+\frac12\right)+O_r(u^2)>0
\]

for all sufficiently small positive \(u\). Therefore
\(D_r>1/\log2\) for each fixed \(r\ge2\).

An interior maximizer satisfies

\[
2rI(t)=t\operatorname{atanh}(t).
\tag{7}
\]

There is a unique solution in \((0,1)\). To see this, set
\(g(t)=2rI(t)-t\operatorname{atanh}(t)\). Then

\[
g'(t)=(2r-1)\operatorname{atanh}(t)-\frac t{1-t^2},
\quad
g''(t)=\frac{(2r-2)-2rt^2}{(1-t^2)^2}.
\]

Thus \(g'\) first increases from zero and then decreases to minus
infinity, with exactly one zero inside the interval. The function \(g\)
first increases from zero and then decreases to minus infinity, again
with exactly one interior zero. Its sign is the derivative sign of the
quotient being maximized.

The script gives the following **numerical**, uncertified approximations;
the proofs of (5)–(7) do not rely on them.

| r | maximizing t | D_r | necessary C_r from (6) |
|---|---:|---:|---:|
| 2 | 0.9906109323 | 1.4517973364 | 0.7258986682 |
| 3 | 0.9995001183 | 1.4432101104 | 0.9621400736 |
| 4 | 0.9999694028 | 1.4427268460 | 1.0820451345 |

The improvement over \(1/\log2=1.4426950409\ldots\) is small. It is a
strict improvement to an obstruction coefficient, not evidence that a
proof of (LM) is close.

## 4. Consequences for the research strategy

The old all-positive example is not the extremal choice even among these
elementary row-pattern constructions. Searches intended to refute (LM)
should include biased Hamming balls, rather than only common
neighborhoods. Nevertheless, the best coefficient in this family is
finite for every fixed \(r\). This family still has exactly logarithmic
growth and therefore does not eliminate (LM).

The implication from (LM) to Paley in
[the frontier](frontier.md#conditional-implication-with-exponents) checks
out: with \(r\varepsilon\ge1\), Hölder gives normalized discrepancy
at most \([C_rp^{-\varepsilon}(1+\log p)]^{1/(2r)}\), and eventually
\(p^{-\varepsilon/(4r)}\). Nothing here establishes the required upper
bound on any new worst-case small-set range.

There is a useful calibration of how strong (LM) already is at fixed
depth. If \(B\) is a clique of size \(n\) in a prime Paley graph, then
\(F_B(a)=n-1\) for every \(a\in B\), so (LM) would imply

\[
n(n-1)^{2r}\le C_r(pn^r+n^{2r}\log p).
\]

For \(n\ge2\), division by \(n^{2r}\) gives
\(n/2^{2r}\le C_r(p/n^r+\log p)\). Splitting according to which term
on the right is larger yields

\[
n=O_r\bigl(p^{1/(r+1)}+\log p\bigr),
\]

where the implicit constant also depends on the hypothesized \(C_r\).
Thus the fourth-moment instance would already give a cubic-root clique
bound. This is a conditional consequence, not an unconditional estimate.

## 5. Reproducibility and limitations

Run `python3 experiments/parallel_classical_2026_09_04.py`.
The result file records exhaustive small-prime checks of (2), concrete
sets witnessing (1), their exact integer moments, and integer checks of
(3). It contains 52,074 row subsets, 33 exact double-counting identities,
99 nonempty-set moment inequalities, and 1,484 hypergeometric checks.
The coefficient table is separately labeled numerical.

No finite search establishes the upper bound (LM), an asymptotic
counterexample to (LM), the full Paley conjecture, or a prize reduction.
The remaining obstacle in this lane is still an upper bound for
worst-case sets, rather than a missing computation in this construction.
