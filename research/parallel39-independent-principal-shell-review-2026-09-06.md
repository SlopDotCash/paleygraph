# Independent review of the principal-prime short-coset construction

Date: 2026-09-06. This review checks the proposed algebraic implication, its norm-cofactor extension, and independently checks the existing order-128 example. Known sublinear cancellation excludes an unbounded bounded-cofactor family in a fixed quartic window; this review makes no asymptotic disproof claim.

## Verdict

The construction is valid with **`n=2N` a power of two and `n>=4`**, `p` prime with `p=1 mod n`, and `g` of exact order `n`. For an integral `f` in `R=Z[zeta_n]` satisfying `Norm(f)=p` and `f(g)=0 mod p`, it produces a nonzero frequency `a` for `H=<g>` with

\[
V(a)\le1,\qquad \eta(a)\ge n-4\pi^2.
\]

The exact frequency is

\[
a=N^{-1}\alpha(g)\pmod p,
\qquad \alpha=p\overline f/f.
\]

If `B=p/f` denotes the adjugate, the numerator is `alpha(g)=bar(f)(g) B(g)`, **not** `B(g)` alone. Use the ordered generator `h=g^(-1)` to match the coefficient pattern from pass 36. This changes the choice of coordinate representatives, not the subgroup or the period being bounded.

For the supplied `p=215535361,n=128` example, all the coefficients of `alpha` are already centered. Thus `V(a)=1` exactly, at `a=38468180`. Exact inequalities give a finite counterexample to the literal bound

\[
M\le2\sqrt{n\log(p/n)}.
\]

## 1. Integrality and the exact nonzero residue

For dyadic `n>=4`, `R` has the integral basis `1,zeta,...,zeta^(N-1)` and defining polynomial `X^N+1`. Let `P_u` be the kernel of evaluation `zeta -> u mod p`, where `u` ranges over the `N` distinct roots of `X^N+1` in `F_p`.

The matrix of multiplication by `f` has determinant `Norm(f)=p`, so the subgroup `fR` has index `p` in the additive group of `R`. Since `f(g)=0`, one has `fR` contained in `P_g`. Both have index `p`, hence

\[
fR=P_g.
\]

This proves `p/f` is integral: `p` belongs to `P_g=fR`. Put `B=p/f` in `R`, so `fB=p`.

The residue of `f` can vanish at no other root `u`: otherwise the same index argument would give `fR=P_u=P_g`, contradicting distinct evaluation kernels. Evaluating `fB=p` shows `B(u)=0` at every `u!=g`.

Moreover `B(g)!=0 mod p`. If it were zero, the degree-below-`N` representative of `B` would vanish at all `N` distinct roots, so `B` would belong to `pR`. Writing `B=pC` would give `fC=1` after canceling `p` in the characteristic-zero domain `R`, contrary to evaluation at `g`.

Conjugation sends `zeta` to `zeta^(-1)`. Thus

\[
\overline f(g)=f(g^{-1})\ne0\pmod p,
\]

because `g^{-1}!=g` for exact order `n>=4`. Therefore `alpha=bar(f)B` is integral and its reduction modulo `p` is nonzero at `g` and zero at every other root.

Equivalently, in split-prime ideal notation,

\[
(\alpha)=P_{g^{-1}}\prod_{u\ne g}P_u.
\]

The elementary argument above establishes the required integrality and residues without relying on an unproved ideal-factorization assertion.

## 2. Archimedean norm and the frequency map

For every complex embedding `sigma`, integer coefficients give `sigma(bar(f))=overline(sigma(f))`. Since `f` is nonzero,

\[
|\sigma(\alpha)|=p\frac{|\overline{\sigma(f)}|}{|\sigma(f)|}=p.
\]

Write `alpha=sum_(j=0)^(N-1) A_j zeta^j` with integer coefficients. The orthogonality of the `N` primitive dyadic roots gives

\[
\sum_jA_j^2=\frac1N\sum_{k\text{ odd}}|\alpha(\zeta^k)|^2=p^2.
\]

Now put `h=g^(-1)` and `a=alpha(g)/N mod p`. The denominator is invertible because `N<p`, and `a!=0` by the preceding residue calculation. The polynomial

\[
P_a(X)=a\sum_{j=0}^{N-1}h^jX^j
\]

takes value `aN=alpha(g)` at `g` and vanishes at the other roots. Since both it and the representative of `alpha` have degree below `N`, equality at all roots proves

\[
A_j\equiv ah^j\pmod p\quad(0\le j<N).
\]

In particular `A_0=a mod p`. Let `r_j` be the centered residue of `A_j`; it is exactly the pass-36 coordinate for this nonzero frequency and generator `h`. Centering minimizes the absolute value among integer representatives, so

\[
\sum_jr_j^2\le\sum_jA_j^2=p^2,
\qquad V(a)\le1.
\]

Because `h^N=-1`, the full subgroup is the disjoint union of these representative pairs. The elementary inequality `cos x>=1-x^2/2` then yields

\[
\eta(a)=2\sum_j\cos(2\pi r_j/p)
\ge 2N-4\pi^2\sum_jr_j^2/p^2
\ge n-4\pi^2.
\]

No claim that arbitrary sign changes preserve the polynomial norm is used: the construction first fixes `alpha`, then centers its coefficients, and finally identifies the resulting subgroup coordinates modulo `p`.

## 3. Necessary edge-case restriction

If the proposed theorem allows `n=2`, it is false at the step asserting a nonzero residue and frequency. Take `R=Z`, `zeta_2=-1`, any odd prime `p`, `g=-1`, and `f=p`. Then `Norm(f)=p` and `f(g)=0 mod p`, but `alpha=p` reduces to zero. Centering its only coefficient gives zero, so it does not produce a nonzero `a`.

The restriction `n>=4` is sufficient and is already the domain used in pass 36's cyclotomic-norm argument. This edge case does not affect the order-128 construction.

## 4. Independent exact check of the finite example

Input: the `p,n,g,adjugate` fields of [parallel31_short_multiples_2026_09_05.json](../results/parallel31_short_multiples_2026_09_05.json), with

\[
p=215535361,\quad n=128,\quad N=64,\quad
g=25525303,\quad f=1+X+X^{19}.
\]

The following were freshly checked using exact integers, separately from the original recursive norm-descent computation:

- Primality by trial division through `floor(sqrt(p))`; `n|(p-1)`; `g^64=-1 mod p`; and `f(g)=0 mod p`.
- `Norm(f)=215535361` by constructing the 64-by-64 negacyclic multiplication matrix and evaluating its determinant with fraction-free elimination, checking exact division at every step.
- Direct negacyclic multiplication of the archived `B` verifies `fB=p`.
- With `alpha=bar(f)B`, direct multiplication verifies `alpha*bar(alpha)=p^2`. Independently summing the squares of its 64 coefficients gives `46455491841400321=p^2`.
- The maximum coefficient magnitude is `92675172`, below `floor(p/2)=107767680`, so every coefficient is already centered.
- All root evaluations other than `g` vanish modulo `p`; `alpha(g)=91074549`. Thus `a=alpha(g)/64=38468180 mod p`. With `h=g^(-1)=165976977`, all 64 coefficients satisfy `A_j=a h^j mod p`.

These establish `V(a)=1` for the actual nonzero frequency. The prime lies in the fixed quartic window

\[
67108864=n^4/4\le p\le n^4=268435456.
\]

The constant-2 obstruction can be proved without numerical trigonometry or logarithms. Positivity of

\[
\int_0^1\frac{x^4(1-x)^4}{1+x^2}\,dx=\frac{22}{7}-\pi
\]

gives `pi<22/7`, hence

\[
\eta(a)\ge128-4\pi^2>\frac{4336}{49}.
\]

Also `e>1+1+1/2+1/6=8/3`, while exact integer arithmetic gives

\[
128\cdot8^{15}-p\cdot3^{15}=1410902777170069>0.
\]

Therefore `log(p/128)<15`. Finally,

\[
\left(\frac{4336}{49}\right)^2-4\cdot128\cdot15
=\frac{361216}{2401}>0.
\]

It follows that `eta(a)>2 sqrt(128 log(p/128))`. This is a rigorous finite counterexample to that literal numerical constant, with no assertion that the displayed frequency is the global maximizer.

## 5. The cofactor extension and known asymptotic exclusion

I also independently checked the extension in [the shell agent's structural note](parallel39-shell-structural-input-2026-09-06.md). Suppose instead that `Norm(f)=d=kp`, where `k>=1` and `p` does not divide `k`. Then the multiplication determinant has `p`-adic valuation one. Its reduction has rank `N-1`: two zero invariant factors modulo `p` would force `p^2` to divide the determinant. Since evaluation at the split roots diagonalizes multiplication, `g` is still the unique root at which `f` vanishes.

The product of all nonidentity Galois conjugates,

\[
B=\prod_{u\text{ odd}\bmod n,\,u\ne1}f(\zeta^u)=d/f,
\]

belongs to `R` and is nonzero upon evaluation at `g`, since none of those factors vanish there. Put `alpha=d bar(f)/f`. The same residue, interpolation, and coefficient-norm calculations now give

\[
\alpha\overline\alpha=d^2,\quad
\sum_jA_j^2=d^2,\quad
V(a)\le k^2,\quad
\eta(a)\ge n-4\pi^2k^2.
\]

The associated negacyclic matrices satisfy `M_alpha M_alpha^T=d^2 I`, because conjugation is the transpose involution in this coefficient basis. Thus the structural note's higher normalized Gram-trace statement is also valid.

For **any established** upper bound `M<=U(n,p)`, the displayed lower bound forces

\[
\boxed{k^2\ge\max\left(0,\frac{n-U(n,p)}{4\pi^2}\right).}
\]

In particular, an upper bound `U=o(n)` implies

\[
k\ge(1-o(1))\frac{\sqrt n}{2\pi}.
\]

Such sublinear cancellation is already known in a fixed quartic window. As a direct primary source, **Bourgain and Konyagin, Theorem 2.1, printed page 78**, state that, for any fixed `delta>0`, subgroups of size at least `p^delta` have `M<=n p^(-gamma)` for some positive `gamma` depending only on `delta`. [Published 2003 paper](https://www.numdam.org/item/10.1016/S1631-073X%2803%2900281-4.pdf). The source statement was freshly checked in this review. In `p<=n^4`, take `delta=1/4`; in any fixed quartic window, any fixed `delta<1/4` applies for sufficiently large parameters. Hence `M=o(n)`.

Consequently, **norm-prime examples, bounded `k`, and more generally `k=o(sqrt(n))` cannot persist at unbounded dyadic orders in that window**. Their exclusion follows from known cancellation, not from the desired square-root estimate. Searching for an unbounded principal quartic-window family is therefore not a viable route to refuting the unspecified-constant target. No effective finite threshold is extracted from the literature constants here, so the order-128 example is consistent with this asymptotic exclusion.

The supported conclusions are the norm-cofactor implication, its necessary cofactor lower bound, and the finite `C=2` counterexample. Neither these nor the flat Gram identities control deviations of every coset at the required scale. The full Paley conjecture and prize goal remain unproved.

No Lean build, external contact, hosted proof submission, or central-file edit was performed in this review.
