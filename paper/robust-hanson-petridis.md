# A robust Hanson–Petridis inequality and constant cancellation in Legendre-symbol sums above $p/2$

*Write-up from the Paley graph research workspace `paleygraph`, 2026-09-29. The results were
obtained and checked in the workspace's "Stepanov wave" (2026-09-26/27) by AI research agents.
No human mathematician has reviewed them.*

**Status.** Theorem A, Corollary B, Theorem C, Theorem D and Corollary E below are proved in
full in this document (Sections 2–8). Their proofs were re-derived independently by a referee
agent [Ref], which found no mathematical error and requested five presentation fixes; all five
are incorporated here. Theorem A, Corollary B and Theorem D (together with the non-vanishing
part of Lemma 4.2) have also been formalized in Lean 4 and compile with only the standard
axioms in two Mathlib versions (Section 10.2). Section 9 collects further results; they were
refereed separately [Ref3], which found no mathematical error and requested the scope and
hypothesis fixes incorporated there. Targeted literature searches
found no earlier statement of these results, but novelty is **not** established
(Section 11). None of this proves the Paley graph conjecture (Section 1.4).

## Abstract

Let $p$ be an odd prime, $\chi$ the Legendre symbol and
$S(A,B)=\sum_{a\in A,\,b\in B}\chi(a+b)$ for $A,B\subseteq\mathbb F_p$. The classical
second-moment bound $|S(A,B)|\le\sqrt{p|A||B|}$ is non-trivial only when $|A||B|>p$. Hanson
and Petridis used Stepanov's method to show that if every sum $a+b$ is a square or zero, then
$|A||B|\le (p-1)/2+|B\cap(-A)|$. We prove a robust form of their inequality. Let $e_b$ be the
number of $a\in A$ with $\chi(a+b)=-1$. The $(e+1)\times(e+1)$ Hankel determinant of the
normalised derivatives of the Hanson–Petridis polynomial has degree exactly
$(e+1)(\frac{p-1}{2}-e)$ whenever $|A|\le(p+1)/2$ and $0\le e\le(|A|-1)/2$. At every $b$ with $e_b\le e$ it vanishes to order at least
$(e+1-e_b)(|A|-\frac{3e+e_b}{2}-[b\in -A])$. This yields three consequences. The first is a
supersaturation inequality: each element of $B$ beyond the Hanson–Petridis bound forces on
average about $e+1$ non-residue sums. The second is a robust Hanson–Petridis inequality: if every $b\in B$
has at most $\eta|A|$ non-residue partners, where $\eta\le 1/8$, then
$|A||B|-|B\cap(-A)|\le(1-\sqrt{2\eta})^{-2}\frac{p-1}{2}+|A|$. The third is a bound for
$|A||B|\ge(\frac12+\kappa)p$, $0<\kappa\le 3/2$:
$$|S(A,B)|\le\Big[1-\big(1-(1+2\kappa)^{-1/2}\big)^2+\frac{\sqrt{(1/2+\kappa)p}+1}{2(p-1)}\Big]|A||B|.$$
The threshold $1/2$ cannot be lowered in a statement covering all set sizes; pairs with
$|A|=2$ are witnesses. Whether it can be lowered when both sets grow is open. For
$p\equiv1\pmod 4$ we obtain two-sided edge-density bounds for induced subgraphs of the Paley
graph on at least $\sqrt{(1/2+\kappa)p}$ vertices. The saving is a constant, not a power of
$p$, and nothing is obtained for $|A||B|\le p/2$. These results do not prove the Paley graph
conjecture.

---

## 1. Introduction

### 1.1 The two-set Paley graph conjecture

Throughout, $p$ is an odd prime and $\chi$ is the Legendre symbol modulo $p$, with
$\chi(0)=0$. For $A,B\subseteq\mathbb F_p$ put
$$S(A,B)=\sum_{a\in A}\sum_{b\in B}\chi(a+b).$$

**Conjecture** (Paley graph conjecture, two-set form; Satake [Sa, Conjecture 7]). For every
$0<\alpha\le 1$ there are $\beta=\beta(\alpha)>0$ and $p(\alpha)$ such that for every prime
$p>p(\alpha)$ and all $S,T\subseteq\mathbb F_p$ with $|S|,|T|>p^{\alpha}$,
$$\Big|\sum_{s\in S,\,t\in T}\chi(s-t)\Big|\le p^{-\beta}|S||T|.$$

Replacing $T$ by $-T$ gives the equivalent form with $s+t$. For $p\equiv1\pmod4$ the
conjecture implies that the clique number of the Paley graph $G_p$ is at most $p^{\varepsilon}$
for every $\varepsilon>0$ once $p$ is large [Sa, Remark 9]. The special case $S=T=\{1,\dots,N\}$ with $N=p^{\varepsilon}$
already implies that the least quadratic non-residue is at most $2N$. Indeed, otherwise every
$s+t\in[2,2N]$ is a residue and the sum equals $N^2$. This is Vinogradov's least-non-residue
conjecture, which is open.

### 1.2 What is known near $|A||B|\asymp p$

*Second moment.* For $B\subseteq\mathbb F_p$,
$\sum_{x\in\mathbb F_p}\big(\sum_{b\in B}\chi(x+b)\big)^2=|B|(p-|B|)$ (Lemma 2.4), and
Cauchy–Schwarz gives
$$|S(A,B)|\le\sqrt{|A||B|(p-|B|)}\le\sqrt{p|A||B|}.$$
Hanson and Petridis quote the second form as Vinogradov's estimate [HP, §1] (from [V, Ch. V,
Ex. 8]). The workspace brief calls it the Chung–Vinogradov bound; we did not consult Chung's
paper. Hanson and Petridis note that it "is still the best known estimate" for general sets
[HP, §1]. The second form is non-trivial only when $|A||B|>p$. The first form is non-trivial
only when $|A||B|>p-|B|$; for $|A|=|B|$ this means $|A||B|>p-\sqrt p+O(1)$.

*Karatsuba.* If $|A|>p^{1/2+\varepsilon}$ and $|B|>p^{\varepsilon}$ then
$|S(A,B)|\le|A||B|p^{-\delta}$ with $\delta=\delta(\varepsilon)>0$ ([Ka]; statement as given
in [FSX, (1.5)–(1.6)], which we used in place of the original). This range needs one of the
sets to be larger than $\sqrt p$ by a power of $p$.

*Hanson–Petridis.* [HP, Theorem 1.2] states the following. Let $p$ be prime, let $d$ be a
proper divisor of $p-1$, and let $Z_d$ be the group of $d$-th roots of unity in $\mathbb F_p$.
If $A,B\subseteq\mathbb F_p$ and $A+B\subseteq Z_d\cup\{0\}$, then $|A||B|\le d+|B\cap(-A)|$.
For $d=(p-1)/2$, $Z_d$ is the set $Q$ of non-zero squares. The corollary $|A|(|A|-1)\le d$ for $A-A\subseteq Z_d\cup\{0\}$
gives the Paley clique bound $\omega(G_p)\le(\sqrt{2p-1}+1)/2$ [HP, Cor. 1.5]. The authors
state that the method works only in prime fields [HP, after Cor. 1.5]. For general sets the
inequality yields very little: by Corollary B with $e=0$ below, it forces only about one
non-residue sum per element of $B$ beyond $(d+r)/|A|$, a saving of order $1/|A|$ in $S$.

So for pairs with $|A|\asymp|B|\asymp\sqrt p$ and $p/2<|A||B|<p-|B|$, none of the bounds
recalled above gives a constant saving.

### 1.3 Results

Notation (fixed for the whole paper): $d=(p-1)/2$; $m=|A|$, $n=|B|$;
$$e_b=\#\{a\in A:\chi(a+b)=-1\},\qquad \delta_b=[\,b\in -A\,],\qquad r=|B\cap(-A)|,\qquad
N_-=\sum_{b\in B}e_b .$$
Always $S(A,B)=mn-r-2N_-$ (Lemma 2.1).

**Theorem A (Hankel-minor Stepanov inequality).** Let $A\subseteq\mathbb F_p$ with
$1\le m\le (p+1)/2$ and let $e$ be an integer with $0\le e\le (m-1)/2$. Then
$$\sum_{\substack{b\in\mathbb F_p\\ e_b\le e}}(e+1-e_b)\Big(m-\frac{3e+e_b}{2}-\delta_b\Big)\;\le\;(e+1)(d-e).
\tag{$\star$}$$
All terms on the left are non-negative. More precisely, a Hankel determinant $H_{e+1}$ built
from the Hanson–Petridis polynomial (Section 4) has degree exactly $(e+1)(d-e)$. At every $b$
with $e_b\le e$ it is divisible by $(x-b)^{T_b}$, where
$T_b=(e+1-e_b)(m-\frac{3e+e_b}{2}-\delta_b)$.

For $e=0$, $(\star)$ reads $m\cdot\#\{b:e_b=0\}-\#\{b\in-A:e_b=0\}\le d$, which is
[HP, Theorem 1.2] for $d=(p-1)/2$. Inequality $(\star)$ is attained with equality at $e=1$
(Remark 4.5).

**Corollary B (supersaturation).** Let $1\le m\le(p+1)/2$, let $B\subseteq\mathbb F_p$ be
arbitrary, and let $0\le e\le(m-1)/2$ be an integer. Then
$$(m-2e)\big((e+1)n-N_-\big)\le(e+1)(d-e+r),$$
equivalently $N_-\ge(e+1)\big[n-\frac{d-e+r}{m-2e}\big]$, and
$$S(A,B)\le mn-r-2(e+1)\Big[n-\frac{d-e+r}{m-2e}\Big].$$

**Theorem C (robust Hanson–Petridis).** Let $A,B\subseteq\mathbb F_p$ and
$0\le\eta\le 1/8$, and suppose that every $b\in B$ has $e_b\le\eta|A|$. Then
$$|A||B|-|B\cap(-A)|\le(1-\sqrt{2\eta})^{-2}\,\frac{p-1}{2}+|A|.$$
The proof rests on the following inequality (Theorem 6.1), which needs no condition on
$\eta$. If $m\le(p+1)/2$ and $e_b\le e'$ for all $b\in B$, then for every integer $e$ with
$e'\le e\le(m-1)/2$,
$$(e+1-e')\Big[\Big(m-\frac{3e+e'}{2}\Big)n-r\Big]\le(e+1)(d-e),\qquad
mn-r\le\frac{(e+1)m\,d}{(e+1-e')(m-2e)}+\frac{2e\,r}{m-2e}.$$

**Theorem D (constant bias above $p/2$).** Let $p\ge 11$, $0<\kappa\le 3/2$,
$u=(1+2\kappa)^{-1/2}$, and let $A,B\subseteq\mathbb F_p$ satisfy $|A||B|\ge(\frac12+\kappa)p$.
Then
$$|S(A,B)|\le\Big[1-(1-u)^2+\frac{\sqrt{(1/2+\kappa)p}+1}{2(p-1)}\Big]\,|A||B|.$$
For $\kappa>3/2$ the case $\kappa=3/2$ applies, a saving of $1/4-O(p^{-1/2})$. For small
$\kappa$, $(1-u)^2=\kappa^2-3\kappa^3+O(\kappa^4)$.

*Comparison with the second-moment bound.* The elementary bound $|S|\le\sqrt{p|A||B|}$ saves
$1-(\frac12+\kappa)^{-1/2}$ when $|A||B|=(\frac12+\kappa)p$. This exceeds the saving
$(1-u)^2$ of Theorem D exactly when $u<2-\sqrt2$, i.e. $\kappa>\frac14+\frac{\sqrt2}2\approx0.957$.
So the new content of Theorem D is the range $p/2<|A||B|<(\frac34+\frac{\sqrt2}2)p\approx1.457p$;
the sharper constant of Section 9.1 extends it to $|A||B|<(1+\frac{\sqrt3}2)p\approx1.866p$
(the crossing $\kappa=(1+\sqrt3)/2$ was found by the second referee [Ref3]).

The error term makes Theorem D vacuous for small $p$. For $p\in\{11,13\}$ the bracket is at
least $1$ for every $\kappa\in(0,3/2]$; at $\kappa=3/2$ it is below $1$ exactly from $p=17$, at
$\kappa=1/2$ from $p=47$, and at $\kappa=1/10$ from $p=2741$ (Remark 7.4).

*Sharpness (Proposition 7.5).* For $p\equiv1\pmod4$, $A=\{0,1\}$ and
$B=\{b:\chi(b)\ne-1,\ \chi(b+1)\ne-1\}$ one has $|A||B|=(p+3)/2$ and $S(A,B)=(p-1)/2$. So the
threshold $1/2$ cannot be lowered in a statement that covers all set sizes. Whether a constant
saving holds for $|A||B|\ge\tau p$ with some $\tau<1/2$ when both $|A|,|B|\to\infty$ is
**open**. For $\kappa\le1$ the same family, padded, shows that no saving larger than
$2\kappa/(1+2\kappa)+O(1/p)$ is possible (Remark 7.6).

**Corollary E (Paley graph).** Let $p\equiv1\pmod4$ with $p\ge13$, let $0<\kappa\le3/2$, and
let $A\subseteq\mathbb F_p$ with $m=|A|\ge\sqrt{(1/2+\kappa)p}$. Let $\rho(A)$ be the edge
density of the subgraph of $G_p$ induced on $A$. Put $c(\kappa)=(1-(1+2\kappa)^{-1/2})^2$,
$\varepsilon_p=\frac{\sqrt{(1/2+\kappa)p}+1}{2(p-1)}$ and $\Phi=1-c(\kappa)+\varepsilon_p$.
Then
$$\frac{c(\kappa)}2-\frac{\varepsilon_p}2-\frac{\Phi}{2(m-1)}\;\le\;\rho(A)\;\le\;
1-\frac{c(\kappa)}2+\frac{\varepsilon_p}2+\frac{\Phi}{2(m-1)} ,$$
and both error terms are $O(p^{-1/2})$.

| here | [R] `stepanov-robust` | [PS] pass summary | referee [Ref] |
|---|---|---|---|
| Lemma 3.1 | Lemma 1.1 | — | §1 |
| Theorem A | Theorem 2.1 | Theorem A | §2 |
| Lemma 4.2 | Lemma 2.2 | — | §3 |
| Corollary B | Corollary 2.3 | Theorem B | §4, fix 1 |
| Theorem C | Theorem 2.4 | Theorem C | §5, fix 2 (§8.2) |
| Theorem D | Theorem 2.5 | Theorem D | §6, fixes 4–5 |
| Corollary E | — | Corollary E | §7 (explicit error terms) |

### 1.4 What these results do not do

They do not prove the Paley graph conjecture, and they do not come close:

1. The saving in Theorem D is a constant, about $\kappa^2$, and not a power $p^{-\delta}$.
2. Theorem D needs $|A||B|>p/2$. For balanced sets this means $|A|,|B|\gtrsim\sqrt{p/2}$,
   while the conjecture concerns $|A|,|B|>p^{\alpha}$ for every $\alpha>0$. The results live
   at the square-root scale only.
3. Below $p/2$ nothing is obtained, and for statements over all set sizes nothing can be
   (Proposition 7.5).
4. The method has a built-in cap. The auxiliary polynomial $H_{e+1}$ has degree
   $(e+1)(d-e)$, while its order at a point of $B$ is at most about $(e+1)|A|$. The count
   $\sum_bT_b\le(e+1)(d-e)$ therefore constrains only configurations with $|A||B|\gtrsim d$.
   Like Hanson–Petridis itself, the argument cannot reach below $|A||B|\asymp p$.

### 1.5 Idea of the proof

Hanson and Petridis choose weights $c_k$ so that
$F(x)=-1+\sum_k c_k(x+a_k)^{D}$, $D=d+m-1$, has degree $d$ and vanishes to order
$m-\delta_b$ at every $b$ for which all $a_k+b$ lie in $Q\cup\{0\}$. A single non-residue $a_k+b$ destroys this: then
$F(b)=-2c_k(a_k+b)^{m-1}\neq0$. The new observation is local. Near a point $b$, $F$ agrees to
order $m-\delta_b$ with $2\sum_{k\in E(b)}c_k(x+a_k)^D$, where $E(b)$ indexes the $e_b$
non-residue partners of $b$. So the normalised derivatives $u_s=F^{(s)}/(D)_s$ split as
$2\rho_s+\varepsilon_s$. The Hankel matrix $[\rho_{i+j}]$ of the first part is a sum of $e_b$
rank-one matrices, identically in $x$, just as the syndromes of a Reed–Solomon code with error
locators $(x+a_k)^{-1}$ are. The second part vanishes at $b$ to order at least $m-\delta_b-s$.
Expanding $\det[u_{i+j}]_{0\le i,j\le e}$ row by row, every surviving term keeps at least
$e+1-e_b$ rows of $\varepsilon$'s. So the determinant vanishes at $b$ to order about
$(e+1-e_b)m$, while its degree is exactly $(e+1)(d-e)$. The exact degree reduces to a binomial
determinant that is not divisible by $p$ because $D\le p-1$; this is where the prime field
enters. The weights $e+1-e_b$ are what the bias application needs, since
$\sum_{b\in B}(e+1-e_b)=(e+1)n-N_-$ is linear in the number of non-residue sums.

### 1.6 Sources

[R] = `research/stepanov-robust-2026-09-26.md` (original proofs),
[Ref] = `research/stepanov-referee2-2026-09-27.md` (independent referee),
[Sh] = `research/stepanov-sharpen-2026-09-27.md` (Section 9 only),
[L] = `research/stepanov-leanrobust-2026-09-27.md` with
`experiments/stepanov_leanrobust_lean/StepanovRobust.lean` (formalization),
[PS] = `research/stepanov-pass-summary-2026-09-27.md`. All paths are relative to the workspace
root. The verifier for this document is `experiments/paper_robust_hp_2026_09_29.py`, with
output in `results/paper_robust_hp_2026_09_29.json`. It recomputes every displayed identity,
inequality and constant at small primes; Section 10 gives details.

---

## 2. Notation, Lagrange weights and elementary facts

$\mathbb F_p$ is the field with $p$ elements, $p$ odd, $d=(p-1)/2$, and $Q$ is the set of
non-zero squares. By Euler's criterion, $y^d=\chi(y)$ for $y\ne0$. We write
$(x)_j=x(x-1)\cdots(x-j+1)$ for the falling factorial, with $(x)_0=1$, and $[P]$ for the
indicator of a statement $P$.

Let $A=\{a_1,\dots,a_m\}\subseteq\mathbb F_p$ with $m\ge1$. Put $D=d+m-1$. Note that
$$m\le\tfrac{p+1}{2}\iff m-1\le d\iff D\le p-1;$$
this is the only role of the size hypothesis on $A$. The **Lagrange weights** are
$c_k=\prod_{l\ne k}(a_k-a_l)^{-1}$ (so $c_1=1$ if $m=1$), and the **Hanson–Petridis
polynomial** is
$$F(x)=-1+\sum_{k=1}^m c_k(x+a_k)^D .$$
For $b\in\mathbb F_p$ put $y_k=b+a_k$ and $E(b)=\{k:\chi(y_k)=-1\}$, so that $e_b=|E(b)|$.
At most one $k$ has $y_k=0$, and this happens exactly when $\delta_b=1$. For $B\subseteq\mathbb F_p$,
the quantities $n,r,N_-$ are as in Section 1.3.

**Lemma 2.1 (counting).** $S(A,B)=mn-r-2N_-$.

*Proof.* Each $b\in B\cap(-A)$ has exactly one $a\in A$ with $a+b=0$. So the $mn$ pairs split
into $r$ pairs with $a+b=0$, $N_-$ pairs with $\chi(a+b)=-1$ and $mn-r-N_-$ pairs with
$\chi(a+b)=1$. $\square$

**Lemma 2.2 (Lagrange weights).** For every $t\in\mathbb F_p$ and $0\le j\le m-1$,
$$\sum_{k=1}^m c_k(t+a_k)^j=[\,j=m-1\,].$$

*Proof.* For $f\in\mathbb F_p[x]$ with $\deg f\le m-1$, Lagrange interpolation gives
$f(x)=\sum_k f(a_k)\prod_{l\ne k}\frac{x-a_l}{a_k-a_l}$. Comparing coefficients of $x^{m-1}$
gives $[x^{m-1}]f=\sum_k c_kf(a_k)$. Apply this to $f(x)=(x+t)^j$, whose coefficient of
$x^{m-1}$ is $\binom{j}{m-1}t^{j-m+1}=[j=m-1]$ for $j\le m-1$. $\square$

**Lemma 2.3 (order criterion).** Let $P\in\mathbb F_p[x]$, $b\in\mathbb F_p$ and
$1\le M\le p$. Then $(x-b)^M\mid P$ if and only if $P^{(j)}(b)=0$ for $0\le j<M$. Also, for
$0\le s\le D\le p-1$ the integer $(D)_s$ is a product of integers in $[1,p-1]$, hence a unit
in $\mathbb F_p$.

*Proof.* Write $P=\sum_j\tau_j(x-b)^j$. Then $P^{(j)}(b)=j!\,\tau_j$, and $j!$ is a unit for
$j<p$. $\square$

**Lemma 2.4 (second moment).** For $c\ne0$, $\sum_{x\in\mathbb F_p}\chi(x)\chi(x+c)=-1$.
Hence for $B\subseteq\mathbb F_p$,
$\sum_{x\in\mathbb F_p}\big(\sum_{b\in B}\chi(x+b)\big)^2=|B|(p-|B|)$, and
$|S(A,B)|\le\sqrt{|A||B|(p-|B|)}$.

*Proof.* For $x\ne0$, $\chi(x)\chi(x+c)=\chi(1+c/x)$, and $1+c/x$ runs over
$\mathbb F_p\setminus\{1\}$ as $x$ runs over $\mathbb F_p^\times$. So the sum is
$\sum_{y\ne1}\chi(y)=-1$. Expanding the square, the $n$ diagonal terms give $n(p-1)$ and the
$n(n-1)$ off-diagonal terms give $-n(n-1)$. Finally
$S(A,B)=\sum_{a\in A}s(a)$ with $s(x)=\sum_{b\in B}\chi(x+b)$, and Cauchy–Schwarz gives
$S^2\le|A|\sum_{x}s(x)^2$. $\square$

---

## 3. The local form of the Hanson–Petridis polynomial

For $E\subseteq\{1,\dots,m\}$ put $R_E(x)=\sum_{k\in E}c_k(x+a_k)^D$.

**Lemma 3.1.** Let $1\le m\le(p+1)/2$.

(a) $\deg F=d$, with leading coefficient $\binom{D}{m-1}\not\equiv0\pmod p$.

(b) If $b\notin-A$ and $0\le j\le m-1$, then
$F^{(j)}(b)=-2(D)_j\sum_{k\in E(b)}c_k\,y_k^{m-1-j}$.

(c) If $b=-a_{k_0}$, the formula of (b) holds for $0\le j\le m-2$, and for $m\ge2$,
$F^{(m-1)}(b)=(D)_{m-1}\big(-c_{k_0}-2\sum_{k\in E(b)}c_k\big)$.

(d) $(x-b)^{m-\delta_b}$ divides $F-2R_{E(b)}$.

*Proof.* (a) The coefficient of $x^{D-l}$ in $F+1$ is $\binom Dl\sum_kc_ka_k^l$. By
Lemma 2.2 with $t=0$, this is $0$ for $l\le m-2$ and $\binom{D}{m-1}$ for $l=m-1$. The degrees
$D-l$ with $l\le m-2$ are exactly those above $d$. Since $d\ge1$, the constant $-1$ does not
affect the coefficient of $x^d$. Finally $\binom{D}{m-1}=D!/((m-1)!\,d!)$ with $D\le p-1$, so
$p\nmid\binom D{m-1}$.

(b) For $j\ge1$, $F^{(j)}=(D)_j\sum_kc_k(x+a_k)^{D-j}$. For $y\ne0$ and $j\le m-1$ we have
$y^{D-j}=y^d\,y^{m-1-j}=\chi(y)\,y^{m-1-j}$. Write $\chi(y_k)=1-2[k\in E(b)]$. Then
$$F^{(j)}(b)=-[j=0]+(D)_j\Big(\sum_kc_ky_k^{m-1-j}-2\sum_{k\in E(b)}c_ky_k^{m-1-j}\Big),$$
where the term $-[j=0]$ is the constant of $F$. By Lemma 2.2 with $t=b$, the first inner sum
is $[j=0]$, and $(D)_0=1$, so the two cancel.

(c) Now $y_{k_0}=0$ and $k_0\notin E(b)$. Since $j\le m-1<D$, the term $k_0$ of
$\sum_kc_ky_k^{D-j}$ vanishes, and the computation in (b) applies to the other terms. The
unweighted sum becomes
$\sum_{k\ne k_0}c_ky_k^{m-1-j}=[j=0]-c_{k_0}\cdot0^{\,m-1-j}$. This equals $[j=0]$ for
$j\le m-2$ and $-c_{k_0}$ for $j=m-1$, which gives (c).

(d) For $0\le j\le m-1$,
$R_{E}^{(j)}(b)=(D)_j\sum_{k\in E}c_k\chi(y_k)y_k^{m-1-j}=-(D)_j\sum_{k\in E}c_ky_k^{m-1-j}$
with $E=E(b)$. Comparing with (b) and (c) gives $(F-2R_{E(b)})^{(j)}(b)=0$ for
$0\le j\le m-1-\delta_b$. Now apply Lemma 2.3 with $M=m-\delta_b\le p$; if $M=0$ there is
nothing to prove. $\square$

**Corollary 3.2 (Hanson–Petridis for $d=(p-1)/2$ [HP, Thm 1.2]).** If
$A+B\subseteq Q\cup\{0\}$ and $B\neq\emptyset$, then $|A||B|\le d+|B\cap(-A)|$.

*Proof.* $|A|\le|A+B|\le d+1$, so $m\le(p+1)/2$. Every $b\in B$ has $E(b)=\emptyset$, so by
Lemma 3.1(d) $(x-b)^{m-\delta_b}\mid F$. Since $F\ne0$ has degree $d$, summing over $b\in B$
gives $mn-r\le d$. $\square$

This is exactly the Hanson–Petridis argument. By Lemma 3.1(b), a point $b\notin-A$ with a
single non-residue partner $a_k$ has $F(b)=-2c_ky_k^{m-1}\ne0$. So $F$ itself says nothing
at such points. The next section repairs this.

---

## 4. The Hankel-minor inequality (proof of Theorem A)

### 4.1 Set-up

For $0\le s\le D$ let
$$u_s(x)=-[s=0]+\sum_{k=1}^mc_k(x+a_k)^{D-s}.$$
Thus $u_0=F$, and $u_s=F^{(s)}/(D)_s$ for $s\ge1$ (a unit by Lemma 2.3). For an integer
$0\le e\le(m-1)/2$ let
$$H_{e+1}(x)=\det\big[u_{i+j}(x)\big]_{0\le i,j\le e}\in\mathbb F_p[x].$$
$H_{e+1}$ depends on $A$ and $e$ only; the point $b$ enters only through the decomposition in
Step 1. Put
$$T_b=\sum_{i=e_b}^{e}(m-e-i-\delta_b)=(e+1-e_b)\Big(m-\frac{3e+e_b}2-\delta_b\Big)\quad(e_b\le e),$$
and $T_b=0$ if $e_b>e$. Every summand is at least $m-1-2e\ge0$.

**Theorem 4.1 (= Theorem A).** Let $1\le m\le(p+1)/2$ and $0\le e\le(m-1)/2$. Then
$\deg H_{e+1}=(e+1)(d-e)$ exactly, and $(x-b)^{T_b}\mid H_{e+1}$ for every
$b\in\mathbb F_p$. Consequently $\sum_{b}T_b\le(e+1)(d-e)$, which is $(\star)$.

### 4.2 Proof

**Step 1 (local splitting).** Fix $b$ with $e_b\le e$ and write $E=E(b)$, $\delta=\delta_b$. By
Lemma 3.1(d), $F=2R_E+(x-b)^{m-\delta}\Psi_b$ for some $\Psi_b\in\mathbb F_p[x]$. For
$0\le s\le2e$ put
$$\rho_s(x)=\sum_{k\in E}c_k(x+a_k)^{D-s},\qquad\varepsilon_s=u_s-2\rho_s .$$
Since $R_E^{(s)}=(D)_s\rho_s$ and $F^{(s)}=(D)_su_s$, we get
$(D)_s\,\varepsilon_s=\big((x-b)^{m-\delta}\Psi_b\big)^{(s)}$. By the Leibniz rule this is
$\sum_{t\le s}\binom st(m-\delta)_t(x-b)^{m-\delta-t}\Psi_b^{(s-t)}$. Here
$s\le2e\le m-1\le m-\delta$, so every term is divisible by $(x-b)^{m-\delta-s}$. As $(D)_s$ is
a unit,
$$(x-b)^{m-\delta-s}\mid\varepsilon_s\qquad(0\le s\le 2e).$$

**Step 2 (rank).** Since $D-2e\ge d\ge0$, for $0\le i,j\le e$
$$\rho_{i+j}=\sum_{k\in E}\big[c_k(x+a_k)^{D-2e}\big](x+a_k)^{e-i}(x+a_k)^{e-j}.$$
Thus $[\rho_{i+j}]_{i,j}=V^{\mathsf T}\Gamma V$, where $V$ is the $e_b\times(e+1)$ matrix of
polynomials $V_{k,j}=(x+a_k)^{e-j}$ and
$\Gamma=\operatorname{diag}\big(c_k(x+a_k)^{D-2e}\big)_{k\in E}$. Over the field
$\mathbb F_p(x)$ this matrix has rank at most $e_b$; it is zero if $e_b=0$. No division by
$x+a_k$ is needed.

**Step 3 (order at $b$).** Row $i$ of $U=[u_{i+j}]$ is $2\cdot(\text{row }i\text{ of }[\rho_{i+j}])$
plus row $i$ of $[\varepsilon_{i+j}]$. The determinant is multilinear in the rows, so
$$H_{e+1}=\sum_{S\subseteq\{0,\dots,e\}}\det M_S,$$
where row $i$ of $M_S$ is the $\varepsilon$-row if $i\in S$ and the $2\rho$-row if $i\notin S$.

* If $|S|<e+1-e_b$, then $M_S$ has more than $e_b$ rows taken from a matrix of rank
  $\le e_b$ over $\mathbb F_p(x)$. These rows are linearly dependent, so $\det M_S=0$ in
  $\mathbb F_p(x)$ and hence as a polynomial.
* If $|S|\ge e+1-e_b$, each entry $\varepsilon_{i+j}$ of a row $i\in S$ is divisible by
  $(x-b)^{m-\delta-i-j}$, hence by $(x-b)^{t_i}$ with $t_i:=m-\delta-i-e$, since $j\le e$.
  Here $t_i\ge m-1-2e\ge0$. Pulling these factors out of the rows in $S$ gives
  $(x-b)^{\sum_{i\in S}t_i}\mid\det M_S$. The sequence $t_i$ is non-negative and decreasing
  in $i$. If $\sigma_1>\sigma_2>\cdots$ are the elements of $S$, then $\sigma_j\le e-j+1$, so
  $t_{\sigma_j}\ge t_{e-j+1}$. Keeping only the first $e+1-e_b$ of them gives
  $\sum_{i\in S}t_i\ge\sum_{i=e_b}^{e}t_i=T_b$.

Hence $(x-b)^{T_b}$ divides every $\det M_S$ and therefore $H_{e+1}$. The closed form of
$T_b$ follows from
$\sum_{i=e_b}^e(m-e-\delta-i)=(e+1-e_b)(m-e-\delta)-(e+1-e_b)\frac{e+e_b}2$.

**Step 4 (degree).** Let $0\le s\le2e$. The coefficient of $x^{D-s-l}$ in
$\sum_kc_k(x+a_k)^{D-s}$ is $\binom{D-s}{l}\sum_kc_ka_k^l$. By Lemma 2.2 this is $0$ for
$l\le m-2$ and $\binom{D-s}{m-1}$ for $l=m-1$. Since $s\le2e\le m-1\le d$ and $d\ge1$, the
constant $-1$ in $u_0$ lies in degree $0<d$. So $\deg u_s\le d-s$, and the coefficient of
$x^{d-s}$ in $u_s$ is $\binom{D-s}{m-1}$ read mod $p$. In the Leibniz expansion
$H_{e+1}=\sum_\sigma\operatorname{sgn}\sigma\prod_iu_{i+\sigma(i)}$, each product has degree
at most $\sum_i(d-i-\sigma(i))=(e+1)d-e(e+1)=(e+1)(d-e)$. Its coefficient of $x^{(e+1)(d-e)}$
is $\prod_i\binom{D-i-\sigma(i)}{m-1}$. Therefore
$$[x^{(e+1)(d-e)}]\,H_{e+1}\equiv\Lambda:=\det\Big[\binom{D-i-j}{m-1}\Big]_{0\le i,j\le e}\pmod p .$$
By Lemma 4.2 below, $p\nmid\Lambda$. So $\deg H_{e+1}=(e+1)(d-e)$ and in particular
$H_{e+1}\ne0$.

**Step 5 (count).** The polynomials $(x-b)^{T_b}$, $b\in\mathbb F_p$, are pairwise coprime and
all divide $H_{e+1}\ne0$. So $\prod_b(x-b)^{T_b}\mid H_{e+1}$ and
$\sum_bT_b\le\deg H_{e+1}=(e+1)(d-e)$. $\blacksquare$

### 4.3 The leading coefficient

**Lemma 4.2.** Let $d\ge1$, $m\ge1$ and $e\ge0$ be integers with $2e\le d$, $2e\le m-1$ and
$D=d+m-1\le p-1$. Then $p\nmid\Lambda=\det\big[\binom{D-i-j}{m-1}\big]_{0\le i,j\le e}$.
Moreover, with $N'=D-2e$ and $c=m-1-e$,
$$\Lambda=(-1)^{e(e+1)/2}\prod_{k=1}^{e}k!\;\prod_{i=1}^{e+1}\frac{(N'+i-1)!}{(c-i+e+1)!\,(N'-c+i-1)!}.$$

We give two proofs of the non-vanishing. The first is from [R, Lemma 2.2] (checked in
[Ref, §3]) and also yields the exact value. The second is due to the Lean formalization
[L, §3.4]; it is self-contained and is the argument that was formalized.

*Proof 1 (Krattenthaler's determinant).* Write $K=m-1$. By Vandermonde's convolution with
$x=e-i\ge0$ and $y=D-e-j\ge0$,
$\binom{D-i-j}{K}=\sum_{l=0}^{e}\binom{e-i}{l}\binom{D-e-j}{K-l}$. So
$[\binom{D-i-j}K]=UW$ with $U_{il}=\binom{e-i}l$ and $W_{lj}=\binom{D-e-j}{K-l}$
($0\le i,l,j\le e$). $U$ vanishes where $i+l>e$ and equals $1$ on the anti-diagonal $i+l=e$.
Hence $\det U$ is the sign of $i\mapsto e-i$, namely $(-1)^{e(e+1)/2}$. Again
$\binom{D-e-j}{K-l}=\sum_t\binom{e-j}t\binom{D-2e}{K-l-t}$, so $W=ZU^{\mathsf T}$ with
$Z_{lt}=\binom{N'}{K-l-t}$. Therefore $\Lambda=\det U\cdot\det Z\cdot\det U=\det Z$.
Reversing the columns of $Z$ ($t\mapsto e-t$, sign $(-1)^{e(e+1)/2}$) turns it into the
Toeplitz matrix $[\binom{N'}{c-l+t}]_{0\le l,t\le e}$ with $c=K-e$. With $n=e+1$ and
indices shifted to $1\le i,j\le n$, this is $\det_{1\le i,j\le n}\binom{N'}{L_i+j}$ with
$L_i=c-i$. Krattenthaler's evaluation [Kr, Theorem 26, eq. (3.12)] at $q=1$ gives
$$\det_{1\le i,j\le n}\binom{A}{L_i+j}=\frac{\prod_{i<j}(L_i-L_j)\,\prod_i(A+i-1)!}{\prod_i(L_i+n)!\,\prod_i(A-L_i-1)!}.$$
Take $A=N'$; then $\prod_{i<j}(L_i-L_j)=\prod_{i<j}(j-i)=\prod_{k=1}^{n-1}k!$. All
factorial arguments are non-negative, because $c\ge e$ (from $2e\le m-1$) and
$N'-c=d-e\ge0$. This gives the displayed value. For non-vanishing, every factorial has
argument in $[0,D-e]$: the numerator has $k!$ with $k\le e$ and $(N'+i-1)!\le(D-e)!$; the
denominator has $(c-i+e+1)!\le(m-1)!$ and $(N'-c+i-1)!\le d!$. Since $D-e\le p-1$, the
identity $\Lambda\cdot(\text{denominator})=\pm(\text{numerator})$ has a product of units mod
$p$ on each side. As $\Lambda\in\mathbb Z$, it follows that $p\nmid\Lambda$. $\square$

*Proof 2 (interpolation; from the formalization [L, §3.4], theorem `hankelLead_ne_zero`).*
Write $K=m-1$, so $D=d+K$ and $e\le K$. For $0\le i,j\le e$ we have $d-i-j\ge0$, and the
integer identity
$$\binom{D-i-j}{K}\,K!\,(d-i)!=(D-i-e)!\;(d-i)_j\;(D-i-j)_{e-j}$$
holds. Indeed, $\binom{D-i-j}KK!\,(d-i-j)!=(D-i-j)!$ because $D-i-j-K=d-i-j$. Multiply by
$(d-i)!/(d-i-j)!=(d-i)_j$, and use $(D-i-j)!=(D-i-e)!\,(D-i-j)_{e-j}$. Hence in
$\mathbb F_p$
$$\binom{D-i-j}K=\alpha_i\,P_j(i),\qquad \alpha_i=\frac{(D-i-e)!}{K!\,(d-i)!},\qquad
P_j(x)=(d-x)_j\,(D-j-x)_{e-j}\in\mathbb F_p[x].$$
All factorials here are of integers at most $D\le p-1$, so $\alpha_i\ne0$. Also
$\deg P_j\le e$. Therefore $\Lambda\equiv\big(\prod_i\alpha_i\big)\det[P_j(i)]_{0\le i,j\le e}$.

Suppose $\det[P_j(i)]=0$ in $\mathbb F_p$. Then there is a non-zero vector
$(\lambda_0,\dots,\lambda_e)$ with $\sum_j\lambda_jP_j(i)=0$ for $i=0,\dots,e$. The
polynomial $G=\sum_j\lambda_jP_j$ has degree at most $e$ and vanishes at the $e+1$ distinct
points $0,1,\dots,e$ of $\mathbb F_p$ (as $e<p$), so $G=0$. Evaluate at $x=d-k$,
$0\le k\le e$:
$$P_j(d-k)=(k)_j\,(K+k-j)_{e-j}.$$
For $j>k$ this is $0$. For $j=k$ it equals $k!\,(K)_{e-k}$, a product of integers in
$[1,p-1]$ (using $e\le K\le D\le p-1$), hence non-zero. So $0=G(d)=\lambda_0P_0(d)$ gives
$\lambda_0=0$. Inductively, $0=G(d-k)=\lambda_kP_k(d-k)$ gives $\lambda_k=0$ for every $k$,
a contradiction. Hence $p\nmid\Lambda$. $\square$

Proof 2 uses only $e\le m-1$, $2e\le d$ and $D\le p-1$. The hypothesis $D\le p-1$ cannot be
dropped for general $d$: $d=2$, $m=3$, $e=1$ gives
$\Lambda=\binom42\binom22-\binom32^2=-3$, which is divisible by $p=3$, where $D=4\ge p$
([L, §2]; verifier §F).

### 4.4 Remarks

**Remark 4.3 (Step 1 without derivatives).** The formalization [L, §3.1] proves Step 1 from
Taylor coefficients. Since
$\varepsilon_s(x+b)=-[s=0]+\sum_kc_k(1-2[k\in E])(x+y_k)^{D-s}$, the coefficient of $x^j$ in
$\varepsilon_s(x+b)$ is
$$-[s=j=0]+\binom{D-s}{j}\sum_kc_k\,(1-2[k\in E])\,y_k^{D-s-j}.$$
Suppose $s+j+\delta_b<m$. If $y_k\ne0$, then $1-2[k\in E]=\chi(y_k)$ and
$\chi(y_k)y_k^{D-s-j}=\chi(y_k)^2y_k^{m-1-s-j}=y_k^{m-1-s-j}$. If $y_k=0$, then
$\delta_b=1$, so $m-1-s-j\ge1$, and both $y_k^{D-s-j}$ and $y_k^{m-1-s-j}$ vanish. So the
coefficient is $-[s=j=0]+\binom{D-s}{j}\sum_kc_ky_k^{m-1-s-j}=-[s=j=0]+\binom{D-s}{j}[s+j=0]=0$
by Lemma 2.2. Hence $(x-b)^{m-\delta_b-s}\mid\varepsilon_s$, with no unit needed. In that
route the hypothesis $D\le p-1$ is used only for $p\nmid\Lambda$ [L, §6].

**Remark 4.4 (where the prime field is used).** Over $\mathbb F_{q}$, $q=p^2$, with its
quadratic character, take $A=B=\mathbb F_p$; the size condition $|A|\le(q+1)/2$ holds. Every
element of $\mathbb F_p$ is a square in $\mathbb F_{q}$, so every $b\in\mathbb F_p$ has
$e_b=0$. But $mn-r=p^2-p>(q-1)/2$, so the analogue of $(\star)$ at $e=0$ is false there. This matches
[HP, after Cor. 1.5]: the argument needs the characteristic to exceed $D$.

**Remark 4.5 (equality at $e=1$; referee fix 3).** $(\star)$ can hold with equality at $e\ge1$.
For $p=13$, $A=\{0,2,3,5\}$, $e=1$ ($d=6$), the points with $e_b\le1$ are
$b=1,7,9,12$ ($e_b=1$, $\delta_b=0$, $T_b=2$) and $b=10,11$ ($e_b=1$, $\delta_b=1$, $T_b=1$).
The left side is $10=2\cdot(6-1)$ [Ref, §2]. In the exhaustive checks the largest ratio
of the two sides at $e\ge1$ is exactly $1$ for $p\in\{5,7,11,13\}$. It is $13/14$, $15/16$,
$19/20$ and $25/26$ at $p=17,19,23,29$ (`results/stepanov_star_exhaustive_2026_09_27.json`;
recomputed for $p\le19$ by the verifier, §C). The sentence in [R] that the maximum at
$e\ge1$ is about $0.92$ is wrong. The Step 3 bound on the order of $H_{e+1}$ at a single
point is also attained: at 976 points with $e\ge1$ among the 578 determinants computed in
§B of the verifier.

---

## 5. Supersaturation (proof of Corollary B)

The hypothesis $m\le(p+1)/2$ is needed for the proof ([Ref, fix 1]: the pass summary's
Theorem B omitted it). No counterexample with $m>(p+1)/2$ was found at $p\le23$, but that
range is unproved.

*Proof of Corollary B.* Let $b\in B$ and put $w_b=(m-2e)(e+1-e_b)-(e+1)\delta_b$. If
$e_b\le e$ then $(3e+e_b)/2\le2e$, so
$$T_b\ge(e+1-e_b)(m-2e)-(e+1-e_b)\delta_b\ge w_b .$$
If $e_b>e$ then $w_b<0=T_b$, because $m-2e\ge1$. Summing over $b\in B$ and using
Theorem A,
$$(m-2e)\big((e+1)n-N_-\big)-(e+1)r=\sum_{b\in B}w_b\le\sum_{b\in\mathbb F_p}T_b\le(e+1)(d-e).$$
This is the first claim; the second is a rearrangement, and the third follows from
Lemma 2.1. $\square$

The sign of the terms with $e_b>e$ is what allows the inequality to be extended from
$\{b:e_b\le e\}$ to all of $B$. In words: once $n>(d-e+r)/(m-2e)$, each further element of
$B$ forces on average $e+1$ further non-residue sums.

---

## 6. Robust Hanson–Petridis (proof of Theorem C)

**Theorem 6.1.** Let $1\le m\le(p+1)/2$ and let $B\subseteq\mathbb F_p$ satisfy $e_b\le e'$
for all $b\in B$. For every integer $e$ with $e'\le e\le(m-1)/2$,
$$(e+1-e')\Big[\Big(m-\frac{3e+e'}2\Big)n-r\Big]\le(e+1)(d-e),\qquad
mn-r\le\frac{(e+1)m\,d}{(e+1-e')(m-2e)}+\frac{2e\,r}{m-2e}.$$

*Proof.* For $b\in B$ we have $e+1-e_b\ge e+1-e'>0$ and
$m-\frac{3e+e_b}2-\delta_b\ge m-\frac{3e+e'}2-\delta_b\ge m-2e-1\ge0$. So
$T_b\ge(e+1-e')(m-\frac{3e+e'}2-\delta_b)$. Summing over $b\in B$ and applying $(\star)$
gives the first display. Since $m-\frac{3e+e'}2\ge m-2e$, it follows that
$(m-2e)n-r\le\frac{(e+1)(d-e)}{e+1-e'}\le\frac{(e+1)d}{e+1-e'}$. Now use the identity
$mn-r=\frac m{m-2e}\big((m-2e)n-r\big)+\frac{2e\,r}{m-2e}$. $\square$

*Proof of Theorem C when $m\le(p+1)/2$.* We may assume $B\ne\emptyset$. Let
$e'=\max_{b\in B}e_b\le\eta m$ and $\eta'=e'/m\le\eta$. Since
$f(\eta)=(1-\sqrt{2\eta})^{-2}$ is increasing, it suffices to prove the bound with $\eta'$.
If $e'=0$, Theorem 6.1 with $e=0$ gives $mn-r\le d$. Otherwise $m\ge1/\eta'\ge8$. Put
$x=\sqrt{\eta'/2}\le1/4$ and $e=\lceil xm\rceil-1$. Then $e+1\ge xm$ and $e<xm\le m/4$, so
$2e\le m-1$. Since $\eta'\le1/8$ is equivalent to $\eta'\le x/2$, we get
$m(x-\eta')\ge mx/2\ge\frac1{2\sqrt{2\eta'}}\ge1$. Hence $e\ge xm-1\ge\eta'm=e'$. Theorem 6.1
applies, and
$$\frac{e+1}{e+1-e'}\le\frac1{1-\eta'/x}=\frac1{1-\sqrt{2\eta'}},\qquad
\frac m{m-2e}<\frac1{1-2x}=\frac1{1-\sqrt{2\eta'}},\qquad
\frac{2e\,r}{m-2e}<m\cdot\frac{2x}{1-2x}\le m,$$
using $r\le m$ and $x\le1/4$. $\square$

**Proposition 6.2 (no size hypothesis needed; [Ref, §8.2], referee fix 2).** Theorem C holds
for $m>(p+1)/2$ as well. This argument is the referee's own. It has not been refereed
separately, but it is checked exhaustively for $p\le19$ here (verifier §C) and for $p\le23$
in [Ref].

*Proof.* Suppose $m>(p+1)/2$, that every $b\in B$ has $e_b\le\eta m\le m/8$, and that
$n\ge2$. Take $b_1\ne b_2$ in $B$ and put $c=b_2-b_1$. Each $b_i$ has at least $7m/8$
partners $a$ with $\chi(a+b_i)\ne-1$, so at least $3m/4$ elements $a\in A$ have both. The
map $a\mapsto y=a+b_1$ is injective into
$\{y:\chi(y)\ne-1,\ \chi(y+c)\ne-1\}$. By Lemma 2.4 (Jacobsthal),
$\sum_{y\ne0,-c}(1+\chi(y))(1+\chi(y+c))=p-3-\chi(c)-\chi(-c)$. So this set has
$$\frac{p-3-\chi(c)-\chi(-c)}4+[\chi(c)=1]+[\chi(-c)=1]\le\frac{p+3}4$$
elements. Hence $3m/4\le(p+3)/4$, i.e. $m\le(p+3)/3\le(p+1)/2$ for $p\ge3$, a contradiction.
So $n\le1$ and $mn-r\le m\le(1-\sqrt{2\eta})^{-2}d+m$. $\square$

**Remark 6.3 (comparison with the second moment).** If $e_b\le e'$ on $B$ and
$2e'\le m-1$, then every $b\in B$ has
$\sum_{a\in A}\chi(a+b)=m-\delta_b-2e_b\ge m-1-2e'\ge0$. Lemma 2.4 with the roles of $A$
and $B$ exchanged then gives $n(m-1-2e')^2\le m(p-m)$. With $e'=\eta m$ and $m$ large this
is $mn\lesssim 2d/(1-2\eta)^2$: a bounded constant, but one that does not tend to $1$ as
$\eta\to0$. Theorem C has $f(\eta)=(1-\sqrt{2\eta})^{-2}=1+2\sqrt{2\eta}+O(\eta)\to1$. It is
smaller than $2/(1-2\eta)^2$ exactly when $\eta<(3-2\sqrt2)/2\approx0.0858$ (with
$s=\sqrt{2\eta}$, equality means $1+s=\sqrt2$).

---

## 7. Constant bias above $p/2$ (proof of Theorem D)

**Lemma 7.1 (core estimate).** Let $1\le m\le(p+1)/2$, $0<\kappa\le3/2$,
$u=(1+2\kappa)^{-1/2}$, and let $B$ satisfy $mn\ge(\frac12+\kappa)p$. Then
$$S(A,B)\le\Big[1-(1-u)^2+u(1-u)\frac md\Big]mn .$$

*Proof* (as in [Ref, §6]). (i) $mn\ge(1+2\kappa)p/2>(1+2\kappa)d=d/u^2$.

(ii) $u\in[\frac12,1)$ and $\theta:=(1-u)/2\in(0,\frac14]$. Let $e=\lfloor\theta m\rfloor$.
Then $e+1>\theta m$ and $m-2e\ge m(1-2\theta)=um>0$. Also $e=0$ if $m=1$, and
$e\le m/4\le(m-1)/2$ if $m\ge2$. By Corollary B and $r\le m$,
$$(m-2e)\big((e+1)n-N_-\big)\le(e+1)(d-e+r)\le(e+1)(d+m).$$

(iii) Put $Y=n-\frac{d+m}{um}$. Then $N_-\ge(e+1)Y$. If $(e+1)n-N_-\le0$ this holds because
$Y\le n$. Otherwise $um((e+1)n-N_-)\le(m-2e)((e+1)n-N_-)\le(e+1)(d+m)$.

(iv) If $Y<0$, then $umn<d+m$. By (i), $umn>d/u$, so $d(1/u-1)<m$, i.e. $(1-u)<um/d$. Then
$(1-u)^2<u(1-u)m/d$, the bracket exceeds $1$, and the claim follows from $S\le mn$. If
$Y\ge0$, then $N_-\ge\theta mY$ and
$$S=mn-r-2N_-\le mn-(1-u)mY=mn\Big[1-(1-u)\Big(1-\frac{d+m}{umn}\Big)\Big].$$
By (i), $\frac{d+m}{umn}<\frac{(d+m)u^2}{ud}=u(1+\frac md)$. Since $1-u>0$, it follows that
$S\le mn[1-(1-u)(1-u-u\frac md)]$, which is the claim. $\square$

**Lemma 7.2 (sub-sampling and symmetries).** (a) For $1\le m_0\le m$,
$S(A,B)=\frac m{m_0}\cdot\operatorname{avg}_{A'\subseteq A,\,|A'|=m_0}S(A',B)$.
(b) $S(A,B)=S(B,A)$, and for a non-residue $\nu$, $S(\nu A,\nu B)=-S(A,B)$, with
$|\nu A|=|A|$ and $|\nu B|=|B|$.

*Proof.* (a) Each $a\in A$ lies in a fraction $\binom{m-1}{m_0-1}/\binom m{m_0}=m_0/m$ of the
$m_0$-subsets, so the average of $S(A',B)=\sum_{a\in A'}\sum_{b}\chi(a+b)$ is
$\frac{m_0}mS(A,B)$. (b) $\chi(\nu(a+b))=-\chi(a+b)$. $\square$

*Proof of Theorem D.* By Lemma 7.2(b) we may assume $m\le n$. Then
$n^2\ge mn\ge(\frac12+\kappa)p$. Let $m_0=\lceil(\frac12+\kappa)p/n\rceil$. Then $m_0\le m$
(as $m\ge(\frac12+\kappa)p/n$ and $m\in\mathbb Z$), $m_0n\ge(\frac12+\kappa)p$, and
$$m_0<\frac{(1/2+\kappa)p}{n}+1\le\sqrt{(1/2+\kappa)p}+1\le\sqrt{2p}+1\le\frac{p+1}2 .$$
The last inequality is $8p\le(p-1)^2$, i.e. $p^2-10p+1\ge0$, which holds for $p\ge11$ and
fails for $p=7$. Lemma 7.1 applies to every $(A',B)$ with $|A'|=m_0$:
$S(A',B)\le[1-(1-u)^2+u(1-u)m_0/d]\,m_0n$. Averaging with Lemma 7.2(a) gives
$S(A,B)\le[1-(1-u)^2+u(1-u)m_0/d]\,mn$. Since $u(1-u)\le\frac14$ and
$m_0/d\le 2(\sqrt{(1/2+\kappa)p}+1)/(p-1)$, we have
$u(1-u)m_0/d\le\frac{\sqrt{(1/2+\kappa)p}+1}{2(p-1)}$. Applying the same bound to
$(\nu A,\nu B)$ bounds $-S(A,B)$. For $\kappa>3/2$ the hypothesis with $\kappa=3/2$ holds.
The expansion of $(1-u)^2$ is a direct series computation (verifier §G). $\blacksquare$

**Remark 7.3.** Theorem D needs no size hypothesis on $A$ or $B$. The condition
$m\le(p+1)/2$ is met automatically by the sub-sampled sets.

**Remark 7.4 (Theorem D is vacuous for small $p$; referee fix 5).** The bracket is below $1$
iff $\varepsilon_p(\kappa):=\frac{\sqrt{(1/2+\kappa)p}+1}{2(p-1)}<(1-u)^2$. For fixed
$\kappa$, $\varepsilon_p$ is decreasing in $p$: the numerator of its derivative is
$-\sqrt{c}\,\frac{p+1}{2\sqrt p}-1<0$ with $c=\frac12+\kappa$. So for each $\kappa$ the bound
is non-trivial exactly from some prime on. That prime is $17$ for $\kappa=3/2$, $47$ for
$\kappa=1/2$, and $2741$ for $\kappa=1/10$ ([Ref, §6]; recomputed with 40-digit arithmetic,
verifier §D).

For $p\in\{11,13\}$ the bracket is at least $1$ for every $\kappa\in(0,3/2]$. Both
$\varepsilon_p(\kappa)$ and $(1-u)^2$ are increasing in $\kappa$. So it suffices that
$\varepsilon_p(\kappa_i)\ge(1-u(\kappa_{i+1}))^2$ on a grid
$0=\kappa_0<\dots<\kappa_{3000}=3/2$. This holds with margin $0.0041$ at $p=13$ and $0.034$
at $p=11$ (40-digit arithmetic). Theorem D is an asymptotic statement with explicit
constants.

**Proposition 7.5 (the threshold $1/2$ over all set sizes; referee fix 4).** Let
$p\equiv1\pmod4$, $A=\{0,1\}$ and $B=\{b\in\mathbb F_p:\chi(b)\ne-1,\ \chi(b+1)\ne-1\}$. Then
$|B|=(p+3)/4$, $r=2$ and $S(A,B)=(p-1)/2$. So $|A||B|=(p+3)/2$ and
$S(A,B)/(|A||B|)=(p-1)/(p+3)$. Consequently, for no $\tau<1/2$ and $c>0$ does
$|S(A,B)|\le(1-c)|A||B|$ hold for all large $p$ and all $A,B$ with $|A||B|\ge\tau p$.

*Proof.* $\chi(-1)=1$, so $b=0$ and $b=-1$ lie in $B$. For $b\notin\{0,-1\}$ we need
$\chi(b)=\chi(b+1)=1$. By Lemma 2.4 the number of such $b$ is
$\frac14\sum_{b\ne0,-1}(1+\chi(b))(1+\chi(b+1))=\frac14\big((p-2)-\chi(-1)-1-1\big)=\frac{p-5}4$.
So $|B|=(p+3)/4$ and $B\cap(-A)=\{0,-1\}$. Every $b\in B$ has $e_b=0$, so
$S=2|B|-2=(p-1)/2$ by Lemma 2.1. $\square$

For $p\equiv3\pmod4$ the same $B$ has $(p+1)/4$ elements, $r=1$ and $S=(p-1)/2$, so the bias
is $(p-1)/(p+1)$ at $|A||B|=(p+1)/2$ (verifier §E). The witness has $|A|=2$. Whether the
threshold can be lowered when both $|A|$ and $|B|$ tend to infinity is **open** [Ref, §6,
Remark (2)]. The Paley graph conjecture would imply that it can; this is a heuristic, not a
result. For $|A|=2$, Corollary B with $e=0$ gives $S(A,B)\le d$ for every $B$. So the family
of Proposition 7.5 is extremal among two-element sets.

**Remark 7.6 (an upper bound for the saving).** Let $p\equiv1\pmod4$ and $A=\{0,1\}$. The
number of $b$ with $e_b=1$ is $(p-1)/2$, and none of them lies in $-A$; each contributes
$\chi(b)+\chi(b+1)=0$ to $S$. For $0<\kappa\le1$, adjoin such points to the set $B$ of
Proposition 7.5 until $|B|=\lceil(\frac12+\kappa)p/2\rceil$. Then
$|A||B|\ge(\frac12+\kappa)p$ and $S=(p-1)/2$. So $S/(|A||B|)=1/(1+2\kappa)+O(1/p)$, and no
bound of the form of Theorem D can have a saving larger than $2\kappa/(1+2\kappa)+O(1/p)$.
Theorem D gives $(1-u)^2\approx\kappa^2$. The true optimal saving lies between these, and
closing the gap is open (see also Section 9.3).

---

## 8. Paley-graph densities (proof of Corollary E)

*Proof of Corollary E.* Let $B=-A$, so that $|A||B|=m^2\ge(\frac12+\kappa)p$ and $p\ge13$. By
Theorem D, $|S(A,-A)|\le\Phi\,m^2$. Since $\chi(-1)=1$ we have $\chi(a-a')=\chi(a'-a)$. The
diagonal contributes $\chi(0)=0$, so with $e(A)$ the number of edges inside $A$,
$$S(A,-A)=\sum_{a\ne a'}\chi(a-a')=2e(A)-2\Big(\binom m2-e(A)\Big)=4e(A)-m(m-1)=m(m-1)\big(2\rho(A)-1\big).$$
Hence $|2\rho(A)-1|\le\Phi\,m/(m-1)=\Phi\,(1+\frac1{m-1})$, which is the claim. Since
$m\ge\sqrt{p/2}$, both $\varepsilon_p$ and $1/(m-1)$ are $O(p^{-1/2})$. $\square$

*Comparison.* From the Hanson–Petridis clique bound one gets only this: every $m$-set
contains at least $m-(\sqrt{2p-1}+1)/2$ non-edges, since deleting one endpoint of each
non-edge leaves a clique. That is an edge density $1-O(1/m)$ at $m\asymp\sqrt p$. The
second-moment bound gives $|2\rho-1|\le\sqrt{p-m}/(m-1)$, which bounds $\rho$ away from $1$
only for $m\ge(1+\kappa')\sqrt p$. Corollary E gives $\rho\in[c(\kappa)/2-o(1),\,1-c(\kappa)/2+o(1)]$
already for $m\ge\sqrt{(1/2+\kappa)p}$. Numerically the bound is far from tight: at
$p=1009$, $\kappa=3/2$, $m=45$, a greedy search finds density $0.716$ against the proven
upper bound $0.895$ (verifier §H). At $p=2017$, $\kappa=3/2$, $m=64$, the referee's local
search finds $0.729$ against $0.889$ [Ref, §7]. For small $\kappa$ and moderate $p$ the bound
is vacuous (Remark 7.4).

---

## 9. Further results

The statements in this section are proved in the workspace notes cited. They take Theorem A
as input and, where stated, Weil's bound. A separate referee agent [Ref3] re-derived each of
them, checked them with independent code (5,555,096 checks, 0 failures) and found no
mathematical error; its scope and hypothesis fixes are incorporated below. The constants
displayed here are also re-checked in §G of our verifier. To avoid a clash with $u=(1+2\kappa)^{-1/2}$ above, we
write $U(w)$ for the function called $u(w)$ in [Sh].

**9.1 A sharper constant** ([Sh, Theorem 1.5]). Put $w=1/(1+2\kappa)$,
$U(w)=(\sqrt{12w-3w^2}-w)/2$, $m_1=\lceil\sqrt{(1/2+\kappa)p}\rceil$ and
$w_1=(p-1+2m_1)/((1+2\kappa)p)$. For $p\ge11$, $0<\kappa\le3/2$ and
$|A||B|\ge(\frac12+\kappa)p$,
$$|S(A,B)|\le U(\min(1,w_1))\,|A||B|\le\Big[1-\eta^\star(\kappa)+2.14\,\frac{\sqrt{2p}+1}{p}\Big]|A||B|,$$
$$\eta^\star(\kappa)=1-U(w)=\tfrac43\kappa^2-\tfrac{40}9\kappa^3+O(\kappa^4).$$
The input is a "mean form" of $(\star)$. With $e=E-1$ and the exact weight
$W_E(x)=(E-x)(m+\frac32-\frac{3E+x}2)$ for $x\le E$ and $W_E(x)=0$ for $x\ge E$, one has
$\sum_{b\in B}W_E(e_b)\le E(d-E+1+r)$. $W_E$ is convex, so Jensen's inequality turns this
into a bound on the mean of $e_b$. A one-variable inequality follows. For $m\le10$ and
$w\in[1/4,11/20]$ it is a computer-assisted rational-interval verification [Sh, Thm 1.3].
Values of $\eta^\star$ against $(1-u)^2$ from Theorem D: $9.838\cdot10^{-3}$ vs
$7.591\cdot10^{-3}$ at $\kappa=0.1$; $0.10436$ vs $0.08579$ at $\kappa=1/2$; $0.28647$ vs
$0.25$ at $\kappa=3/2$. The ratio tends to $4/3$ as $\kappa\to0$. For
$\kappa>(1+\sqrt3)/2\approx1.366$ (equivalently $w<2-\sqrt3$) the second-moment saving
$1-(\frac12+\kappa)^{-1/2}$ is larger than $\eta^\star$ (e.g. $0.2929$ against $0.2865$ at
$\kappa=3/2$) [Ref3].

**9.2 This is the limit of $(\star)$** ([Sh, Prop. 2.1]). Consider balanced pairs with
$m\to\infty$, $m/d\to0$. There, no non-negative combination of the inequalities $(\star)$
applied to the fixed set $A$ (or to $B$), over all $e$ and with the exact weights, certifies a
saving larger than $\eta^\star(\kappa)+o(1)$. Combinations that also apply $(\star)$ to subsets
$A'\subseteq A$ are not covered; that part is heuristic [Sh]. (The feasible profiles have
$|A||B|\le(\frac12+\kappa)(p-1)$, just below the threshold; a continuity argument in $w$ closes
the gap [Ref3].) Explicit feasible profiles of the corresponding linear programme,
checked constraint by constraint in exact integers, have every $e_b\approx t^\ast m$ with
$t^\ast=(1-U(w))/2\approx\frac23\kappa^2$. So a saving linear in $\kappa$ needs an input
beyond $(\star)$. The claim "combining $(\star)$ over several $e$ gives a linear saving" is
refuted by these profiles.

**9.3 Upper bounds for the true saving** ([Sh, Cor. 3.3], given Weil). As $p\to\infty$ the
optimal saving over $|A||B|\ge(\frac12+\kappa)p$ is at most
$1-\sup_m\beta_m(1/(2mw))$, where $\beta_m(q)$ is the mean of $1-2j/m$ over the lowest mass
$q$ of the Binomial$(m,\frac12)$ distribution. The two-element bound $2\kappa/(1+2\kappa)$ of
Remark 7.6 is the best bounded-$|A|$ upper bound exactly for $\kappa\le11/48$. Beyond that,
five-element sets do better; for example the saving is at most $29/80=0.3625$ at
$\kappa=1/2$.

**9.4 Balanced sets** ([Sh, §4], given Weil). Let $\tau_k$ be the infimum of $\tau$ such that
a constant saving holds for $|A||B|\ge(\tau+\kappa)p$, $\min(|A|,|B|)\ge k$ and large $p$.
Then $\tau_k\le1/2$ by Theorem D, $\tau_1=\tau_2=1/2$, and $\tau_k\ge k/2^k$ [Sh, Prop. 4.1].
Conversely, for pairs outside a *balanced window*
$M(k,\kappa)<|A|\le|B|<12\sqrt p$, with $M=\lceil12/\tau\rceil$ and $\tau=k/2^k+\kappa$, the
threshold $k/2^k$ suffices, with saving $\min(\kappa/(\tau M),\,0.159)$, provided $k\ge2$,
$0<\kappa\le1$ and $p\ge\max(25,(2M^2/\kappa)^2)$ [Sh, Thm 4.3; Ref3]. Inside the window only the threshold $1/2$ of Theorem D is known. A constant
saving for balanced sets of size $c\sqrt p$ with $c<1/\sqrt2$ is open. It would imply that
the bipartite graph $\chi(a+b)=1$ contains no $K_{c\sqrt p,c\sqrt p}$, a constant-factor
improvement of Hanson–Petridis for balanced bicliques.

**9.5 Characters of order $k$** ([R, Prop. 2.9], [Sh, §5]). Let $k\mid p-1$ and
$d_k=(p-1)/k$. Redefine $e_b=\#\{a\in A:a+b\ne0,\ (a+b)^{d_k}\ne1\}$. Then $(\star)$ holds with
$d$ replaced by $d_k$, for $m+d_k-1\le p-1$ and $2e\le\min(m-1,d_k)$. For a character $\psi$
of order $k$ and $|A||B|\ge(1+\lambda)(p-1)/k$ with $0<\lambda\le3$ and $p\ge(1+\lambda)k+1$,
this gives
$|S_\psi(A,B)|\le\sqrt{1-2\theta_1(1-\theta_1)(1-\cos(2\pi/k))}\,|A||B|$. Here
$m_1=\lceil\sqrt{(1+\lambda)d_k}\rceil$, $w_1=(d_k+m_1)/((1+\lambda)d_k)$ and
$\theta_1=(1-U(\min(1,w_1)))/2$ [Sh, Cor. 5.2], so
$\theta_1\to\theta(\lambda)=(1-U(1/(1+\lambda)))/2=\lambda^2/6+O(\lambda^3)$. The threshold
$(p-1)/k$ is sharp: take $A=\{0\}$ and $B$ the subgroup of index $k$.

**9.6 Why the Hankel step is needed** ([R, §1.2, §4]). Averaging Hanson–Petridis over
subsets $T\subseteq A$ gives only $(m-e)N_e-r_e\le K(m,e)\,d$, where
$N_e=\#\{b:e_b\le e\}$, $r_e=\#\{b\in-A:e_b\le e\}$ and
$K(m,e)=\min_{1\le t\le m-e}\binom mt/\binom{m-1-e}{t-1}\ge e\cdot\exp(1-e/m)$. This
constant grows with $e$ instead of tending to $1$. An explicit
abstract set system shows that "Hanson–Petridis for every subset" together with set sizes and
the second moment does not imply the robust inequality: there,
$m\cdot\#\{\text{coverage}\ge m-1\}=2(m/(m-1))^2d$.

---

## 10. Verification

### 10.1 Exact computations

| check | scope | count | failures | files |
|---|---|---|---|---|
| worker verifier [R] | every step of [R, §§1–2]; exhaustive for $p\le13$ and $p=17$ with $m\le6$; structured and random $A$ up to $p=211$ | 1,208,814 | 0 | `experiments/stepanov_robust_2026_09_26.py` → `results/stepanov_robust_2026_09_26.json` |
| independent exhaustive C check of $(\star)$ | every $A\ni0$ with $\lvert A\rvert\le(p+1)/2$, every admissible $e$, all primes $5\le p\le29$ | 1,074,401,579 cases | 0 | `experiments/stepanov_star_exhaustive_2026_09_27.c` → `results/stepanov_star_exhaustive_2026_09_27.json` |
| Theorem D on adversarial rectangles | the 165 most biased rectangles found by the adversary search | 165 | 0 | same JSON |
| referee [Ref] (code independent of the worker's) | own proof of every step; full $H_{e+1}$ in 2,509 cases with exact root multiplicities; $p\nmid\Lambda$ for every $p\le1500$, $m\le(p+1)/2$, $e\le12$; exhaustive over $A\ni0$ for $p\le23$; local search up to $p=2003$; Corollary E at 40 primes | 257,753,636 | 0 | `experiments/stepanov_referee2_2026_09_27.py` → `results/stepanov_referee2_2026_09_27.json` |
| second referee [Ref3] for Section 9 and Proposition 2.9 (independent code) | own re-derivations; $(\star)_k$ for every $A\ni0$ over 15 $(p,k)$ pairs with $p\le19$; full $H_{e+1}$ in 413 order-$k$ cases; Theorems 1.4–1.5 of [Sh] over all pairs up to affine maps for $p\le23$; order-$k$ bias exhaustively at $p=7,13$ | 5,555,096 | 0 | `experiments/stepanov_referee3_2026_09_29.py` → `results/stepanov_referee3_2026_09_29.json` |
| sharpening verifier [Sh] (Section 9 only) | as described in [Sh] | 5,911,513 | 0 | `experiments/stepanov_sharpen_2026_09_27.py` → `results/stepanov_sharpen_2026_09_27.json` |
| this paper's verifier | see below | 24,349,465 | 0 | `experiments/paper_robust_hp_2026_09_29.py` → `results/paper_robust_hp_2026_09_29.json` |

The workspace's orchestrator also read Steps 1–5, Lemma 4.2 and Theorems C and D line by line
and found no error [PS, §2].

The verifier for this document uses the standard library plus numpy, sympy and mpmath (a sympy dependency), and runs in
about 18 s of CPU time. It recomputes the following.

* **§A.** Lemma 2.2 and Lemma 3.1(a)–(d), for every $A\ni0$ with $|A|\le(p+1)/2$ at
  $p\le13$, and for random $A$ at $17\le p\le43$.
* **§B.** 578 full determinants $H_{e+1}$ for $11\le p\le151$ and $e\le4$, using intervals,
  $Q\cup\{0\}$ (the extreme case $D=p-1$), random sets and greedy sum-cliques. For each it
  checks $u_s=F^{(s)}/(D)_s$, the vanishing top coefficients and the leading coefficients
  $\binom{D-s}{m-1}$, the exact degree $(e+1)(d-e)$, the leading coefficient $\Lambda$ (by
  separate modular elimination), and root multiplicities by synthetic division at every $b$
  (44,056 order checks). It also checks Step 1 directly (75,911 checks).
* **§C.** $(\star)$; Corollary B for all $B$ simultaneously (the worst $B$ is
  $\{b:w_b>0\}$); both displays of Theorem 6.1 on level sets; the $\eta$-form of Theorem C
  in exact rational arithmetic; Proposition 6.2 ($n\le1$ when $m>(p+1)/2$); the second-moment
  inequality of Remark 6.3. All of this is exhaustive over $A\ni0$ for $p\le19$, with random
  and structured $A$ at $p=101,211,409$ for Theorem 6.1 and the $\eta$-form. It also checks
  the count in Proposition 6.2 for every $c\ne0$ and $p\le97$, and the equality witness of
  Remark 4.5.
* **§D.** Theorem D and Lemma 7.1 for every $A\ni0$ and every $n$ at $p\le19$, using that
  the extremal $B$ for given $A,n$ consists of the $n$ largest or smallest row sums (float,
  margin $10^{-9}$, 200 values of $\kappa$ per pair). Also the family $|A|=2$ and random small
  $A$ up to $p=10009$ (100 values of $\kappa$ per pair), and the non-vacuity thresholds of
  Remark 7.4.
* **§E.** Proposition 7.5 for all primes $5\le p\le3000$, and Remark 7.6 at $p=10009$ and
  $p=100049$.
* **§F.** The product formula of Lemma 4.2 against exact integer determinants for $D\le36$;
  $p\nmid\Lambda$ for every $p<300$, $m\le(p+1)/2$, $e\le5$, and for general $d$ at $p\le41$;
  each step of Proof 2; the counterexample with $D\ge p$.
* **§G.** The series for $(1-u)^2$ and $\eta^\star$, the crossing point $(3-2\sqrt2)/2$, and
  the constants of Section 9.
* **§H.** Corollary E exhaustively for small $m$ at $p\in\{13,17,29,37,41\}$, and for greedy
  dense and sparse sets at $p\in\{101,401,1009\}$.
* **§I.** Lemma 2.4.

The checks are not vacuous: the bounds of Theorem A and of Step 3 are attained (Remark 4.5),
and Proposition 7.5 is an equality case of Corollary B.

### 10.2 Lean formalization

`experiments/stepanov_leanrobust_lean/StepanovRobust.lean` (1,384 lines, `import Mathlib`) is
described in [L]. It states and proves the following.

* `hankel_stepanov_paley`: $(\star)$ for $p\ne2$, $|A|\le(p+1)/2$, $2e+1\le|A|$.
* `hankel_det_natDegree_eq`: the exact degree $(e+1)(d-e)$.
* `hankel_order_paley`: the divisibility $(x-b)^{T_b}\mid H_{e+1}$.
* `hankelLead_ne_zero`: $p\nmid\Lambda$ by Proof 2, for any $d$ with $2e\le d$, $2e+1\le m$,
  $d+m-1<p$.
* `bias_inequality` and `bias_bound`: Corollary B.
* `subset_bound`: Lemma 7.1.
* `constant_bias`: Theorem D, verbatim.

Not formalized: Theorem C, the index-$k$ version (Section 9.5), the identity
$u_s=F^{(s)}/(D)_s$ (not needed there), and the exact value of $\Lambda$.

Compilation status (checked by the orchestrator on 2026-09-29; [L, §4]):

* The file compiles with exit code 0 and no errors in two environments, by calling the
  toolchain's `lean` binary with `LEAN_PATH` set to prebuilt package directories
  (read-only): Lean v4.30.0-rc2 with Mathlib `5450b53`, and Lean v4.29.1 with Mathlib
  `5e932f9`.
* In both environments all 17 `#print axioms` lines read exactly
  `propext, Classical.choice, Quot.sound`; in particular `sorryAx` does not occur, so
  `hankel_stepanov_paley`, `bias_inequality`, `bias_bound` and `constant_bias` are
  machine-checked.
* One resource setting was needed: under Mathlib `5e932f9` the proof of `subset_bound`
  exceeds the default heartbeat limit, so it carries `set_option maxHeartbeats 1000000`;
  this changes no logic.
* `experiments/stepanov_leanrobust_2026_09_27.py` re-runs both compilations and 202,148
  arithmetic checks (0 failures); its results are in
  `results/stepanov_leanrobust_2026_09_27.json`.
* Not done: compilation against the prove2me Mathlib pin `c5ea003`, and server-side
  verification on prove2me.

---

## 11. Related work and novelty

* **Hanson–Petridis** [HP]: exact containment $A+B\subseteq Z_d\cup\{0\}$. They describe
  Vinogradov's estimate as the best known bound for general sets, and state that their method
  is restricted to prime fields. There is no robust form in [HP].
* **Kalmynin**, arXiv:2504.10202 [Kal]: $d$-critical pairs, where $A+B\subseteq\mu_d\cup\{0\}$
  and $|A||B|=d+|(-A)\cap B|$, and the resolution of Sárközy's conjecture for quadratic
  residues. Only exact containment and equality are treated.
* **Rudnev–Tyrrell**, arXiv:2607.24270 [RT]: proper subgroups are not sumsets except in the
  case $|A|=|B|=2$, $|H|=4$. Only exact decompositions are treated.
* The following were checked from abstracts or local texts, and all concern exact
  decompositions or containments:
  * Yip, arXiv:2501.16620;
  * Yip–Yoo, arXiv:2608.02568;
  * Cochrane, arXiv:2607.28559;
  * Kim–Yip–Yoo, arXiv:2602.20919 and arXiv:2607.24370.
* **Shkredov**, *Sumsets in quadratic residues* (Acta Arith. 164 (2014); arXiv:1305.4093):
  Theorem 3.2 treats approximate equality $A+A\approx R$ by second-moment arguments. It is not
  a bias bound for arbitrary pairs in the window.
* **Heath-Brown–Konyagin** (Q. J. Math. 2000): Stepanov's method for several shifts, with
  exact membership in a subgroup.
* **Volostnov** (arXiv:1712.09355) and **Schoen–Shkredov** (arXiv:2004.01885):
  character-sum bounds near $\sqrt p$ under small-doubling hypotheses.
* **Yip's PhD thesis** (UBC, 2024): not fetched (**UNVERIFIED**).
* The rank observation of Step 2 is the familiar fact that a Hankel matrix of power sums
  with $e_b$ terms has rank at most $e_b$, the principle behind syndrome decoding. We did not
  find it combined with Stepanov's method.

Sources: [R, §3] (searched 2026-09-26; that search hit a session limit before two items could
be fetched) and [Ref, §10] (searched 2026-09-27).

**Novelty.** Targeted searches found no prior statement of $(\star)$, of a robust
Hanson–Petridis inequality, or of any constant-factor bias bound for arbitrary sets in the
window $p/2<|A||B|<p$. This is a negative search result by AI agents with web tools, not a
proof of novelty, and novelty is therefore **not established**. A specialist's reading is
needed.

---

## 12. Open problems

1. **Threshold for growing sets.** Is there $\tau<1/2$ such that $|S(A,B)|\le(1-c)|A||B|$
   whenever $|A||B|\ge\tau p$ and $\min(|A|,|B|)\to\infty$? The case of balanced sets of size
   $c\sqrt p$ with $c<1/\sqrt2$ is the heart of this (Section 9.4).
2. **The optimal saving.** Between $\eta^\star(\kappa)\approx\frac43\kappa^2$ (Section 9.1;
   $(1-u)^2\approx\kappa^2$ refereed) and $2\kappa/(1+2\kappa)$ (Remark 7.6). By Section 9.2,
   $(\star)$ alone cannot give a saving linear in $\kappa$.
3. **Degree budget.** The degree $(e+1)(d-e)$ caps the method at $|A||B|\asymp p$. Any
   progress below the square-root scale, or any power saving, needs an auxiliary polynomial
   with a smaller degree budget. Such an argument must fail over $\mathbb F_{p^2}$
   (Remark 4.4).
4. **Review and formal verification.** A human specialist's reading of Section 4; formal
   verification of Theorem C; a stored compile log and server-side checking of [L].
5. **Characters of order $k$.** A referee of Section 9.5.

---

## References

* [HP] B. Hanson and G. Petridis, *Refined estimates concerning sumsets contained in the roots
  of unity*, arXiv:1905.09134v3 (2020); Theorem 1.2, Corollary 1.5 and §2. Local copy:
  `sources/sigma-hanson-petridis-1905.09134.txt`. Published in Proc. London Math. Soc. (2021)
  according to [R, §3]; we did not check the journal version.
* [Sa] S. Satake, *On the restricted isometry property of the Paley matrix*, arXiv:2011.02907v2
  (2020); Conjecture 7 and Remark 9. Local copy: `sources/satake-2011.02907.txt`.
* [Ka] A. A. Karatsuba, *Distribution of values of Dirichlet characters on additive
  sequences*, Soviet Math. Dokl. 44 (1992), 145–148. Not consulted directly; the statement
  used is that of [FSX].
* [FSX] É. Fouvry, I. E. Shparlinski and P. Xi, *Estimates for trilinear and quadrilinear
  character sums*, Rev. Mat. Iberoam. 41 (2025), 1925–1956; (1.4)–(1.6). Local copy:
  `sources/bilinear-2025.txt`.
* [V] I. M. Vinogradov, *Elements of Number Theory*, Dover (1954), Ch. V, Ex. 8, as cited in
  [HP].
* [Kr] C. Krattenthaler, *Advanced determinant calculus*, Sém. Lothar. Combin. 42 (1999),
  Art. B42q; arXiv:math/9902004, Theorem 26, eq. (3.12). Local copy:
  `sources/stepanov-robust-krattenthaler-math9902004.txt`.
* [Kal] A. Kalmynin, *On additive irreducibility of multiplicative subgroups*,
  arXiv:2504.10202v2 (2025). Local copies: `sources/sigma-kalmynin-2504.10202v2.txt`,
  `sources/stepanov-robust-kalmynin-2504.10202.txt`.
* [RT] M. Rudnev and F. Tyrrell, *Multiplicative subgroups of prime fields are not sumsets*,
  arXiv:2607.24270v2 (2026). Local copy: `sources/stepanov-robust-rudnev-tyrrell-2607.24270.txt`.
* [Sh14] I. D. Shkredov, *Sumsets in quadratic residues*, Acta Arith. 164 (2014), 221–243;
  arXiv:1305.4093.
* [HBK] D. R. Heath-Brown and S. V. Konyagin, *New bounds for Gauss sums derived from $k$th
  powers, and for Heilbronn's exponential sum*, Q. J. Math. (2000), 221–235.
* [Vol] A. S. Volostnov, *On some double sums with multiplicative characters*,
  arXiv:1712.09355 (2017). Local copy: `sources/sigma-volostnov-shkredov-1712.09355.txt`.
* Workspace notes: [R] `research/stepanov-robust-2026-09-26.md`;
  [Ref] `research/stepanov-referee2-2026-09-27.md`;
  [Ref3] `research/stepanov-referee3-2026-09-29.md` (second referee: Section 9 and Proposition 2.9 of [R]).
  [Sh] `research/stepanov-sharpen-2026-09-27.md`;
  [L] `research/stepanov-leanrobust-2026-09-27.md`;
  [PS] `research/stepanov-pass-summary-2026-09-27.md`;
  brief `research/sigma-brief-2026-09-05.md`.
