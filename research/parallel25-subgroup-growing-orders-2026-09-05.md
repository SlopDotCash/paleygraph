# Growing subgroup orders: repeated-coordinate upper bound and the remaining six-distinct family

**Status: a uniform upper bound is proved for the repeated-coordinate part
of the opposite-free sixth energy, aggregated over every product-ratio
fiber. It is (O(n^{69/20}(1+\log n)^{1/5})), using a published subgroup
additive-energy theorem. No improved upper bound for the full unbalanced
aggregate, full sixth-energy exponent, square-root character cancellation,
Paley conjecture, or official-prize bridge is proved.**

This pass also proves an exact rank interpretation of the cyclic norm
restriction for arbitrary dyadic orders. Actual subgroup relations at
orders 64, 128, and 256 show that all rationally independent cyclic shifts
of one relation can contribute only one factor of the characteristic to
its norm. These are limitations of that proposed amplification argument,
not counterexamples to the target conjecture. The conclusions below are
mathematical derivations and finite machine checks, not human or formal
certification. No claim of literature novelty is made.

## 1. Definitions and the existing intrinsic ledger

Let (H\le\mathbb F_p^*) have dyadic order (n\ge4), where (p) is
prime and (n^4/4\le p\le n^4). Thus (p>6), (n^2<p), and
(-1\in H). Set

\[
 r_2(s)=\#\{(a,b)\in H^2:a+b=s\},\qquad
 r_4(s)=\#\{(a,b,c,d)\in H^4:a+b+c+d=s\},
\]
\[
 E_2(H)=\sum_s r_2(s)^2,\qquad
 E_3(H)=\#\{h\in H^6:h_1+\cdots+h_6=0\}.
\]

Our (E_3) equals the equal-sums count for two ordered triples, because
(-H=H). It is **not** the third moment of difference multiplicities.

Call a six-word intrinsic when its multiset is a union of three opposite
pairs, and opposite-free when it contains no opposite pair at all. The
intrinsic fourth- and sixth-word counts are

\[
 T_4=3n^2-3n,\qquad T_6=15n^3-45n^2+40n.
\]

For a fixed split into two triples, replace the second triple by its
negative and let (B_\rho) count equal-sum triple pairs with product
ratio (\rho\in H). The exact identities from
[pass 23](parallel23-subgroup-upper-2026-09-05.md) are

\[
 B_\rho=nN_\rho,\quad
 N_\rho=\sum_z w(z)w(z/\rho),\quad
 w(z)=\#\{(a,d)\in H^2:ad-a-d=z\}.
\]

The intrinsic contribution to (B_\rho), proved in
[pass 24](parallel24-subgroup-unbalanced-2026-09-05.md), is exactly

\[
 I_\rho=\begin{cases}
 6n^3-9n^2+4n,&\rho=1,\\
 18n(n-2),&\rho\in H^2\setminus\{1\},\\
 0,&\rho\notin H^2.
 \end{cases}                                                   \tag{1}
\]

Here (H^2=\{h^2:h\in H\}). Thus (\sum_\rho I_\rho=T_6), and
the residual (nN_\rho-I_\rho) counts nonintrinsic words, including
those that still contain an opposite pair. These conventions remain
unchanged throughout this pass.

## 2. A uniform upper bound for every repeated-coordinate residual

For an ordered word (h=(h_1,\ldots,h_6)), let (m_a(h)) be the
multiplicity of (a\in H), and write

\[
 \nu(h)=\sum_{a\in H}\binom{m_a(h)}2.
\]

This counts the unordered positional pairs with equal entries. Fix one
of the 15 positional pairs and set its common value to (a\in H).
The remaining four entries must sum to (-2a). Multiplication by
(a^{-1}) is a bijection of (H), and negation is another one, so their
number is (r_4(-2)=r_4(2)), independent of (a). Hence

\[
 \sum_{h\in H^6:\,\sum h_i=0}\nu(h)=15n r_4(2).                 \tag{2}
\]

For the same fixed positional pair, normalize its repeated value to
1 and count only intrinsic words. The remaining multiset must be

- (\{-1,-1,v,-v\}), with (\{v,-v\}\ne\{1,-1\}): there are
  (n/2-1) choices of the opposite pair and 12 orderings;
- (\{1,-1,-1,-1\}): there are 4 orderings.

These cases are disjoint and exhaustive, giving (6n-8). Therefore

\[
 \sum_{h\text{ intrinsic}}\nu(h)=15n(6n-8),\qquad
 \sum_{h\text{ nonintrinsic},\,\sum h_i=0}\nu(h)
   =15n\bigl[r_4(2)-6n+8\bigr].                                \tag{3}
\]

Let (R_{6,\mathrm{rep}}) count opposite-free zero-sum words with at
least one repeated coordinate. Each is nonintrinsic and has
(\nu(h)\ge1). Also (r_4(2)=\sum_s r_2(s)r_2(2-s)\le E_2)
by Cauchy–Schwarz. We obtain the explicit uniform bound

\[
 \boxed{0\le R_{6,\mathrm{rep}}
 \le15n\bigl[r_4(2)-6n+8\bigr]\le15nE_2(H).}                  \tag{4}
\]

The middle expression is an exact weighted count of all nonintrinsic
repetitions; it is generally larger than the opposite-free subset.
In particular, it is not being asserted equal to that subset. Formula
(3) also gives the lower bound (r_4(2)\ge6n-8).

Let (B^{\mathrm{rep}}_\rho) be the portion of (B_\rho) whose
corresponding six-word is opposite-free and repeated. The fixed split
partitions this subset by its ratio, so (4) says

\[
 \sum_{\rho\in H}\frac{B^{\mathrm{rep}}_\rho}{n}
 \le15\bigl[r_4(2)-6n+8\bigr].                               \tag{5}
\]

Thus this controls a specified part of the actual residual aggregate
over all ratios, including the unbalanced ones; there is no extra factor
of (n) from summing a separate bound for each ratio.

Theorem 3 of [Murphy–Rudnev–Shkredov–Shteinikov, *On the few products,
many sums problem*, arXiv:1712.00410v1](https://arxiv.org/html/1712.00410)
states, for multiplicative subgroups of size at most (\sqrt p),

\[
 E_2(H)\ll n^{49/20}\log^{1/5}n.
\]

The subgroup hypothesis and exact exponent were checked against the
primary text. Its (\mathsf E(A)) is our (E_2). Its separately
defined (\mathsf E_3(A)=\sum_x r_{A-A}(x)^3) is a different
quantity; its cubic-energy estimates are not used here. The quartic
window satisfies the theorem's size hypothesis. Combining it with (4)
gives, uniformly over growing dyadic orders,

\[
 \boxed{R_{6,\mathrm{rep}}
   \ll n^{69/20}(1+\log n)^{1/5}.}                            \tag{6}
\]

This is a sub-fourth-power bound for this restricted class. It does not
bound the remaining six-distinct class.

## 3. Quantitative localization to six distinct, fully unbalanced words

Let (J_6) count nonintrinsic zero-sum words containing an opposite
pair. Let (R_{6,\mathrm{dist,bal}}) count opposite-free words with
six distinct entries that have opposite products across at least one
split into triples. Finally let (D_6) count opposite-free words with
six distinct entries whose products are unbalanced at all ten splits.
There is the exact disjoint partition

\[
 E_3=T_6+J_6+R_{6,\mathrm{rep}}+R_{6,\mathrm{dist,bal}}+D_6.    \tag{7}
\]

For completeness, the earlier bound (J_6\le15n(E_2-T_4)) follows by
marking an opposite positional pair, choosing its first value, and
noting that the remaining four-word is nonintrinsic. Summing the 15
marked possibilities only overcounts.

As in pass 23, put

\[
 X=E^\times((H-1)\setminus\{0\})-
       \bigl[2(n-1)^2-(n-1)\bigr]\ge0.
\]

The previous balanced-word estimate gives
(R_{6,\mathrm{dist,bal}}\le R_{6,\mathrm{bal}}\le10nX).
Its zero-product convention is
(E^\times(H-1)=6n^2-9n+4+X). Substituting (4) into (7) gives the
explicit restriction

\[
 \boxed{0\le E_3-T_6-D_6
 \le15n(E_2-T_4)+15n\bigl[r_4(2)-6n+8\bigr]+10nX.}           \tag{8}
\]

[Shkredov, *On tripling constant of multiplicative subgroups*,
arXiv:1504.04522v1, §4, Theorem 6](https://arxiv.org/html/1504.04522v1)
gives (E^\times(H-1)\ll n^2(1+\log n)) when (n^2<p), by
taking both subgroups to be (H) and both nonzero shifts to be (-1).
The theorem and these hypotheses were rechecked in the primary text.
Together with (6), this yields

\[
 E_3=T_6+D_6+
 O\!\left(n^{69/20}(1+\log n)^{1/5}+n^3(1+\log n)\right),    \tag{9}
\]

with a nonnegative error before applying the (O)-notation. Thus any
excess above this error scale must be carried by the six-distinct,
fully unbalanced family. No upper bound improving the known scale for
(D_6) is obtained. Equation (9) is a localization result with an
independently bounded error, not a bound on the unknown (E_3) obtained
by renaming it.

## 4. What multiple cyclic relations do and do not force

Write (n=2^k), (d=n/2), and choose a generator (g\in H). Encode
an opposite-free six-word by its exponent polynomial modulo
(X^d+1\), and divide its integer coefficients by their common content:

\[
 f(X)=\sum_{j=0}^{d-1}a_jX^j,\qquad Q=\sum_j a_j^2.
\]

No opposite entries means that coefficient cancellation does not occur
in this reduction. Division by content preserves the zero at (g),
since the content is at most 6 and (p>6). One supported value would
give (6a=0), which is impossible; consequently (Q\le26), with
the largest possible undivided pattern ((5,1)). Six distinct entries
give (Q=6).

The polynomial (X^d+1=\Phi_n(X)) is irreducible over (\mathbb Q):
its translate by 1 is Eisenstein at 2, since (d) is a power of 2.
Thus the nonzero polynomial (f), of degree less than (d), has
nonzero integer norm

\[
 N=\left|\operatorname{Norm}_{\mathbb Q(\zeta_n)/\mathbb Q}
                      f(\zeta_n)\right|.
\]

As in pass 24, orthogonality gives
(d^{-1}\sum_{u\bmod n\text{ odd}}|f(\zeta_n^u)|^2=Q), and
the arithmetic-geometric mean inequality gives (0<N\le Q^{d/2}).

Let (M_f) be multiplication by (f) on the integer module
(\mathbb Z[X]/(X^d+1)), in the usual monomial basis, and define

\[
 z_p(f)=\#\{u\bmod n:\ u\text{ odd},\ f(g^u)=0\pmod p\}.
\]

Because (p\equiv1\pmod n), the reduction of (X^d+1) splits
into (d) distinct roots (g^u), with (u) odd. Evaluation at
these roots diagonalizes multiplication by (f). It follows exactly
that

\[
 \operatorname{rank}_{\mathbb Q}M_f=d,\qquad
 \operatorname{rank}_{\mathbb F_p}M_f=d-z_p(f).                \tag{10}
\]

For an integer nonsingular matrix, a rank defect (z) modulo (p)
forces (p^z) to divide its determinant. One direct proof performs
invertible row and column operations over the localization
(\mathbb Z_{(p)}) until the mod-(p) independent block is the
identity; the remaining (z\times z) block has all entries divisible
by (p), and its determinant supplies (z) factors. Since
(|\det M_f|=N), (10) proves

\[
 \boxed{z_p(f)\le v_p(N),\qquad
 z_p(f)\log p\le\frac n4\log Q.}                            \tag{11}
\]

In the quartic window this is an (O(n/\log n)) restriction on the
number of conjugate zeros of any one such relation polynomial. For
six distinct entries, use (Q=6). This is a uniform structural
restriction, but it does not bound how many different relation
polynomials can vanish at the same characteristic.

Multiplicative closure produces the (d) reduced relations
(X^jf\pmod{X^d+1}), (0\le j<d), by scaling the original word
by (g^j). They are the columns of (M_f), and hence are rationally
independent. All vanish at (g). Rational independence does not imply
that they supply different mod-(p) root constraints: the mod-(p)
rank can be (d-1), and (v_p(N)) can be 1. Multiplying the norms
of all the shifts multiplies the logarithmic norm upper bound by the
same number as the guaranteed prime factors, so this alone gives no
stronger inequality. A useful bound would require additional
arithmetic independence, not simply closure of one relation's orbit.

## 5. Actual arithmetic witnesses and exact finite checks

The new [verifier](../experiments/parallel25_subgroup_growing_orders_2026_09_05.py)
and [result](../results/parallel25_subgroup_growing_orders_2026_09_05.json)
use four actual prime subgroups, with trial-division primality checks,
exact generators, and explicit multiplicative closure checks. Exhaustive
normalized six-word enumeration is compared with all actual
product-ratio fibers computed separately from (w). It checks (2),
(3), (4), (7), and (8), and retains every repeated multiplicity pattern.
Zero entries in the next table mean that the corresponding counted
class is empty.

| (p) | (n) | (R_{6,\mathrm{rep}}) | Upper bound in (4) | (D_6) | (E_3-T_6-D_6) | Upper bound in (8) |
|---:|---:|---:|---:|---:|---:|---:|
| 6,700,417 | 64 | 367,680 | 1,477,440 | 0 | 1,066,560 | 2,287,680 |
| 7,204,033 | 64 | 168,960 | 276,480 | 0 | 1,735,680 | 2,188,800 |
| 67,403,009 | 128 | 0 | 184,320 | 276,480 | 6,266,880 | 7,004,160 |
| 1,073,748,737 | 256 | 0 | 0 | 368,640 | 0 | 0 |

The final row has (E_2=T_4=195840), (X=0), and
(r_4(2)=1528=6n-8). Accordingly the entire nonintrinsic excess is
(E_3-T_6=D_6=368640). This is a concrete remaining six-distinct,
fully unbalanced class, not a hypothetical obstruction outside actual
subgroups. Its size does not violate an (O(n^3)) target.

The following norm witnesses use the displayed generators. Exponents
list the six entries (g^{e_i}); all their zero-sum words are
opposite-free and unbalanced at every split.

| (p,n) | (g) | Exponents (e_i) | (Q) | Integer norm (N) |
|---|---:|---|---:|---:|
| 6,700,417; 64 | 1,688,191 | 0,14,21,0,0,0 | 18 | 18447047002047976322 |
| 7,204,033; 64 | 4,105,673 | 0,19,47,0,10,34 | 8 | 3164818145296 |
| 67,403,009; 128 | 64,701,253 | 0,116,26,91,38,14 | 6 | 103665827842 |
| 1,073,748,737; 256 | 1,064,280,392 | 0,26,56,129,217,116 | 6 | 4419283227438373820201815271409083396 |

In every row, (z_p(f)=v_p(N)=1), the sole primitive zero exponent
is (u=1), and the ranks over (\mathbb Q,\mathbb F_p) are
((d,d-1)). The verifier computes the norm by repeated exact quadratic
norm descent, computes mod-(p) matrix rank separately, and evaluates
all primitive roots and all cyclic-shift columns. The nonzero norm
certifies the rational rank via the multiplication determinant identity.
The last two rows have six distinct entries. These finite witnesses
establish the stated failure of cyclic-shift amplification at growing
orders beyond the prior finite classification; they do not establish
an infinite family with those properties.

The completed run reports **PASS**, with 120,585 normalized multiset
matches, 512 independent product-ratio fiber checks, 256 cyclic-shift
vanishing checks, and four exact checks each of the displayed incidence,
intrinsic subtraction, rank-defect, and norm-divisibility assertions.
The imported energy theorems and the general algebraic arguments are
not certified by these finite tests. No estimate for the unresolved
six-distinct aggregate follows from the tests.

## 6. Pinned inputs and outputs

SHA-256 values at the completed run:

| Artifact | SHA-256 |
|---|---|
| `experiments/parallel25_subgroup_growing_orders_2026_09_05.py` | `5b5462ca986bbf08e2f581d898cef78525d220c9f2e1cad2a25163d4c795f83d` |
| `results/parallel25_subgroup_growing_orders_2026_09_05.json` | `000370af584720ee83e8bbabeed282024a0a630a625a129b274d0bd98fa2c35a` |
| `research/parallel24-subgroup-unbalanced-2026-09-05.md` | `ffe63245f28880ff53ace7210bf54f3baafca77c110b5aec0ac42cd3bf99bf1b` |
| `research/parallel23-subgroup-upper-2026-09-05.md` | `1b97e491bc1c0b2cfe391388b0340d89d1f2a7168822db116513b63f4148bd41` |
| `research/parallel21-subgroup-next-input-2026-09-05.md` | `84005053f3954ee2d36eccf85b2786f1fe00d6bd956c76b23b1b40e82b991c22` |
| `results/parallel23_subgroup_upper_2026_09_05.json` | `d675ecc47be55cde21ba176fbd5d3bfdf9c0da35af79c093079d31e47cb82770` |
| `sources/sigma-subgroup-2026-09-05/mrss-1712.00410.pdf` | `9154e144adaecd1aeb7ce87265b548328efd5ec811bb48d211043f41b8493477` |
| `sources/sigma-subgroup-2026-09-05/mrss-1712.00410.txt` | `b58780e7733f9a82373ab98a92d6b60e3c8ed6445dfc79b5291b8ab188acff1c` |
| `sources/mixed-periods-2026-09-04/shkredov-1504.04522.html` | `be13a230cc3cfdb780fd8dca0606d35a6f561221cb3ba67a860be6f759a1f459` |

Prior notes, experiments, and central artifacts were not edited.
