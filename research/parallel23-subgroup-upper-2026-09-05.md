# Product-ratio fibers and a partial upper bound for opposite-free sixth energy

**Status: an actual subgroup upper bound is proved for the product-balanced
part of the opposite-free sixth energy. The unbalanced remainder is not
controlled at the required order. No new total-energy or character-sum
exponent, full Paley conjecture, or official-prize bridge is proved.**

The new bound is

\[
 R_{6,\mathrm{bal}}(H)\le 10n\mathcal X(H)
 \ll n^3(1+\log n),
 \qquad n^2<p.                                      \tag{1}
\]

Here “balanced” means that at least one partition of the six positions
into two triples has opposite products. The proof uses multiplicative
closure essentially. It also identifies every product-ratio fiber with
a dilation correlation of a shifted-product distribution. Applying the
available shifted-energy theorem to all fibers separately gives only
the fourth-power scale for the total sixth energy.

A genuine circular dyadic subgroup in the requested quartic window has
\(\mathcal X=0\) but \(R_6>0\). Thus circularity removes the balanced
part completely while leaving actual unbalanced relations. This example
does not contradict an \(O(n^3)\) bound or the target conjecture.

## Definitions and imported inputs

Let \(H\le\mathbb F_p^*\) have dyadic order \(n\), with \(p>3\).
In the requested range assume \(n\ge4\) and
\(n^4/4\le p\le n^4\); in particular \(n^2<p\). Since \(-1\in H\),
the ordered additive energy

\[
 E_3(H)=\#\{(x,y)\in H^3\times H^3:\sum_i x_i=\sum_i y_i\}
\]

also counts ordered zero-sum six-tuples. Let \(R_6(H)\) count those
zero-sum tuples containing no pair of opposite entries.

Use the previously proved collision quantity

\[
 R=(H-1)\setminus\{0\},\qquad
 \mathcal X(H)=E^\times(R)-[2(n-1)^2-(n-1)]
             =\sum_{i,j}(C_{ij})_3\ge0,              \tag{2}
\]

where \(C_{ij}=\#\{z\in H_i:1+z\in H_j\}\), \(H_i\) are the
multiplicative cosets, and \((c)_3=c(c-1)(c-2)\). The exact collision
bijection is in [the mixed-period note](mixed-periods-and-shifted-energy.md).
In particular, \(\mathcal X=0\) is equivalent to \(C_{ij}\le2\) for
every cell, or circularity of the affine copies of \(H\).

The sole analytic input for (1) is [Shkredov, *On tripling constant of
multiplicative subgroups*, arXiv:1504.04522v1, §4, Theorem 6](https://arxiv.org/html/1504.04522v1).
For subgroups \(\Gamma,\Pi\) with \(|\Gamma||\Pi|<p\) and nonzero
shifts \(a,b\), that theorem gives

\[
 E^\times(\Gamma+a,\Pi+b)
 \ll |\Gamma||\Pi|\log\min(|\Gamma|,|\Pi|)
       +|\Gamma|^2+|\Pi|^2.
\]

Its energy convention includes zero. Taking both subgroups to be \(H\)
and both shifts to be \(-1\) gives

\[
 E^\times(H-1)\ll n^2(1+\log n).                     \tag{3}
\]

The primary HTML statement and its convention were inspected during
this pass; the previously archived HTML is pinned below. This is an
imported theorem, not a formalized result of this workspace.

## Exact product-ratio fibers

For \(\rho\in H\), put

\[
 B_\rho=\#\{(x,y)\in H^3\times H^3:
           \sum_i x_i=\sum_i y_i,
           \prod_i x_i=\rho\prod_i y_i\}.
\]

Define the shifted-product distribution and its translate by

\[
 r(v)=\#\{(a,d)\in H^2:(a-1)(d-1)=v\},\qquad
 w(z)=r(z+1)=\#\{(a,d)\in H^2:ad-a-d=z\}.
\]

Then the following identities hold exactly:

\[
 \boxed{B_\rho=nN_\rho,\qquad
 N_\rho=\sum_{z\in\mathbb F_p}w(z)w(z/\rho),\qquad
 E_3(H)=n\sum_{\rho\in H}N_\rho.}                    \tag{4}
\]

To prove the first identity, start with a pair of triples counted by
\(B_\rho\) and set

\[
 \lambda=\frac{x_1}{y_2y_3}\in H,\qquad
 a=\lambda y_2,\quad d=\lambda y_3,\quad
 b=\lambda x_2/\rho,\quad c=\lambda x_3/\rho.
\]

Multiplicative closure puts all four new variables in \(H\). The
product-ratio equation gives the scaled triples

\[
 \lambda x=(ad,\rho b,\rho c),\qquad
 \lambda y=(\rho bc,a,d).
\]

Their sums are equal exactly when

\[
 ad-a-d=\rho(bc-b-c).                                 \tag{5}
\]

Conversely, any solution \((a,b,c,d)\in H^4\) of (5), together with
any \(\lambda\in H\), gives the two displayed triples divided by
\(\lambda\). Their product ratio is \(\rho\), and recovering
\(x_1/(y_2y_3)\) recovers precisely \(\lambda\). The constructions
are inverse. Counting (5) gives \(N_\rho\), proving (4).

At \(\rho=1\), translation of the index of summation gives

\[
 B_1=nE^\times(H-1).                                  \tag{6}
\]

For every \(\rho\), multiplication of the index by \(\rho\) is a
permutation of \(\mathbb F_p\), so Cauchy–Schwarz yields

\[
 N_\rho\le N_1=E^\times(H-1).
\]

Consequently (3) and (4) give only

\[
 E_3(H)\le n^2E^\times(H-1)\ll n^4(1+\log n).         \tag{7}
\]

The factor lost is explicit: there are \(n\) correlations. A total
\(E_3=O(n^3(1+\log n)^A)\) bound would require
\(\sum_{\rho\in H}N_\rho=O(n^2(1+\log n)^A)\), whereas (3)
bounds each correlation separately by that order. Identity (4) alone
supplies no aggregate saving.

## Removing permutation pairs and isolating the bounded part

The shifted set \(H-1\) has \(n\) elements, one of them zero. Hence
there are \(2n-1\) ordered pairs with zero product, and

\[
 E^\times(H-1)=E^\times(R)+(2n-1)^2
             =6n^2-9n+4+\mathcal X(H).                \tag{8}
\]

The count of ordered pairs of triples that are permutations of each
other is

\[
 D_3(n)=36\binom n3+9n(n-1)+n
       =6n^3-9n^2+4n.
\]

The three terms correspond respectively to three distinct entries,
exactly two equal entries, and all three equal. Combining (6) and (8)
therefore proves

\[
 \boxed{B_1-D_3(n)=n\mathcal X(H).}                   \tag{9}
\]

This is the exact count of nonpermutation equal-sum, equal-product
triple pairs. It is not the exact count of opposite-free six-tuples.
If two equal-sum, equal-product triples share an entry, canceling that
entry leaves two pairs with the same sum and product. The quadratic
polynomial with those elementary symmetric coefficients recovers each
unordered pair, so the full triples are permutations. Thus a
nonpermutation pair has disjoint cross-support. This argument does not
exclude opposite entries inside either individual triple.

For a fixed partition \(I\mid I^c\) of the six positions into triples,
call a zero-sum word balanced at this partition when

\[
 \prod_{i\in I}h_i=-\prod_{i\in I^c}h_i.              \tag{10}
\]

Keep the first triple and negate the second. Its two sums and two
products now agree. If the original word is opposite-free, these
triples cannot be permutations, since even a shared value across them
would exhibit an opposite pair in the original word. Therefore (9)
gives the subset inequality

\[
 R_{6,\mathrm{bal},I}(H)\le n\mathcal X(H).            \tag{11}
\]

There are \(\binom63/2=10\) unordered partitions. Taking their union,
using \(\mathcal X\le E^\times(H-1)\), and then (3) proves (1).
Overlaps between the ten partition classes only strengthen the upper
bound. Define the remaining count by

\[
 R_{6,\mathrm{unbal}}=R_6-R_{6,\mathrm{bal}}.
\]

No \(O(n^3\operatorname{polylog}n)\) estimate for this remaining
count follows here. The strictness of (11) also occurs in the target
window: at \((p,n)=(7204033,64)\), \(n\mathcal X=43776\), but
the fixed-partition opposite-free count is only \(25344\).

## Rooted coset walks and what circularity actually bounds

The independent coset calculation gives another exact way to locate
the missing estimate. Index \(H_0=H\), let \(e_0\) be its coordinate
vector, and put \(\kappa=Ce_0\). For \(x\in H_j\), the ordered
two-sum count is \(r_2(x)=\kappa_j\), while \(r_2(0)=n\).
The standard changes of variable give
\(C_{ij}=C_{ji}=C_{-i,j-i}\) and row sums \(n-\mathbf1_{i=0}\).

For fixed \(x\in H_j\), the number of \(h\in H\) with
\(x-h\in H_s\) is \(C_{-j,s-j}=C_{j,s}\), using \(z=-h/x\)
followed by inversion. Convolution thus gives

\[
 r_3(x)=n\mathbf1_{j=0}+(C^2)_{j0},\qquad
 r_3(0)=n\kappa_0.
\]

Squaring and summing over zero and all nonzero cosets proves

\[
 \boxed{E_3(H)=n^3+2n^2\sum_s\kappa_s^2+n^2\kappa_0^2
                    +n(C^4)_{00}.}                 \tag{12}
\]

Suppose now that the actual subgroup is circular, so \(C_{ij}\le2\).
Every \(\kappa_j\) is even except at the coset \(2H\), where it is
odd: interchange the two summands and account for the unique equal
pair when the target belongs to \(2H\). Thus the odd entry equals one,
all other entries are zero or two, and

\[
 \sum_j\kappa_j=n-1,\qquad \sum_j\kappa_j^2=2n-3.
\]

In the dyadic case one also has \(\kappa_0=0\). Here is a short proof
of this extra fact, rather than an assumption about small examples.
The set \(S=\{x\in H:1-x\in H\}\), of size \(\kappa_0\le2\),
is invariant under both \(x\mapsto1-x\) and \(x\mapsto1/x\).
If it contains \(-1\) or \(1/2\), these maps force all three distinct
elements \(-1,2,1/2\) into it, impossible for \(p>3\). Otherwise
both involutions have no fixed point on \(S\). If \(S\) is nonempty
it must therefore have two elements and the involutions must coincide.
This implies \(x^2-x+1=0\), hence \(x\) has order six when \(p>3\),
contradicting dyadic order. So \(S\) is empty.

Write \(v=C\kappa\). Positivity and the cell bound imply

\[
 \max_j v_j\le2(n-1),\qquad
 \sum_jv_j=n(n-1)-\kappa_0=n(n-1),
\]

and hence \(\|v\|_2^2\le2n(n-1)^2\). Substitution into (12) yields
the explicit circular-subgroup estimate

\[
 E_3(H)\le2n^4+n^3-4n^2.                             \tag{13}
\]

Circularity gives \(E_2=3n^2-3n\). The previously proved
[opposite-pair decomposition](parallel21-subgroup-next-input-2026-09-05.md)
then gives \(E_3=T_6+R_6\), where
\(T_6=15n^3-45n^2+40n\). Thus

\[
 R_6(H)\le2n^4-14n^3+41n^2-40n.                     \tag{14}
\]

These are actual bounds on the checkable circular subgroup class,
with explicit constants and no logarithm. Their exponent remains four.
Neither bounded cells nor (12) has supplied the needed third-power
estimate for the rooted fourth walk.

## Exact finite checks and an unbalanced circular example

The [new verifier](../experiments/parallel23_subgroup_upper_2026_09_05.py)
completed with status `PASS`; its full output is in the
[recorded results](../results/parallel23_subgroup_upper_2026_09_05.json).
All calculations use integer arithmetic. The verifier constructs actual
subgroups from an element of exact dyadic order, checks primality and
closure, and counts zero-sum six-tuples by normalizing the first entry
to one. It enumerates a multiset pair and a multiset triple, matches
their sums, restores their permutation multiplicities, and multiplies
by \(n\) for the free scaling coordinate. Product-balance and absence
of opposite pairs are tested separately.

It independently computes the shifted-product distribution, all 560
product-ratio correlations across the nine groups, and their agreement
with the normalized six-tuple counts. Direct triple sum/product buckets
give an additional check for five small groups. Four complete coset
matrices check (12), 397 nonzero-coset triple-count identities, and the
collision identity (2). Larger circular cases use the exact \(\mathcal
X=0\) certificate from (2), rather than constructing the full matrix.

The seven cases below are all in \(n^4/4\le p\le n^4\). Two additional
small checks, \((17,8)\) and \((97,8)\), lie outside that window and
are identified as such in the result file.

| p | n | \(\mathcal X\) | \(R_6\) | balanced at fixed split | balanced at some split | unbalanced |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1049 | 8 | 0 | 0 | 0 | 0 | 0 |
| 2017 | 8 | 0 | 0 | 0 | 0 | 0 |
| 17393 | 16 | 0 | 0 | 0 | 0 | 0 |
| 6700417 | 64 | 114 | 367680 | 7296 | 38400 | 329280 |
| 7204033 | 64 | 684 | 353280 | 25344 | 238080 | 115200 |
| 67403009 | 128 | 720 | 829440 | 55296 | 552960 | 276480 |
| 1073748737 | 256 | 0 | 368640 | 0 | 0 | 368640 |

For the last subgroup,

\[
 E_2=195840=T_4,\qquad E_3=249088000,
 \qquad T_6=248719360.
\]

An explicit normalized opposite-free zero-sum six-tuple is

\[
 (1,914267366,972187974,9468345,10646993,240926795).
\]

All entries lie in this actual subgroup; their integer sum is
\(2147497474=2p\). The verifier checks that no two entries are
opposites modulo \(p\), and that all ten partitions fail (10).
This rules out replacing the unbalanced remainder by zero under
circularity. It makes no asymptotic claim about the size of that
remainder.

## Input hashes and scope of verification

| File | SHA-256 |
| --- | --- |
| `experiments/parallel23_subgroup_upper_2026_09_05.py` | `ea2023a2cf0b78434efebb9913249c6dc74aec9247da2ad5b8d8a2000251414f` |
| `results/parallel23_subgroup_upper_2026_09_05.json` | `d675ecc47be55cde21ba176fbd5d3bfdf9c0da35af79c093079d31e47cb82770` |
| `research/parallel21-subgroup-next-input-2026-09-05.md` | `84005053f3954ee2d36eccf85b2786f1fe00d6bd956c76b23b1b40e82b991c22` |
| `research/parallel22-subgroup-independent-review-2026-09-05.md` | `3915e01c73041265eac37be84b9f9756ffc57fa72ba257bf546ab2e6b23fe693` |
| `research/mixed-periods-and-shifted-energy.md` | `de7d7202e55253955b62f73c973d2bdee42b48123eb344c875d5f0527a21f4aa` |
| `sources/mixed-periods-2026-09-04/shkredov-1504.04522.html` | `be13a230cc3cfdb780fd8dca0606d35a6f561221cb3ba67a860be6f759a1f459` |

The completed finite checks support the stated identities and examples;
they do not establish a uniform bound for the unbalanced remainder.
This pass contains a mathematical derivation, an imported analytic
theorem, and finite exact checks. It is not independent human review or
formal certification, and makes no literature novelty claim.
