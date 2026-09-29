# Twenty-first pass: the top squarefree aggregate and a weighted obstruction

Status: elementary deductions with exact finite checks; independent review
is outstanding. No uniform character-moment bound or Paley proof is claimed.
The weighted construction below is not a prime-field character kernel.

## 1. Only a one-sided top aggregate is needed

Let p be an odd prime, χ its quadratic character with χ(0)=0,
C⊂F_p, and n=|C|. Put

    F_C(x)=Σ_(c∈C) χ(x−c),       M_(2r)(C)=Σ_x F_C(x)^(2r),
    e_j(x)=Σ_(Q⊂C, |Q|=j) ∏_(c∈Q) χ(x−c),
    T_j(C)=Σ_x e_j(x)=Σ_(Q⊂C, |Q|=j) K(Q),
    K(Q)=Σ_x χ(∏_(c∈Q)(x−c)).

An empty sum over j-subsets is zero. For every integer r≥1,

    T_(2r) ≥ −(8r)^r p n^r/(2r)!,                         (1)
    |T_(2r)| ≤ [M_(2r)+(8r)^r p n^r]/(2r)!,              (2)
    M_(2r) ≤ 2^r(2r)! T_(2r)+2(16r)^r p n^r.             (3)

Consequently, at fixed r and β≥0, a uniform one-sided upper estimate
T_(2r)≤D_r p^(1+β)n^r is equivalent, up to constants depending on r,
to the corresponding moment estimate M_(2r)≤C_r p^(1+β)n^r.
The lower side already has the required order. This statement needs
no additive-relation hypothesis and works on any prescribed size slice.

### Proof, including the zero character values

For a row ε_i∈{−1,0,1}, let N count its nonzero entries, and let
F=Σ_i ε_i. Its elementary-symmetric generating polynomial is

    E(t)=∏_i(1+tε_i)=(1+t)^((N+F)/2)(1−t)^((N−F)/2).

Differentiation gives (1−t²)E′=(F−Nt)E. Therefore H_j=j!e_j satisfies

    H_0=1, H_1=F,
    H_(j+1)=F H_j−j(N−j+1)H_(j−1).                       (4)

For N≥2r, H_(2r)(F;N) is the characteristic polynomial of the
real symmetric tridiagonal matrix of size 2r with zero diagonal and
off-diagonals √(j(N−j+1)), j=1,…,2r−1. The determinant recurrence
is exactly (4). Conjugating by diag(1,−1,1,−1,…) negates the matrix,
so its real eigenvalues occur as r pairs ±λ_i. Its absolute row sums
are at most R=2√((2r−1)N). Thus |λ_i|≤R and R²≤8rN, and

    H_(2r)(F;N)=∏_(i=1)^r (F²−λ_i²).

If F²≥R², the product is nonnegative and at most F^(2r).
Otherwise every factor has absolute value at most R². Hence

    H_(2r)≥−(8rN)^r,       |H_(2r)|≤F^(2r)+(8rN)^r.

If F²≥2R², every factor is at least F²/2, so F^(2r)≤2^r H_(2r).
If F²<2R², the preceding lower bound gives
F^(2r)−2^r H_(2r)≤2(2R²)^r. Both cases imply

    F^(2r)≤2^r H_(2r)+2(16rN)^r.

For N<2r, the actual elementary symmetric coefficient is zero and
|F|≤N, so F^(2r)≤(2r)^r N^r; all three inequalities still hold.
For N=0 they hold directly. In a character row N=n or n−1.
Sum the inequalities, use N≤n, and divide H_(2r) by (2r)! to obtain
(1)–(3). This argument does not assume that every row has n signs.

### Consequence for the previous restricted criterion

The [twentieth pass](parallel20-inversion-moments-2026-09-05.md)
shows, for fixed r≥3 and 2≤h≤⌊(r+1)/2⌋, that a uniform moment
estimate on B_h sets of size k=⌊p^(1/(r+1))⌋ transfers to all sets
of that same size, for sufficiently large primes. B_h means that
equal h-term sums, with repetitions allowed, have identical multisets.

In its sufficient criterion SS-B*, one can replace the assumed moment
estimate by just the one-sided top aggregate estimate

    Σ_(Q⊂C, |Q|=2r) K(Q) ≤ D_r p^(1+β) k^r             (SF)

on the same B_h sets. A sequence of fixed orders r→∞ with β→0
would imply the full classical conjecture through the unchanged
single-size sampling transfer. No such sequence of estimates is proved.
The equivalence above is between moment and aggregate hypotheses;
it is not a converse from Paley to SS-B*.

### Exact sixth-moment ledger

The recurrence also gives

    H_4(X;N)=X⁴−(6N−8)X²+3N(N−2),
    H_6(X;N)=X⁶−(15N−40)X⁴
              +(45N²−210N+184)X²−15N(N−2)(N−4).

For n fully nonzero signs, reduction modulo ε_i²=1 gives

    F⁶=a_0+a_2 e_2+a_4 e_4+720 e_6,
    a_0=15n³−30n²+16n,
    a_2=90n²−300n+272,       a_4=360n−960.

These coefficients follow either by classifying multiplicities of six
indices, or by eliminating H_2,H_4,H_6 using (4). At a row with
one zero, average over replacing that zero by ±1. The right side
does not change, while the left side increases by
15F⁴+15F²+1. Since Σ_x χ(x−c)χ(x−d)=−1 for distinct c,d,

    M_6=p a_0−a_2 binom(n,2)+a_4 T_4+720 T_6
          −15Σ_(c∈C)F_C(c)⁴−15Σ_(c∈C)F_C(c)²−n.     (5)

In particular, at n=⌊p^(1/4)⌋ the ordinary individual degree-four
Weil estimate bounds the a_4 T_4 term by O(p n³). The exceptional
row corrections in (5) are favorable for an upper bound. The
uncontrolled positive aggregate at order six is the substantive term.

## 2. A weighted model passes low-order tests and fails the next moment

For every sufficiently large prime p, put k=⌊p^(1/4)⌋≥256. There
is a finite rationally weighted measure on {−1,0,1}^k, of total
mass p, for which:

- each coordinate has mean zero, zero mass 1, and squared norm p−1;
- the column Gram matrix is exactly pI−J;
- all distinct odd-product sums vanish, and the degree-two sum is −1;
- each distinct degree-four sum has absolute value at most 2√p;
- each distinct degree-six sum has absolute value at most 5√p;
- its fourth moment is at most 3p k²;
- its sixth moment divided by p k³ is at least k/2, tending to infinity.

Thus exact first/second moments, a Gaussian fourth-moment bound, and
the usual individual Weil-shaped bounds through degree six do not,
by themselves, imply a Gaussian sixth-moment bound on this size slice
in the class of weighted sign kernels. The column labels can also
be a genuine additive Sidon subset of F_p.

**Scope:** rows carry rational weights, not p equal unit weights;
each zero-column fibre has mass one but can contain many rows. There
is no claimed identification ε_c(x)=χ(x−c), no asserted field-translation
law, and no full p-column Paley identity. Labeling the columns by a
Sidon set does not supply any of these missing properties. This is
an obstruction to deductions from the listed relaxed constraints,
not a counterexample to a prime-field moment conjecture or Paley.

### Construction and exact second moment

Let P=p−k, L=⌊√p⌋, W=P−L, and e=(k−1) mod 2∈{0,1}.
We have k²≤L≤k²+2k. Split the row mass into:

1. Extreme rows, total mass L: all signs +1 or all −1, equally.
2. Bulk rows, total mass W: no zeros, with the sum distribution below.
3. Zero rows, total mass k: choose the zero coordinate uniformly;
   choose the remaining sign sum +e or −e equally, then choose a
   sign vector with that sum uniformly.

In every sum class, use the uniform distribution on its sign vectors.
If the sum is zero the two copies are simply combined. Define

    v=(Pk−ke²−Lk²)/W
      = k−[Lk(k−1)+ke²]/W.

Then k−2≤v<k. Indeed W≥k⁴−k²−3k>0 and, with
A=Lk(k−1)+ke²≤k⁴+k³−2k²+k,

    2W−A ≥ k⁴−k³−7k >0       (k≥3).

Let t be the largest nonnegative integer of parity k with t²≤v,
let b=t+2, and λ=(v−t²)/(b²−t²). In the bulk choose the sums
±t with total conditional probability 1−λ and ±b with total
conditional probability λ. Their signs are equiprobable. These are
available sign sums, with b≤k, and 0≤λ≤1. All weights are rational.

Permutation and global sign symmetries are built into the measure.
Every coordinate has zero mass one. The total second moment is

    M_2=Lk²+Wv+ke²=Pk=pk−k².

If g is any off-diagonal column inner product, symmetry and the
preceding identity give k(p−1)+k(k−1)g=pk−k², hence g=−1.
Global sign reversal gives zero means and zero odd-product sums.

### Formula for every distinct-product correlation

For N nonzero coordinates with specified sum F, let

    κ_d(F,N)=[z^d](1+z)^((N+F)/2)(1−z)^((N−F)/2)/binom(N,d)
            = H_d(F;N)/(N)_d,                             (6)

where (N)_d=N(N−1)…(N−d+1), and d≤N. Averaging over all
sign vectors with that sum proves (6) by permutation symmetry.
For an even d≤6, the total correlation of any specified d-set is

    C_d=L+W[(1−λ)κ_d(t,k)+λκ_d(b,k)]
                +(k−d)κ_d(e,k−1).                       (7)

The last coefficient is k times the probability that the zero is
outside the chosen d-set. If d=k there is no such contribution.
For d=2 the formula gives exactly −1, consistently with the Gram
argument. Formula (6) is used only where its denominator is nonzero.

### Bulk moment estimates

Let Y=F² in a bulk row. It has mean v and

    Var(Y)=(v−t²)(b²−v)≤(b²−t²)²/4
          =4(t+1)²≤8k.

Here t≤√k and k≥6 justify the last inequality. It follows that

    (k−2)²≤E(F⁴)≤k²+8k.

For k≥8, Y≤(√k+2)²≤3k, and for k≥24 this yields

    E(F⁶)≤3k(k²+8k)≤4k³.

The total fourth moment therefore satisfies, for k≥9,

    M_4=Lk⁴+W E(F⁴)+ke⁴
        ≤pk²+p(k²+8k)+k≤3pk².                            (8)

The first term is bounded by pk² because L≤√p and k²≤√p.

### Sharp enough degree-four bound

The quadratic v²−(6k−8)v decreases for v∈[k−2,k]. Thus

    −2k²+2k ≤ E H_4(F;k) ≤ −2k²+18k−12 <0.

Its absolute value is at most 2k². Also (k)_4≥k⁴(1−6/k) and
√p<(k+1)². For k≥39,

    (1+1/k)²/(1−6/k)≤5/4,

equivalently k²−38k−4≥0. The bulk term in C_4 therefore lies
between −(5/2)√p and zero. With N=k−1, the zero-row numerator
H_4(e;N) equals 3N(N−2) if e=0 and 3(N−1)(N−3) if e=1.
It is nonnegative and at most 3k². Since (k−1)_4≥k⁴/2 for
k≥20, its contribution is at most 6/k<1. Finally
√p−1≤L≤√p. Combining these estimates gives

    −2√p≤C_4≤2√p.                                       (9)

### Degree-six bound and sixth-moment obstruction

The displayed H_6 polynomial, the bulk moment bounds, and v≤k give

    |E H_6(F;k)|
       ≤4k³+15k(k²+8k)+45k²·k+15k³
       =79k³+120k²≤100k³.

For the zero rows, |H_6(e;k−1)|≤1+15k+45k²+15k³≤100k³.
The factorial products satisfy (k)_6≥k⁶/2 for k≥30 and
(k−1)_6≥k⁶/2 for k≥42. Formula (7) now gives

    |C_6|≤√p+200p/k³+200/k²
          ≤√p+400k+1≤3√p+1≤5√p,                        (10)

using p<(k+1)⁴≤2k⁴, k≥256, and √p≥k². Together with
C_2=−1 and the vanishing odd sums, this verifies all the standard
individual bounds |C_d|≤(d−1)√p through d=6.

But the extreme rows alone contribute

    M_6≥Lk⁶≥k⁸,        p k³<2k⁷,
    M_6/(p k³)≥k/2 →∞.                                  (11)

This construction works for every prime with k≥256; it assumes no
theorem about primes in short intervals. As p ranges to infinity
through primes, k also tends to infinity.

### Additive labels and what they do not encode

The [nineteenth-pass](parallel19-relation-free-reduction-2026-09-05.md)
greedy extension bound is f_2(s)=s³+s²+s. Since f_2(k−1)<k³<p,
a Sidon k-set exists in F_p. Relabeling the k columns by this set
does not change any statistic above. It does not connect the row
signs to field differences. Therefore the construction does not
disprove the Sidon character estimate at the r=3 SS-B* slice.

For finite reproducible label fixtures, one can use an odd prime
q≥k and c_t=t+2q(t² mod q), 0≤t<k, when 2max c_t<p. Equality
of two integer pair sums first gives equality of t-sums by reduction
modulo 2q (no carry), and then equality of their square sums modulo q.
In odd characteristic their sums and products determine the unordered
pair. The bound on 2max c_t precludes modular wraparound. The verifier
also checks every unordered pair sum directly; no prime-gap claim
is inferred from its selected q or p values.

## 3. Sources, verification and the remaining obligation

The comparison with actual character kernels uses only the standard
Weil bound for a product of d distinct linear factors. Its formulation
is checked in Lemma 2.1 of
[McDonald–Sahay–Wyman, arXiv:2210.03789v2](https://arxiv.org/html/2210.03789v2).
That is a classical input, not a new result. No theorem about VC
dimension is imported. The polynomial reduction and weighted model
are derived above using elementary algebra and a symmetric matrix.
No new PDF was used or visually reviewed.

The [exact verifier](../experiments/parallel21_squarefree_moments_2026_09_05.py)
compares the recurrence with binomial expansion, checks the inequalities,
independently evaluates actual small-field character products and (5),
and explicitly enumerates small rationally weighted sign measures.
Large fixtures use exact compressed formulas, exact primality checks,
Sidon pair-sum checks, and rational comparisons. They do not evaluate
actual character kernels over their large fields. The prose arguments
cover arbitrary allowed parameters; finite experiments do not establish
the missing asymptotic estimate. No formal or independent proof review
has been obtained for this pass.

The next sufficient input is the positive upper estimate (SF) on
actual character translates. The listed low-order constraints leave
room for spikes in a larger model class. A successful argument must
use additional structure, such as the actual translation and
factorization relations. All full classical, subgroup, spectral and
official prize targets remain unproved.

