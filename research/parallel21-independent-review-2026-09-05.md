# Separate-agent proof review of pass 21

Result: I found no mathematical error in the all-parameter Krawtchouk
inequalities, the sixth-moment identity, or the weighted construction
under its stated hypothesis k≥256. The construction's limitation to
rationally weighted sign kernels is essential. No bound on the missing
positive character aggregate, Paley conjecture, or prize target follows
without an additional input.

This is a separate agent's algebraic and source-code review. It is not
independent human review, external expert certification, or formal proof
verification. I read the proof and verifier, recalculated the relevant
inequalities and combinatorial coefficients, and inspected the cited
local HTML statement. I did not rerun the verifier, independently
enumerate the finite fixtures, or audit its generated results file.

## Inputs and scope

The following SHA-256 values identify the bytes reviewed:

| Input | SHA-256 |
| --- | --- |
| [Pass-21 proof](parallel21-squarefree-moments-2026-09-05.md) | `42bfee31ccd53abaa9081128a3dfa250c5316570c407740c439b7175759f0210` |
| [Pass-21 verifier](../experiments/parallel21_squarefree_moments_2026_09_05.py) | `59094f846116bb56cd5975a2ba41986b282b3f09e6ccf44ff04d2b8fc2047269` |
| [Pass-20 transfer](parallel20-inversion-moments-2026-09-05.md) | `31ec520dcb38decced35eb111128d65383f4a4f6fb2eca545d8bc5406b55922f` |
| [Archived comparison source](../sources/mcdonald-sahay-wyman-2210.03789v2.html) | `3ae1bc2ca60ef412b271e3a03056df70fe1465fda26a89727e1929cfe9e1f622` |

The pass-20 document was read to check the scope of the imported
fixed-size transfer. The earlier chain of reductions, pass-19 extension
theorem, and full prize connection were not independently re-audited in
this bounded review. No new literature search, PDF review, or formal
proof build was performed.

## Krawtchouk comparison

The generating-polynomial differentiation gives

    H_(j+1)=F H_j−j(N−j+1)H_(j−1).

When N≥2r the positive off-diagonal squares of the stated symmetric
tridiagonal matrix give exactly this determinant recurrence. Alternating
diagonal conjugation pairs its real eigenvalues as ±λ. Every matrix
entry on an off-diagonal is bounded by √((2r−1)N), so the row-sum
bound R=2√((2r−1)N) is valid. In particular,

    H_(2r)=∏(F²−λ_i²),       |λ_i|≤R.

For F²≤R², each factor has absolute value at most R². For F²≥R²,
all factors lie between zero and F². These two cases establish the
claimed lower and absolute bounds. In the small-F case for the final
comparison, the sharper intermediate bound H_(2r)≥−R^(2r) gives

    F^(2r)−2^r H_(2r)≤2(2R²)^r≤2(16rN)^r.

For F²≥2R², each factor is at least F²/2. Thus the same comparison
holds globally. These arguments use real F in the matrix branch;
they need no distributional assumption.

For N<2r the symmetric coefficient vanishes at every actual sign sum,
and F^(2r)≤N^(2r)≤(2r)^r N^r. This also covers the small-support
cases that do not belong to the matrix branch. N=0 is treated
separately without division. Actual character rows contain either n
or n−1 nonzero entries, so summing with N≤n proves all three stated
inequalities, including empty C.

The equivalence is between two families of upper-bound hypotheses,
with constants depending on the fixed order. Absorbing the additional
p n^r term requires the stated β≥0. It does not control the positive
aggregate. The fixed-r quantifiers agree with the cited transfer;
letting r increase means choosing a sequence of fixed orders, not
silently claiming uniformity of the constants in r.

## Sixth-moment ledger

Independent multiplicity counting gives the displayed coefficients.
For a fixed pair of odd-multiplicity indices, the contributions are

    12+20+120(n−2)+30(n−2)+180 binom(n−2,2)
      =90n²−300n+272.

For four such indices the coefficient is 480+360(n−4)=360n−960.
The constant contribution is n+15n(n−1)+90 binom(n,3), giving
15n³−30n²+16n. Six distinct indices contribute 720.

Replacing one zero by an equiprobable sign leaves every multilinear
term unchanged in expectation and adds exactly 15F⁴+15F²+1 to F⁶.
Consequently the corrections and their signs in equation (5) are
correct. At n=⌊p^(1/4)⌋, an individual quartic bound gives
|a_4 T_4|=O(√p n^5)=O(p n³); the corrections are favorable for
an upper bound. This calculation does not supply the needed T_6 upper
bound.

## Weighted-model constants and quantifiers

All model assertions below are checked at the stated k≥256, even
where the proof records weaker sufficient thresholds.

- W>0 and k−2≤v<k follow from the displayed lower bound for W and
  upper bound for A. Their subtraction gives exactly k⁴−k³−7k>0.
  Parity-compatible t exists, b=t+2≤√k+2≤k, and the two mixing
  probabilities are nonnegative rational numbers. Combining duplicate
  zero-sum classes preserves their total mass.
- The measure's total mass is L+W+k=p. Every coordinate has zero
  mass one. Permutation symmetry, total second moment pk−k², and
  diagonal norm p−1 imply off-diagonal Gram entry −1. Sign reversal
  cancels every odd distinct product.
- The correlation formula includes the correct zero-row factor k−d.
  The d=k case has no zero-row term and must bypass the denominator
  (k−1)_k=0; the verifier does so. For the large-model estimates
  d≤6<k, all denominators are positive.
- The two-point variance is exactly (v−t²)(b²−v). The bounds
  4(t+1)²≤8k and F²≤3k give E F⁴≤k²+8k and E F⁶≤4k³.
  The final fourth-moment inequality has surplus
  p k(k−8)−k≥0 for k≥9, so the displayed 3p k² bound is valid.
- For H_4 the lower numerator is −2k²+2k and the upper numerator
  is −2k²+18k−12. The bulk term is between −(5/2)√p and zero.
  The nonnegative zero term is below one. Thus C_4 is between
  −(3/2)√p−1 and √p+1, which lies inside [−2√p,2√p].
- In H_6, all polynomial coefficient signs used to bound absolute
  values are valid for k≥256. The bulk bound is
  79k³+120k²≤100k³; the zero numerator is also ≤100k³ in absolute
  value. The product estimate ∏(1−j/k)≥1−Σj/k proves the four
  falling-factorial lower bounds at their stated thresholds. Hence
  |C_6|≤√p+400k+1≤3√p+1≤5√p follows.
- L≥k² and p<(k+1)⁴≤2k⁴ give M_6≥k⁸ and
  M_6/(p k³)≥k/2. The construction uses each prime's own k, and
  therefore does not require primes between consecutive fourth powers.
  The verifier's five chosen primes are finite fixtures only.

The no-carry argument for the explicit Sidon label fixtures is valid:
t_i+t_j<2q, reduction modulo 2q recovers this integer sum, and
the square-residue digits then recover the square sum modulo the odd
prime q. Sum and square sum determine the unordered pair modulo q.
The checked bound 2max c_t<p prevents wraparound. This proves the
fixture argument under its stated q and p conditions; it is not an
independent all-prime proof of the imported greedy existence result.

The cited local HTML's Lemma 2.1 supplies the comparison bound
(d−1)√p for the quadratic character of a product of d distinct
linear factors. Such a product is not a square. No use of the source's
VC-dimension results is needed.

## Verifier coverage and limits

The code has several useful independent representations: binomial
expansion versus recurrence; direct small-field product characters
versus row polynomials; and explicit small weighted row enumeration
versus compressed correlation formulas. I found no discrepancy between
the code's formulas and the proof. Exact `Fraction` comparisons and
integer primality cutoffs avoid floating-point acceptance.

The finite pointwise checks test only parity-compatible sign sums,
although the matrix argument itself covers real F in its branch.
Actual field checks cover the listed primes only, sets of size at most
six, and moments through order eight. Large fixtures use compressed
model formulas, not actual field characters. The constant ledger
checks only k=256,…,4096; the all-k claims rely on the algebra above.
No finite test coverage establishes an asymptotic estimate.

The construction does not retain uniform sampling over p rows, a
single unit row in every zero fibre, actual translation/factorization
relations, or a full p-column Paley matrix. The fourth-moment assertion
as written concerns the sum of the model's k columns; the review does
not upgrade it to a uniform hypothesis over other cardinalities or
arbitrary coefficient vectors. Sidon column labels add an additive
property of the labels but impose no relation on row signs.

Consequently this is a valid counterexample to a sixth-moment
deduction from the explicitly listed constraints in the larger weighted
model class. It neither refutes an actual prime-field moment conjecture
nor proves that those constraints are insufficient after extra Paley
identities are imposed. The existing scope paragraphs state this
distinction adequately. No mathematical correction is requested.
