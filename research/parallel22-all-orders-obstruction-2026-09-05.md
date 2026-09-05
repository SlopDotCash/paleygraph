# Twenty-second pass: all individual orders and all lower moments are insufficient in a weighted model

Status: elementary proof for a relaxed class, not a character-sum bound
or Paley counterexample. Every statement below concerns an explicitly
weighted sign measure. No identification with field translates is made.

## 1. Statement and scope

Fix an integer r≥3. For a sufficiently large prime p, set
k=⌊p^(1/(r+1))⌋ and assume k≥2(r+1). There is a finite positive
rational measure μ on {−1,0,1}^k with total mass p satisfying:

1. Every column has mean zero, zero mass one and squared norm p−1.
   The column Gram matrix is exactly pI−J.
2. For every subset D of distinct columns, its correlation is
   C_D=0 for |D| odd, C_D=−1 for |D|=2, and C_D=k² for
   even |D|≥4. Thus all individual distinct-product bounds
   |C_D|≤(|D|−1)√p hold, simultaneously for 1≤|D|≤k.
3. For every real coefficient vector a and every integer 1≤s≤r−1,

       ∫|Σ_i a_i ε_i|^(2s) dμ
          ≤[(3/2)(2s−1)!!+1] p (Σ_i a_i²)^s.          (1)

   For unsigned subset indicators, the sharper constant
   (2s−1)!!+1 works, uniformly over all subsets of the k columns.
4. For the sum of all k columns,

       ∫(Σ_i ε_i)^(2r) dμ / (p k^r) ≥ k/2 →∞.         (2)

The labels of the k columns may be chosen to form a B_h subset of
F_p for every fixed 2≤h≤⌊(r+1)/2⌋. This is only a property of the
labels: the row measure does not obey ε_c(x)=χ(x−c).

Therefore the listed exact first/second statistics, every individual
distinct-product bound, and even Gaussian-order bounds for all lower
even moments of all coefficient vectors do not imply the next
Gaussian-order moment bound in this weighted model class.

This does not disprove the SS-B* hypothesis for actual character
translates. In particular, μ is not uniform sampling over p unit rows;
each zero-column fibre contains many rows; and there is no full
p-column translation identity or multiplicative-character factorization
law. The all-coefficient assertion (1) is an extra property verified
for this model, not a claim about actual Paley matrices.

## 2. The measure

Write L=k², W=p−k−L, and a_0=(L+1)/W. For a full sign vector
ε∈{−1,1}^k put

    F(ε)=Σ_i ε_i,       e_2(ε)=Σ_(i<j) ε_i ε_j=(F²−k)/2.

The measure is the sum of three positive components:

- total mass L on the two extreme vectors ±(1,…,1), equally;
- total mass W on the full sign cube with density
  ρ(ε)=1−a_0 e_2(ε) relative to uniform cube measure;
- for each column i, mass one on the row vectors with ε_i=0 and
  all other signs independent uniform ±1.

All the weights are rational. Since E_cube e_2=0, the cube
component has precisely mass W and the whole measure has mass p.

The positivity estimate actually holds for every integer k≥2 and
p≥k⁴. First,

    W≥k⁴−k²−k>0,
    W−(k²+1)k(k−1)≥k²(k−2)≥0,
    W−(k²+1)k≥k(k−2)(k²+k+1)≥0.

Because −k/2≤e_2≤binom(k,2), these imply

    a_0 binom(k,2)≤1/2,       a_0 k≤1,
    1/2≤ρ(ε)≤3/2.                                    (3)

The factorization in the third line is exact:
k⁴−k³−k²−2k=k(k−2)(k²+k+1). The prime and size assumptions
in Section 1 ensure p≥k^(r+1)≥k⁴.

## 3. All individual correlations at once

For a distinct subset D, write ε_D=∏_(i∈D)ε_i. Uniform cube
orthogonality says E ε_D=0 unless D is empty, and

    E ε_D e_2 = 1 if |D|=2, and 0 otherwise.

Indeed each term is E ε_(D△{i,j}), which is nonzero precisely
when D={i,j}. For nonempty D, every zero-row component contributes
zero: if its zero lies in D the product is zero, and otherwise at
least one independent mean-zero sign remains.

The extreme component contributes L when |D| is even and zero
when it is odd. The cube contribution for |D|=2 is
−Wa_0=−(L+1). Consequently C_2=−1, even C_d=L for d≥4,
and odd C_d=0, exactly as claimed. Since L=k²≤√p, the ordinary
Weil-shaped individual bounds hold in every available degree.

Sign reversal makes every column mean zero. Only its own zero-row
component vanishes, of total mass one, so its squared norm is p−1.
Together with C_2=−1 this proves the Gram statement.

The comparison bound is the same classical distinct-root Weil bound
already checked in [pass 21](parallel21-squarefree-moments-2026-09-05.md).
No new external result is needed to construct or analyze the model.

## 4. All lower moments, including arbitrary real coefficients

Let G=Σ_i b_i ξ_i with independent uniform signs ξ_i. For integer s≥1,

    E G^(2s)≤D_s (Σ_i b_i²)^s,       D_s=(2s−1)!!.     (4)

To prove this directly, expand G^(2s). Only index words in which
every multiplicity is even survive. Their monomials are nonnegative,
even for real b_i. Pair the 2s positions in all D_s ways and assign
an index to each pair. Every surviving word is counted at least once.
The sum over a fixed pairing is (Σ_i b_i²)^s. This proves (4).

For a coefficient vector b on all k columns, the extreme component
has moment at most L k^s(Σ_i b_i²)^s by Cauchy–Schwarz.
The cube component is at most (3/2)W D_s(Σ_i b_i²)^s by (3)
and (4). Each zero-row component has moment at most
D_s(Σ_i b_i²)^s, so their sum is at most kD_s times that quantity.
As W+k=p−L≤p and kD_s≤(3/2)kD_s, we get

    ∫|Σ_i b_i ε_i|^(2s) dμ
      ≤[L k^s+(3/2)pD_s](Σ_i b_i²)^s.

If s≤r−1 then Lk^s=k^(s+2)≤k^(r+1)≤p, proving (1).
The constants depend on the fixed moment order; there is no growing
s estimate with an absolute constant.

### Sharper unsigned subset bound and an exact moment formula

For m≤k, let U_s(m)=E(ξ_1+…+ξ_m)^(2s), with U_0(m)=1.
Let

    V_s(m)=E (ξ_1+…+ξ_m)^(2s) e_2(ξ_1,…,ξ_m)
          =[U_(s+1)(m)−m U_s(m)]/2.

We have V_s(m)≥0. For Y=(ξ_1+…+ξ_m)² and an independent
copy Y′, monotonicity gives

    2 Cov(Y^s,Y)=E[(Y^s−Y′^s)(Y−Y′)]≥0.

Since E Y=m, this is precisely the numerator defining V_s(m).
The same conclusion follows by counting the words whose odd
multiplicity set is a specified pair.

For any fixed subset B of m columns, direct averaging gives

    M_(2s)(B)=L m^(2s)+(p−L−m)U_s(m)
                −(L+1)V_s(m)+m U_s(m−1).              (5)

The zero outside B contributes U_s(m); the zero inside B contributes
U_s(m−1). Marginalizing the cube density drops every pair not fully
inside B. Formula (5) is valid for m≥1; the empty subset is zero
at positive moments.

Since V_s(m)≥0 and U_s(m)≤D_s m^s by (4),

    M_(2s)(B)≤L m^(2s)+(p−L)D_s m^s
              ≤(D_s+1)p m^s       (1≤s≤r−1).          (6)

Thus at r=3 this model has M_4(B)≤4p|B|² for every subset B.
It does not retain the sharper constant 3 of pass 21's different
model. Formula (5) at s=1 gives exactly M_2(B)=pm−m².

## 5. The next moment diverges on the precise size slice

The extreme rows alone give M_(2r)≥L k^(2r)=k^(2r+2).
Let q=r+1. Since k=⌊p^(1/q)⌋,

    k^q≤p<(k+1)^q.

For k≥2q, the binomial expansion and binom(q,j)≤q^j give

    (1+1/k)^q≤Σ_(j=0)^q (q/k)^j<2.

Consequently p<2k^(r+1), so

    M_(2r)/(pk^r)≥k^(r+2)/p>k/2.

Every fixed r is chosen before p tends to infinity; k then tends
to infinity through its actual values for primes. No theorem about
primes in short intervals is assumed.

## 6. Relation-free labels do not repair the missing row structure

The [pass-19 extension bound](parallel19-relation-free-reduction-2026-09-05.md)
uses f_h(s)=s+Σ_(j=1)^h s^(2h−j), for p>h. If
2h≤r+1 and k≥2(r+1), then k>h+1 and

    f_h(k−1)<(h+1)k^(2h−1)<k^(2h)≤k^(r+1)≤p.

Greedy extension therefore produces a B_h k-set in F_p.
Choose it as the column labels. All row statistics remain unchanged,
and no relation between those labels and ε is created. This covers
the allowed relation orders in SS-B* without asserting a character
kernel realization.

## 7. Verification and remaining obligation

The new exact verifier will check positivity, explicit small weighted
row measures, all distinct products in those small measures, formula
(5), lower moments for varied coefficient vectors, and compressed
examples on actual prime-size slices. Its bounded checks support the
formulas, while the all-parameter assertions rely on the proof above.
Large compressed models are not actual field character evaluations.

The missing input remains an actual uniform character-sensitive upper
bound, such as the squarefree aggregate of pass 21 or the quartic-star
energy in the [separate classical review](parallel21-positive-upper-review-2026-09-05.md).
Those identities use additional full-field structure absent here.
The full classical, subgroup, spectral and official prize targets
remain unproved. Independent review and formal verification of this
new deduction are not yet complete.

