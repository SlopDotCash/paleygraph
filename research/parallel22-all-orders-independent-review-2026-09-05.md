# Independent review of the all-orders weighted obstruction

**Result: the mathematical claims in the reviewed version pass this
independent audit. No mathematical correction is requested.** This is
an elementary obstruction within a weighted sign-model class. It does
not supply a character-sum bound, a counterexample to SS-B*, or a Paley
proof. This review is not a formal proof-assistant certification.

Reviewed source: [parallel22-all-orders-obstruction-2026-09-05.md](parallel22-all-orders-obstruction-2026-09-05.md).

Reviewed SHA-256:

    b959300e5c9c886e2936c306a2b0846f01d47e77b795bf1c36b5325877ee926c

The snapshot was read in full. A subsequent byte comparison against the
working source found no changes before this review was written. The
reviewed version's last section still describes the root verifier and
independent review as unfinished; those are workflow-status statements,
not defects in the mathematical argument. A later source revision must
be compared with this pinned version before extending this verdict.

## Positivity and total mass

The construction has L=k², W=p−k−k², and full-cube density
ρ=1−((k²+1)/W)e₂. For k≥2 and p≥k⁴, W>0. The two exact
lower bounds used in the proof are correct:

    W−(k²+1)k(k−1) ≥ k³−2k²,
    W−(k²+1)k ≥ k⁴−k³−k²−2k
                         =k(k−2)(k²+k+1).

The first is deliberately stronger than the inequality needed for
nonnegativity. After division by W it gives
((k²+1)/W)binom(k,2)≤1/2. The second gives
((k²+1)/W)k≤1. Combining these with −k/2≤e₂≤binom(k,2)
proves 1/2≤ρ≤3/2, including k=2. In particular, the cube component
is positive; its negative degree-two perturbation is not an unnoticed
signed measure. Its average density is one because E_cube e₂=0.
The component masses add to W+L+k=p.

## Exact correlations and the comparison in every order

The distinct-coordinate products are the Walsh basis of the full sign
cube. The perturbing polynomial e₂ contains only basis elements of
cardinality two. Hence its contribution to a nonempty D is exactly
−(L+1) when |D|=2 and zero otherwise. This computation is valid for
all |D|≤k; there is no truncation at a fixed degree.

For each zero-row component, either D contains its zero coordinate or
D contains at least one remaining independent mean-zero sign. Its
contribution is therefore zero in both cases. The two extreme rows
contribute L precisely for even |D|. This gives odd C_D=0,
C_D=−1 at degree two, and C_D=k² at every even degree at least four.
Since p≥k⁴, these satisfy the stated Weil-shaped comparisons in all
available degrees. At degree one the right side is zero and the mean
really is zero. The empty product is separately the total mass p.

Exactly one zero-row fibre, of mass one, vanishes in any specified
column. This gives squared norm p−1. With the off-diagonal correlation
−1, the Gram matrix is pI−J. In particular the exact arbitrary-vector
second moment is

    ∫(Σ_i b_i ε_i)²dμ = pΣ_i b_i²−(Σ_i b_i)²,

consistent with, and stronger than, the stated s=1 upper bound.

## Arbitrary real coefficients and lower even moments

The pairing argument for the independent cube is valid for real,
possibly negative, coefficients. After averaging, only words with
even index multiplicities survive. Every surviving monomial is
nonnegative; a word with multiplicities 2m_i admits
∏_i(2m_i−1)!!≥1 pairings. Summing over all (2s−1)!! pairings therefore
bounds the cube moment by (2s−1)!!(Σ_i b_i²)^s without an illicit
absolute-value or sign assumption.

For the extreme component, Cauchy–Schwarz gives
L|Σ_i b_i|^(2s)≤Lk^s(Σ_i b_i²)^s. For the cube and zero-row components,
ρ≤3/2 and the pairing bound give the claimed upper contributions.
Their combined coefficient is at most (3/2)p(2s−1)!! because
W+k=p−L≤p. Finally

    Lk^s=k^(s+2)≤k^(r+1)≤p  for 1≤s≤r−1.

Thus the displayed all-coefficient bound holds simultaneously over all
real coefficient vectors, with a constant depending only on the fixed
moment order. The proof does not assume that the coefficients are
indicators, nonnegative, normalized, or drawn from a finite list.

## Unsigned formula and its sharper constant

The formula

    M_(2s)(B)=Lm^(2s)+(p−L−m)U_s(m)
                      −(L+1)V_s(m)+mU_s(m−1)

correctly separates zeros inside and outside B. Marginalizing e₂ drops
all pairs not contained in B, including pairs with exactly one endpoint
outside B. The coefficient p−L−m is nonnegative because m≤k and
p−L−k=W>0.

The covariance argument establishes V_s(m)≥0: with Y=(Σξ_i)²,
its defining numerator is Cov(Y^s,Y), and the independent-copy identity
has a nonnegative integrand. Dropping the negative V term is therefore
valid for an upper bound. Applying the pairing bound to U_s(m) and
U_s(m−1) yields the claimed coefficient p−L. The extreme contribution
is at most pm^s because Lm^s≤Lk^s≤p. These steps prove the sharper
constant (2s−1)!!+1 uniformly over all subsets. The empty subset is
handled separately at positive moments; no use of U_s(−1) is needed.
At s=1 the formula simplifies to pm−m² exactly.

## Prime-size and B_h quantifiers

For every fixed r≥3, the actual cardinality k=⌊p^(1/(r+1))⌋ tends
to infinity as the prime p tends to infinity. No requirement that a
prime occur for every preselected k is made. The condition k≥2(r+1)
implies (k+1)^(r+1)<2k^(r+1) by the finite geometric-series estimate in
the proof. Positive extreme mass alone then gives a next-moment ratio
strictly greater than k/2. A lower bound of k/2, as in the statement,
is consequently valid. Positivity also follows on this slice because
r+1≥4. The threshold ensures all degree-2r products are available,
although the next moment does not rely on this observation.

For every allowed h, p>h and k>h+1. The cited extension function obeys

    f_h(k−1)<(h+1)k^(2h−1)<k^(2h)≤k^(r+1)≤p.

This guarantees the successive extensions needed to reach k labels.
If a simultaneous labeling for all allowed h is desired, one may choose
h=⌊(r+1)/2⌋: a nonempty B_h set is B_j for every j≤h, since any
j-term relation can be padded on both sides by h−j copies of a fixed
label. This is a clarification, not a correction to the stated
existence claims. Labels do not impose a relation between those field
elements and the row signs.

## Independent finite checks

A separate [exact audit script](../experiments/parallel22_all_orders_independent_audit_2026_09_05.py)
was written without importing the root verifier. It represents every
weight with common denominator 2^k. The full-cube numerator is
W−(L+1)e₂; each zero-row numerator is 2; the extreme rows receive the
additional numerator L·2^(k−1). All evaluations use Python integers or
exact rational numbers.

The [saved results](../results/parallel22_all_orders_independent_audit_2026_09_05.json)
contain:

- 297 positivity parameter cases and 15,444 sign-sum density checks;
- 9 explicitly enumerated weighted models, containing 37,116 rows;
- all 501 nonempty distinct-subset correlations for the models k=2,…,8;
- 57 individual zero-mass checks;
- 105 lower-moment checks for varied real-coefficient specializations;
- 148 direct checks of the exact unsigned moment formula and its bound;
- 3 next-moment growth checks at admissible thresholds;
- 84 size-ratio checks and 630 B_h-label threshold checks.

The larger explicitly enumerated fixtures use the integer parameters
(r,k,p)=(3,8,4096),(4,10,100000),(5,12,2985984). These are deliberately
integer-parameter checks of the stronger positivity/moment construction,
not claims that those p values are prime. The prime-size quantifiers
were audited through the general argument above. No large-field
character kernel was evaluated in this independent model audit.

These finite checks supplement the algebraic review. They do not
establish statements for every real vector or every prime by themselves;
those assertions rely on the reviewed proofs. No prior artifact or
central project document was edited by this reviewer.
