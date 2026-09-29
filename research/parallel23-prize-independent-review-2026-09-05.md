# Independent review of the pinned prize bridge

**Verdict: no mathematical correction is requested for the reviewed
version.** The restricted root/list formula, conditional phase estimate,
pigeonhole obstruction, and remainder-count interpretation check out.
The additional exact-pinned ArkLib definitions confirm that the final
counts describe the archived Lambda and IsMCA quantities. No bound for
the production-size maxima, prize certificate, or Paley theorem follows.
This review is neither a Lean build nor formal or human referee approval.

## Pinned inputs and source scope

The [bridge note](parallel23-prize-bridge-2026-09-05.md), its
[verifier](../experiments/parallel23_prize_bridge_2026_09_05.py), and
[recorded results](../results/parallel23_prize_bridge_2026_09_05.json)
were read directly. The seven input hashes in the recorded results
match the current files. The verifier was not rerun wholesale.

| Input | SHA-256 |
|---|---|
| research/parallel23-prize-bridge-2026-09-05.md | ab9ccb9ffca81ce34411d2ddc9b86e18a01c4be142787aa98185226e85a20b6b |
| experiments/parallel23_prize_bridge_2026_09_05.py | 3c2d820fdd87c6c9760848e4d65c43dd1f47f370c62a7bb3478ec1f65951e24a |
| results/parallel23_prize_bridge_2026_09_05.json | ab5f8bbb4a94574f33a8ee155a5d3a8b3a010cb3b17acd9ead5782a131d7f641 |
| sources/parallel23-prize-dependencies/ListDecodability.lean | d2b9f145b82056b26ee46c8c42c1463936379c33a3fa31802decd4eb85711912 |
| sources/parallel23-prize-dependencies/ProximityGap_ProximityGenerators.lean | af2bcba4034bee9f33042f5378d0f0a0315a46cabd2de14a7d6357d4405187a7 |
| sources/official-prize-2026-09-04/ProximityPrize/Benchmark/IRSProfile.lean | 75d40b94814bdc0ad00b107a2cf6a995c79f45a4b8c1ec1ce2f18ec612e81126 |
| sources/official-prize-2026-09-04/ProximityPrize/Benchmark/TargetLower.lean | 7d11762452db502cc244e2522eb83f01d3d4b8318ec1771d0c86041181bdc474 |
| sources/official-prize-2026-09-04/dependencies/ArkLib/ProofSystem/ToyProblem/Impl/IRS.lean | 9a06a01f2cd6529208fb91ee87777467d1a2b45364a6da6bc7ee12ea7e4f6357 |
| sources/official-prize-2026-09-04/dependencies/CompPoly/Fields/KoalaBear/Basic.lean | a9f5b2ba2348a53dc15f6c942f63929d01ebdec4cda088dd4d93c3d8bc19f32f |
| sources/official-prize-2026-09-04/dependencies/CompPoly/Fields/KoalaBear/Ext6.lean | adba31e20d793d6316aab2589a7f88de91e4860e2d6d358237eb911076090b6e |
| results/parallel23_prize_additional_sources_2026_09_05.json | 7cb13a0866ac3dc65d23757dffb2afde510d5c0d66df1f3a489dfc2f903b6d39 |

The [additional-source record](../results/parallel23_prize_additional_sources_2026_09_05.json)
identifies ArkLib commit e65197892890b8fd9b0dc05b8980273cf1d595cc,
matching the original archived dependency. All four files in that
record were checked against their byte lengths and SHA-256 values.
The semantic review uses ListDecodability and ProximityGenerators,
with the archived IRS implementation and benchmark profile. It does
not substitute current-branch definitions or claim an audit of every
transitive Lean dependency. No external retrieval was performed by
this reviewer, and no previous artifact was edited.

## Exact root sets and coefficient conditions

With 1≤k<s=k+a≤n<p and W monic of degree s, any degree-<k
polynomial F agreeing with W on at least s domain points has W−F
monic of degree s. The root bound forces exactly s distinct roots,
so W−F is precisely the monic polynomial Q_T of that size-s root
set. Conversely matching the coefficients of degrees k through s−1
makes W−Q_T have degree less than k. The leading degree-s terms
cancel automatically. Different T give different polynomials F and
different evaluation codewords because k≤n and the domain has n
distinct points. This proves the exact list equality, not merely
an upper bound counting explanations.

Because T lies in F_p, every coefficient of Q_T is in F_p. Thus a
top coefficient of W outside F_p really makes the list empty.
Extension-valued coefficients below degree k are unrestricted and
are absorbed by F. The two F_25 examples in the verifier respect
this distinction; the representation uses θ²=2, with 2 a nonsquare
in F_5.

For Q_T=X^s+c_1X^(s−1)+…, Newton's identities read

    v_j+c_1v_(j−1)+…+c_(j−1)v_1+jc_j=0.

This verifies the sign and indexing of the displayed recursion.
The hypotheses imply a=s−k≤n−1<p, so every j=1,…,a is invertible
in F_p. Therefore the triangular map from the first a coefficients
to the first a power sums is bijective. Without a<p, the reverse
implication could fail; here it is supplied by the stated size
hypotheses and has not been omitted.

## The generating-function estimate is correctly conditional

For b≠0, write w_x=e_p(P_b(x)). The logarithm of
∏_x(1+zw_x) has coefficient (−1)^(j−1)Σ_xw_x^j/j at degree j.
For j≤s<p, the coefficient vector of jP_b is nonzero, so each such
power sum is bounded by M_a(D). Coefficientwise absolute-value
majorization of the exponential therefore gives

    |[z^s]∏_x(1+zw_x)|
       ≤[z^s]exp(UΣ_(j≥1)z^j/j)
       =binom(U+s−1,s).

Only logarithm terms through degree s matter. Complementing subsets
replaces the coefficient by the product of all w_x, of modulus one,
times the (n−s)-coefficient for w_x^(-1); these are the phases of
−P_b. Thus replacing s by t=min(s,n−s) is valid. Fourier inversion
has principal coefficient binom(n,s)/p^a and exactly p^a−1 other
frequencies, proving the stated prefactor and error term.

The argument does not allow M_1 to replace M_a for a>1. The explicit
F_17 example checks the proposed numerical substitution exactly:
54·289−12870=2736 exceeds 288. This refutes that substitution, while
leaving a conjectured linear-period estimate untouched.

The verifier's additional finite count bounds are direct checks of
those inequalities, not an independent computation of every M_a.
For its particular full-domain quadratic example in F_13, the usual
quadratic Gauss bound gives M_2≤√13+1<5, consistent with U=5.
The general theorem remains conditional on the displayed M_a input;
no production-size bound for that input is asserted.

## Polynomial-phase pigeonholing and official arithmetic

Evaluation of zero-constant polynomials of degree at most a is
injective on n distinct points when a<n: a nonzero difference
polynomial cannot have n roots. The p^a coefficient vectors therefore
produce p^a distinct evaluation vectors. Partitioning each coordinate
into R intervals gives at most R^n boxes. Under R^n<p^a, two vectors
in one box differ by an actual nonzero polynomial with zero constant
coefficient, and every coordinate difference has an integer
representative d_x satisfying |d_x|<p/R.

The elementary inequality cos θ≥1−θ²/2, followed by π²<10,
yields the strict bound 1−20/R² on each real exponential. The
coefficient vector is not a prohibited nonzero constant phase; its
constant coefficient remains zero after subtraction.

If a/n tends to a positive constant and p tends to infinity, the
choice R=⌊p^(a/(2n))⌋ tends to infinity and satisfies R^n<p^a.
Combining the lower estimate with M_a≤n proves M_a/n→1. To contrast
this with O(√(n log p)), one also needs log p=o(n); the stated
quartic subgroup window has exactly that property. The note's
application in that window is therefore valid. This does not provide
a large list: Fourier inversion sums many complex coefficients whose
aggregate can still cancel.

The archived field parameters agree with the arithmetic used:
KoalaBear.fieldSize is 2^31−2^24+1=2130706433, and the archived
Ext6 file states card_ext6=fieldSize^6. The profile has total message
dimension 2^20 and eight rows, giving scalar polynomial dimension
k=2^17. It has n=2^18 domain columns. Doubling the interleaving gives
sixteen scalar rows with the **same** degree bound k, by reindexing
Fin 2 × Fin 8. It does not mean using total dimension 2^20 with
sixteen rows and thereby halving the scalar degree bound.

The numerical assertions were independently checked with integers:

    floor(2130706433^6 / 2^128)=274980728111395087,
    30·26215−3·262144=18,
    131072+26215=157287,
    262144−157287=104857,
    (11/16)·262144=180224.

Since p>2^30, the middle exponent inequality really proves
p^26215>8^262144. The phase lower bound at this fixed radius is
therefore rigorous even though no production-size polynomial witness
is constructed. It does not assert the same numerical lower bound
at other radii.

For 1≤j≤131071, j/n lies strictly between zero and the archived
minimum relative distance 131073/262144, so it is admissible. The
range also ensures s=n−j≥k+1 for the witness-shrinking argument
below. The separate spot-check score clause of ProtocolClaim is
still required; the combined-error budget alone does not establish
that clause or a prize-winning score.

## Remainder lists agree with the archived Lambda

[ListDecodability, line 139](../sources/parallel23-prize-dependencies/ListDecodability.lean)
defines Lambda as the supremum over received words of the cardinality
of closeCodewordsRel. The latter is the intersection of the code
with the relative Hamming ball. At radius j/n this is exactly
agreement in at least s=n−j columns, with agreement of every scalar
row in each accepted column.

Every received word has a unique degree-<n polynomial interpolant
in each scalar row. On a size-s subset T, agreement with a degree-<k
polynomial is equivalent to the degree-<s remainder modulo Q_T having
no coefficients in degrees k,…,s−1. If that condition holds, the
remainder itself is the unique explanation. Taking a **set of distinct
remainder tuples** therefore gives exactly the point list. Counting
root subsets instead would overcount explanations whose agreement
set is larger than s. The note correctly avoids that error.

All domains and alphabets here are finite. Thus the point-list
cardinalities and their maximum are finite, Lambda is not top, and
the toNat appearing in the archived certifiedGammaError expression
is the actual integer maximum L_s. The reindexing from two blocks of
eight rows to sixteen rows preserves column Hamming distance and
therefore preserves the maximized list size.

## Ratio counts agree with the archived IsMCA and mcaError

[ProximityGenerators, lines 109–115](../sources/parallel23-prize-dependencies/ProximityGap_ProximityGenerators.lean)
defines IsMCA by existence of a set T with sufficient cardinality,
explainability of the combined word on T, and nonexplainability of
at least one input word there. Its AffineLineGenerator, at line 280,
is exactly γ↦(1,γ). The combination is consequently W+γZ, including
γ=0.

Remainder is linear, so folded explainability is
ρ_T(W)+γρ_T(Z)=0. If ρ_T(Z)=0 as well, this equation forces
ρ_T(W)=0. Both original words would then be explainable on T,
contradicting the MCA predicate. Conversely, a solution with
ρ_T(Z)≠0 supplies a folded explanation and a nonexplainable input
Z, hence an actual MCA witness. A nonzero coordinate of ρ_T(Z)
determines γ uniquely. Different witnesses may give the same γ,
which is why the set of distinct ratios, rather than a subset count,
is necessary.

The archived predicate allows |T|≥s, whereas the proposed ratio set
uses |T|=s. Their equivalence needs the argument given in the note:
from a larger witness choose any k points and interpolate each row
of Z. Nonexplainability supplies a point where a row disagrees with
that interpolant. Those k+1 points already prevent an explanation;
extend them to s points within the witness. Nonexplainability of Z
persists and the folded explanation restricts. The hypothesis
s≥k+1 is exactly what makes this reduction valid. No monotonicity
of nonexplainability under arbitrary restriction is being assumed.

ProximityGenerators lines 152–154 define mcaError as the supremum,
over all input families, of this event's probability under a uniform
scalar. Finite alphabets make the supremum a maximum. Thus the
quantity is exactly B_s/q, where B_s is the maximum distinct-ratio
count described in the note. The archived IRS equality then gives

    Gamma(j/n)=(B_s+L_s)/q.

Therefore its inequality against 2^(−128) is equivalent to the
stated **joint** integer budget. It is neither two independent
budgets of that size nor equality with the exact winning-set error.
The reviewed text preserves both distinctions.

## Verification scope and final limitations

The recorded verifier supplies 841 Newton identities, 10 monic list
classifications, two extension-valued examples, 24 distinct-remainder
list comparisons, 72 MCA ratio comparisons, and an explicit small
polynomial-phase box collision. Its polynomial remainder, Newton
indexing, field arithmetic, and distinct-set counting were reviewed.
The source parses successfully. The recorded hashes match, so those
results are consistent with the reviewed code. They are historical
checks, not a new full execution by this reviewer.

A separate bounded check was run in memory on the constant code over
F_5 with four columns. For twelve deterministically sampled words at
each of widths one and two, every ordered pair was checked at s=2
and s=3. All 576 cases gave equality between: MCA witnesses of any
size at least s; witnesses of size exactly s; and nonzero-direction
ratio conditions. This independently checks the shrinking issue
that the companion verifier's fixed-size comparisons do not directly
exercise. It is finite evidence supporting the proof, not an
exhaustive production-field calculation.

No issue requires changing the mathematical statements in this
version. The exact final count is still unbounded, and the conditional
monic family does not dominate arbitrary received words. The phase
obstruction supplies no list or MCA lower bound. No current prize
status, production-size maximum, submission, Lean build, or complete
Paley proof is claimed. Only this new review file was written.
