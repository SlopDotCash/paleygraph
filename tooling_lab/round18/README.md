# Round18: counting an entire arithmetic section with digit constraints

[Round17](../round17/README.md) showed that exact norm and digit energy do not determine whether a compressed word represents an actual centered coset. This round replaces the large digit cube by an exact language for the generator-2 family, counts the entire section, and examines every vector with minimum digit support.

Signed-digit numeration and transfer matrices are established mathematics. The [prior-art audit](prior_art.md) identifies close published and local antecedents. These are new executable tools and finite results in this lab; global historical novelty is unestablished. Neither prize is being proved.

## What the tool retains

For Q=2^N+1=k*p and f=2-X^-1, write D=f*F_a/p. The exact nonzero membership rule is:

1. Nonzero digits alternate in sign, including across intervening zeros.
2. Their number is odd, enforcing the antiperiodic endpoint.
3. S(D)=sum_j2^(N-1-j)*D_j is divisible by k.

The centered scalar is S(D)/k; zero has its own state. The [derivation](DERIVATION.md) proves necessity and sufficiency. The [recognizer and counter](digit_language.py) keep sign state, residue modulo k and digit support. They compute exact distributions without enumerating p scalars or3^N digit words. This is a specialized sufficient interface for realizability, rather than an additional norm statistic.

## Complete arithmetic support counts

| Prime p | Subgroup order n | Cofactor k | All scalars checked, including zero | Minimum nonzero digit support | Maximum digit support |
|---:|---:|---:|---:|---:|---:|
|17|8|1|17|1|3|
|257|16|1|257|1|7|
|65537|32|1|65537|1|15|
|6700417|64|641|6700417|5|27|

At the largest input, the full alternating language has4294967297 words including zero. The cofactor selects exactly6700417. The producer uses862064 transitions and at most36025 populated states at a layer. Runtime in the log is uncontrolled; no benchmark speedup is claimed. Complexity depends on k: O(k*N^2) integer additions, with integer bit lengths and storage additional costs. A huge cofactor can defeat this representation.

The [separate review](language_review.json) uses binary sign words and positional residue weights, instead of signed-digit transitions. It agrees on every support count in all644 residue profiles across the four cases. It also directly enumerates all6766228 actual scalars by centered modular doubling, checking215464060 scalar-coordinate transitions. The largest result thus covers every scalar. Norms are separately checked only for the512 minimum-support vectors below, not for all6.7 million scalars.

Tiny controls exhaust all6642 ternary words at N4 and N8 using direct recurrence inversion. The review also checks130 saved generator-2 records,48 modified words, and three words that each satisfy two of the main conditions and fail the third.

## Sparsest actual digits still have substantially different norms

Every nonzero actual D at p6700417 has at least five nonzero digits. Exactly512 attain five. These counts refer to scalar encodings; each multiplicative subgroup coset has64 such encodings, related by signed rotation of the digits. A complete support search examines412736 candidate signed words of support1,3 or5 and retrieves exactly the same512 scalars as the full independent census.

| Norm defect tau at minimum digit support | Scalar encodings |
|---:|---:|
|449|128|
|1217|64|
|8513|64|
|15937|64|
|24001|64|
|40193|64|
|84481|64|

Every listed norm is independently computed by an integer resultant, with Norm(D)=641*tau. All these words have digit energy five, while their norm defects vary by more than a factor188. Even exact minimum support does not determine norm in this actual section. This is not a lower bound for every vector in the ideal, a global minimum-norm theorem, or an asymptotic bound. It says nothing by itself about the coefficient energy of the original F.

## Reusing other saved encodings

The [orientation adapter](orientation_adapter.py) recognizes when X->X^e changes the saved relation into a signed monomial times2-X^-1. It converts digits with one signed permutation. The [separate polynomial-substitution review](orientation_review.json) verifies all451 supported records:

- p257,g249: e11, relation multiplier X^5;256 records.
- p65537,g65529: e19, multiplier-1;65 records.
- The two p6700417,g2 walks: identity;65 records each.

At p2013265921,N64, generator2 has the wrong order, so the adapter declines the33 saved records. The general round17 inverse checker still verifies them. A deliberately doubled relation is rejected as a nonassociate.

## Next missing interface

For general short f, round17 gives an exact lattice intersected with a parallelotope, but the inverse can be dense and the cofactor enormous. Next compare small portfolios of short-relation digit bounds with exact membership, and learn which constraints they omit. In the generator-2 laboratory, compare bounded digit windows with constant sign memory before attempting a general relation. Preserve the arithmetic congruence and endpoint; dropping either is already falsified by controls.

```sh
/opt/miniconda3/bin/python3 tooling_lab/round18/digit_language.py
/opt/miniconda3/bin/python3 tooling_lab/round18/language_review.py
/opt/miniconda3/bin/python3 tooling_lab/round18/orientation_adapter.py
/opt/miniconda3/bin/python3 tooling_lab/round18/orientation_review.py
/opt/miniconda3/bin/python3 tooling_lab/round18/verify_round.py
```

The [manifest](manifest.json) binds sources, reports and inputs and checks preservation of round17 and round16 artifacts. Reviews are root-run separate implementations, not external peer review or Lean certification.
