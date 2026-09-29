# Seventh parallel pass: stronger subgroup cancellation and one crossing class

**Assessment: mathematical progress, with the full goal still unproved.**
This pass improves a subgroup maximum estimate on an asymptotically
full proportion of the eligible quartic-window primes. It also proves
one more fixed-length necklace class. Neither result supplies uniform
square-root cancellation, the classical two-set conjecture, or a
reduction to the official Reed–Solomon prize.

Started September 4; final local verification completed September 5, 2026.
The three parallel workers covered subgroup moments, necklace contractions,
and the classical arbitrary-set approach. Their account usage limit
interrupted the workers. Root preserved their drafts, completed the
subgroup verifier, independently audited the necklace sources and
boundaries, and ran all four final verifiers locally.

## 1. The main improvement

For dyadic N, let
P_N={p prime:N⁴/4≤p≤N⁴, p≡1 mod N}.
On the same class K_N from the fifth pass, with
|K_N|/|P_N|=1−O(1/log N), the new result is

\[
 \boxed{\max_{a\ne0}\left|\sum_{h\in H_N}e_p(ah)\right|
       \le(17+\log N)^{1/9}N^{8/9}.}
\]

The previous exponent on this class was 17/18. The desired exponent
is 1/2+o(1), at every eligible prime, so there are still two gaps:
the strength of cancellation and the exceptional primes.

The [multiplier proof](parallel7-multiplier-2026-09-04.md) is elementary
apart from the previously proved sixth-energy prime class and its
imported prime-counting input. For any multiplicative subgroup H of
size n≥2, write θ=η/n and let α_r be the inverse Fourier transform
of |θ|^r. This function may be signed. Jensen is used only against
the genuine s-fold sum probability measure μ_s. Bilinear
orthogonality then yields

\[
 |\theta(a)|^{rs}\le\frac2n+
 \frac{\sqrt{p\,\widetilde{\mathcal E}_r\widetilde E_s}}{n^{r+s}}
 \quad(r\ge2,\ s\ge1\text{ an integer}).
\]

The proof also records the exact origin atoms and negative constant
on the nonzero coordinates. With r=s=3, p≤n⁴ and E₃≤Bn³, this gives
M≤(B+2)^(1/9)n^(8/9). No parity, zero-triple exclusion, or extra
fourth-energy assumption is required.

Higher-energy feedback now gives centered E₄ at exponent 43/9,
improving 44/9. An additional hypothesis centered E₄≤Dn⁴ would
give M≤2^(1/4)D^(1/8)n^(7/8); that hypothesis is not proved here.

The [independent subgroup audit](parallel7-subgroup-2026-09-04.md)
accepted the multiplier argument and established two limits:

- Actual quartic-window examples have negative α₃ coordinates.
  Therefore α₃ cannot be used as a probability measure for another
  unrestricted Jensen step.
- Across every fixed real r≥1 and integer s≥1, substituting the
  present energy estimates into this gate gives its largest exponent
  saving, 1/9, at r=s=3. This is a limitation of these substitutions,
  not a lower bound on actual sums or a barrier to other methods.

## 2. A new necklace contraction

The [crossing proof](parallel7-necklace-2026-09-04.md) gives

\[
 |N(ABABCC)|\le56p^{7/2},\qquad
 |N(w)|\le58p^{7/2}
\]

for its full 36-word orbit under rotations, reversal and label
permutation, for every prime p≡1 mod 4.

Two rank-four irreducible middle convolutions have local type
χ⊕1³ at the crossing. After the final quadratic twist one has
χ³⊕1, which excludes an invariant in their tensor product.
The proof restores the infinity correction, finite zero masks,
exceptional fibers and graph-normalization terms explicitly.
Katz's rank, irreducibility, local-monodromy and weight statements
were checked in the archived primary manuscript.

Length-six coverage advances from 651 to **687 of 729** words.
The remaining classes are AABCCB (18), ABACBC (18), ABCABC (6).
This count is not a percentage of progress toward Paley. Even all
length-six words would leave the growing-depth signed aggregate
and its spectral consequences unresolved.

## 3. What the classical lane established

The [classical analysis](parallel7-classical-2026-09-04.md) derives an
exact fourth-moment conditional mean given a four-point core and an
exact local perturbation estimate. It also constructs an artificial
nonnegative affine-invariant quartic on n=floor(p^(1/3)) subsets with
the same mean and variance as the true fourth moment, bounded
quartic coefficients, and very large exceptional values.

This proves that those structural and statistical facts alone do
not control every set. The artificial function does not have the
full character-sum structure: it is not a counterexample to Paley
or to the desired uniform moment estimate. No new worst-case
classical character-sum upper bound was obtained.

## 4. Verification

All four verifiers completed successfully on their final recorded
inputs. Their checks supplement ordinary proofs; they do not prove
asymptotic statements or formally verify imported theorems.

| Lane | Recorded checks |
|---|---|
| [Multiplier](../experiments/parallel7_multiplier_2026_09_04.py) | 8 subgroup cases, 24 norm checks, 261 Fourier-inversion enclosures, 924 strict gate certificates and 24 flat-spectrum algebraic cases, covering 3,816 frequency/order pairs |
| [Subgroup audit](../experiments/parallel7_subgroup_2026_09_04.py) | 2 negative cosets, 84 positive-measure Jensen checks, 15,936 rational order pairs, 3,162 derivative checks, exact order-64 zero-triple witness |
| [Necklace](../experiments/parallel7_necklace_2026_09_04.py) | 2,700 cross-ratio entries, 2,795 affine Legendre entries, 285 zero-stalk checks, 190 compact trace vectors, 5,590 pointwise traces, 95 inner expansions, 180 orbit values, 48 extension-field traces and 12 degree-four polynomials |
| [Classical](../experiments/parallel7_classical_2026_09_04.py) | 1,750 conditional identities over 23,720 completions, 54,285 local set pairs, 4 exact affine-hypergraph fields |

The outputs are [multiplier](../results/parallel7_multiplier_2026_09_04.json),
[subgroup](../results/parallel7_subgroup_2026_09_04.json),
[necklace](../results/parallel7_necklace_2026_09_04.json), and
[classical](../results/parallel7_classical_2026_09_04.json).
The [integration audit](../results/parallel7_pass_audit_2026_09_04.json)
pins input, artifact and source hashes, local links and review scope.

## 5. Remaining work

The useful next subgroup question is a stronger centered eighth-energy
budget or an estimate using arithmetic information that the current
energy substitution discards. The necklace question is the remaining
three classes and a method that controls the full signed aggregate
as length grows. The classical question is uniform control of
exceptional sets using their character-sum structure.

The full restricted operator beyond the enlarged kernel span is also
uncontrolled. The exact connection to the official prize remains a
separate unproved obligation. No novelty, best-known-bound, Lean
formalization, prize submission or completion claim is made.
