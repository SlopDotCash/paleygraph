# Third parallel pass: the arbitrary-gap family is proved

**There is a concrete advance on the necklace subproblem; the Paley
conjectures and the general prize bridge remain unproved.** The preceding
pass handled only two separations between exceptional labels. This pass
handles every separation, with an explicit constant uniform in the length.
That closes one stated proof obligation, while leaving the weighted
aggregate required for clique bounds open.

Three agents worked on necklaces, subgroup mixed energy, and classical
fourth moments. The primary agent worked on a conic representation of
the two-anchor graph, audited the new analytic inputs, and integrated the
results. Independent agent reviews checked the conic calculation and the
classical counterexample. A subgroup endpoint inequality was corrected;
the conic review prompted an explicit centering-matrix check. Final
artifacts incorporate these changes.

## The new uniform family

For A={0}, B={1}, C={0,1}, every prime p=1 mod 4, and every pair of
positive gaps j,m, the [necklace proof](parallel3-necklace-2026-09-04.md)
establishes

\[
 \left|N\left(BA^{j-1}CA^{m-1}\right)\right|
 \le(\min(j,m)+1)p^{(j+m+1)/2}+1.
\]

The proof identifies a rank-j hypergeometric kernel exactly, including
its singular value. Tensoring it with the rank-two kernel produces a
cohomology space of dimension j+1 when j differs from two; the equal-rank
exception is corrected explicitly using the preceding pass. The argument
uses [Katz's hypergeometric theorem](https://web.math.princeton.edu/~nmk/g2hyper62finalcorrected.pdf)
and the curve trace and weight machinery in
[Katz's monograph](https://web.math.princeton.edu/~nmk/Katz-GKM.pdf).
The relevant primary pages were archived, rendered, and independently
inspected. These are imported established theorems, not conjectural inputs.

Label permutations cost at most kp^(k/2), where k=j+m. Thus the
normalized bound tends to zero uniformly for k=o(sqrt(p)), including
logarithmic length. Words with multiple occurrences of both exceptional
labels and their weighted aggregate remain outside this theorem.

For finite bookkeeping, the proved two-label families and this family
cover at least 279 of the 729 words of length six over {A,B,C}: 189 use
at most two labels, and 90 have label multiplicities (4,1,1). This count
measures coverage at one length; it is not a fraction of Paley proved.

## The other three lanes

| Lane | Verified result | Limit |
|---|---|---|
| [Subgroup arithmetic](parallel3-subgroup-2026-09-04.md) | The primitive polynomial discriminant is an explicit order index squared times a power of two. Mixed energy satisfies B<=k^2+4k v_p(I_(2k)). Quartic examples show that unramifiedness, nonzero derivatives, and Galois symmetry do not force collision-free fibers. | The pointwise bound v_p(I_(2k))=O(k log k) is unproved and can be stronger than the energy requirement. Appropriate smaller tower levels and higher centered moments still need control. |
| [Classical moments](parallel3-classical-2026-09-04.md) | The elliptic trace has an exact Fourier transform in terms of Kloosterman sums. An interval construction refutes an unrestricted Fourier-L1 sufficient estimate, even after deleting any fixed number of cross-ratios. | The prime construction has no polynomial upper bound in the interval length. It does not refute a bound restricted to n>=p^epsilon for fixed epsilon, or refute the moment hypothesis or Paley. |
| [Two-anchor matrix](parallel3-conic-2026-09-04.md) | A conic parametrization gives an exact lifted kernel f(t/u)f(tu), a boundary block, and a Mellin matrix with Jacobi-product entries. Six automorphisms give exact smaller spectral blocks. | The boundary coupling has scale sqrt(p), and the Mellin matrix has nonzero off-diagonal entries. Neither representation yields an improved uniform spectral bound. |

The main clique obligation remains the
[signed aggregate estimate from the second pass](parallel2-spectral-transfer-2026-09-04.md).
Adding individual bounds by the triangle inequality still loses an
exponential factor. The new family supplies a useful analytic tool and
a larger proved class, but does not supply that cancellation. The uniform
thin-subgroup target, full arbitrary-two-set conjecture, and proposed
connection to the official proximity challenges remain separate gaps.

## Verification and reproducibility

- [Necklace verifier](../experiments/parallel3_necklace_2026_09_04.py),
  [results](../results/parallel3_necklace_2026_09_04.json):
  320 gap identities and bounds, 192 zero-retaining averages,
  46 exact norm certificates with 864 positive principal minors, and
  48 literal cyclotomic hypergeometric fibers. Tests include ranks at
  least the characteristic, as well as 10 literal necklace sums.
- [Subgroup verifier](../experiments/parallel3_subgroup_2026_09_04.py),
  [results](../results/parallel3_subgroup_2026_09_04.json):
  five integer index/discriminant certificates and three
  modular certificates, including exact index valuations at the two
  quartic resonances. No new prime scan or factorization was used.
- [Classical verifier](../experiments/parallel3_classical_2026_09_04.py),
  [results](../results/parallel3_classical_2026_09_04.json):
  106 exact Fourier transforms, 900 symmetry checks, 54
  directional-pairing and affine-invariance checks, and 20 integer
  configuration records, of which 12 were independently counted.
- [Conic verifier](../experiments/parallel3_conic_2026_09_04.py),
  [results](../results/parallel3_conic_2026_09_04.json):
  93312 kernel entries across 14 primes, 2800 exact cyclotomic
  Fourier entries, 84 quotient/centering checks, and 84 automorphisms.

These finite checks support the
displayed algebra and examples; they do not replace the ordinary proofs
or certify the imported general theorems. No new Lean certificate,
external submission, or claim of literature novelty is made.
