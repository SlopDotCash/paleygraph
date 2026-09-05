# Exact observables, moment trades, and retained incidence

This tool asks which measurements lose information about an actual arithmetic target. It produces exact pairs of prime-field inputs that agree on chosen measurements but disagree on the target, and then tests additional measurements against those pairs. It does not attempt a Paley or prize proof.

The target is `M6(C) = sum_x (sum_(c in C) chi(x-c))^6`. All character values and moments are computed as integers. The tested measurements are `M2`, `M4`, additive energy `E2`, the two boundary statistics `B2=sum_(c in C) F_C(c)^2`, `B4=sum_(c in C) F_C(c)^4`, and the **unordered** deck of quartic correlations `K(Q)=sum_x prod_(c in Q) chi(x-c)` for all four-subsets of C.

## Iteration 1: actual kernels with indistinguishable measurements

`fiber_audit.py` partitions actual sets by their exact measurement vectors. For each group, or *fiber*, it stores target minima/maxima and an extremal witness pair. It also reports whether the fiber straddles the explicitly chosen baseline `15 p |C|^3`. That Gaussian-shaped baseline is not assumed to be a theorem.

At p=61, consider:

```
C = {0,1,2,3,7,9}
D = {0,1,2,9,10,17}.
```

Both have `M2=330`, `M4=5766`, `B2=14`, `B4=86`, `E2=98`, and quartic deck

```
[-11,-7,-3,-3,-3,1,1,1,1,5,9,9,9,13,13].
```

But `M6(C)=140190` and `M6(D)=163230`. A second exact pair is stored at p=257. These are actual character evaluations, unlike the earlier weighted sign-model obstruction in pass 22.

**What this establishes:** M6 is not determined by these measurements on this finite domain. It does **not** show that no useful upper bound can depend on those measurements. A valid bound must cover the largest target in its fiber; the finite fiber maximum is itself a finite upper bound. The examples do not violate the Gaussian baseline, and they are not Paley counterexamples.

The exhaustive suites cover all sets of size six containing 0 and 1 at p=13,17,29,41,61. Every affine orbit has a representative in this slice, and all the selected measurements and even moments are affine invariant. Frequencies on this slice are **not** uniform frequencies on affine orbits or on all subsets.

## Iteration 2: the primitive missing motion

`moment_trades.py` explains the collision by an integer histogram move. For six columns, away from the six zero fibers, `|F_C(x)|` lies in `{0,2,4,6}`. Write h for its four histogram counts. Fixing its mass, second moment and fourth moment leaves the primitive integer-kernel direction

```
v = (-10, 15, -6, 1).
```

Indeed, `sum v_i x_i^(2j)=0` for `j=0,1,2`, while `sum v_i x_i^6=23040`. The three-by-four constraint matrix has rank three, and the last coordinate of v is one, so every integral kernel element is an integer multiple of v. The boundary absolute values are `{1,3,5}`; their mass, B2 and B4 uniquely determine their histogram by an invertible Vandermonde matrix.

Thus every pair with those fixed lower/boundary measurements has an interior histogram change `t v`, and its M6 changes by `23040 t`. Nonnegativity supplies an exact integer interval for t. The tool compares this *histogram relaxation* with the genuinely realizable character-set fiber.

| p | Normalized sets | Actual lower/boundary fibers | Fibers with a strict gap below the relaxed upper endpoint | Largest gap |
|---:|---:|---:|---:|---:|
| 29 | 17,550 | 64 | 2 | 23,040 |
| 41 | 82,251 | 135 | 34 | 23,040 |
| 61 | 455,126 | 257 | 205 | 46,080 |

These exhaustive finite comparisons identify a concrete missing object: **which integer moment-preserving moves can actually be realized by character translates?** The integer move itself is elementary known mathematics; the arithmetic realization tests are the specialized experiment.

## Iterations 3 and 4: recover information already present in the quartics

For six columns, a quartic subset is indexed by its omitted pair. Attach its K-value to that edge of a weighted complete graph W on the six columns. The unordered deck forgets this incidence. The first refinement retains the multiset of the six sorted incident edge lists, a familiar graph invariant.

In the full p=61 census, adding this incidence signature reduces 63 target-ambiguous fibers to one. The remaining fiber has 60 normalized representatives. Its two extremal sets are

```
C = {0,1,2,4,38,52},  M6=101790
D = {0,1,3,8,10,21}, M6=124830.
```

Their rooted edge signatures coincide, but their weighted graphs are nonisomorphic. The next refinement retains the triangle contraction `tr(W^3)`, computed from the **same already evaluated quartics**, with no additional character sums. Its values are `-1992` and `4152` for that pair. The refined measurements have no remaining M6 ambiguity in the entire normalized p=61 census. At p=257 they also resolve the three sampled ambiguous fibers.

This is a finite feature-refinement result, not an identity expressing M6 through triangle contractions, a universal determination theorem, or a proved asymptotic bound. Rich fingerprints can simply identify each input. The tool therefore records fiber counts and retains counterexamples rather than reporting predictive accuracy as a theorem. Graph incidence and triangle traces are established tools; no invention of either is claimed.

## Scale and the full-problem gap

Seeded tests additionally cover 8,192 six-sets at p=257; 4,096 six-sets at p=1297; and 1,024 eight-sets at p=4099. The last two lie at `p >= n^4`, the first size slice relevant to the existing sixth-moment reduction. **Every full quartic fingerprint in those two samples is a singleton.** The absence of collisions there is therefore uninformative about sufficiency, and it must not be presented as successful transfer of a theorem. Strong exact collisions at p=61 and p=257 lie outside that size slice.

The next useful experiment is a constrained search for critical-scale sets preserving a *small* incidence-based statistic, with new primes reserved for evaluation. The mathematical opportunity is an estimate on quartic-incidence contractions or on arithmetic realizability of histogram moves. Nothing here supplies such an estimate. Computing all quartics costs `O(p binom(n,4))` per set; the validation also computes sixth correlations. These combinatorial costs prevent direct large-n use. Sparse incidence sampling, cached correlations and integer-fiber search are unimplemented extensions.

## Run and verify

Use a Python environment with NumPy and SymPy. The environment used here is `/opt/miniconda3/bin/python3`.

```sh
python3 fiber_audit.py --suite toy
python3 fiber_audit.py --suite scale
python3 fiber_audit.py --suite holdout
python3 moment_trades.py
python3 incidence_refinement.py
python3 ../novelty/replay_observable_certificates.py
```

The discovery code verifies the exact second-moment and boundary-corrected sixth-moment identities, directly replays every exported witness through a separate scalar implementation, checks affine invariance, tests monotonicity of refinement, and includes M6 itself as a deliberately excluded leakage control. `T6` is also excluded because the existing identity makes it equivalent to the target after the other terms are fixed. Conservative integer-overflow and invalid-sampling guards restrict the API to its supported exact domain.

The independent standard-library verifier imports none of the discovery code. It replays all 60 stored first-stage witness endpoints and independently enumerates the complete p=29 realization test. Results, source hashes, literal inputs, and scope are saved in JSON. Numerical timings vary; exact mathematical outputs are deterministic.

## Closest prior art

- [Clarke et al., counterexample-guided abstraction refinement](https://www.cs.cmu.edu/~emc/papers/Conference%20Papers/Counterexample-guided%20Abstraction%20Refinement.pdf): the refinement strategy is established.
- [Kemper–Lopatin–Reimers, separating invariants over finite fields](https://arxiv.org/abs/2011.07408): separation by invariants is established.
- [Lasserre, bounding support from marginal moments](https://arxiv.org/abs/1011.0138): moment information and support bounds have an established theory.
- [Diaconis–Sturmfels, algebraic algorithms for conditional distributions](https://www.math.ucdavis.edu/~deloera/MISC/LA-BIBLIO/trunk/Diaconis/Diaconis8.pdf): integer fibers and moment-preserving moves belong to established algebraic statistics.

The potential original contribution is the particular arithmetic diagnostic, its realization-gap certificates and the experimentally selected quartic-incidence refinement. Global mathematical novelty remains unestablished; the [literature audit](../novelty/literature_ledger.md) makes the search scope explicit.
