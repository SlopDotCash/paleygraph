# Critical-size audit and exact exchange decomposition

This round tests whether round one's successful toy separator carries useful information at `p approximately n^4`. Its first result is negative: adding the raw and centered quartic-incidence triangle to a small linear predictor does not help on the three critical-size random holdouts. Its second result is a new diagnostic implementation: exact Johnson-scheme projections identify which interaction orders comprise the missing sixth character statistic, without enumerating all replacement sets.

These answer different questions. The first tests a particular target-free feature representation. The second deliberately evaluates and decomposes the target; it is **not a new predictor or an estimate that avoids the target**. Neither provides a worst-case bound.

## Frozen arithmetic quantities

For `C` an n-element subset of Fp, define

```
K(Q) = sum_x product_(c in Q) chi(x-c)
T6(C) = sum_(Q subset C, |Q|=6) K(Q).
```

The generator computes the elementary coefficients of the character rows through degree six in `O(6np)` integer operations, and checks the existing exact identity relating T6 to M6, T4, M2 and boundary terms. It independently checks selected sets with round one's scalar implementation. All characteristic arithmetic is exact; the regression is numerical.

Using only the already computed quartics, form a symmetric graph with zero diagonal and

```
W_ij = sum_(Q subset C minus {i,j}, |Q|=4) K(Q).
```

At n=6 each edge contains one quartic, exactly as in round one. At larger n it contains a sum. Two low-cost contractions are `tr(W^3)` and `tr((PWP)^3)`, with `P=I-J/n`. The implementation stores `n^3 tr((PWP)^3)` as an integer. Three controls independently shuffle the quartic values before attaching them to their column subsets. This preserves the unordered quartic deck and removes its actual incidence assignment.

The baseline has intercept, normalized M4, B2, B4 and additive energy. The candidate adds the two triangle contractions. The matched-size control adds the corresponding pair from one fixed shuffle. The target is `T6/sqrt(p*binomial(n,6))`. Normalizations, seeds, all input sets and both initial and revised cohort analyses are saved.

## Result: the toy prediction advantage does not survive

The revised analysis canonicalizes every set under all affine maps determined by ordered pairs of its own points. It retains one representative per affine orbit, then uses a fixed hash to separate training and evaluation. The p=61 sample lost 143 duplicate-orbit records. The critical-size samples had none. Three-anchor stress inputs are also excluded if their orbit appeared in training.

| p | n | Training / evaluation orbits | Baseline holdout R² | Add true incidence | Add shuffled incidence |
|---:|---:|---:|---:|---:|---:|
| 61 | 6 | 958 / 947 | 0.001708 | 0.019456 | 0.000846 |
| 1297 | 6 | 4083 / 4109 | -0.001856 | -0.002074 | -0.001972 |
| 2437 | 7 | 2024 / 2072 | -0.000698 | -0.002547 | -0.001700 |
| 4129 | 8 | 1023 / 1025 | -0.000763 | -0.000893 | -0.003707 |

Here R² is measured against the evaluation cohort's own mean. Negative values mean worse mean-square prediction than that constant benchmark. These are finite fixed-seed measurements, without confidence intervals or a model-selection significance claim. The conclusion is limited to this linear model and these contractions; it does not show that the full incidence graph lacks information or rule out a nonlinear use of it.

Transferring the p=61 candidate coefficients directly to the three larger primes gives R² of approximately -0.038, -0.163 and -0.525. Thus the toy gain is a poor basis for an extrapolated estimate. Three-anchor stress cases provide a separate selected distribution, not an exhaustive exceptional-set search. No uniform tail control follows from random-cohort prediction.

`cohort_review.py` independently reconstructs all quartics, T6 and M6 by scalar character products on twelve selected extreme inputs. It directly multiplies the centered graph matrices to verify the contraction formula. All checks pass. `cohort_review.json` is the revised result; `results.json` preserves the initial literal-set split so the correction is visible.

## A diagnostic that asks where the missing information lives

`exchange_spectrum.py` treats all n-subsets of Fp as a Johnson graph: one step replaces one selected point by one unselected point. This changes the actual arithmetic input while staying on the same finite-set domain. Let P be the uniform one-swap operator. Its known eigenvalues on harmonic level d are

```
lambda_d = 1 - d(p-d+1)/(n(p-n)).
```

T6 is a multilinear polynomial of degree six in the set-membership indicators, so it belongs to levels 0 through 6. The exact level-d projector is

```
T6_d = product_(e=0..6, e!=d) (P-lambda_e)/(lambda_d-lambda_e) applied to T6.
```

Ordinarily, applying P repeatedly would enumerate a huge number of nearby sets. Here the character-row generating functions evaluate every needed distance average directly. For a fixed row x, let `e_r(in)` and `e_r(out)` be elementary coefficients of its signs on C and its complement. Their full product is

```
product_(c in Fp) (1+t*chi(x-c)) = (1-t²)^((p-1)/2).
```

Division by the short inside polynomial gives the outside coefficients. If D is a uniformly chosen n-set obtained by replacing j points of C, then

```
E[T6(D) | |D minus C|=j]
  = sum_(r=0..6) [binom(n-j,r)/binom(n,r)]
                    [binom(j,6-r)/binom(p-n,6-r)]
                    sum_x e_r(in)*e_(6-r)(out).
```

This is elementary inclusion probability, not a probabilistic approximation. Character rows with the same numbers of positive and negative entries are grouped. Only seven coefficients and seven distance averages are needed. A short exact birth-death recurrence converts distance averages into powers of P and hence all seven projections. Python rational arithmetic avoids unstable subtraction of close floating-point quantities.

The tool uses established slice harmonic analysis. Its specific contribution is the arithmetic generating-function backend and its connection to the earlier observable-loss experiment. The projections are a way to diagnose which interaction orders an argument must address, not a claim to have invented Johnson schemes or Hoeffding decompositions.

## Observed harmonic content

The p=13,n=6 and p=17,n=6 runs exhaust all normalized sets containing 0 and 1. This gives exact uniform-subset averages for these affine-invariant quantities: the affine group is sharply two-transitive, so the number of affine maps taking two selected points to 0 and 1 is constant. All off-diagonal entries of the exact harmonic Gram matrix vanish in both exhaustive runs.

| p | n | Inputs | Highest-level squared energy / centered second moment |
|---:|---:|---|---:|
| 13 | 6 | Exhaustive, 330 | 0.511380 |
| 17 | 6 | Exhaustive, 1365 | 0.838498 |
| 61 | 6 | 512 seeded samples | 0.998980 |
| 1297 | 6 | 512 seeded samples | 0.999864 |
| 2437 | 7 | 512 seeded samples | 0.999971 |
| 4129 | 8 | 512 seeded samples | 1.000026 |

The last four rows are sampled ratios, not exact proportions. The last ratio slightly exceeds one because cross terms need not vanish in a finite sample. At the three critical-size cases the sum of squared lower nonconstant components divided by the centered second moment is about `5e-6` to `9e-6` in these samples. This is evidence that this target's ordinary fluctuations sit overwhelmingly in its highest interaction level on the tested distribution. Rare exceptional sets can behave differently.

### Revision: exact population energy, then compression beyond field enumeration

An [independent implementation](review_exchange.md) now replaces the sampled energy calculation by exact all-input second moments. It uses bivariate character generating functions and Eberlein eigenvalues, without importing the original projector algebra. At p=1297,2437,4129 the exact lower-degree shares are respectively `9.020826e-6`, `7.628071e-6`, and `5.302080e-6`. At p=61 the exact share is0.005227816868, demonstrating the sampling error in the earlier row.

The [grouped backend](exchange_grouped.md) compresses the calculation again: only six character-pair types occur. Raising their short generating polynomials to the known multiplicities eliminates field-element enumeration. Falling-factorial containment probabilities eliminate factorials of full field sizes. At fixed degree, the coefficient work is independent of p and n apart from integer bit complexity and the prototype's primality test.

This gives an exact all-input lower-degree share of `3.458075708225e-10` at p=6,700,417,n=50, in a fraction of a second. All exact fractions and six earlier cross-checks are saved. This is a new usable diagnostic capability in the lab: determine the complete second-moment distribution among interaction levels at a field size where enumerating input sets is impossible. It leaves rare-set behavior unresolved.

Orthogonality gives a precise conceptual limit: in the full uniform L2 space, a component on level six is orthogonal to every polynomial of membership degree at most five. That statement concerns polynomial degree on the **whole p-coordinate slice**, not the number of scalar features. A nonlinear function of lower-order measurements can have higher polynomial degree. In particular, a quartic triangle can contain level-six content; the projection result does not itself prove that the tested triangle predictor must fail.

The next useful experiment is to measure a proposed representation's overlap with the exact level-six component, while independently testing exceptional arithmetic families. An L2 statement, even an exact one, will still leave the conjecture's worst-case requirement open.

`stress_exchange.py` applies the frozen diagnostic to 28 selected arithmetic inputs: multiplicative subgroups, arithmetic and geometric progressions, and extremes from the saved random/three-anchor cohorts. Lower components remain small in absolute value on these inputs, but can matter when the total target is small: the size-eight subgroup at p=4129 has T6=-12 and lower nonconstant sum approximately3.213. The tool therefore preserves signed component values rather than declaring every selected input "99.99% high order" from a population-level measurement.

The same script revisits both exact p=61 collision pairs from round one. **All components of degrees0 through5 agree exactly within each pair; their degree-six components differ by32.** This localizes the old ambiguity to a specific irreducible interaction level, instead of merely reporting a difference in M6. It is an exact certificate on those pairs, not a classification of every ambiguity fiber.

## Verification and reproduction

The exchange generator separately enumerates all 1716 six-subsets of F13, directly computes their target values, and checks three complete distance profiles. It also checks 21 exact eigen-equations by summing each projector over every one-swap neighbor. The all-input Gram orthogonality checks and saved exact fractions provide additional evidence beyond numerical agreement.

From this directory:

```sh
python3 critical_audit.py
python3 cohort_review.py
python3 exchange_spectrum.py
```

The generator uses NumPy for character rows and Python integers/rationals for large contractions and projection coefficients. All sources and results remain in the lab; the earlier proof records are read-only.

## Sources and prior-art boundary

- [Filmus and Mossel, Harmonicity and invariance on slices of the Boolean cube](https://cris.technion.ac.il/en/publications/harmonicity-and-invariance-on-slices-of-the-boolean-cube/) supplies the established harmonic framework; the authors' [publication list](https://yuvalfilmus.cs.technion.ac.il/publications/papers/) records the slice/Johnson-scheme setting and earlier work.
- [Suzuki's Johnson graph lecture](https://icu-hsuzuki.github.io/t-algebra/lec7.html), Lemma7.1, records the adjacency eigenvalues used by the normalized swap operator. The prototype also checks the resulting eigen-equations directly on tiny arithmetic inputs.
- [Kunisky, Moore and Wein, Tensor cumulants for statistical inference on invariant distributions](https://arxiv.org/abs/2404.18735) is prior art for organized invariant polynomial contractions. No theorem from that paper is used to certify these particular arithmetic contractions.
- Local `research/parallel22-all-orders-independent-review-2026-09-05.md` already discusses Walsh information loss; pass22/pass23 notes already distinguish moment reformulations from estimates. The current implementation does not turn those reformulations into an upper bound.

Queries on September5 included `harmonic multilinear polynomials slice Johnson scheme Hoeffding decomposition Filmus` and `Johnson graph Bernoulli Laplace eigenvalues`. Local searches covered Walsh, Hoeffding, ANOVA, noise operators, multislices, and Möbius inversion in the inspected Paley research and novelty notes. These searches are bounded. The exact workflow may be a useful local contribution; historical originality remains unestablished.
