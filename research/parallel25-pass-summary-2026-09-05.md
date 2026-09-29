# Twenty-fifth research assessment

**The full Paley conjecture and the Proximity Prize remain unproved.**
This pass makes a useful exact reduction on the pinned prize problem
and bounds a further part of the subgroup energy. It does not improve
the worst-case cancellation exponent. There is no justified percentage
complete or estimate of the work left: the remaining uniform estimates
can still contain the central difficulty.

## What is closer

The [prize projection argument](parallel25-prize-projection-2026-09-05.md)
removes the row count from the pinned combination-round certificate.
For any linear code over F_q, its maximum number of MCA-bad scalars
is exactly invariant under row interleaving. List maxima are also
invariant when the scalar list bound L satisfies binom(L+1,2)<q.
A random nonzero row projection proves both statements. The strict
projection probability and integer rounding retain the MCA endpoint
q−1; arbitrary received words are included.

For the previously archived official profile, q=2130706433^6 and
R=floor(q/2^128)=274980728111395087 satisfy binom(R+1,2)<q.
Consequently, at each allowed radius, its count inequality is exactly

    B8 + L16 ≤ R  if and only if  B1 + L1 ≤ R.

All quantities here use the same extension field. Neither scalar maximum
is bounded in this pass. The separate spot-check obligation remains.
This is a reduction of the archived combination-round condition, not
a current prize-status claim or a proved equivalence between Paley and
the grand prize. The argument has an author audit and finite verification;
its separate-author review was interrupted and remains outstanding.

The [subgroup argument](parallel25-subgroup-growing-orders-2026-09-05.md)
now works uniformly in growing dyadic orders in the quartic window.
It bounds all opposite-free six-term relations having repeated entries,
aggregated over every product-ratio fiber, by

    R6,rep ≤ 15n[r4(2)−6n+8] ≤ 15n E2(H)
           ≪ n^(69/20)(1+log n)^(1/5).

Combining this with earlier bounds leaves

    E3(H) = T6 + D6
            + O(n^(69/20)(1+log n)^(1/5) + n^3(1+log n)),

where the error before applying O-notation is nonnegative and D6 counts
six distinct, opposite-free entries unbalanced at all ten triple splits.
The exponent 69/20 is still above the cubic target scale. D6 is
uncontrolled, so even bounding D6 alone would not automatically provide
the full square-root target. Actual order-256 data have no other
nonintrinsic contribution and D6=368640. This identifies a surviving
family; it is not a counterexample to any asymptotic conjecture.

## What the other approaches resolve

The [signed inversion argument](parallel25-signed-inversion-2026-09-05.md)
derives exact joint identities for both sign halves and their averages.
At the thin slice their averaged difference is negligible at the target
scale, but their common size still contains the unknown original sixth
moment. Repeated fractional-linear changes preserve an augmented signed
moment. They do not independently randomize the transported signs or
supply a contraction strong enough to close the uniform bound.

The [spectral argument](parallel25-spectral-full-operator-2026-09-05.md)
decomposes the full two-anchor operator under its S3 symmetries and
retains all coupling and boundary terms. Starting from the uniform
vector and iterating at any depth can visit only a sector of asymptotic
dimension m/6. About five-sixths of the space remains outside this
procedure. On the other sectors the desired edge reduces to a specific
upper bound on an elliptic conjugacy average, improving its elementary
constant from 1 to 2/3. That saving is unproved. On the trivial sector
the nonvanishing rank-two border must also be controlled. The earlier
logarithmic depth allowance remains insufficient for the full edge.

The subgroup norm calculation also rules out one proposed amplification:
all rationally independent cyclic shifts of an actual relation can
share a single prime factor in its norm. Their rational independence
does not imply independent congruence restrictions.

## Verification and parallel execution

Root worked on the prize reduction while three workers pursued signed
inversion, the full spectral operator, and growing subgroup orders.
These four lanes cover the current distinct approaches; they are not
a claim to exhaust every possible mathematical method. All three
workers then ended with usage-limit errors. Their substantive files were
saved. Root finished the spectral verifier run and independently reviewed
the three worker proofs using a newly written verifier.

All four author verifiers and the independent root verifier pass.
The [root review](parallel25-root-review-2026-09-05.md) describes the
checked algebra, source hypotheses, and finite cases. In particular:

- Prize: 401427 exhaustive tiny-code word pairs, 25252 nonzero list
  projections, exact extension-field affine-line checks, and exact
  official budget arithmetic. These are finite checks of the written
  reductions, not computations of the official maxima.
- Signed inversion: 666 input sets, 15140 actual sign-half identities,
  and full missing-row and composed-weight checks.
- Subgroup: four actual quartic cases at orders 64,128,256, 120585
  normalized multiset matches, and independently recomputed norm evidence.
- Spectral: 24 primes in the author verifier; root separately checks
  the full matrix, sectors, and both leakage coefficients at four primes,
  including p=269 beyond the author's range.

No human referee or new formal Lean proof is claimed. The root-authored
prize reduction still needs separate-author review. Worker usage limits
are an execution interruption; the full research goal remains active.
No model change, reset, credit purchase, external submission, automation,
or memory update was performed.

The [final artifact audit](../results/parallel25_pass_audit_2026_09_05.json)
checks current hashes, prior proof preservation, syntax, links, and source
archives. The source manifest retains all 32 earlier entries and adds
two primary HTML archives. MRSS Theorem 3 supplies the subgroup E2
bound; the Gopalan–Guruswami–Raghavendra abstract is only provenance for
classical interleaving work. No theorem is imported from that abstract.

## Inputs still needed

1. A uniform signed moment or exceptional-input estimate sufficient for
   the full classical two-set conjecture.
2. Stronger subgroup energy control, including six distinct fully
   unbalanced relations and exceptional primes, sufficient for uniform
   square-root cancellation.
3. The elliptic operator saving in all sectors, retaining the border,
   or a sufficient longer-depth bound for the actual full spectrum.
4. Bounds on the scalar MCA and list maxima over the full official
   extension field meeting the integer budget, plus the spot check.

The useful next mathematical step is a bound on one of these actual
quantities. The identities and finite experiments in this pass do not
already supply such a bound.
