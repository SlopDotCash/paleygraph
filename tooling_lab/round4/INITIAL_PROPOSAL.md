# Next bounded prototype: exact moments with arithmetic marks retained

The completed round3 counterfactuals rule out interpreting unmarked global L2 data as a complete distributional description. The next experiment should preserve an explicit arithmetic mark instead of merely recomputing that global signature more accurately. Pointed/coloured Johnson ideas already occur in the local July11 record; this is an implementation proposal, with historical novelty unestablished.

## Proposed object and computational mechanism

Fix a small ordered set M of field columns, initially three points. Compute exact mean and second moment of T6 over uniformly chosen n-sets C **containing M**. Retain the mark's arithmetic identities in the output. This is a conditional finite diagnostic, not a worst-case bound or a predictor that avoids the target.

Partition row indices x by the sign vector `(S_(x,m):m in M)`. There are at most `3^|M|` cells, usually far fewer. For each ordered pair of cells U,V, retain counts of row pairs with mutual sign0,+1,-1. These coloured edge counts retain information absent from the unmarked pair-type table.

For each cell pair and mutual sign, the full paired-row type histogram is known from the conference identities. Subtract the explicitly known pairs at the marked columns to get the unmarked type multiplicities. The contribution from each marked column is fixed; only the other n−|M| selected columns are sampled.

For one row pair with entries a_y,b_y, write

```
G_y(t,u) = (1+a_y*t)(1+b_y*u).
```

For a fixed C, `[t^6 u^6] product_(y in C) G_y` is the product of that row pair's sixth elementary coefficients. Marked columns contribute their G factors deterministically. For unmarked columns, expand

```
product_(y outside M) [1+v*(G_y(t,u)-1)].
```

Replace each v^k coefficient by `(n-|M|)_k/(q-|M|)_k`, the exact inclusion probability for its k distinct columns. Truncation at degrees6 in t,u and12 in v suffices. Group equal `(a,b)` types and use multinomial powers, as in round2. Sum over coloured cell-pair multiplicities. A corresponding one-row formula supplies the conditional mean.

The expensive new input is the coloured cell-pair count, not the bounded-degree coefficient calculation. For a Paley matrix it can be obtained through exact character convolutions of cell indicators:

```
sum_(x in U,y in V) chi(x-y).
```

Together with cell sizes and diagonal intersections, this determines the +/− counts. One convolution per cell can supply all cell-pair counts after aggregation. Use direct integer sums for the first prototype; use the existing exact NTT backend for a scaled version, with explicit transform-length and coefficient-recovery guards. Do not silently reuse the present NTT prime beyond its supported length.

## Acceptance and rejection tests

1. Verify the formula by direct enumeration on small prime fields and on the Paley49/Peisert49 pair. Condition on the same actual marked columns, not an averaged orbit chosen after observing the target.
2. One or two fixed marks should leave affine-invariant Paley mean/second-moment averages unchanged. This is a required null control, not a new positive result. A three-point mark is the first intended nontrivial test.
3. Verify every coloured count directly before replacing it with convolution. Preserve sign, order, zero entries, and whether the mark lies in a special subfield or multiplicative configuration.
4. Test whether the marked statistic distinguishes any of the complete twin histogram features that the unmarked compiler loses. Report failure if it does not; do not infer improved worst-case information from greater complexity alone.
5. Scale a fixed mark to p=1297 and p=65537 with exact arithmetic. State cell count, integer bit cost, transform count, and whether memory grows with the field or the number of input sets.

The meaningful outcome would be a new exact arithmetic-conditioned computation that was previously infeasible, plus evidence identifying which marks expose useful variation. It would still leave the uniform exceptional-set problem open. A calculation that merely restates T6 pointwise, drops the mark during averaging, or violates field/domain constraints does not meet the objective.

## Independent remaining directions

- Coding: certify whole interpolation-track batches or subtract covered basis families symbolically; retain all codewords, supports and overlap. Measure the remaining fixed-rate cost rather than extrapolating the current small-rate speedups.
- Spectral: a new phase feature must use alignment with the fixed localized operator and survive S3-preserving controls. Global chirp-invariant summaries alone fail the sufficiency test already saved in round3. Do not disguise the Rayleigh quotient itself as a new explanatory feature.
