# Pass 47: anchored rows and arithmetic profile recognition

The Paley graph conjecture and the Proximity Prize remain unproved.
This pass gives stronger tests for candidate subgroup data, but no new
uniform estimate. The best weighted-triangle power remains86/15, still
1/15 above the intermediate17/3 target for this route. Even that target
would not by itself prove the full goals.

## Results

Known matrix entries give

    h_i>=a_i(a_i-1)+a_(-i)(a_(-i)-1),  i!=0,
    h_0>=a_0(a_0-1).

They exclude both earlier uncentered interval families directly, without
an excess estimate. The pass46 row-capacity inequality remains valid;
its total-X consequences alone had not excluded the long-filler example.

Centering the high interval repairs the anchored-row defect. A choice of
the single odd label also enforces the necessary product congruence.
This repaired abstract family retains ||T||_(3/2)>=n^(61/15), obeys the
compared scalar bounds and every total-X subset inequality at large n,
but has no supplied local excess allocations, incidence matrix, or field
realization.

For an actual subgroup of even order n, with quotient labels q^i, root
derived the exact field identities

    sum_i a_i q^(ki)=n sum_(j=0,...,k) binom(nk,nj) mod p,
    k>=1.

If a proposed nonnegative integer profile has mass n-1, these identities
for k=1,...,n-1 characterize the actual profile. Newton's recurrence
recovers its monic root polynomial, and distinct q^i have unique root
multiplicities. This is a recognition theorem, not a bound on those
multiplicities, energy, the weighted triangle, or the spectral edge.

The full derivations and model scope are in the
[mathematical note](parallel47-anchored-arithmetic-2026-09-06.md).

## Exact checks

The [checker](../experiments/parallel47_anchored_arithmetic.py) passed:

- Short exact primality certificates for m=1094909953 and
  p=1121187791873=1024m+1, in the quartic interval1024^4<p<2*1024^4.
- All2,048 positions of the repaired family in the specified R range
  fail the first power-sum condition at this one n,p. The first three
  geometric evaluations also agree with term-by-term sums.
- The actual order1024 subgroup has a0=0, maximum coefficient2 and
  A2=2045. Its product identity, first three moments, and all1,024
  support-relevant anchored-index checks pass. Other nonzero indices
  have zero right side in the anchored inequality.
- All53 power sums and all three full Newton reconstructions agree with
  direct root products at (p,n)=(97,8),(353,16),(278177,32).
- The first repaired abstract profile passes257 total-X threshold
  checks and the listed anchored, product, mass and norm checks.

The finite sieve is exhaustive only over its stated family, range and
one field. In particular, the actual subgroup's a0=0 already distinguishes
it from these high-centered profiles. No uniform rejection, density
statement, or performance advantage over enumerating H is asserted.
The [results](../results/parallel47_anchored_arithmetic_2026_09_06.json)
contain the exact residues, hashes, coefficients and counts.

## Review and next action

Root completed ordinary proofs and exact checks. The previously observed
parallel lanes had errored at the usage limit; no separate-agent, Lean,
external review or novelty claim is made. No Lean process was polled,
started or changed, and no Prove2Me submission was made in this pass.
Earlier noncentral artifacts and source hashes are checked in the
[pass audit](../results/parallel47_pass_audit_2026_09_06.json).

The next mathematical task is to obtain a uniform coefficient or excess
bound from the arithmetic constraints or the full matrix multiplication
law. Recognition alone leaves that task open. Uniform W86/15,
energy49/20, period71/72, the quartic absolute exception exponent5/3,
and the full goals remain unchanged. The goal remains active.
