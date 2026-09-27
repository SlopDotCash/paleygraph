# Pass 49: an exact content valuation and its prime support

The Paley graph conjecture and Proximity Prize remain unproved. This
pass clarifies the arithmetic input from pass 48; it gives no improvement
to the uniform bound at the target primes.

## Proved reductions

Let Ccrit_n be the gcd of all coefficients of the existing polynomial
Res(R_n,R_n'+T R_n''). It divides the earlier two-resultant gcd G_n.
For every odd prime, Ccrit_n is divisible by that prime exactly when
R_n has a root of multiplicity at least three. Using coefficients removes
the small-characteristic restriction needed by the earlier evaluation
test.

Root derives the exact formula

    v_p(Ccrit_n)=sum_i min(v_p(R_n'(lambda_i)),
                          v_p(R_n''(lambda_i))).

It dominates the sum of the triple masses over all p-adic precisions.
The difference is an explicit nonnegative correction from cancellation
among reciprocal root differences. A separate integer-polynomial
example has valuation 4 but total triple mass 3, so equality cannot
be assumed as a general algebraic rule. Another example has content
16 and two-resultant gcd 400, showing that the refinement can be strict.
These two examples are not subgroup data.

An ordinary norm argument also proves

    p divides Ccrit_n, p odd  implies  p=1 mod n.

Outside this splitting class, either the field of definition doubles
and the mixed representation is unique, or two norm conditions restrict
it to the roots of a quadratic. Primitive fibers there have size at most
two. The target primes already satisfy p=1 mod n, so this support
restriction does not exclude them or improve their estimate.

The sufficient tower criterion may now use v_p(Ccrit_s) in place of
v_p(G_s). Its required uniform bound remains unproved. The full
[derivation](parallel49-critical-content-2026-09-06.md) keeps the
higher-precision and cancellation costs explicit.

## Exact checks and next action

The [checker](../experiments/parallel49_critical_content.py) and
[results](../results/parallel49_critical_content_2026_09_06.json) verify:

- Five full resultant polynomials through order 128, using 82 exact
  resultant evaluations including checks beyond the interpolation points.
- Seventy p-adic cases, with verified precision up to 16, matching the
  exact derivative-minimum formula and its correction.
- Twenty-one nonsplitting fields, including extensions of degrees four
  and eight, with 10,920 direct balanced pairs and 9,984 norm-quadratic
  checks.
- Both auxiliary polynomial examples with exact integer coefficients.

In all five checked subgroup orders, Ccrit_n equals G_n. In all seventy
checked splitting cases, higher triple masses and derivative corrections
vanish. These finite facts give neither a numerical improvement here
nor a uniform vanishing theorem.

Root completed ordinary proofs and exact checks. No separate-agent,
Lean, external review or novelty claim is made. The previously observed
usage-limited lanes were not restarted. No Lean process was polled,
started or changed, and no Prove2Me submission was made. The source and
artifact checks are recorded in the
[audit](../results/parallel49_pass_audit_2026_09_06.json).

The next mathematical task is a uniform estimate at splitting primes,
on the actual triple aggregate or a sufficient arithmetic certificate.
Uniform W86/15, energy49/20, period71/72, absolute prime exception5/3,
and all full targets remain unchanged. The goal remains active.
