# Pass45: full weighted matrix constraints

The full Paley conjecture and Proximity Prize remain unproved. This pass
extends the exact matrix algebra and checks two source comparisons.
It does not improve the current uniform triangle, energy, or period bound.

For arbitrary real weights u,v, the
[weighted multiplication law](parallel45-weighted-matrix-algebra-2026-09-06.md)
is

    L(u)L(v)=n<u,v>I-nuv^T+L(L(u)v).

It gives a rank-two commutator and exact cubic and quartic trace formulas.
The augmented field-convolution matrices commute; the compressed matrices
generally do not. For u=a^2, the cubic trace retains exactly the weighted
triangle W. This is a representation of the remaining quantity, not a
bound on it.

An exact orthogonal projection onto the span of the shifted incidence
matrices isolates the mean, centered correlation, and unused squared norm.
Inserting only the current scalar moment bounds gives an error of order
n^7 in the quartic range, weaker than the existing n^(86/15) triangle
estimate. The full multiplication constraints remain available for more
targeted work; the failure of that particular norm insertion is not an
impossibility theorem for this approach.

The [source comparison](parallel45-shifted-product-scope-2026-09-06.md)
shows that Warren's published shifted-product consequence is compatible
with the surviving concentration parameters. A specific cubic product
inequality would extend the mass cap to low-doubling levels, but that
inequality is unproved. The 2026 algebraic-group result checked here is
over C; a uniform finite-field transfer is not supplied.

Finite checks cover23,000 weighted multiplication entries,40 trace
identities,20 rational orthogonal projections, and20 direct field
convolutions. All passed. Root derived and reviewed the ordinary arguments;
no separate-agent, Lean, or external peer review is claimed. The three
previous agent lanes were not restarted after their usage-limit failures.

The preceding performance turn yielded new live process evidence and
prepared unvalidated Lean branch candidates. This mathematical pass did
not start, change, or poll Lean and did not submit to Prove2Me. It verified
the saved pass44 state before proceeding. The current bounds remain
W86/15, energy49/20, period71/72, and quartic absolute exception count5/3.
The next useful task is to exploit the full weighted multiplication law
against actual high-level incidence concentration, or prove a valid
shifted-product input for those intersections. The full goal stays active.
