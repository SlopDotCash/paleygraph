# Coincident edge cosets in the weighted triangle sum

For an odd prime and a subgroup H=-H of order n, put
f(x)=r_(H-H)(x) and Fk*=sum_(x!=0) f(x)^k. These are moments of
difference multiplicities, not the project's equal-sum energies T_k.
Consider the nondegenerate weighted triangle sum

    W = sum_(alpha,beta!=0, alpha!=beta)
          f(alpha)^2 f(beta)^2 f(alpha-beta)^2.

Let P be its part where beta/alpha belongs to H. On G=F_p^*/H,
the function f is constant on each coset. Fixing d in a coset D,
the maps h -> (d/(1-h), hd/(1-h)) and (a,b) -> b/a prove

    #{h in H\{1}: 1-h in D} = f(D).

Thus, with all quotient sums over G,

    P = n sum_(A,D) f(A)^4 f(D) f(DA)^2
      <= n (sum_A f(A)^4)(sum_D f(D)^3)
      = F4* F3*/n.

The inequality is Holder with exponents 3 and 3/2, followed by the
bijection D -> DA. Permuting the triangle vertices preserves the
weights and permutes the three edge cosets, since -1 belongs to H.
Each of the three pairwise coincidence regions therefore has weight P.

There is also an exact inclusion-exclusion identity. Put kappa=f(1).
All three edge cosets coincide precisely when beta=h alpha with both
h and 1-h in H. There are kappa such h, independent of alpha, and
the weight is f(alpha)^6. The triple-coincidence weight is kappa F6*.
Any intersection of two pairwise coincidence events is this same event.
Consequently

    W_coin = 3P - 2 kappa F6* <= 3 F3* F4*/n.                 (1)

The scoped published inputs F3* << n^3 log n and F4* << n^(11/3)
give W_coin << n^(17/3) log n when n << p^(2/3). See the
[independent review](parallel39-independent-coset-coincidence-review-2026-09-06.md)
for the exact source page, zero-term subtraction, and conditional
energy implication. The contribution W_distinct, with three distinct
edge cosets, still requires a new estimate. General Young gives only
W <= (F3*)^2 << n^6(log n)^2.

Distinctness alone does not restore the incidence hypothesis audited
in pass38. For p=97, H of order8 inside S of order24, both positive
normalized pairs have three distinct edge cosets and multiplicity3.
Dropping that multiplicity still changes the incidence count48 to16.
This S is a general invariant set, not a claimed difference level set.

The [exact finite checker](../experiments/parallel39_edge_cosets.py)
partitions W directly in seven cases, up to p=6700417,n=64, and checks
the quotient identity, (1), the Holder bound and Young bound. Its
[results](../results/parallel39_edge_cosets_2026_09_06.json) are finite
implementation checks. The universal identity follows from the proof
above; no full energy improvement, uniform period improvement, Lean
verification, or literature novelty is claimed.
