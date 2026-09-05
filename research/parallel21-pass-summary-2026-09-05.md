# Twenty-first pass: a sharper moment target and a weighted obstruction

We have isolated the missing character estimate more precisely.
For every fixed moment order 2r, only a one-sided upper bound on
the aggregate of the squarefree 2r-shift correlations is needed.
Its negative side is already bounded at the Gaussian scale.
**The required upper bound remains unproved.**

The [derivation](parallel21-squarefree-moments-2026-09-05.md) uses
the elementary symmetric polynomial of each row's signs. Its
three-term recurrence is the characteristic-polynomial recurrence
of a symmetric tridiagonal matrix. The resulting real-root bounds give

    T_(2r) ≥ −(8r)^r p n^r/(2r)!,
    M_(2r) ≤ 2^r(2r)! T_(2r)+2(16r)^r p n^r,

where T_(2r) sums the complete character products over all distinct
2r-subsets of the input set. A moment bound also bounds |T_(2r)|.
Thus the one-sided aggregate bound and moment bound are equivalent
up to fixed-order constants. This applies to the restricted size
slice established in pass 20. The proof retains the character's
zero rows, and an exact sixth-moment identity makes their corrections
explicit.

A separate rationally weighted sign model shows the limit of several
available constraints. For every sufficiently large prime, at exactly
k=⌊p^(1/4)⌋, the model has the correct means, zero masses and Gram
matrix pI−J. It satisfies the usual individual Weil-shaped bounds
through degree six and M_4≤3pk², yet M_6/(pk³)≥k/2→∞. Its
columns can be labeled by an actual additive Sidon set.

The qualification is essential: these are weighted sign rows, not
uniform prime-field character translates. A zero-column fibre has
mass one and can contain many rows. Sidon labels do not impose the
missing field-difference identities. The construction rules out a
deduction from the listed relaxed constraints alone; it does not
refute the actual character-moment conjecture or Paley.

The [verifier](../experiments/parallel21_squarefree_moments_2026_09_05.py)
and [results](../results/parallel21_squarefree_moments_2026_09_05.json)
passed exact arithmetic checks:

- 67,080 recurrence/binomial comparisons and 201,240 pointwise inequalities;
- 22,996 actual small-field aggregate identities and 5,749 exact
  sixth-moment ledgers;
- direct enumeration of 1,250 weighted rows, checking 190 Gram entries
  and 486 distinct-subset correlations;
- five larger compressed weighted fixtures with k=256,257,384,512,1024,
  plus 796,097 explicit Sidon pair-sum checks.

The large fixtures do not evaluate actual character kernels over
their fields. Finite checks support the algebra and implementation;
the written proofs cover the general deductions, and neither proves
the missing asymptotic character estimate.

The [final audit](../results/parallel21_pass_audit_2026_09_05.json)
pins the artifacts and preserves all five pass-20 artifacts and all
25 preceding source-manifest entries. One primary HTML source was
added to verify the standard distinct-root Weil comparison bound:
[McDonald–Sahay–Wyman, Lemma 2.1](https://arxiv.org/html/2210.03789v2).
No VC-dimension theorem is imported. No PDF was newly used or
visually reviewed.

The next task is to prove the positive aggregate estimate using
actual character-translate structure. There is no new Paley
cancellation exponent, clique bound, subgroup square-root bound or
official prize implication from this pass.

A [separate-agent proof review](parallel21-independent-review-2026-09-05.md)
found no mathematical correction in the main inequalities, sixth-moment
ledger or weighted-model constants. It inspected the verifier without
rerunning it. This completed review supersedes the initial pending-review
annotations in the pinned proof snapshot. Independent human review and
formal verification remain outstanding.

The full goal remains active and unachieved. Root completed the
original deductions and exact verification. Three fresh bounded lanes
are now investigating the classical positive aggregate, the subgroup
exceptional-prime input and the spectral estimate. The reviewer moved
to the spectral lane after finishing its report. Their further results
will be assessed separately as they complete.

