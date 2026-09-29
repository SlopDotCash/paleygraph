# The algebraic projection core now checked in Lean

**The local Lean checks pass with only propext, Classical.choice, and
Quot.sound. They contain no unproved imported hypotheses or sorryAx.**
The checked statements below do not formalize codes, MCA witnesses,
the full prize reduction, or either conjecture. The private server
submission currently queued is the earlier finite incidence lemma,
not this whole file. No literature novelty is claimed for these
elementary linear-algebra and counting facts.

## 1. Formal statements and evidence

Let F be a finite field, V a finite F-vector space, and q=|F|.

For a nonzero linear functional f:V→F, the
[kernel and survival proof](../prove2me/Check_paley_functional_kernel_card.lean)
checks

    |V| = |ker(f)|q,
    for some integer t≥1, |V|=qt and #{v : f(v)≠0}=(q−1)t.

Surjectivity follows by scaling a vector on which f is nonzero. The
quotient-by-kernel equivalence gives the first cardinal identity;
partitioning V into zero and nonzero fibers gives the second. The
first identity in the Lean file is stated more generally with Nat.card;
the finite hypothesis is explicit for the positive survival count.

The [combined proof](../prove2me/Check_paley_nonzero_functional_count.lean)
then checks two further statements. Let S be a finite index set and
f_b:V→F a nonzero linear functional for every b∈S. Repeated functionals
at different indices are allowed and counted separately.

1. If K<q and every nonzero a∈V satisfies
   #{b∈S:f_b(a)≠0}≤K, then |S|≤K.
2. If S is nonempty and |S|≤q, then there exists a nonzero a∈V
   such that f_b(a)≠0 for every b∈S.

For the first statement, every functional has the same positive
kernel cardinality t, since q times this cardinality equals |V|.
Apply the finite incidence lemma to P=V\{0}: its size is qt−1,
and every b survives on exactly (q−1)t projections.

For the second statement, if no common vector existed then each
nonzero a would detect at most |S|−1 functionals. The first statement
with K=|S|−1<q would force |S|≤|S|−1, a contradiction. The
nonempty hypothesis avoids demanding a nonzero vector in the zero
vector space when there are no functionals. No restriction that the
functionals themselves be distinct is needed.

The [kernel result record](../results/prove2me_functional_counts_local_2026_09_05.json)
and [core result record](../results/prove2me_projection_core_local_2026_09_05.json)
pin the proof bytes and the reported axiom lists. The files were
compiled under Lean 4.30.0 and Mathlib c5ea00351c28e24afc9f0f84379aa41082b1188f.
The proof imports only Mathlib and Lean modules; it imports no platform
theorem or target statement containing a placeholder.

## 2. Simultaneous preservation in the written MCA argument

This gives a shorter, stronger form of the averaging step in
[pass 25, Section 2](parallel25-prize-projection-2026-09-05.md).
For a fixed interleaved pair, let S⊆F be its bad scalars. If S is
nonempty, choose one allowed agreement witness T_gamma for each
gamma∈S. The restricted direction is nonzero in the corresponding
quotient by the restricted scalar code. A linear functional on that
quotient, nonzero on one row, gives a nonzero linear functional
l_gamma on the row-projection space F^r.

There are at most q such functionals, one per bad scalar. The common
projection theorem supplies one nonzero row projection a with
l_gamma(a)≠0 for every gamma∈S. At each witness, the projected
folded word is still in the restricted code, and the projected direction
is still outside it. Thus every original bad scalar remains bad for
this single projected pair. Extra bad scalars created by projection do
not invalidate the inclusion.

If S is empty the required inclusion is automatic. Scalar embedding
supplies the reverse inequality between maxima, recovering B_r=B_1.
This argument works over the same scalar field and does not descend
the official extension field to its prime subfield.

**The code-theoretic application in this section is a written argument,
not part of the Lean declarations above.** The next formalization step
is to construct the quotient witnesses, prove survival under the common
projection, and then define and compare the actual MCA maxima. The
full root-authored prize reduction still lacks a separate-author review.

Similarly, for a fixed list of M codewords, one nonzero column of each
pairwise difference supplies a nonzero functional. A common projection
preserves distinctness whenever binom(M,2)≤q. This permits equality in
the sufficient union bound; the earlier strict inequality was sufficient
and remains correct. The archived official budget already satisfies the
strict inequality, so its numerical certificate does not change.

## 3. What remains

The local formalization checks the algebra supporting a reduction.
It gives no upper bound for the official scalar MCA/list maxima, no
new worst-case Paley exponent, and no proof of the subgroup or spectral
targets. Server queue state is recorded in
[the project setup](../PROVE2ME.md); queued compilation is never
reported as an accepted server proof.
