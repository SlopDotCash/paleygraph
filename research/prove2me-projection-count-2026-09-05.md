# Finite counting lemma for the projection reduction

This note isolates the finite incidence lemma behind Section 2 of
the [prize projection argument](parallel25-prize-projection-2026-09-05.md).
Its formal proof does not formalize finite-field hyperplane counts,
quotient witnesses, code maxima, or the full Paley or prize statements.

## 1. Statement

Let S and P be finite sets, with a relation R between S and P. Let
q,t,K be natural numbers satisfying q≥2, t≥1, and K<q. Suppose

    |P|+1 = qt,
    #{b in P : R(a,b)} ≥ (q−1)t for every a in S,
    #{a in S : R(a,b)} ≤ K for every b in P.

Then |S|≤K. Empty S and K=0 are included, as are t=1 and K=q−1.
This is an elementary finite double-counting theorem. No claim of
literature novelty is made.

## 2. Proof

Counting the related pairs in the two orders yields

    |S|(q−1)t ≤ |P|K = (qt−1)K.

If |S|≥K+1, the lower bound would exceed the upper bound, because

    (K+1)(q−1)t − (qt−1)K
      = (q−1)t − K(t−1)
      ≥ (q−1)t − (q−1)(t−1)
      = q−1 > 0.

This contradiction gives the result. The strict positive gap is the
reason that averaging loses no integer count, even at K=q−1.

## 3. Relation to the unformalized code theorem

For r-row projections over a field of order q, take t=q^(r−1),
P the nonzero row projections, and S the bad scalars of a fixed pair.
Let R express survival as a bad scalar after projection. The written
code argument proves the three incidence hypotheses using a quotient
witness and a hyperplane count. The lemma then bounds |S| by the scalar
maximum K when K<q; K=q is separately trivial. Embedding supplies the
opposite inequality between maxima.

The formal theorem here takes those incidence hypotheses explicitly.
It does not assume the desired conclusion |S|≤K, and it does not
establish the code-theoretic hypotheses by itself. Its server status and
exact artifact hashes are recorded separately when verification completes.
