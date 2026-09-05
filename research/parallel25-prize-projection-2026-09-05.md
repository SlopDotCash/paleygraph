# Removing the interleaving loss from the pinned certificate

**Status: an exact reduction of the pinned combination-round certificate
to scalar-code quantities, with no extra assumption on received words.
The scalar quantities and the full prize goal remain unbounded here.**
The separate spot-check obligation is unchanged. This does not prove a
Paley-to-prize implication or replace the full goal by a smaller target.

The new argument within this project uses random linear projections of
the alphabet. It does not require that the received words lie in a
proper subfield or a low-dimensional alphabet subspace. Preservation
of list-decoding radius under interleaving is classical; see the
[primary Gopalan–Guruswami–Raghavendra abstract](https://arxiv.org/abs/0811.4395).
Only that provenance was checked, not a theorem from the paper. All
finite inequalities used below have self-contained proofs. No claim
of literature novelty is made.

## 1. Definitions and a projection probability

Let F be a finite field of order q, and C≤F^D a linear code. For an
integer r≥1, let C^r be its r-row interleaving, with column Hamming
distance. Write L_r(δ) for its largest list at radius δ, over every
received word in (F^r)^D. Write B_r(δ) for the largest number of
MCA-bad scalars for a pair W,Z of such words. A scalar γ is bad when
some allowed agreement set T has

    (W+γZ)|_T∈(C|_T)^r, but (W|_T,Z|_T)∉(C|_T)^(2r).

These are counts, not probabilities; the MCA error is B_r/q.
For fixed γ and T satisfying the first condition, the second condition
is equivalent to Z|_T∉(C|_T)^r, since W=(W+γZ)−γZ.

Choose a uniformly from F^r\{0}. Its projection on alphabet vectors
is π_a(v)=a·v. It maps C^r into C. For any nonzero v∈F^r,

    Pr[a·v=0]=α_r=(q^(r−1)−1)/(q^r−1)<1/q.         (1)

This includes α_1=0. Uniformity over nonzero a matters for the exact
integer conclusion below; including the zero map loses that endpoint.

## 2. MCA is exactly invariant under row interleaving

**Theorem.** For every r and every radius, with the same allowed
agreement subsets in the definition,

    B_r(δ)=B_1(δ).                                   (2)

Proof: fix W,Z and its bad-scalar set S. For every γ∈S, fix one
witness T_γ. The r restrictions of Z give a nonzero vector in the
rth power of the quotient F^(T_γ)/(C|_(T_γ)). Choose a linear
functional on this quotient which is nonzero on at least one entry.
It produces a nonzero vector v_γ∈F^r. If a·v_γ≠0, then
π_a(Z)|_(T_γ) is still outside C|_(T_γ). The projected folded
explanation still belongs to C|_(T_γ), so γ remains MCA-bad for
the single fixed projected pair (π_aW,π_aZ).

Thus each original bad scalar survives projection with probability at
least 1−α_r. Linearity of expectation gives

    (1−α_r)|S|≤E_a |Bad(π_aW,π_aZ)|≤B_1.             (3)

Let K=B_1. If K=q, the desired upper bound is trivial. If K≤q−1,

    |S|≤K/(1−α_r)
       =K+K(q^(r−1)−1)/[q^(r−1)(q−1)]<K+1.

Integrality yields |S|≤K. Taking the maximum over W,Z proves the
upper bound. Embedding scalar words in one row gives B_r≥B_1.
This proves (2), including the endpoint K=q−1, γ=0, and r=1.
The chosen quotient functionals may depend on the witness; the random
row projection a is the same across all scalars in the expectation.

No list-size estimate or Reed–Solomon-specific hypothesis is used in
this theorem. In particular, it applies directly over the full
extension field of the pinned profile, to arbitrary received pairs.

## 3. Lists over a large alphabet have no interleaving loss

**Theorem.** If L_1≤L and binom(L+1,2)<q, then L_r≤L for every r.
In particular, if binom(L_1+1,2)<q, then L_r=L_1.

Suppose a received r-row word has at least L+1 distinct explaining
codewords F_1,…,F_(L+1). Every projection π_aF_i is a scalar
codeword within the same radius of π_aW: an agreeing column stays
agreeing. For each pair i≠j, the nonzero matrix F_i−F_j has a
nonzero column. By (1), the probability that their entire projected
codewords coincide is at most α_r<1/q. A union bound gives

    Pr[some projected pair coincides]
       ≤binom(L+1,2)α_r<1.

Some projection therefore gives L+1 distinct codewords in a scalar
list, a contradiction. The reverse inequality L_r≥L_1 follows by
embedding. This is a statement about distinct codewords; it never
counts agreement subsets as separate explanations.

A second useful bound follows without the square-root-sized list
condition. If L_1≤L<q, then

    L_r≤floor[L(q−1)/(q−L)].                         (4)

For a particular list of size M, choose a uniformly from all F^r,
including zero. If its projected multiplicities are m_j, their sum
is M and there are at most L distinct images, so Σ_j m_j²≥M²/L.
Every distinct pair collides with probability at most 1/q, giving

    M²/L≤EΣ_jm_j²≤M+M(M−1)/q.

Cancel M>0 and rearrange to get (4). If M=0 it is immediate.
The bound contains no factor depending on the number of rows.

## 4. Exact scalar form of the pinned official certificate

Keep the archived September-4 parameters and contract versions from
[pass 23](parallel23-prize-bridge-2026-09-05.md):

    p=2130706433, q=p^6,
    D=μ_(262144)⊂F_p⊂F_q, n=262144, k=131072,
    C=RS[F_q,D,k], δ=j/n, 1≤j≤131071,
    R=floor(q/2^128)=274980728111395087.

The combination certificate is B_8(δ)+L_16(δ)≤R. The exact
integer computation gives

    q=93571093019388561295270373781649880353786165192103559169,
    R(R+1)<2q.

Consequently the entire combined condition is equivalent to

    B_8+L_16≤R  if and only if  B_1+L_1≤R.           (5)

One direction follows by scalar embedding. Conversely, B_1+L_1≤R
implies L_1≤R, so Section 3 gives L_16=L_1; Section 2 always gives
B_8=B_1. This proves (5). The list-preservation hypothesis is not
assumed independently: it follows from the proposed scalar certificate
itself and the verified integer budget.

Thus the eight- and sixteen-row maxima in the earlier exact remainder
description can both be replaced by single-row maxima at the same
extension field and radius, without a numerical loss in this
certificate. Neither single-row maximum has been bounded here.
The separate spot-check inequality and all official admissibility
requirements remain necessary; equation (5) proves only the equivalence
of the combination-round count inequality already isolated in pass 23.
No current website, leaderboard, or source-version claim is made.

## 5. Base-field list transfer and the distinct MCA obstacle

Because D⊂F_p, expanding each F_q coefficient in an F_p basis
identifies C^r exactly with the 6r-row interleaving of
C_0=RS[F_p,D,k], with the same column metric. If its scalar base-field
list maximum is bounded by A<p, equation (4) gives

    L_r(RS[F_q,D,k],δ)≤floor[A(p−1)/(p−A)]           (6)

for every r. If binom(A+1,2)<p, the stronger theorem gives exact
equality to the base-field list maximum. This treats every
extension-valued center, not just base-field centers. For instance,
a hypothetical base-field list bound A≤p−18 would give

    L_16≤252217214619120071<R.

That particular scalar bound is not proved here, and the leftover
budget must still cover the MCA term. The example merely quantifies
the sufficient input rather than asserting that it holds.

The same coordinate expansion does not make multiplication by an
arbitrary γ∈F_q into multiplication by an F_p scalar. The exact
MCA invariance in Section 2 applies at fixed scalar field; it does not
identify B_1 over F_q with B_1 over F_p. The earlier
[subfield-valued-word descent](parallel-prize-subfield-2026-09-04.md)
covered a restricted received-word class and cannot be applied to all
extension-valued words.

One exact statement for arbitrary words is available. Let K be the
base-field scalar MCA maximum. If K<p, then for every extension-valued
pair W,Z and every affine F_p-line ℓ=α+βF_p with β≠0,

    |Bad_(F_q)(W,Z)∩ℓ|≤K.                            (7)

To prove this, write γ=α+βt, t∈F_p. Expand W+αZ and βZ in an
F_p basis, obtaining a multirow base-field pair. Code membership on
each witness is preserved exactly by this expansion. Thus the bad
parameters t are exactly those of that base-field interleaved pair.
Section 2 bounds them by K (indeed the assertion also holds trivially
when K=p). This uses one fixed expansion for the whole line, not a
different projection for each t.

For K=0 there are no bad extension scalars, and for K=1 there is at
most one, since any two distinct points lie on an affine F_p-line.
For 2≤K<p, choosing a point of a nonempty bad set and partitioning
all its other points by lines through that point gives only

    |Bad_(F_q)(W,Z)|≤1+(K−1)(q−1)/(p−1).             (8)

This elementary line-count bound does not establish the required
security numerator. Additional structure of the actual remainder
ratios is still needed for an extension-field MCA estimate. No claim
that (8) is sharp for an actual Reed–Solomon bad set is made.

## 6. Verification and the remaining proof task

The [verifier](../experiments/parallel25_prize_projection_2026_09_05.py)
and [results](../results/parallel25_prize_projection_2026_09_05.json)
test actual finite-code lists, projected-codeword collisions, bad-scalar
sets and witness survival. They also verify the official integer
condition, the scalar/interleaved equivalence on small exact cases,
and extension-field affine-line restrictions. Coverage and its limits
are recorded with the current input hashes. The run checks 401,427
received-word pairs across its exhaustive scalar and two-row cases,
including the endpoint B_1=q−1 for the actual RS code over F_5 with
n=4,k=2 and agreement at least three. Its binary length-three case
is a linear repetition code with abstract coordinate labels, not an
RS evaluation at three distinct elements of F_2. Two-row maxima are
exhausted only for the two small repetition codes; the F_5 two-row
case uses 80 sampled pairs and every nonzero row projection.

List checks use 25,252 nonzero projections in three cases over F_5
and F_29. The scalar list maximum six there is certified by the
six information subsets and a quadratic received word attaining six;
it is not a claim of exhaustive received-word enumeration. The F_25
checks use 50 received pairs, all 25 challenges and allowed witnesses,
and all 30 affine F_5-lines. Finite examples check the implementation
and endpoints; the general statements follow from the proofs above.
No formal verification is claimed.

The next required estimate is now scalar and fully quantified:

    max_(f,g∈F_q^D)|Bad_C(f,g;δ)|
       +max_(w∈F_q^D)|List_C(w;δ)|≤R,

at the pinned field, domain, degree and radius, together with the
separate contract obligations. The reduction removes interleaving as
an additional loss. It supplies no general upper estimate for these
two maxima and does not establish the sought subgroup cancellation
or full classical Paley conjecture.
