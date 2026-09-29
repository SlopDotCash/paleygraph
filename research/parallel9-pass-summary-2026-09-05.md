# Ninth pass: individual necklaces at every fixed degree

**Progress: the finite word catalog has been replaced by a general argument.**
The [written proof](parallel9-all-degrees-2026-09-05.md) obtains

\[
 |\operatorname{tr}(D_{Z_1}S\cdots D_{Z_k}S)|
 \le 3a(2a+2)^{k-2}p^{(k+1)/2},\qquad k\ge3,
\]

for p≡1 mod 4, a<p distinct anchors, and arbitrary nonempty labels
Z_i drawn from those anchors. The order-one trace is zero and the
order-two trace has magnitude at most a²p. The estimate is uniform
in the anchor positions. Its proof has been checked locally against
the cited primary theorems; independent mathematical review and
formal verification have not been obtained.

For every fixed a and k, this gains a factor √p over the elementary
trace bound. It proves the fixed-degree necklace estimate stated as
[Conjecture 1.14 in Kunisky's paper](https://arxiv.org/html/2303.16475v1#S1.SS5).
The paper's Theorem 1.17 then supplies weak empirical spectral
convergence at every fixed degree. This is a substantial extension
of the preceding all-word result at length six; it does not prove
the minimum-eigenvalue conjecture or the full Paley conjecture.

The principal sheaf starts at rank one. Each nonempty label twist,
followed by quadratic middle convolution, increases rank by at
least one; the rank increments are nondecreasing. That statement
holds for every number of anchors. The final rank-one pairing
therefore has no invariant. The raw convolution also has boundary
corrections. An exact perverse-sheaf induction shows that these
retain a strict weight gap, while their total constituent mass
grows by at most 2a+2 per step. The proof keeps all finite masks,
new invariant stalks and infinity corrections.

The [verifier](../experiments/parallel9_all_degrees_2026_09_05.py)
checks 114,510 rank transitions and inverse local types, including
1,691,778 block identities. It also checks 1,170 original-matrix
endpoint identities, 980,530 trace-sign entries and 410 literal
necklace sums. The [results](../results/parallel9_all_degrees_2026_09_05.json)
and [artifact audit](../results/parallel9_pass_audit_2026_09_05.json)
pin final proof, script and primary-source hashes. The small-field
bound checks can be vacuous; their purpose is to test exact formulas,
not establish the asymptotic theorem. The uniform proof uses the
symbolic rank induction and the imported cohomology theorems.

The remaining spectral obstacle is now explicit. The coefficient
grows exponentially with length, and taking absolute values across
all words makes the normalized aggregate bound grow like
√p((2a+2)√(2^a−1))^k. This is insufficient for the
[extreme-eigenvalue criterion](parallel2-spectral-transfer-2026-09-04.md).
We need cancellation between words, a different aggregate argument,
or another way of controlling the extreme eigenvalues.

The subgroup bound remains M≤(17+log N)^(1/9)N^(8/9) on the
existing density-one class of quartic splitting primes. Its uniform
square-root target, exceptional primes, stronger centered energy
budget, full restricted operator, classical arbitrary-set bound and
exact official prize bridge remain open. No new result in those
directions is claimed in this pass.

Root completed this extension locally. A current agent-status check
confirms that all three parallel workers remain stopped by the
account usage limit. No reset, purchase, external submission or
formal certification was made. The full research goal stays active.
