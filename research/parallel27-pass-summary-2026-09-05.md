# Twenty-seventh research assessment

**The full Paley conjecture and the Proximity Prize remain unproved.**
This pass improves exact access to the subgroup remainder that the
previous estimates leave open. It proves a new-to-this-project disjointness
fact and derives a faster exact counting formula. It does not prove a
uniform upper bound on that remainder or improve a cancellation exponent.

## Completed mathematical and computational work

A six-term word with distinct nonzero entries and no opposite pair can
be product-balanced at at most one of its ten three-versus-three splits.
Two such splits would force an equal or opposite pair. The algebraic
core passes Lean with only standard axioms. This makes the sum of the
ten distinct-balanced counts exact, without a union-bound overcount.

Marking a repeated coordinate pair and dividing by its number of possible
marks also gives an exact repeated-word count. Combining these two facts
with the existing coset formula for E3 yields an exact sparse formula for
D6, the six-distinct, opposite-free, fully unbalanced remainder. With the
subgroup supplied, the formula costs O(n^(49/20)polylog n) field and
comparison operations under the previously imported energy estimates,
and O(n²) storage. This is a cost bound for computing D6, not a bound on
its value. No claim of literature novelty is made.

The [proof and formulas](parallel27-six-distinct-2026-09-05.md) include
all normalization factors, the zero coset, the exact repeated-coordinate
weights, and a 720n divisibility check from free scaling orbits. The
initial implementation used the wrong power to label cosets; direct
checks caught it. The corrected final formula uses z^n, whose kernel is H.

The [sparse results](../results/parallel27_six_distinct_2026_09_05.json)
cover 17 actual subgroups. Six new quartic-window samples at n=512 and
1024 have E3=T6+D6: all nonintrinsic six-words lie in the remaining
distinct unbalanced family. Their D6 counts range from zero to 5,898,240.
They are below the cubic scale in these selected cases. This does not
bound their order classes uniformly, identify the worst prime, or justify
an asymptotic extrapolation.

A separate C++ program directly enumerates all normalized six-word
classes in each of the same 17 fields. It checks each candidate's
opposite pairs, repetitions, and ten product splits instead of using
the sparse formula. Every category agrees exactly. The
[verification record](../results/parallel27_independent_check_2026_09_05.json)
also pins the three Lean statements. Both programs are root-authored:
independent implementation is not separate-author or human review.
The full combinatorial formula remains an ordinary mathematical proof.

## Remaining quantitative gap

No new uniform estimate on D6 follows from its exact isolation. The
existing upper bounds on repeated and balanced words cannot control this
remaining family. Even a cubic sixth-energy estimate would be only a
partial period bound; it would not establish uniform square-root
cancellation at logarithmic moment depth. The classical signed-moment,
full spectral, and scalar MCA/list estimates are unchanged.

The primary-source check retained the same applicable MRSS additive
energy and Shkredov shifted multiplicative-energy inputs. A limited
literature search identified no additional bound to import into this
argument. This is not a claim that the literature is exhausted or that
the recorded estimates are optimal. The
[source scope record](../results/parallel27_source_scope_2026_09_05.json)
states exactly what was used.

The previous goal turn is classified as progress: verified formalization
and corrected attribution. This turn adds a proved disjointness fact,
exact formula, and independently checked larger finite cases. The full
objective remains active and unachieved. The private Prove2Me incidence
statement job finished compiling and now has theorem ID
d025929b-75ed-4c1c-b661-8e491017b153. Its locally checked proof was
submitted once, with submission ID 27ff60e7-7d6d-48da-8c34-7b123c6c0a3f;
server verification returned ACCEPTED and the theorem readback is Proved
in the pinned environment. The
[current readback](../results/prove2me_projection_server_2026_09_05.json)
records both. There are now five accepted elementary proofs, including
this finite incidence lemma. Neither the code-specific MCA reduction
nor either full conjecture has received a new hosted proof verdict.
No duplicate statement job was created. No worker completion, reset, credit
purchase, public release, mission launch, or memory write is claimed.

Earlier proof artifacts and all 36 source-manifest entries are retained.
The [final audit](../results/parallel27_pass_audit_2026_09_05.json) pins the
current results and records the limitations above. The next mathematical
step still requires a bound on actual arithmetic quantities, not further
subtraction of known contributions from E3.
