# Formal projection preservation and its existing source

The [standalone Lean proof](../prove2me/Check_paley_mca_projection.lean)
now proves the code-specific MCA projection argument and exact equivalence
of uniform bounds on bad-scalar counts. This is independent verification
of a result already in the pinned ArkLib source. It is not a new MCA
invariance theorem or a bound on the remaining scalar problem.

## Provenance correction

The previously archived
[Errors.lean](../sources/parallel23-prize-dependencies/ProximityGap_Errors.lean)
contains `mcaError_interleaved_le` and `mcaError_interleaved_eq` near its end.
Root read their complete proofs and checked the
[primary file at the same commit](https://github.com/Verified-zkEVM/ArkLib/blob/e65197892890b8fd9b0dc05b8980273cf1d595cc/ArkLib/Data/CodingTheory/ProximityGap/Errors.lean#L1091).
For a module code, a positive number of rows, and radius strictly between
zero and one, ArkLib already proves exact equality of affine-line MCA
under row interleaving. Its alphabet may itself be a vector space.

The source selects one agreement witness for each bad scalar. It constructs
a proper subspace of row combinations for each such witness and avoids all
of them at once. There are at most |F| bad scalars, and at most |F| proper
subspaces cannot cover a nontrivial finite F-vector space. One projection
therefore preserves every originally bad scalar. Taking suprema gives the
upper inequality; the existing constant-row embedding gives the reverse.

Pass 25 independently derived the same MCA conclusion and disclaimed
literature novelty, but its assessment did not identify this exact prior
theorem. Treat its MCA portion as a rediscovery. The current central
assessment corrects that omission; historical proof and result bytes are
preserved. This source discovery does not invalidate the reduction.

Root also inspected the archived ListDecodability file and the two newly
archived files in the [source ledger](../results/parallel26_primary_sources_2026_09_05.json).
No matching finite collision-count list transfer was found in those
inspected files. This limited search establishes neither novelty nor
absence elsewhere. The pass-25 list argument remains an author-checked
ordinary proof with finite tests, without a new separate-author review or
Lean formalization.

## Exact standalone statement

Let F be a finite field, I a finite nonempty row set, V an F-vector space,
and (U_y) any indexed family of subspaces of V. Define a scalar γ to be
bad for row families W,Z if some y satisfies

    W_i + γ Z_i ∈ U_y for every i, and Z_i ∉ U_y for some i.

For scalar w,z, badness means w+γz ∈ U_y and z ∉ U_y for some y.
The theorem `projection_preserves_bad` constructs a nonzero a in F^I
such that every scalar bad for W,Z is bad for the projected pair

    w = Σ_i a_i W_i,   z = Σ_i a_i Z_i.

For each bad γ choose y and i with Z_i outside U_y. A separating linear
functional annihilates U_y but not Z_i. Its composition with row projection
is a nonzero linear functional on F^I. The previously checked common
projection lemma avoids the kernels of all these at most |F| forms.
Linear combinations preserve membership of the folded rows in U_y, while
the selected nonzero functional value prevents membership of projected z.
The empty bad set is handled separately using the nonempty row set.

`rowsBad_iff_input_failure` proves that, conditional on the folded rows
belonging to U_y, failure of some Z_i to belong is equivalent to failure
of some W_i or Z_i to belong. This uses subtraction in the subspace and
includes γ=0. The scalar analogue is also checked.

For every natural number K, `uniform_bad_count_bound_iff` proves

    (every row pair has at most K bad scalars)
      iff (every scalar pair has at most K bad scalars).

The forward implication repeats the scalar pair in every row. The reverse
injects the original bad-scalar subtype into that of the projected pair.
No restriction K<|F| is imposed on this final theorem.

## Codeword agreement and official definitions

For a linear code C contained in F^D and T contained in D, define U_T as
the words agreeing on T with some codeword in C. The file proves directly
that U_T is a subspace. For any family A of allowed coordinate sets, the
code badness predicate is

    ∃ T∈A: every folded row agrees on T with a codeword,
             and some input row has no such codeword.

The file proves its identification with the abstract predicate and then
proves `code_projection_preserves_bad`, `code_uniform_bad_count_bound_iff`,
and `codeRowsBad_iff_input_failure`. No finiteness assumption on D or A is
needed for these standalone statements; the row set and F are finite.

For the finite official domain, take A to contain precisely the sets with
|T| ≥ |D|(1−δ). The source's `projectedCodeSubmod` membership is existence
of a codeword with the given restricted word, exactly the EqOn witness
used here. The source lemma `projectedCodeSubmod_moduleInterleavedCode_iff`
identifies membership with membership of every row after transposition.
The affine-line generator has coefficients (1,γ), and failure of either
input is equivalent to failure of Z under the folded-membership condition.
Thus the definitions match as a mathematical interpretation.

This last comparison is a source audit, not a Lean theorem importing
ArkLib. The compiled file uses Mathlib only. The exact source theorem
`mcaError_interleaved_eq` was inspected, but its whole dependency closure
was not rebuilt or axiom-audited in this pass. Do not confuse these two
forms of evidence.

## Verification and remaining scope

The [Lean result](../results/parallel26_mca_projection_2026_09_05.json)
and [compiler log](../results/parallel26_mca_lean_2026_09_05.log) pin the
checked proof. Seven principal statements print only standard axioms:
`propext`, `Classical.choice`, and `Quot.sound`, or a subset. The file has
no `sorry`, `admit`, custom axiom, unsafe declaration, or imported open
project theorem. The checked environment is Lean 4.30.0 and Mathlib
c5ea00351c28e24afc9f0f84379aa41082b1188f. Reproduce using

```sh
python3 experiments/parallel26_lean_check_2026_09_05.py
```

The larger standalone result has no Prove2Me server verdict. The earlier
private incidence-lemma statement job remains pending in the
[server readback](../results/prove2me_projection_server_2026_09_05.json).
No duplicate submission was made.

This proves a reduction for arbitrary linear codes, not a numerical
upper bound on scalar MCA or scalar list size. It neither descends the
official extension field to its prime subfield nor establishes a Paley
implication. Both full conjectures remain unproved.
