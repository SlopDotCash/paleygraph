# Prove2Me connection for this project

This project reuses `/Users/shawwalters/prove2me_workspace` for Lean and
the Prove2Me API. The authenticated server check on 2026-09-05 confirms
the supported environment below. The service updated from version 0.9.6
to 0.9.7 during setup; the skill and four relevant API references were
refreshed from the official site before proceeding:

- Lean `leanprover/lean4:v4.30.0`.
- Mathlib `c5ea00351c28e24afc9f0f84379aa41082b1188f`.
- Existing Paley proposal `b5e7121a-eb3c-48f1-a495-29fffd24e71d`, private
  and in Draft status.
- Main theorem `paley_two_set_conjecture`, ID
  `973dac7a-fac5-4e8e-961a-936f59818faa`, still Open.
- The four earlier elementary proofs still have server verdict ACCEPTED.

The [connection record](results/prove2me_connection_2026_09_05.json)
contains the current readbacks, supported environments, and IDs. It
contains no API key or bearer token. Authentication uses the existing
shared `credentials.json` outside this project, with mode 600 and a
verified gitignore rule. Access tokens are kept only in process memory;
requests refuse redirects and use only `https://prove2.me/api/v1`.

From this project, check the connection with:

```sh
python3 scripts/prove2me.py
```

Replace an expired credential using the hidden terminal prompt:

```sh
python3 scripts/prove2me.py --login
```

Read the shared workspace's `AGENTS.md`, `PROXIMITY.md`, `SKILL.md`,
and relevant references before editing or submitting Lean files.
Keep new theorem statements in `Theorems/`, proofs in `Solutions/`,
and reusable definitions in `Definitions/`. Use the supported environment
explicitly for new private submissions. Compile locally and inspect
the final theorem's axioms before uploading; a statement containing
`sorry` or a reduction importing open children is not a completed proof.

The mission's launch is a separate operation. The Prove2Me skill says,
"only your human can self-audit and launch your mission proposal."
The project connection and private theorem verification work while the
proposal remains a draft. No launch is needed to continue this research.

The current mathematical frontier is in
[the pass-27 assessment](research/parallel27-pass-summary-2026-09-05.md).
Prove2Me verifies submitted Lean statements; it does not turn the
remaining unproved cancellation estimates into assumptions we may use.

The new [projection counting lemma](research/prove2me-projection-count-2026-09-05.md)
passes local Lean verification with only `propext`, `Classical.choice`,
and `Quot.sound`. Its private statement job is
`843b8456-f9ad-4158-809d-8378028cc958`. The
[server record](results/prove2me_projection_server_2026_09_05.json)
contains the latest queue status; a queued statement is not a verified
server proof. The idempotent commands below reuse the saved job and
submission IDs:

```sh
python3 scripts/prove2me_projection.py poll
# Once the private statement has status PUBLISHED:
python3 scripts/prove2me_projection.py verify
python3 scripts/prove2me_projection.py poll
```

Two supporting finite-field cardinality identities also pass local Lean
checks in [the saved proof](prove2me/Check_paley_functional_kernel_card.lean).
The [combined algebraic core](research/prove2me-projection-core-2026-09-05.md)
now passes local Lean checks as well, including one common nonzero
projection for any nonempty family of at most |F| nonzero functionals.
The [code-specific projection and uniform count-bound equivalence](research/parallel26-mca-formalization-2026-09-05.md)
now also pass Lean, with seven principal statements checked. The pinned
ArkLib source already proves exact MCA interleaving invariance, so this
is independent verification of an existing result. It supplies no scalar
MCA upper bound and has no hosted verdict. The earlier private statement
job has now published theorem `d025929b-75ed-4c1c-b661-8e491017b153`.
Its proof was submitted once as `27ff60e7-7d6d-48da-8c34-7b123c6c0a3f`
and received an **ACCEPTED** server verdict. The theorem readback is
`Proved` in the pinned environment. Together with the four earlier
elementary proofs, this makes five accepted proofs. The hosted statement
is the finite incidence lemma; it is not the code-specific MCA reduction
or either full conjecture. Polling reuses these IDs.
The [setup audit](results/prove2me_setup_audit_2026_09_05.json) pins the
historical connection and counting-lemma artifacts while distinguishing
local proof status from the hosted verdict. The
[pass-26 audit](results/parallel26_pass_audit_2026_09_05.json) records the
code formalization and its then-pending queue readback. The
[pass-27 audit](results/parallel27_pass_audit_2026_09_05.json) records the
published statement, subsequent proof submission, and later research.
