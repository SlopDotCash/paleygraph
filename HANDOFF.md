# Handoff notes for the next agent (2026-09-05)

**State of the problem: the prime-field Paley graph conjecture is open. Nothing
in this repository proves it. Do not inherit any claim without re-running its
verifier.**

The conjecture (Satake 2020, Conjecture 7): for every `0<ε<1` there are
`δ>0`, `p₀` such that for all primes `p>p₀` and all `A,B ⊆ F_p` with
`|A|,|B|>p^ε`, `|Σ_{a∈A,b∈B} χ(a+b)| ≤ p^{−δ}|A||B|`.

## 1. Read these first, in this order

1. `README.md` — status and the index of every pass.
2. `research/sigma-pass-summary-2026-09-05.md` — the sigma pass (four waves,
   sixteen notes). Sections 6 and 9 list the open obligations; section 9 states
   the obstruction in three equivalent forms.
3. `research/sigma-barriers-2026-09-05.md` — the barrier map and the `F_{p²}`
   test: any step of a proposed proof that holds verbatim over `F_{p²}` cannot
   pass `ε = 1/2`, because the subfield `F_p` is a `√q`-clique there.
4. `research/frontier.md` — the concurrent session's frontier (passes 1–29).

## 2. What is actually established (all with exact verifiers)

- Machine-checked in Lean 4 (Mathlib `c5ea003`, axioms only
  `propext, Classical.choice, Quot.sound`), and ACCEPTED by prove2me's server:
  shift orthogonality `Σ_x χ(x)χ(x+c) = −1`; the second moment
  `Σ_x(Σ_{b∈B}χ(x−b))² = |B|(p−|B|)`; the Chung bound; the interval-to-
  least-non-residue implication. Files: `~/prove2me_workspace/Theorems/Thm_paley_*.lean`,
  `Solutions/Sol_paley_*.lean` (outside this repo; see `PROVE2ME.md`).
- Exact reformulations (proved): the conjecture is equivalent to domination
  of the dilation `2k`-th moment over `t ≠ 0` by its square tuples
  (`sigma-tuple`); the subgroup case `B = H` is equivalent to Bourgain's open
  shifted-subgroup bound (`sigma-structured`); the extractor statement is
  equivalent to the conjecture (`sigma-extractor`); every termwise-Weil
  dilation-moment bound is invariant under `A ↦ uA` and cannot single out
  `t = 1` (`sigma-crux`).
- Literature facts (cited with sources): Chang's `4/9` theorem needs small
  additive doubling; no published exponent below `1/2` exists for two arbitrary
  sets (Fouvry–Shparlinski–Xi 2024); the shifted-subgroup bound at
  `|H| = p^ε` is Bourgain's Problem 5 in Chang's survey; best clique bound is
  Hanson–Petridis `(√(2p−1)+1)/2`; Sárközy's conjecture was proved by Kalmynin
  (arXiv:2504.10202); no 2024–2026 source claims a proof (`sigma-lit2026`).
- Refuted with exact witnesses (do not re-propose): Conjecture SI and SI*
  (biased dilate sets are not logarithmic), Conjecture T(k) for every
  constant (`A = B = Q` gives ratio `(2k−1)!!√p`), the shift-chain rescue by
  energy, the `k!` tuple constant (correct is `(2k−1)!!`), container-to-subset
  inheritance for subgroups, coset bootstrapping, "cancellation across cosets"
  as an easier problem (it is a relabelling), four structural hypotheses about
  biased rectangles, and the orbit-average hypothesis for extractors.
- Referee verdict on the concurrent session's sheaf-theoretic passes 9–10:
  no counterexample in 1.4 million checks, classified plausible with
  citation-level gaps; passes 11–29 were not refereed.

## 3. The obstruction, in one paragraph

Every bound written down in this workspace, except the sum-product inputs,
is field-agnostic and therefore cannot pass `ε = 1/2`. In the dilation
picture, moment bounds depend on `(A,B)` only through dilation-invariant
ratio-class sums. In the tuple picture, the conjecture needs a saving of
`(2k−1)p^{1/2+2kδ'}` over the triangle bound on the non-square Weil sums,
which are Frobenius traces of hyperelliptic curves branched on subsets of
the ratio set (`sigma-frobenius`). In the subgroup picture it is
cancellation in `Σ_{ψ∈Ψ} ψ(u)J(ψ,χ)` over character subgroups of order
`p^{1−ε}` (`sigma-dual`). A proof needs an input that fails over `F_{p²}`.

## 4. Candidates still open (highest leverage first)

1. T′(k): `|T_ns| ≤ C_{k,η} p^{1/2+η}(mn)^{2k}/min(m,n)^{k−1/2}` (sharp
   exponent, survives every family tested, implies the conjecture for
   `δ < ε/2 − 1/(2k)`); for subgroup pairs it is a Gaussian-moment statement
   for `T_H` across cosets (`sigma-stress`).
2. The Jacobi-sum statement J(ε,δ) (`sigma-dual` §4) — the same open input
   as Bourgain's Problem 5.
3. Single-prime `E_3(μ_N) ≤ N^{3+o(1)}` at every quartic prime, which would
   make the subgroup exponent `8/9` unconditional (`sigma-subgroup`).
4. Closing the citation-level gaps of passes 9–10 from text sources of Weil II
   and BBD (`sigma-referee` §G1–G3).

## 5. prove2me

- Workspace: `~/prove2me_workspace` (Lean v4.30.0, Mathlib `c5ea003`).
  `credentials.json` lives there (mode 600, not in this repo); never print it.
- Everything created is PRIVATE: eight theorems (`paley_*`), four with
  ACCEPTED proofs, and a draft mission proposal
  `b5e7121a-eb3c-48f1-a495-29fffd24e71d`. Ids and the request ledger are in
  `research/sigma-p2m-2026-09-05.md`; the launch commands are in
  `~/prove2me_workspace/proposals/README.md`. Launching or making anything
  public is the human's decision.

## 6. Reproducing

Every claim has a verifier under `experiments/` writing to `results/`;
`README.md` lists the commands. All sigma verifiers run in under ten minutes
each except where the note says otherwise. Python 3 with numpy, scipy and
sympy suffices; Lean checks need the prove2me workspace.

The `sources/` directory is not in this repository except for
`sources/manifest.json` (checksums): the PDFs and page images are
copyrighted papers. Re-fetch them from the arXiv ids in the notes.

## 7. Conventions and pitfalls

- One prefix per pass (`parallelN-*` for the concurrent session, `sigma-*`
  for this one). Re-read `README.md` from disk before editing it; two agents
  were writing it concurrently on 2026-09-05.
- Honesty labels in every note: PROVED / CITED / CONDITIONAL / HEURISTIC /
  REFUTED / OPEN. A claim without a complete argument is not PROVED.
- Test every inequality by exact integer computation before writing it down;
  three constants in this pass were wrong until a verifier caught them.
- Sub-agents were killed once by an account usage limit. Have workers create
  their note and verifier in their first actions and update incrementally.
- Sub-agents cannot be messaged in this harness; relaunch with a pointer to
  their scratch artifacts.
- A worker's `pkill -f` on its own verifier name will also kill the
  orchestrator's reproduction run; reproduce after the worker finishes.
