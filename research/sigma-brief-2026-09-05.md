# Sigma pass: shared brief for parallel workers (2026-09-05)

Read this before doing anything else. Then read `README.md` and
`research/frontier.md` in `/Users/shawwalters/Desktop/paleygraph`.

## The target

The prime-field Paley graph conjecture (two-set form, Satake 2020, Conjecture 7):
for every `0 < ε < 1` there exist `δ > 0` and `p₀` such that for every prime
`p > p₀` and all `A, B ⊆ F_p` with `|A|, |B| > p^ε`,

    | Σ_{a∈A} Σ_{b∈B} χ_p(a+b) | ≤ p^{-δ} |A| |B|,

where `χ_p` is the Legendre symbol with `χ_p(0)=0`. The `a-b` and `a+b`
forms are equivalent (replace `B` by `-B`).

The conjecture is open. Ten prior passes in this workspace have produced
finite certificates, moment identities, a spectral (Kunisky-type) programme,
and a subgroup Gauss-sum line. None proves the conjecture. Do not assume any
prior "source-dependent" claim in `research/` is correct; treat them as
unreviewed unless you check them.

## Non-negotiable rules

1. **Honesty over ambition.** Never write "proved", "theorem", or "we show"
   for anything you have not written a complete, gap-free argument for.
   Label everything as one of: `PROVED` (complete argument in your note),
   `CITED` (with exact reference and statement), `CONDITIONAL` (state the
   hypothesis precisely), `HEURISTIC`, `REFUTED` (with an explicit witness),
   or `OPEN`.
2. **Exact checks.** Every identity, inequality, or constant you use must be
   tested by exact integer/rational computation at small primes in a Python
   script. If a claimed bound fails a numerical test, the claim is refuted,
   not the test. Record witnesses.
3. **Cite precisely.** When you rely on the literature, fetch the source
   (arXiv HTML/abstract or a locally saved PDF/text in `sources/`) and quote
   the statement number and hypotheses. Do not cite from memory alone; if
   you cannot fetch, say `UNVERIFIED CITATION`.
4. **No editing of shared files.** Do not edit `README.md`,
   `research/frontier.md`, or any file that another pass created. Create
   only files with your own prefix (below). Do not run `git`.
5. **Tooling.** Python 3 with numpy, scipy, sympy is available
   (`/opt/miniconda3/bin/python3`). Web search and fetch tools are available
   through `ToolSearch` (load `WebSearch`, `WebFetch`). Lean 4 toolchains
   exist under elan (`~/.elan`), including v4.30.0. Do not use `/tmp`; use
   the scratchpad path given in your prompt if you need temporary files.
6. **Do not spawn further sub-agents.**

## Files you must produce

Let `<dir>` be your direction slug from your prompt. **Another Claude
session is concurrently writing files with the `parallel11` prefix in this
same workspace. Never touch those. Your prefix is `sigma`.**

- `research/sigma-<dir>-2026-09-05.md` — your write-up. First line
  after the title must be a bold **Status:** sentence stating what is
  PROVED, what is REFUTED and what is OPEN. Then precise statements,
  complete proofs, exact citations, witnesses, and an explicit list of
  remaining obligations. Write it like a paper section, not a diary.
- `experiments/sigma_<dir>_2026_09_05.py` — a runnable verifier using
  only the standard library plus numpy/scipy/sympy. It must finish in under
  ten minutes on this machine and write
  `results/sigma_<dir>_2026_09_05.json` with the counts of checks
  performed and every witness.
- Your final message to the orchestrator: a summary of at most 400 words
  with the status line, the strongest fully proved statement (stated
  exactly), the strongest refuted hypothesis (with witness), and the open
  obligations. No file dumps.

## Background facts you may use freely (all elementary, all PROVED in
`research/moments-and-obstructions.md` or standard)

- For `c ≠ 0`: `Σ_x χ(x)χ(x+c) = -1`. Hence for `B ⊆ F_p`,
  `Σ_x |Σ_{b∈B} χ(x-b)|² = |B|(p-|B|)`, and (Chung/Vinogradov)
  `|Σ_{a∈A,b∈B} χ(a-b)| ≤ √(|A||B|(p-|B|))`. This is nontrivial only when
  `|A||B| > p`: the "square-root barrier".
- Weil: for a polynomial `f` over `F_p` of degree `d` that is not a constant
  times a square, `|Σ_x χ(f(x))| ≤ (d-1)√p`.
- Karatsuba-type amplification (Hölder in one variable plus Weil) gives
  nontrivial bounds when `|A| > p^{1/2+ε}` and `|B| > p^ε`; it cannot pass
  below `|A| ≈ √p` because the Hölder step sums over all of `F_p`.
- The conjecture for `A = B = {1,…,N}` with `N = p^ε` implies the least
  quadratic non-residue is at most `2N`, i.e. it implies Vinogradov's
  conjecture, which is open. Any general proof must therefore beat the
  Burgess barrier for short interval character sums as a special case.
- Graham–Ringrose: there are infinitely many `p` with cliques of size
  `≫ log p · log log log p`, so no `O(log p)` clique bound holds.

Root (the orchestrator) will integrate all notes into a pass summary and
update `README.md` and `frontier.md`. Your job is one direction, done with
full rigor and full honesty.
