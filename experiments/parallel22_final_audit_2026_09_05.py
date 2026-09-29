#!/usr/bin/env python3
"""Reconcile the parallel pass, reviewed snapshots and current artifacts."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import ast
import json
import re

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/"results/parallel22_pass_audit_2026_09_05.json"
H = lambda p: sha256((ROOT/p).read_bytes()).hexdigest()
J = lambda p: json.loads((ROOT/p).read_text())
state = J("results/parallel22_prior_state_2026_09_05.json")
prior = state["prior_sha256"]
preserved = []
for p, h in prior.items():
    if p == "sources/manifest.json":
        continue
    assert H(p) == h, (p, "prior artifact changed")
    preserved.append(p)
old_bytes = state["prior_manifest_text"].encode()
assert sha256(old_bytes).hexdigest() == prior["sources/manifest.json"]
old_manifest = json.loads(old_bytes)
manifest = J("sources/manifest.json")
assert len(old_manifest) == 26 and len(manifest) == 27
assert manifest[:-1] == old_manifest
assert manifest[-1]["name"] == "nist-dlmf-20-11-gauss-2026-09-05"
assert "archive_error" in manifest[-1]
assert not any(k in manifest[-1] for k in ("path", "sha256", "bytes"))

input_results = [
    "results/parallel22_all_orders_obstruction_2026_09_05.json",
    "results/parallel21_positive_upper_review_2026_09_05.json",
    "results/parallel21_subgroup_next_input_2026_09_05.json",
    "results/parallel21_spectral_flat_direction_2026_09_05.json",
    "results/parallel21_squarefree_moments_2026_09_05.json",
]
input_checks = []
for result in input_results:
    for p, h in J(result)["input_sha256"].items():
        assert H(p) == h, (result, p)
        input_checks.append({"result": result, "path": p, "sha256": h})

snapshot_path = "results/parallel22_spectral_reviewed_snapshot_2026_09_05.json"
snapshot = J(snapshot_path)
for kind in ("note", "result"):
    assert sha256(snapshot[kind+"_text"].encode()).hexdigest() == snapshot[kind+"_sha256"]
historical = {snapshot["note_path"]: snapshot["note_sha256"],
              snapshot["result_path"]: snapshot["result_sha256"]}
old_spectral = json.loads(snapshot["result_text"])
new_spectral = J(snapshot["result_path"])
changed_result_keys = sorted(k for k in set(old_spectral) | set(new_spectral)
                             if old_spectral.get(k) != new_spectral.get(k))
assert changed_result_keys == ["input_sha256", "post_verification_metadata"]

original = snapshot["note_text"]
needle = "onto nonzero quadratic-residue frequencies. For u=1_C/√m, define"
replacement = """onto nonzero quadratic-residue frequencies. This frequency labeling uses
the classical positive quadratic Gauss evaluation. Its missing sign
normalization is supplied, from NIST DLMF 20.11.1–2, in
[the independent review, Section 4](parallel22-spectral-independent-review-2026-09-05.md).
The Gauss-magnitude argument below suffices for the Jacobi norm but
cannot determine that frequency label by itself.
For u=1_C/√m, define"""
needle2 = "All needed identities are derived explicitly here. This lane has\nnot received separate human or formal proof review."
replacement2 = """The Jacobi and anchor identities are derived explicitly here. The
Fourier-frequency identification additionally uses the classical positive
quadratic Gauss evaluation, as supplied in the linked separate-agent
review. Its review found no formula error and identified this omitted
source dependency. Separate human and formal proof review remain open."""
assert original.count(needle) == original.count(needle2) == 1
assert original.replace(needle, replacement).replace(needle2, replacement2) == (ROOT/snapshot["note_path"]).read_text()

review_paths = [
    "research/parallel22-all-orders-independent-review-2026-09-05.md",
    "research/parallel22-quartic-star-independent-review-2026-09-05.md",
    "research/parallel22-subgroup-independent-review-2026-09-05.md",
    "research/parallel22-spectral-independent-review-2026-09-05.md",
]
review_checks = []
for review in review_paths[1:]:
    text = (ROOT/review).read_text()
    pattern = r"\]\(([^)]+)\)\s*\|\s*" + chr(96) + r"([0-9a-f]{64})" + chr(96)
    linked = re.findall(pattern, text)
    pairs = [(str((Path(review).parent/p).as_posix()), h) for p, h in linked]
    pairs += re.findall(r"\|\s*((?:research|experiments|results)/[^\s|]+)\s*\|\s*([0-9a-f]{64})\s*\|", text)
    assert len(pairs) == 5, (review, pairs)
    for p, h in pairs:
        p = str((ROOT/p).resolve().relative_to(ROOT))
        source = "current file"
        if H(p) != h:
            assert historical.get(p) == h, (review, p, "unexplained changed input")
            source = snapshot_path
        review_checks.append({"review": review, "path": p, "reviewed_sha256": h, "resolved_against": source})
model_proof = "research/parallel22-all-orders-obstruction-2026-09-05.md"
model_hash = re.search(r"Reviewed SHA-256:\s+([a-f0-9]{64})", (ROOT/review_paths[0]).read_text()).group(1)
assert H(model_proof) == model_hash
independent_result = J("results/parallel22_all_orders_independent_audit_2026_09_05.json")
assert independent_result["reviewed_proof_sha256"] == model_hash
review_checks.append({"review": review_paths[0], "path": model_proof,
                     "reviewed_sha256": model_hash, "resolved_against": "current file"})

source_checks = []
for entry in manifest:
    if entry["name"] in ("thorner-zaman-2108.10878v2", "mcdonald-sahay-wyman-2210.03789v2"):
        assert H(entry["path"]) == entry["sha256"]
        source_checks.append({"path": entry["path"], "sha256": entry["sha256"]})
assert len(source_checks) == 2

artifacts = set()
for folder in ("research", "experiments", "results"):
    artifacts.update(str(p.relative_to(ROOT)) for p in (ROOT/folder).glob("parallel22*") if p != OUT)
artifacts.update([
    "research/parallel21-positive-upper-review-2026-09-05.md",
    "research/parallel21-subgroup-next-input-2026-09-05.md",
    "research/parallel21-spectral-next-input-2026-09-05.md",
    "experiments/parallel21_positive_upper_review_2026_09_05.py",
    "experiments/parallel21_subgroup_next_input_2026_09_05.py",
    "experiments/parallel21_spectral_flat_direction_2026_09_05.py",
    "results/parallel21_positive_upper_review_2026_09_05.json",
    "results/parallel21_subgroup_next_input_2026_09_05.json",
    "results/parallel21_subgroup_next_input_audit_2026_09_05.json",
    "results/parallel21_spectral_flat_direction_2026_09_05.json",
])
syntax_checks = []
for p in sorted(artifacts):
    if p.endswith(".py"):
        ast.parse((ROOT/p).read_text())
        syntax_checks.append(p)

central = ["README.md", "research/frontier.md", "research/source-audit.md",
           "research/checkpoint-2026-09-04.md", "sources/manifest.json"]
links, excluded = [], []
docs = [p for p in artifacts if p.endswith(".md")] + central[:-1]
for doc in sorted(docs):
    for line in (ROOT/doc).read_text().splitlines():
        for target in re.findall(r"\[[^\]\n]*\]\(([^)\n]+)\)", line):
            if target.startswith(("https://", "http://", "#", "mailto:", "codex://")):
                continue
            if target == "u" and "RawHyp" in line and doc == "research/source-audit.md":
                excluded.append({"source":doc, "target":target, "reason":"Historical mathematical function evaluation."})
                continue
            if target == "Σ_i b_i²" and "[L k^s+(3/2)pD_s](Σ_i b_i²)" in line and doc == model_proof:
                excluded.append({"source":doc, "target":target, "reason":"Explicit multiplication of mathematical bracketed factors."})
                continue
            name = target.split("#")[0]
            if not name:
                continue
            p = (ROOT/Path(doc).parent/name).resolve()
            assert p.exists() or p == OUT, (doc, target)
            links.append({"source":doc, "target":target})
counts = {p:J(p).get("counts", J(p).get("checks")) for p in input_results[:-1]}
counts["results/parallel22_all_orders_independent_audit_2026_09_05.json"] = independent_result["counts"]
audit = {
    "status": "Passed combined artifact and proof-scope audit. Four limited deductions reviewed; full Paley/subgroup/spectral/prize goals remain unproved.",
    "audited_at_utc": datetime.now(timezone.utc).isoformat(),
    "previous_turn_classification": "progress",
    "input_hash_checks": input_checks,
    "separate_agent_review_input_checks": review_checks,
    "spectral_normalization_repair": {
        "review": review_paths[-1], "source": manifest[-1], "reviewed_snapshot": snapshot_path,
        "only_two_prose_blocks_amended": True, "unchanged_verifier_rerun_completed": True,
        "all_nonmetadata_results_unchanged": True, "changed_result_keys": changed_result_keys,
        "scope": "NIST's classical Gauss formula supplies the missing sign. Finite surd checks do not independently prove this normalization.",
    },
    "source_input_hash_checks": source_checks,
    "preserved_pass21_artifacts": preserved,
    "source_manifest_prior_entries_preserved": 26, "source_manifest_entries": 27,
    "syntax_checks": syntax_checks, "verification_counts_by_program": counts,
    "local_link_checks": links, "mathematical_notation_excluded_from_link_check": excluded,
    "artifact_sha256": {p:H(p) for p in sorted(artifacts)},
    "central_file_sha256": {p:H(p) for p in central},
    "worker_status": "All three research lanes and four required separate-agent reviews completed. One spectral review attempt hit model capacity; a different available worker completed it. Latest live listing has no running child agent.",
    "review_limitations": [
        "Separate-agent review is not human refereeing or formal verification.",
        "Both weighted model implementations are not actual character kernels.",
        "Classical quartic-star inequalities reformulate an unproved positive upper bound.",
        "The subgroup payoff requires unproved uniform fourth/sixth relation bounds.",
        "The spectral theorem controls the uniform vector only.",
    ],
    "goal_status": "active and unachieved",
    "open_obligations": [
        "A uniform actual character upper estimate sufficient for full classical Paley.",
        "Uniform square-root subgroup cancellation, including exceptional primes.",
        "The all-vector spectral edge or another route to the required clique bounds.",
        "A verified exact quantitative bridge to the official Reed-Solomon prize.",
        "Human mathematical review and formal verification where claimed.",
    ],
}
OUT.write_text(json.dumps(audit, indent=2, ensure_ascii=False)+"\n")
for link in links:
    assert (ROOT/Path(link["source"]).parent/link["target"].split("#")[0]).resolve().exists()
print(json.dumps({"audit":str(OUT), "input_hash_checks":len(input_checks),
                  "review_input_checks":len(review_checks), "artifacts":len(artifacts),
                  "syntax_checks":len(syntax_checks), "local_links":len(links),
                  "prior_artifacts_preserved":len(preserved), "manifest_entries":len(manifest),
                  "goal_status":audit["goal_status"]}, indent=2))
