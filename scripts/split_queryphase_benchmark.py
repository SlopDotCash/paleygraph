#!/usr/bin/env python3
"""Prepare isolated diagnostic candidates; does not run or validate Lean.

This transformation is deliberately pinned to the observed v19 source. The
original statements and proof branches are preserved, modulo indentation and
the explicit branch hypothesis. Generated files are unvalidated experiments.
"""
import argparse
import hashlib
import json
from pathlib import Path


EXPECTED_SHA256 = "138b2b0eb454e2b0010b59fe92462ef41f21496a93d68c624a493182cf740abb"
NAME = "query_phase_step_preserves_fold"


def sha(text):
    return hashlib.sha256(text.encode()).hexdigest()


def prepare(source, output):
    original = source.read_text()
    if sha(original) != EXPECTED_SHA256:
        raise ValueError("Source differs from reviewed v19; re-review before transforming it")
    option = "set_option maxHeartbeats 4000000 in\n"
    start = original.index(option)
    declaration = original.index("lemma " + NAME + "\n", start)
    end_header = original.index("            omega)) := by\n", declaration)
    end_header += len("            omega)) := by\n")
    header = original[declaration:end_header]
    tail_start = original.index("\n\nend FinalQueryRoundIOR", end_header)
    proof = original[end_header:tail_start]
    marker = "  by_cases h_k_pos : k.val > 0\n"
    assert proof.count(marker) == 1
    common, cases = proof.split(marker)
    pos_marker = "  · -- Case k > 0: The guard is present.\n"
    zero_marker = "  · -- Case k = 0: No guard.\n"
    assert cases.startswith(pos_marker) and cases.count(zero_marker) == 1
    positive, zero = cases[len(pos_marker):].split(zero_marker)

    def unindent(branch):
        assert all(not line.strip() or line.startswith("    ") or line.lstrip().startswith("--")
                   for line in branch.splitlines())
        result = "\n".join(line[2:] if line.strip() else line for line in branch.split("\n"))
        # Check that no non-whitespace proof content changed.
        assert "".join(branch.split()) == "".join(result.split())
        return result

    binder_end = "))))) :\n"
    assert header.count(binder_end) == 1

    def helper(suffix, hypothesis, branch):
        signature = header.replace("lemma " + NAME, "private lemma " + NAME + suffix, 1)
        signature = signature.replace(binder_end, ")))))\n    (h_k_pos : " + hypothesis + ") :\n", 1)
        return option + signature + common + unindent(branch).rstrip() + "\n\n"

    positive_decl = helper("_positive", "k.val > 0", positive)
    zero_decl = helper("_zero", "¬ k.val > 0", zero)
    wrapper = option + header + marker
    arguments = (" (𝔽q := 𝔽q) (β := β) (γ_repetitions := γ_repetitions)\n"
                 "      v c_k s' stmtIn oStmtIn h_relIn h_c_k_correct_of_k_pos\n"
                 "      challenges h_s'_mem h_k_pos\n")
    for suffix in ["_positive", "_zero"]:
        wrapper += "  · exact " + NAME + suffix + arguments
    assert wrapper.startswith(option + header)

    prefix = original[:start].replace(
        "set_option linter.style.longFile 3100\n",
        "set_option linter.style.longFile 3100\n"
        "set_option profiler true\nset_option profiler.threshold 100\n", 1)
    tail = original[tail_start:]
    notice = ("/- Diagnostic candidate: NOT Lean-validated or benchmarked.\n"
              "Original proof branches preserved; source SHA pinned by generator. -/\n")
    bodies = {"prefix": "", "positive": positive_decl, "zero": zero_decl,
              "combined": positive_decl + zero_decl + wrapper}
    output.mkdir(parents=True, exist_ok=False)
    files = {}
    for kind, body in bodies.items():
        text = prefix + notice + body + tail
        path = output / ("queryphase-" + kind + ".lean")
        path.write_text(text)
        files[kind] = {"path": str(path.resolve()), "sha256": sha(text)}
    manifest = {"original_source": str(source.resolve()), "original_sha256": sha(original),
                "files": files, "source_branches_preserved": True,
                "original_wrapper_signature_preserved": True,
                "lean_validated": False, "benchmarked": False,
                "existing_job_changed": False}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    print(json.dumps(prepare(args.source, args.output), indent=2))
