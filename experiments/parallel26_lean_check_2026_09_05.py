#!/usr/bin/env python3
"""Check the standalone MCA reduction in the existing pinned Lean environment."""
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import re
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = Path.home() / 'prove2me_workspace'
LAKE = Path.home() / '.local/bin/lake'
PROOF = ROOT / 'prove2me/Check_paley_mca_projection.lean'
OUT = ROOT / 'results/parallel26_mca_projection_2026_09_05.json'
LOG = ROOT / 'results/parallel26_mca_lean_2026_09_05.log'
EXPECTED_REV = 'c5ea00351c28e24afc9f0f84379aa41082b1188f'
EXPECTED_THEOREMS = {
    'projection_preserves_bad', 'rowsBad_iff_input_failure',
    'scalarBad_iff_input_failure', 'uniform_bad_count_bound_iff',
    'code_projection_preserves_bad', 'code_uniform_bad_count_bound_iff',
    'codeRowsBad_iff_input_failure',
}

def run(args):
    return subprocess.run(args, cwd=WORKSPACE, text=True,
                          capture_output=True, check=True).stdout.strip()

def main():
    proof = PROOF.read_bytes()
    before = sha256(proof).hexdigest()
    assert not re.search(r'\b(sorry|admit|axiom|unsafe)\b|import Theorems\.', proof.decode())
    assert proof == (WORKSPACE / 'Solutions' / PROOF.name).read_bytes()
    toolchain = (WORKSPACE / 'lean-toolchain').read_text().strip()
    assert toolchain == 'leanprover/lean4:v4.30.0'
    rev = run(['git', '-C', str(WORKSPACE / '.lake/packages/mathlib'), 'rev-parse', 'HEAD'])
    assert rev == EXPECTED_REV
    version = run([str(LAKE), 'env', 'lean', '--version'])
    started = time.monotonic()
    proc = subprocess.run([str(LAKE), 'env', 'lean', str(PROOF)],
                          cwd=WORKSPACE, text=True, capture_output=True)
    output = proc.stdout + proc.stderr
    LOG.write_text(output)
    assert proc.returncode == 0, 'Lean failed; inspect the saved log.'
    assert before == sha256(PROOF.read_bytes()).hexdigest()
    assert 'sorryAx' not in output
    axioms = {}
    for name, values in re.findall(r"'([^']+)' depends on axioms: \[([^\]]*)\]", output):
        axioms[name] = [s.strip() for s in values.split(',') if s.strip()]
    assert set(axioms) == {'PaleyMCAProjection.' + s for s in EXPECTED_THEOREMS}
    assert all(set(v) <= {'propext', 'Classical.choice', 'Quot.sound'} for v in axioms.values())
    record = {
        'status': 'passed', 'checked_at_utc': datetime.now(timezone.utc).isoformat(),
        'wall_seconds': time.monotonic() - started, 'exit_code': proc.returncode,
        'proof_sha256': before, 'log_sha256': sha256(LOG.read_bytes()).hexdigest(),
        'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
        'toolchain': toolchain, 'lean_version': version, 'mathlib_rev': rev,
        'axioms_by_theorem': axioms,
        'scope': 'Standalone code-agreement projection and uniform count-bound equivalence. '
                 'Independent verification of a result already present in pinned ArkLib; '
                 'no direct ArkLib import, new scalar bound, or full conjecture proof.',
    }
    OUT.write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record, indent=2))

if __name__ == '__main__':
    main()
