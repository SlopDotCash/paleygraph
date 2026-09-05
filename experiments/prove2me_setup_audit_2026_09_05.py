#!/usr/bin/env python3
"""Audit connection and exact submitted proof artifacts, without reading secrets."""
import ast
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import re
import stat
import subprocess

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = Path.home() / 'prove2me_workspace'
OUT = ROOT / 'results/prove2me_setup_audit_2026_09_05.json'
def J(p): return json.loads((ROOT / p).read_text())
def H(p): return sha256((ROOT / p).read_bytes()).hexdigest()

connection = J('results/prove2me_connection_2026_09_05.json')
assert connection['authenticated']
assert connection['goal']['status'] == 'Open'
assert connection['proposal']['visibility'] == 'private'
assert connection['proposal']['status'] == 'Draft'
assert len(connection['verified_existing_submissions']) == 4
assert all(v['status'] == 'ACCEPTED' for v in connection['verified_existing_submissions'].values())
assert connection['environment'] == 'c5ea00351c28e24afc9f0f84379aa41082b1188f'
assert connection['environment'] in json.dumps(connection['supported_environments'])
assert stat.S_IMODE((WORKSPACE / 'credentials.json').stat().st_mode) == 0o600
ignored = subprocess.run(['git', 'check-ignore', '-q', 'credentials.json'],
    cwd=WORKSPACE, capture_output=True)
assert ignored.returncode == 0

payload = J('prove2me/projection-survival-problem.json')
local = J('results/prove2me_projection_local_2026_09_05.json')
server = J('results/prove2me_projection_server_2026_09_05.json')
assert payload['private'] is True and payload['env'] == connection['environment']
proof_path = 'prove2me/Sol_paley_projection_survival_count.lean'
proof = (ROOT / proof_path).read_text()
assert H(proof_path) == local['proof_sha256'] == server['proof_sha256']
assert H('prove2me/projection-survival-problem.json') == local['payload_sha256'] == server['payload_sha256']
assert H('research/prove2me-projection-count-2026-09-05.md') == local['source_sha256']
assert local['local_proof_exit_code'] == local['target_compilation_exit_code'] == 0
assert local['local_proof_axioms'] == ['propext', 'Classical.choice', 'Quot.sound']
assert not local['sorryAx'] and not re.search(r'\b(sorry|admit|axiom)\b', proof)
pre, body = proof.split('theorem solution', 1)
statement = 'theorem paley_projection_survival_count' + body.split(':= by', 1)[0] + ':= by sorry'
assert statement == payload['problems'][0]['formal_statement']
assert pre.strip() == payload['problems'][0]['preamble']
shared_proof = WORKSPACE / 'Solutions/Sol_paley_projection_survival_count.lean'
shared_target = WORKSPACE / 'Theorems/Thm_paley_projection_survival_count.lean'
assert shared_proof.read_text() == proof
assert shared_target.read_text() == pre + statement + '\n'
extra_proofs = {
    'prove2me/Check_paley_functional_kernel_card.lean': 'results/prove2me_functional_counts_local_2026_09_05.json',
    'prove2me/Check_paley_nonzero_functional_count.lean': 'results/prove2me_projection_core_local_2026_09_05.json',
}
for path, record in extra_proofs.items():
    checked = J(record)
    assert checked['exit_code'] == 0 and H(path) == checked['proof_sha256']
    assert all(v == ['propext', 'Classical.choice', 'Quot.sound'] for v in checked['theorems'].values())
    assert (ROOT / path).read_bytes() == (WORKSPACE / 'Solutions' / Path(path).name).read_bytes()
    assert not re.search(r'\b(sorry|admit|axiom)\b|import Theorems\.', (ROOT / path).read_text())
job = server.get('publish_job', {})
if job.get('status') == 'PUBLISHED':
    assert job['visibility'] == 'private' and job['formal_statement'] == statement
    assert job['theorem_id'] == server['theorem_id']
verdict = server.get('verdict', {}).get('status')
if verdict == 'ACCEPTED':
    assert server['verdict']['theorem_id'] == server['theorem_id']
    assert server['theorem_readback']['status'] == 'Proved'
    assert server['theorem_readback']['mathlib_rev'] == connection['environment']

artifacts = ['PROVE2ME.md', 'scripts/prove2me.py', 'scripts/prove2me_projection.py',
    'research/prove2me-projection-count-2026-09-05.md',
    'results/prove2me_connection_2026_09_05.json',
    'results/prove2me_projection_local_2026_09_05.json',
    'results/prove2me_projection_server_2026_09_05.json',
    'prove2me/projection-survival-problem.json', proof_path,
    'prove2me/Sol_paley_projection_survival_count.md',
    'prove2me/service-skill-current.md',
    'experiments/prove2me_setup_audit_2026_09_05.py']
artifacts += [str(p.relative_to(ROOT)) for p in (ROOT / 'prove2me/service-references').glob('*.md')]
artifacts += list(extra_proofs) + list(extra_proofs.values())
artifacts.append('research/prove2me-projection-core-2026-09-05.md')
for p in artifacts:
    text = (ROOT / p).read_text()
    assert not re.search(r'p2m_[A-Za-z0-9_-]+\.[A-Za-z0-9_-]+\.', text), ('credential found', p)
    if p.endswith('.py'): ast.parse(text)
links = []
for p in ['PROVE2ME.md', 'research/prove2me-projection-count-2026-09-05.md',
          'research/prove2me-projection-core-2026-09-05.md']:
    for target in re.findall(r'\[[^\]\n]*\]\(([^)\n]+)\)', (ROOT / p).read_text()):
        if target.startswith(('https://', 'http://', '#')): continue
        resolved = (ROOT / Path(p).parent / target.split('#')[0]).resolve()
        assert resolved.exists() or resolved == OUT, (p, target)
        links.append({'source': p, 'target': target})
audit = {
    'status': 'setup verified', 'audited_at_utc': datetime.now(timezone.utc).isoformat(),
    'platform_version': connection['platform_version'],
    'credential_permissions': '0600', 'credential_gitignored': True,
    'credentials_read_by_audit': False,
    'existing_accepted_proofs': 4, 'main_goal': 'Open',
    'mission_status': 'private Draft',
    'new_local_proof': 'passed with only standard axioms',
    'additional_local_core_theorems': 4,
    'new_publish_status': job.get('status'), 'new_server_verdict': verdict,
    'local_link_checks': links,
    'artifact_sha256': {p: H(p) for p in sorted(artifacts)},
    'scope': 'Hosted payload is the finite incidence lemma. Four further local theorems cover kernel counts and a common projection; the code-specific reduction and both conjectures remain unproved here.',
}
OUT.write_text(json.dumps(audit, indent=2) + '\n')
print(json.dumps({k: v for k, v in audit.items() if k != 'artifact_sha256'}, indent=2))
