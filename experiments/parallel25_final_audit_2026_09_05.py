#!/usr/bin/env python3
"""Reconcile pass25 results and distinguish completed from interrupted reviews."""
import ast
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/parallel25_pass_audit_2026_09_05.json'
def H(p): return sha256((ROOT / p).read_bytes()).hexdigest()
def J(p): return json.loads((ROOT / p).read_text())

central = ['README.md', 'research/frontier.md', 'research/source-audit.md',
           'research/checkpoint-2026-09-04.md', 'sources/manifest.json']
prior = J('results/parallel25_prior_state_2026_09_05.json')
preserved = []
for p, expected in prior['prior_sha256'].items():
    if p not in central:
        assert H(p) == expected, ('prior artifact changed', p)
        preserved.append(p)
assert len(preserved) == 22
old = json.loads(prior['prior_manifest_text'])
assert sha256(prior['prior_manifest_text'].encode()).hexdigest() == prior['prior_sha256']['sources/manifest.json']
manifest = J('sources/manifest.json')
assert len(old) == 32 and manifest[:32] == old and len(manifest) == 34
ledgers = ['results/parallel25_prize_projection_source_2026_09_05.json',
           'results/parallel25_subgroup_source_2026_09_05.json']
assert manifest[32:] == [J(p) for p in ledgers]

lanes = {
    'prize': 'results/parallel25_prize_projection_2026_09_05.json',
    'signed': 'results/parallel25_signed_inversion_2026_09_05.json',
    'subgroup': 'results/parallel25_subgroup_growing_orders_2026_09_05.json',
    'spectral': 'results/parallel25_spectral_full_operator_2026_09_05.json',
    'root_review': 'results/parallel25_root_review_2026_09_05.json',
}
input_checks = []
for lane, p in lanes.items():
    d = J(p)
    assert d['status'].lower().startswith('pass') and d['counts']
    for source, expected in d['input_sha256'].items():
        assert H(source) == expected, (lane, 'input changed', source)
        input_checks.append({'lane': lane, 'path': source, 'sha256': expected})

prize = J(lanes['prize'])
assert prize['counts']['exhaustive_MCA_word_pairs'] == 401427
assert prize['counts']['projected_list_agreement_preservation'] == 25252
official = prize['official']
assert official['q'] == 2130706433 ** 6
assert official['budget'] == official['q'] // 2**128 == 274980728111395087
assert official['budget'] * (official['budget'] + 1) < 2 * official['q']
assert J(lanes['signed'])['counts']['actual_sign_half_moment_identities'] == 15140
assert J(lanes['subgroup'])['counts']['normalized_six_multiset_matches'] == 120585
assert J(lanes['spectral'])['counts']['actual_S3_action'] == 24
review = J(lanes['root_review'])
assert review['counts']['independent_128_dim_Bareiss_norm'] == 1
assert review['counts']['independent_direct_leakage_surd'] == 4
assert any(f['p'] == 269 and f['outside_worker_range'] for f in review['spectral'])

# Resolve both the main registry and original source-package manifests.
registry = {}
for entry in manifest:
    if entry.get('path') and entry.get('sha256'):
        registry[entry['path']] = ('sources/manifest.json', entry)
packages = ['sources/official-prize-2026-09-04/manifest.json',
            'sources/sigma-subgroup-2026-09-05/manifest.json',
            'sources/mixed-periods-2026-09-04/manifest.json']
def walk(obj, parent, source):
    if isinstance(obj, list):
        for value in obj: walk(value, parent, source)
    elif isinstance(obj, dict):
        name = obj.get('path', obj.get('file'))
        if name and obj.get('sha256'):
            registry[str(parent / name)] = (source, obj)
        for value in obj.values():
            if isinstance(value, (list, dict)): walk(value, parent, source)
for p in packages:
    walk(J(p), Path(p).parent, p)
source_checks = []
for p in sorted({r['path'] for r in input_checks if r['path'].startswith('sources/')}):
    assert p in registry, ('source missing registry', p)
    reg, entry = registry[p]
    assert H(p) == entry['sha256']
    assert (ROOT / p).stat().st_size == entry['bytes']
    source_checks.append({'path': p, 'sha256': H(p), 'bytes': entry['bytes'],
                          'manifest': reg})

artifacts = set()
for folder in ('research', 'experiments', 'results'):
    artifacts.update(str(p.relative_to(ROOT)) for p in (ROOT / folder).glob('parallel25*')
                     if p.is_file() and p != OUT)
syntax = []
for p in sorted(artifacts):
    if p.endswith('.py'):
        ast.parse((ROOT / p).read_text()); syntax.append(p)
links, exclusions = [], []
for doc in sorted(p for p in artifacts if p.endswith('.md')) + central[:-1]:
    for line in (ROOT / doc).read_text().splitlines():
        for target in re.findall(r'\[[^\]\n]*\]\(([^)\n]+)\)', line):
            if target.startswith(('http://', 'https://', '#', 'mailto:', 'codex://')):
                continue
            if target == 'u' and 'RawHyp' in line and doc == 'research/source-audit.md':
                exclusions.append({'source': doc, 'target': target, 'reason': 'Historical mathematical function evaluation.'})
                continue
            resolved = (ROOT / Path(doc).parent / target.split('#')[0]).resolve()
            assert resolved.exists() or resolved == OUT, (doc, target)
            links.append({'source': doc, 'target': target})

audit = {
    'status': 'Passed bounded-result and artifact audit; full goal remains unproved.',
    'audited_at_utc': datetime.now(timezone.utc).isoformat(),
    'previous_turn_classification': 'progress',
    'input_hash_checks': input_checks,
    'source_hash_checks': source_checks,
    'source_package_manifest_sha256': {p: H(p) for p in packages},
    'prior_artifacts_preserved': preserved,
    'source_manifest_previous_entries_preserved': 32,
    'source_manifest_added_entries': 2,
    'source_manifest_total_entries': 34,
    'verification_counts_by_lane': {k: J(v)['counts'] for k, v in lanes.items()},
    'syntax_checks': syntax, 'local_link_checks': links,
    'mathematical_notation_excluded': exclusions,
    'artifact_sha256': {p: H(p) for p in sorted(artifacts)},
    'central_file_sha256': {p: H(p) for p in central},
    'review_status': {
        'signed': 'Root, distinct from author; independent verifier passed.',
        'subgroup': 'Root, distinct from author; primary HTML and independent verifier checked.',
        'spectral': 'Root, distinct from author; completed author run and independent verifier.',
        'prize': 'Author audit and finite verifier only; separate-author review interrupted.',
    },
    'worker_state_at_reconciliation': 'All three errored with usage limits; no worker review completion claimed.',
    'goal_status': 'active and unachieved',
    'scope_limits': [
        'Exact scalar certificate equivalence leaves both scalar maxima and the spot check unbounded.',
        'MCA interleaving invariance keeps the same scalar field; it is not descent from F_q to F_p.',
        'Subgroup repeated-entry bound is above cubic scale and the six-distinct remainder is uncontrolled.',
        'Spectral block decomposition preserves full dimension and does not prove the operator saving.',
        'Signed inversion averages retain the unknown original moment.',
        'No worst-case cancellation exponent, full proof, human referee claim, or new formal theorem in these four notes.',
    ],
}
OUT.write_text(json.dumps(audit, indent=2) + '\n')
print(json.dumps({'status': 'passed', 'input_hashes': len(input_checks),
      'prior_artifacts_preserved': len(preserved), 'artifacts': len(artifacts),
      'source_files': len(source_checks), 'syntax_checks': len(syntax),
      'local_links': len(links), 'manifest_entries': len(manifest),
      'goal_status': audit['goal_status']}, indent=2))
