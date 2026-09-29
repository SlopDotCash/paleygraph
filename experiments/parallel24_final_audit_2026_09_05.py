#!/usr/bin/env python3
"""Reconcile the four bounded pass24 lanes, reviews, and prior artifacts."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import ast
import json
import re

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/parallel24_pass_audit_2026_09_05.json'
def H(p): return sha256((ROOT/p).read_bytes()).hexdigest()
def J(p): return json.loads((ROOT/p).read_text())

central = ['README.md', 'research/frontier.md', 'research/source-audit.md',
           'research/checkpoint-2026-09-04.md', 'sources/manifest.json']
prior = J('results/parallel24_prior_state_2026_09_05.json')
preserved = []
for p, expected in prior['prior_sha256'].items():
    if p in central:
        continue
    assert H(p) == expected, ('prior artifact changed', p)
    preserved.append(p)
old_manifest = json.loads(prior['prior_manifest_text'])
assert sha256(prior['prior_manifest_text'].encode()).hexdigest() == prior['prior_sha256']['sources/manifest.json']
assert H('sources/manifest.json') == prior['prior_sha256']['sources/manifest.json']
assert J('sources/manifest.json') == old_manifest and len(old_manifest) == 32

lanes = {
    'inversion': 'results/parallel24_inversion_orbit_upper_2026_09_05.json',
    'classical': 'results/parallel24_classical_exceptions_2026_09_05.json',
    'spectral': 'results/parallel24_spectral_operator_2026_09_05.json',
    'subgroup': 'results/parallel24_subgroup_unbalanced_2026_09_05.json',
}
input_checks = []
for lane, path in list(lanes.items())+[('classical_root_review', 'results/parallel24_classical_review_2026_09_05.json')]:
    d = J(path)
    assert 'status' in d and d.get('counts', d.get('checks'))
    for p, expected in d['input_sha256'].items():
        assert H(p) == expected, (lane, 'input hash mismatch', p)
        input_checks.append({'lane': lane, 'path': p, 'sha256': expected})

reviews = []
review_checks = []
for lane in lanes:
    p = 'research/parallel24-'+lane+'-independent-review-2026-09-05.md'
    doc = (ROOT/p).read_text()
    found = []
    for line in doc.splitlines():
        cols = line.split('|')
        if len(cols) < 4:
            continue
        match = re.fullmatch(r'\s*`?([0-9a-f]{64})`?\s*', cols[2])
        if not match:
            continue
        name = cols[1].strip().strip('`')
        link = re.search(r'\]\(([^)]+)\)', name)
        if link:
            name = str((ROOT/Path(p).parent/link.group(1)).resolve().relative_to(ROOT))
        assert H(name) == match.group(1), (p, 'review hash mismatch', name)
        found.append(name)
        review_checks.append({'review': p, 'path': name, 'sha256': match.group(1)})
    assert len(found) >= 3 and lanes[lane] in found, (p, 'incomplete review ledger')
    reviews.append(p)

inv = J(lanes['inversion'])
assert inv['cases'] == 573
assert inv['checks']['signed_int64_accumulation_safe'] == 573
assert inv['checks']['ordered_tuple_expansion'] == 27
assert inv['checks']['bounded_good_representative'] == 1149
assert inv['checks']['unsigned_inversion_coordinates'] == inv['checks']['signed_inversion_coordinates'] == 3570540
cl = J(lanes['classical'])
assert cl['counts']['complete_swap_drift_identities'] == 1745
assert cl['counts']['thin_slice_volume_and_drift_certificates'] == 7
sp = J(lanes['spectral'])
assert sp['counts']['actual_prime_neighborhood'] == 88
assert sp['counts']['exact_next_diagonal_coefficient'] == 85
assert sp['counts']['zero_first_residual_fields'] == 3
sg = J(lanes['subgroup'])
assert sg['counts']['eligible_actual_prime'] == 690
assert sg['counts']['independent_integer_norm_determinants'] == 1770
assert sg['counts']['actual_product_ratio_correlation'] == 10088

artifacts = set()
for folder in ('research', 'experiments', 'results'):
    artifacts.update(str(p.relative_to(ROOT)) for p in (ROOT/folder).glob('parallel24*')
                     if p.is_file() and p != OUT)
syntax = []
for p in sorted(artifacts):
    if p.endswith('.py'):
        ast.parse((ROOT/p).read_text())
        syntax.append(p)

links, exclusions = [], []
for doc in sorted(p for p in artifacts if p.endswith('.md'))+central[:-1]:
    for line in (ROOT/doc).read_text().splitlines():
        for target in re.findall(r'\[[^\]\n]*\]\(([^)\n]+)\)', line):
            if target.startswith(('http://', 'https://', '#', 'mailto:', 'codex://')):
                continue
            if target == 'u' and 'RawHyp' in line and doc == 'research/source-audit.md':
                exclusions.append({'source': doc, 'target': target, 'reason': 'Historical mathematical function evaluation.'})
                continue
            resolved = (ROOT/Path(doc).parent/target.split('#')[0]).resolve()
            assert resolved.exists() or resolved == OUT, (doc, target)
            links.append({'source': doc, 'target': target})

source_checks = []
source_paths = sorted({row['path'] for row in input_checks if row['path'].startswith('sources/')})
for p in source_paths:
    entry = next((e for e in old_manifest if e.get('path') == p), None)
    if entry:
        assert H(p) == entry['sha256']
        assert (ROOT/p).stat().st_size == entry['bytes']
    else:
        assert p.endswith('.txt'), ('unregistered non-text source', p)
    source_checks.append({'path': p, 'sha256': H(p), 'bytes': (ROOT/p).stat().st_size,
                          'manifest_entry': bool(entry), 'scope': 'Existing pinned primary archive or associated text extraction.'})

audit = {
    'status': 'Passed bounded-result and artifact audit; full goal remains unproved.',
    'audited_at_utc': datetime.now(timezone.utc).isoformat(),
    'previous_turn_classification': 'progress',
    'input_hash_checks': input_checks,
    'review_input_hash_checks': review_checks,
    'source_hash_checks': source_checks,
    'prior_artifacts_preserved': preserved,
    'source_manifest_unchanged_entries': 32,
    'verification_counts_by_lane': {k: J(v).get('counts', J(v).get('checks')) for k, v in lanes.items()},
    'independent_root_review_counts': J('results/parallel24_classical_review_2026_09_05.json')['counts'],
    'review_files': reviews,
    'syntax_checks': syntax,
    'local_link_checks': links,
    'mathematical_notation_excluded': exclusions,
    'artifact_sha256': {p: H(p) for p in sorted(artifacts)},
    'central_file_sha256': {p: H(p) for p in central},
    'review_scope': 'Every author has a distinct reviewing agent; root reviewed the classical lane. Agent review is not human refereeing or formal verification.',
    'goal_status': 'active and unachieved',
    'open_obligations': [
        'Transported signed moments or uniform exceptional-input bounds sufficient for full classical Paley.',
        'Unbalanced subgroup relations at growing orders and uniform square-root cancellation including exceptional primes.',
        'The remaining actual spectral operator, growing depth, and other sectors with coupling retained.',
        'Arbitrary-word remainder and scalar-ratio bounds meeting the pinned official certificate and separate spot-check obligation.',
    ],
    'scope_limits': [
        'Unsigned good representative in every inversion family does not bound the original moment.',
        'Local persistence/sign-split calculations do not exclude all possible arithmetic approaches.',
        'The spectral upper estimate is on a specified two-dimensional domain, with its leakage included.',
        'The subgroup classification covers only the finite order classes 4,8,16.',
        'No worst-case cancellation exponent, full goal proof, or official prize certificate is established.',
        'No source version or official prize state refreshed; all prior manifest entries unchanged.',
    ],
}
OUT.write_text(json.dumps(audit, indent=2)+'\n')
print(json.dumps({'audit': str(OUT), 'input_hash_checks': len(input_checks),
                  'review_hash_checks': len(review_checks), 'prior_artifacts_preserved': len(preserved),
                  'artifacts': len(artifacts), 'syntax_checks': len(syntax),
                  'local_links': len(links), 'source_files': len(source_checks),
                  'goal_status': audit['goal_status']}, indent=2))
