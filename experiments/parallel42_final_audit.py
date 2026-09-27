#!/usr/bin/env python3
"""Verify pass42 artifact scope and preserve the completed prior pass."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import re

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def main():
    prior = json.loads((ROOT / 'results/parallel42_prior_state_2026_09_06.json').read_text())
    central = ['README.md', 'research/frontier.md', 'research/source-audit.md',
               'research/checkpoint-2026-09-04.md']
    unchanged = {p: h for p, h in prior['files'].items() if p not in central}
    for path, expected in unchanged.items():
        assert digest(path) == expected, path
    scope = json.loads((ROOT / 'results/parallel42_source_scope_2026_09_06.json').read_text())
    for source in scope['sources']:
        assert digest(source['path']) == source['sha256']
        if 'text_path' in source:
            assert digest(source['text_path']) == source['text_sha256']
    assert digest(scope['archive_manifest_path']) == scope['archive_manifest_sha256']
    assert digest('sources/manifest.json') == scope['main_manifest_sha256']
    norm = json.loads((ROOT / 'results/parallel42_norm_budget_2026_09_06.json').read_text())
    affine = json.loads((ROOT / 'results/parallel42_classical_affine_2026_09_06.json').read_text())
    assert norm['all_passed'] and affine['all_passed']
    assert sum(c['independent_determinants'] for c in norm['cases']) == 446
    assert sum(c['extra_split_primes_checked'] for c in norm['cases']) == 101
    assert [len(c['all_split_exceptional_primes']) for c in norm['cases']] == [1, 2, 18]
    witness = norm['zero_diagonal_witness']
    assert witness['positive_a_histogram'] == {'1': 1, '2': 7}
    assert witness['X'] == witness['X_dist'] == 72
    assert len(witness['rich_cells']) == 12
    assert not witness['quartic_window_n4_over4_to_n4']
    assert affine['center_tuples'] == 2266
    artifacts = sorted(str(p.relative_to(ROOT)) for folder in ['experiments', 'research', 'results']
                       for p in (ROOT / folder).glob('parallel42*')
                       if p.is_file() and p.name != 'parallel42_pass_audit_2026_09_06.json')
    local_links = 0
    for path in [p for p in artifacts if p.endswith('.md')]:
        contents = (ROOT / path).read_text()
        for target in re.findall(r'\]\(([^)]+)\)', contents):
            if '://' in target or target.startswith('#'):
                continue
            target = target.split('#')[0]
            assert ((ROOT / path).parent / target).exists(), (path, target)
            local_links += 1
    for path in central:
        assert 'parallel42' in (ROOT / path).read_text()[:4000]
    audit = dict(
        observed_at_utc=datetime.now(timezone.utc).isoformat(),
        classification='progress: compatible norm budget, absolute prime-exception count, zero-diagonal witness, and exact affine recovery comparison',
        previous_turn_classification=prior['previous_turn_classification'],
        goal_achieved=False,
        weighted_budget='sum_(p=1 mod n) X_p log p <= (M/2) log([8n^2(n-1)^2-2n^4]/M)',
        quartic_absolute_exception_count='O_c(n^(63/31)/log n) for X_p>n^(61/31)',
        uniform_X_bound_proved=False,
        energy_exponent='49/20 unchanged', period_exponent='71/72 unchanged',
        finite_scope=dict(complete_conductors=[4, 8, 16], exceptional_counts=[1, 2, 18],
                          independent_determinants=446, extra_split_primes=101,
                          reflection_center_tuples=2266, zero_diagonal_witness_p=353,
                          zero_diagonal_witness_in_quartic_window=False),
        review_scope=scope['review_scope'],
        parallel_lanes='Main findings returned; all three subsequently failed at account usage limit.',
        no_lean_build=True, no_process_changes=True, no_prove2me_submission=True,
        prior_noncentral_files_unchanged=len(unchanged), local_links_checked=local_links,
        artifacts_sha256={p: digest(p) for p in artifacts},
        central_changes={p: dict(before=prior['files'][p], after=digest(p)) for p in central},
        source_manifest_sha256=digest('sources/manifest.json'),
        next_action='Control individual exceptional primes or exploit additional structure of actual joint incidence and difference levels; retain full original goal.')
    dest = ROOT / 'results/parallel42_pass_audit_2026_09_06.json'
    dest.write_text(json.dumps(audit, indent=2) + '\n')
    print(json.dumps(dict(all_passed=True, artifacts=len(artifacts), local_links=local_links,
                          prior_noncentral_unchanged=len(unchanged), goal_achieved=False,
                          audit=str(dest))))


if __name__ == '__main__':
    main()
