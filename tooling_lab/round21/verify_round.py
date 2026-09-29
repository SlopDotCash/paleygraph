#!/usr/bin/env python3
"""Bind erasure/residual certificates and preserve the prior frozen rounds."""
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
BINDINGS = []


def digest(path): return sha256(path.read_bytes()).hexdigest()


def bind(path, expected):
    assert path.is_file() and digest(path) == expected, str(path)
    BINDINGS.append(str(path.resolve()))


def scan(owner, value):
    if isinstance(value, list):
        for v in value: scan(owner, v)
    elif isinstance(value, dict):
        for key, v in value.items():
            if key.endswith('_sha256') and isinstance(v, dict):
                for name, expected in v.items(): bind(owner.parent/name, expected)
            else: scan(owner, v)


def main():
    preserved = {}
    for name, filename in [('round20', 'manifest.json'), ('round19', 'manifest.json'), ('round18', 'manifest.json'), ('round17', 'checkpoint.json'), ('round16', 'manifest.json')]:
        path = HERE.parent/name/filename; old = json.loads(path.read_text())
        for rel, entry in old['artifacts'].items(): bind(path.parent/rel, entry['sha256'])
        preserved[name] = len(old['artifacts'])
    reports = {}
    for path in HERE.glob('*.json'):
        if path.name == 'manifest.json': continue
        data = json.loads(path.read_text()); scan(path, data); reports[path.stem] = data
    for name in ['erasure_summary', 'translation_audit', 'pullback_summary', 'heldout_summary', 'prefix_certificate', 'plot_metadata']:
        assert reports[name]['status'] == 'produced'
    for name in ['erasure_review', 'translation_review', 'pullback_review', 'prefix_review', 'query_review']:
        assert reports[name]['status'] == 'passed'
    initial = reports['erasure_review']; assert len(initial['cases']) == 30
    assert sum(c['queries_replayed'] for c in initial['cases']) == 2621
    assert sum(c['scalar_candidates_checked'] for c in initial['cases']) == 234893
    assert sum(c['full_small_codebook_oracles'] for c in initial['cases']) == 1825
    assert sum(c['cyclic_controls_checked'] for c in initial['cases']) == 512
    assert sum(c['initial_budget_controls_checked'] for c in initial['cases']) == 6
    assert len(initial['corrupt_certificates_rejected']) == 6
    translation = reports['translation_review']; assert len(translation['cases']) == 5
    assert sum(c['pair_comparisons'] for c in translation['cases']) == 364029
    assert sum(c['difference_certificates'] for c in translation['cases']) == 27441
    assert sum(c['common_suffix_incidents'] for c in translation['cases']) == 218071
    assert sum(c['orbits_checked_in_boxes'] for c in translation['cases']) == 409158
    assert sum(c['witnesses_checked'] for c in translation['cases']) == 37
    assert len(translation['corrupt_translation_certificates_rejected']) == 4
    improved = reports['pullback_review']; assert len(improved['cases']) == 12
    assert sum(c['dual_functional_identities_checked'] for c in improved['cases']) == 242
    assert sum(c['queries_replayed'] for c in improved['cases']) == 892
    assert sum(c['scalar_candidates_checked'] for c in improved['cases']) == 57010
    assert sum(c['cyclic_controls_checked'] for c in improved['cases']) == 576
    assert len(improved['corrupt_pullback_certificates_rejected']) == 4
    caps = {c['erasures']: c for c in reports['pullback_summary']['cases']}
    assert caps[20]['universal_candidate_box_cap'] == 1 and caps[20]['universal_unique_completion']
    assert caps[24]['universal_candidate_box_cap'] == 4
    assert caps[32]['baseline_universal_cap'] == 110612791296 and caps[32]['universal_candidate_box_cap'] == 41472
    assert not caps[32]['universal_unique_completion']
    heldout = reports['heldout_summary']['cases']
    assert heldout[0]['statuses'] == {'complete': 32} and heldout[0]['cyclic_controls'] == 64
    assert all(c['statuses'] == {'budget_exceeded': 32} for c in heldout[1:])
    prefix = reports['prefix_review']
    assert prefix['scalar_words_checked'] == prefix['certified_width_lower_bound'] == 100000
    assert prefix['pairwise_distinctions_implied'] == 4999950000
    controls = reports['query_review']['controls']; assert len(controls) == 9
    assert sum(c['status'] == 'rejected' for c in controls) == 5
    assert reports['erasure_output']['count'] == 1 and reports['erasure_output']['completions'][0]['scalar'] == 1234567
    sources = list(HERE.glob('*.py'))
    for path in sources: compile(path.read_text(), str(path), 'exec')
    for path in [*HERE.glob('*.md'), HERE.parent/'NEXT_ITERATION.md']:
        prose = re.sub(r'```.*?```|`[^`\n]*`', '', path.read_text(), flags=re.S)
        for target in re.findall(r'\]\(([^)]+)\)', prose):
            if target.startswith(('http://', 'https://', '#')): continue
            resolved = path.parent/target.split('#')[0].strip('<>')
            assert resolved.exists() or resolved == HERE/'manifest.json', (path, target)
    artifacts = {str(path.relative_to(HERE)): {'sha256': digest(path), 'bytes': path.stat().st_size}
                 for path in sorted(HERE.rglob('*')) if path.is_file() and '__pycache__' not in path.parts
                 and path.name not in ('manifest.json', 'verification.log')}
    out = {'status': 'verified', 'created_utc': datetime.now(timezone.utc).isoformat(),
           'scope': 'Exact translated residual boxes, scalar-quotient erasure decoding, improved centered pullback bounds, cyclic completion guarantees, held-out geometry limits and a finite100000-state bound. General geometry and prize aggregates remain open.',
           'previous_artifacts_preserved': preserved, 'source_and_input_bindings': len(BINDINGS),
           'python_sources_compiled': len(sources), 'markdown_local_links_passed': True,
           'translation_comparisons_checked': 364029, 'initial_query_replays': 2621,
           'improved_query_replays': 892, 'cyclic_controls_checked': 1088,
           'prefix_scalar_words_checked': 100000, 'certified_width_lower_bound': 100000,
           'corrupted_math_certificates_rejected': 14, 'command_line_rejection_controls': 5,
           'artifacts': artifacts}
    (HERE/'manifest.json').write_text(json.dumps(out, indent=2)+'\n')
    print(f'Round21 verified: {len(artifacts)} artifacts, {len(BINDINGS)} bindings; preserved {preserved}.', flush=True)


if __name__ == '__main__': main()
