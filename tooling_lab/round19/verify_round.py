#!/usr/bin/env python3
"""Bind relation-interface certificates and preserve the completed rounds."""
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import re
import sympy as sp

HERE = Path(__file__).resolve().parent
BINDINGS = []


def digest(path): return sha256(path.read_bytes()).hexdigest()


def bind(path, value):
    assert path.is_file() and digest(path) == value, str(path)
    BINDINGS.append(str(path.resolve()))


def scan(owner, value):
    if isinstance(value, list):
        for item in value: scan(owner, item)
    elif isinstance(value, dict):
        for key, item in value.items():
            if key.endswith('_sha256') and isinstance(item, dict):
                for name, expected in item.items(): bind(owner.parent/name, expected)
            else: scan(owner, item)


def main():
    preserved = {}
    for name, manifest in [('round18', 'manifest.json'), ('round17', 'checkpoint.json'), ('round16', 'manifest.json')]:
        path = HERE.parent/name/manifest; old = json.loads(path.read_text())
        for rel, entry in old['artifacts'].items(): bind(path.parent/rel, entry['sha256'])
        preserved[name] = len(old['artifacts'])
    reports = {}
    for path in HERE.glob('*.json'):
        if path.name == 'manifest.json': continue
        data = json.loads(path.read_text()); scan(path, data); reports[path.stem] = data
    for name in ('window_hierarchy', 'relation_barrier', 'projection_codec'): assert reports[name]['status'] == 'produced'
    for name in ('window_review', 'barrier_review', 'separation_bounds', 'codec_review', 'codec_family_controls'):
        assert reports[name]['status'] == 'passed'
    windows = reports['window_review']; assert sum(c['cube_words_checked'] for c in windows['cases']) == 6642
    assert len(windows['sharpness_controls']) == 6
    assert [r['windows_only'] for r in windows['cases'][-1]['profiles']] == [6561, 1153, 561, 385, 321, 289, 273, 257]
    assert [r['with_odd_support'] for r in windows['cases'][-1]['profiles']] == [3281, 577, 305, 257, 257, 257, 257, 257]
    barrier = reports['barrier_review']; assert sum(c['independent_resultants'] for c in barrier['cases']) == 25
    assert [c['certified_all_relation_l1_budget'] for c in barrier['cases']] == [5, 65, 16385, 26174, 11607941]
    assert barrier['finite_portfolio_relations_checked'] == 16 and barrier['finite_portfolio_aliases_checked'] == 128
    assert barrier['finite_portfolio_survivors'] == [127]*16
    separators = reports['separation_bounds']; assert len(separators['cases']) == 5
    assert separators['cases'][-1]['minimum_separating_l1_lower_bound'] == 11607942
    assert separators['cases'][-1]['minimum_separating_l1_upper_bound'] == 23215882
    codec = reports['codec_review']; assert codec['small_cube_words_checked'] == 6561 and codec['small_actual_words'] == 257
    assert sum(c['groups']['positives']['records'] for c in codec['cases']) == 484
    assert sum(c['groups']['probes']['records'] for c in codec['cases']) == 450
    assert all(c['mode'] == 'field_projection' for c in codec['cases'])
    assert codec['cases'][-1]['coefficient_products_per_reencoding'] == 576
    assert codec['barrier_pairs_checked'] == 5 and sum(c['records_checked'] for c in codec['degeneracy_controls']) == 204
    assert codec['malformed_controls_checked'] == 3 and len(codec['corrupt_artifacts_rejected']) == 12
    family = reports['codec_family_controls']; assert sum(c['relations_checked'] for c in family['cases']) == 80
    assert sum(c['digit_words_checked'] for c in family['cases']) == 117200
    assert family['cases'][-1]['radius_two_nonzero_relations'] == 0
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
           'scope': 'Exact window hierarchy, universal bounded-L1 relation obstruction and separator bracket, general projection/re-encoding codec and degenerate fallback. Membership is resolved within the stated algebraic model; efficient general enumeration and norm aggregates remain open.',
           'previous_artifacts_preserved': preserved, 'source_and_input_bindings': len(BINDINGS),
           'python_sources_compiled': len(sources), 'sympy_version': sp.__version__, 'markdown_local_links_passed': True,
           'window_words_checked': 6642, 'all_relation_budget_certificates': 5,
           'codec_small_census_checked': 6561, 'codec_family_relations_checked': 80,
           'codec_family_words_checked': 117200, 'artifacts': artifacts}
    (HERE/'manifest.json').write_text(json.dumps(out, indent=2)+'\n')
    print(f'Round19 verified: {len(artifacts)} artifacts, {len(BINDINGS)} bindings; preserved {preserved}.', flush=True)


if __name__ == '__main__': main()
