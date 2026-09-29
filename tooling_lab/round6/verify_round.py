#!/usr/bin/env python3
"""Integrate completed evidence and bind current bytes; not a full rerun."""
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
LAB = HERE.parent
BINDINGS = []


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def read(name):
    return json.loads((HERE / name).read_text())


def bind(path, expected):
    assert path.is_file() and digest(path) == expected, str(path)
    BINDINGS.append(str(path.resolve().relative_to(LAB)))


def resolve(owner, name):
    if owner.parent.name == 'pencil_tracks_review' and name == 'pencil_tracks.py':
        return HERE / 'pencil_tracks' / name
    for path in (owner.parent / name, LAB / name, HERE / name):
        if path.is_file():
            return path
    raise AssertionError((str(owner), name))


def source_file(owner):
    return HERE / {
        'plot_metadata.json': 'plot_results.py',
        'trace_envelopes/results.json': 'trace_envelopes/trace_envelopes.py',
        'trace_envelopes_review/review_results.json': 'trace_envelopes_review/audit_certificates.py',
        'critical_third/scaling_results.json': 'critical_third/derive_scaling.py',
        'critical_third/inventory_free_results.json': 'critical_third/inventory_free.py',
    }[str(owner.relative_to(HERE))]


def bindings(owner, obj):
    if isinstance(obj, list):
        for item in obj:
            bindings(owner, item)
    elif isinstance(obj, dict):
        for key, value in obj.items():
            if key.endswith('_sha256') and isinstance(value, dict):
                for name, expected in value.items():
                    bind(resolve(owner, name), expected)
            elif 'sha256' in key and isinstance(value, str):
                if key == 'sha256':
                    target = resolve(owner, obj['path'])
                elif key == 'source_sha256':
                    target = source_file(owner)
                elif key == 'inventory_sha256':
                    target = resolve(owner, obj['inventory'])
                elif key == 'results_sha256':
                    folder = owner.parent.name.removesuffix('_review')
                    target = HERE / folder / 'results.json'
                elif key in ('reviewed_results_sha256', 'reviewed_source_sha256'):
                    target = HERE / 'trace_envelopes' / ('results.json' if 'results' in key else 'trace_envelopes.py')
                elif key == 'q_polynomials_sha256':
                    target = HERE / 'critical_third/q_polynomials.json'
                elif key == 'helper_sha256':
                    target = HERE / 'critical_third/derive_scaling.py'
                elif key == 'candidate_sha256':
                    target = HERE / 'pencil_tracks/pencil_tracks.py'
                elif key == 'reviewer_sha256':
                    target = owner.parent / 'review_tracks.py'
                elif key == 'field_oracle_sha256':
                    target = LAB / 'round4/novelty/review_support_batches.py'
                else:
                    raise AssertionError((str(owner), key))
                bind(target, value)
            else:
                bindings(owner, value)


def main():
    previous = json.loads((LAB / 'round5/manifest.json').read_text())
    for name, record in previous['artifacts'].items():
        bind(LAB / 'round5' / name, record['sha256'])
    nested = 0
    for folder in ('critical_third', 'overlap_preflight', 'trace_envelopes_review'):
        for record in read(folder + '/manifest.json')['files']:
            bind(HERE / folder / record['path'], record['sha256'])
            nested += 1
    for path in HERE.rglob('*.json'):
        if path.name != 'manifest.json':
            bindings(path, json.loads(path.read_text()))
    sources = list(HERE.rglob('*.py'))
    for path in sources:
        compile(path.read_text(), str(path), 'exec')
    for name in ('trace_envelopes/verification.json', 'trace_envelopes_review/review_results.json',
                 'critical_third_review/results.json', 'critical_third_review/enclosure_results.json',
                 'pencil_tracks/verification.json', 'pencil_tracks_review/review_tracks.json',
                 'pencil_tracks_review/review_discovery.json', 'cover_barriers/verification.json',
                 'cover_barriers_review/results.json', 'overlap_preflight/verification.json'):
        expected = {'pencil_tracks_review/review_tracks.json': 'whole-field track compiler independently verified',
                    'pencil_tracks_review/review_discovery.json': 'discovery repair and failure claims independently verified'}.get(name, 'passed')
        assert read(name)['status'] == expected, name
    trace = read('trace_envelopes_review/review_results.json')
    assert (trace['pointwise_dual_inequalities'], trace['primal_equalities'], trace['corrupted_certificates_rejected']) == (5208, 104, 56)
    critical = read('critical_third_review/results.json')
    assert critical['raw_union_equalities'] == 361 and critical['second_moment_union_equalities'] == 39
    assert len(read('critical_third_review/enclosure_results.json')['large_cases']) == 2
    barrier = read('cover_barriers_review/results.json')
    assert barrier['saved_ledger_sample_base_rows_read'] == 893464
    assert barrier['joint_family_obstructions'] == 16
    overlap = read('overlap_preflight/results.json')['results']
    assert [r['after_redundancy_deletion'] for r in overlap] == [118, 49]
    assert all(r['complete_output_equal_to_frozen_s_subset_oracle'] for r in overlap)
    for row in read('pencil_tracks/results.json')['large']:
        cert = row['result']['certificate']
        assert cert['complete_symbolic_node_count'] in (131074, 3)
    broken = []
    for path in HERE.rglob('*.md'):
        prose = re.sub(r'```.*?```|`[^`\n]*`', '', path.read_text(), flags=re.S)
        for target in re.findall(r'\]\(([^)]+)\)', prose):
            if target.startswith(('http://', 'https://', '#', 'mailto:')):
                continue
            linked = path.parent / re.sub(r':\d+$', '', target.split('#')[0].strip('<>'))
            if linked != HERE / 'manifest.json' and not linked.exists():
                broken.append((str(path), target))
    assert not broken, broken
    artifacts = {str(path.relative_to(HERE)): {'sha256': digest(path), 'bytes': path.stat().st_size}
                 for path in sorted(HERE.rglob('*')) if path.is_file() and '__pycache__' not in path.parts
                 and path.name != 'verification.log' and path != HERE / 'manifest.json'}
    out = {'created_utc': datetime.now(timezone.utc).isoformat(), 'status': 'passed',
           'scope': 'Saved-evidence integration; completed reviews bound to current sources and inputs. No full census rerun, human review, prize proof, or historical novelty certification.',
           'source_and_input_bindings': len(BINDINGS), 'previous_round5_artifacts_preserved': len(previous['artifacts']),
           'nested_snapshot_files': nested, 'python_sources_compiled': len(sources),
           'markdown_local_links_passed': True, 'artifacts': artifacts}
    (HERE / 'manifest.json').write_text(json.dumps(out, indent=2) + '\n')
    print(f'Round6 verification passed: {len(artifacts)} artifacts, {len(BINDINGS)} source/input bindings; all {len(previous["artifacts"])} round5 artifacts preserved.', flush=True)


if __name__ == '__main__':
    main()
