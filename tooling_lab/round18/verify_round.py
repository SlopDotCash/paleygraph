#!/usr/bin/env python3
"""Bind the digit-language experiment and preserve earlier rounds."""
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import re
import numpy as np
import sympy as sp

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
    previous = HERE.parent/'round17/checkpoint.json'; old = json.loads(previous.read_text())
    assert old['status'] == 'carry_compression_realizability_checkpoint_passed'
    for name, metadata in old['artifacts'].items(): bind(previous.parent/name, metadata['sha256'])
    older_path = HERE.parent/'round16/manifest.json'; older = json.loads(older_path.read_text())
    for name, metadata in older['artifacts'].items(): bind(older_path.parent/name, metadata['sha256'])
    reports = {}
    for path in HERE.glob('*.json'):
        if path.name == 'manifest.json': continue
        value = json.loads(path.read_text()); scan(path, value); reports[path.stem] = value
    assert reports['digit_language']['status'] == reports['orientation_adapter']['status'] == 'produced'
    assert reports['language_review']['status'] == reports['orientation_review']['status'] == 'passed'
    review = reports['language_review']; cases = review['cases']
    assert sum(c['actual_scalar_words_checked'] for c in cases) == 6766228
    assert sum(c['scalar_coordinate_transitions'] for c in cases) == 215464060
    assert sum(c['all_residue_profiles_checked'] for c in cases) == 644
    assert sum(c['words_checked'] for c in review['small_controls']) == 6642
    assert (cases[-1]['minimum_nonzero_support'], cases[-1]['maximum_support']) == (5, 27)
    assert cases[-1]['actual_support_distribution'][5] == review['low_support_resultants'] == 512
    assert review['saved_g2_records_checked'] == 130 and review['mutations_checked'] == 48
    assert len(review['necessary_condition_controls']) == 3
    assert sum(c['vectors_checked'] for c in reports['orientation_review']['cases']) == 451
    assert reports['orientation_review']['cases'][-1]['status'] == 'generator_2_has_wrong_order'
    assert reports['orientation_review']['nonassociate_control_checked']
    assert reports['digit_language']['low_support']['norm_defect_counts_by_support'] == {
        '5': [[449, 128], [1217, 64], [8513, 64], [15937, 64], [24001, 64], [40193, 64], [84481, 64]]}
    sources = list(HERE.glob('*.py'))
    for path in sources: compile(path.read_text(), str(path), 'exec')
    for path in [*HERE.glob('*.md'), HERE.parent/'NEXT_ITERATION.md']:
        prose = re.sub(r'```.*?```|`[^`\n]*`', '', path.read_text(), flags=re.S)
        for target in re.findall(r'\]\(([^)]+)\)', prose):
            if target.startswith(('http://', 'https://', '#')): continue
            linked = path.parent/target.split('#')[0].strip('<>')
            assert linked.exists() or linked == HERE/'manifest.json', (str(path), target)
    artifacts = {str(path.relative_to(HERE)): {'sha256': digest(path), 'bytes': path.stat().st_size}
                 for path in sorted(HERE.rglob('*')) if path.is_file() and '__pycache__' not in path.parts
                 and path.name not in ('manifest.json', 'verification.log')}
    output = {'status': 'verified', 'created_utc': datetime.now(timezone.utc).isoformat(),
              'scope': 'Exact specialized digit-language counts, all-scalar support review, minimum-support norms and supported orientation adapters. No general short-relation automaton or uniform arithmetic estimate.',
              'previous_round17_artifacts_preserved': len(old['artifacts']),
              'previous_round16_artifacts_preserved': len(older['artifacts']),
              'source_and_input_bindings': len(BINDINGS), 'python_sources_compiled': len(sources),
              'numpy_version': np.__version__, 'sympy_version': sp.__version__,
              'markdown_local_links_passed': True,
              'previous_checkpoint_sha256': {'../round17/checkpoint.json': digest(previous)},
              'actual_scalars_checked': 6766228, 'scalar_coordinate_transitions_checked': 215464060,
              'minimum_support_norms_checked': 512, 'orientation_records_checked': 451,
              'artifacts': artifacts}
    (HERE/'manifest.json').write_text(json.dumps(output, indent=2)+'\n')
    print(f"Round18 verified: {len(artifacts)} artifacts, {len(BINDINGS)} bindings; all {len(old['artifacts'])} round17 and {len(older['artifacts'])} round16 artifacts preserved.", flush=True)


if __name__ == '__main__': main()
