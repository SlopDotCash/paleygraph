#!/usr/bin/env python3
"""Freeze the completion/state certificates and preserve preceding artifacts."""
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
    for name, filename in [('round19', 'manifest.json'), ('round18', 'manifest.json'), ('round17', 'checkpoint.json'), ('round16', 'manifest.json')]:
        path = HERE.parent/name/filename; old = json.loads(path.read_text())
        for rel, entry in old['artifacts'].items(): bind(path.parent/rel, entry['sha256'])
        preserved[name] = len(old['artifacts'])
    reports = {}
    for path in HERE.glob('*.json'):
        if path.name == 'manifest.json': continue
        data = json.loads(path.read_text()); scan(path, data); reports[path.stem] = data
    for name in ['completion_summary', 'splice_summary', 'residual_refinement', 'plot_metadata']:
        assert reports[name]['status'] == 'produced'
    for name in ['completion_review', 'splice_review', 'refinement_review', 'query_review']:
        assert reports[name]['status'] == 'passed'
    completion = reports['completion_review']; assert len(completion['cases']) == 13
    assert sum(c['direct_scalar_words'] for c in completion['cases']) == 197520
    assert sum(c['binary_language_pairs'] for c in completion['cases']) == 37471
    assert sum(c['query_counts'] for c in completion['cases']) == 234
    assert sum(c['listed_completions'] for c in completion['cases']) == 1804
    assert len(completion['corrupt_artifacts_rejected']) == 5
    summary = {c['name']: c for c in reports['completion_summary']['cases']}
    assert summary['p65537_saved']['max_width'] == 1597 and summary['p65537_saved']['nodes'] == 6650
    assert summary['p65537_saved_reordered']['max_width'] == 5 and summary['p65537_saved_reordered']['nodes'] == 74
    large = summary['symbolic_N32_k641']
    assert (large['word_count'], large['nodes'], large['edges'], large['max_width']) == (6700417, 36381, 67668, 2565)
    for N in [8, 16, 64, 128]:
        c = summary[f'symbolic_N{N}_k1']
        assert c['nodes'] == 5*N-6 and c['edges'] == 11*N-17 and c['word_count'] == (1 << N)+1
    splice = reports['splice_review']; assert len(splice['cases']) == 8
    assert sum(c['cross_splices'] for c in splice['cases']) == 524288
    assert sum(c['scalar_encodings'] for c in splice['cases']) == 2048
    assert len(splice['corrupt_artifacts_rejected']) == 5
    refinement = reports['refinement_review']
    assert sum(c['adaptive_suffixes_checked'] for c in refinement['cases']) == 9
    assert sum(c['adaptive_arithmetic_words'] for c in refinement['cases']) == 2304
    assert len(refinement['corrupt_artifacts_rejected']) == 2
    assert refinement['cases'][1]['final_supplied_prefix_classes'] == 205
    assert len(reports['query_review']['cases']) == 16 and len(reports['query_review']['invalid_queries_rejected']) == 4
    assert reports['query_example']['count'] == 209394 and reports['query_example']['truncated']
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
           'scope': 'Complete fixed-order digit diagrams, exact partial completions, finite state lower bounds and certified adaptive suffix refinement. Efficient general residual equivalence and norm aggregates remain open.',
           'previous_artifacts_preserved': preserved, 'source_and_input_bindings': len(BINDINGS),
           'python_sources_compiled': len(sources), 'markdown_local_links_passed': True,
           'complete_scalar_words_checked': 197520, 'binary_language_pairs_checked': 37471,
           'cross_splices_checked': 524288, 'adaptive_suffix_words_checked': 2304,
           'corrupt_artifacts_rejected': 12, 'invalid_query_controls': 4,
           'artifacts': artifacts}
    (HERE/'manifest.json').write_text(json.dumps(out, indent=2)+'\n')
    print(f'Round20 verified: {len(artifacts)} artifacts, {len(BINDINGS)} bindings; preserved {preserved}.', flush=True)


if __name__ == '__main__': main()
