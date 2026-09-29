#!/usr/bin/env python3
"""Bind completed exact evidence to current source/input bytes and cross-check it.

This fast integration audit reads saved independent reviews; it does not claim
to re-enumerate their large censuses. Lane READMEs give full replay commands.
"""
from datetime import datetime, timezone
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import re
import sys
import warnings

HERE = Path(__file__).resolve().parent
BINDINGS = []


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def read(name):
    return json.loads((HERE/name).read_text())


def bind(path, expected):
    assert path.is_file(), path
    assert digest(path) == expected, path
    BINDINGS.append(str(path.relative_to(HERE.parent)))


def bind_map(folder, values):
    for name, expected in values.items():
        bind(HERE/folder/name, expected)


def main():
    sources = list(HERE.rglob('*.py'))
    with warnings.catch_warnings():
        # A historical archived docstring contains a literal backslash.
        # Syntax errors still fail; do not edit the frozen source snapshot.
        warnings.simplefilter('ignore', SyntaxWarning)
        for path in sources:
            compile(path.read_text(), str(path), 'exec')
    nested_count = 0
    for folder in ('marked_counts','marked_counts_ablation','coefficient_structure'):
        manifest = read(folder+'/manifest.json')
        for row in manifest['files']:
            bind(HERE/folder/row['path'], row['sha256'])
            nested_count += 1
    for folder, names in {
        'marked_moments': ['results.json','contraction_envelopes.json'],
        'support_batches': ['results.json','verification_results.json','iteration_results.json'],
        'scalar_fibers': ['results.json','degree_six_results.json'],
        'novelty': ['direct_conditional_oracle.json'],
    }.items():
        for name in names:
            bind_map(folder, read(folder+'/'+name)['source_sha256'])
    bind_map('support_batches', read('support_batches/initial_results.json')['archive_source_sha256'])

    review = read('novelty/conditional_compiler_review.json')
    assert review['status'] == 'passed' and review['saved_oracle_records'] == 2454
    bind(HERE/'marked_moments/conditional_moments.py', review['compiler_sha256'])
    bind(HERE/'novelty/direct_conditional_oracle.json', review['oracle_sha256'])
    bind(HERE/'novelty/review_conditional_compiler.py', review['reviewer_sha256'])
    review = read('novelty/contraction_review.json')
    assert review['status'] == 'passed' and review['complete_subset_moment_checks'] == 18
    bind(HERE/'marked_moments/contraction_reduction.py', review['candidate_sha256'])
    bind(HERE/'novelty/review_contraction_reduction.py', review['reviewer_sha256'])
    review = read('novelty/envelope_review.json')
    assert review['status'] == 'passed' and review['rational_rounding_checks'] == 424
    for name, expected in review['source_sha256'].items():
        bind(HERE/('novelty' if name.startswith('review_') else 'marked_moments')/name, expected)
    review = read('novelty/support_batches_review.json')
    assert review['status'] == 'passed' and review['complete_polynomial_scalar_comparisons'] == 261
    bind_map('support_batches', review['candidate_sha256'])
    bind(HERE/'novelty/review_support_batches.py', review['reviewer_sha256'])
    review = read('novelty/scalar_fibers_review.json')
    assert review['status'] == 'passed' and review['complete_polynomial_scalar_comparisons'] == 1700
    bind(HERE/'scalar_fibers/scalar_fibers.py', review['candidate_sha256'])
    bind(HERE/'scalar_fibers/results.json', review['results_sha256'])
    bind(HERE/'novelty/review_scalar_fibers.py', review['reviewer_sha256'])
    bind(HERE/'novelty/review_support_batches.py', review['independent_field_oracle_sha256'])

    # The optimized backend and full relation table must describe the same
    # cells and integer contraction; moment readbacks must obey the identity.
    count = 0
    for path in sorted((HERE/'marked_counts_ablation').glob('contraction_p*_marks_*.json')):
        row = json.loads(path.read_text())
        provenance = row['provenance']
        bind(HERE/provenance['input_file'], provenance['input_sha256'])
        for key, name in [('cpp_sha256','single_contraction.cpp'),('header_sha256','ntt_exact.hpp'),
                          ('binary_sha256','single_contraction'),('verifier_sha256','verify_contraction.py')]:
            bind(HERE/'marked_counts_ablation'/name, provenance[key])
        bind(HERE/'marked_moments/contraction_reduction.py', provenance['reduction_source_sha256'])
        bind(HERE/'marked_moments/conditional_moments.py', provenance['conditional_compiler_source_sha256'])
        full = read(provenance['input_file'])
        if 'pair' in full:
            full = next(x['record'] for x in full['pair'] if x['marks'] == row['marks'])
        assert row['cells'] == full['cells']
        Q = 0
        for rel in full['relations']:
            u, v = rel['left_pattern'], rel['right_pattern']
            if 0 in u or 0 in v:
                continue
            h2 = u[0]*u[1]+u[0]*u[2]+u[1]*u[2]
            h3 = v[0]*v[1]*v[2]
            Q += rel['count']*rel['sign']*h2*h3
        assert Q == row['Q']
        assert row['metadata']['transforms'] == 3
        for result in row['reduced_moments']:
            assert Fraction(*result['second_moment']) == Fraction(*result['cell_only_term']) + Fraction(*result['Q_coefficient'])*Q
            assert Fraction(*result['variance']) == Fraction(*result['second_moment'])-Fraction(*result['mean'])**2
        count += 1
    assert count == 13
    summary = read('marked_counts_ablation/summary.json')
    assert summary['total_direct_completions'] == 6363136
    bind(HERE/'marked_counts_ablation'/summary['contraction_validation_file'], summary['contraction_validation_sha256'])

    envelopes = read('marked_moments/contraction_envelopes.json')['cases']
    for row in envelopes:
        assert row['Q']**2 <= Fraction(*row['Q_bound_squared'])
        radius2 = Fraction(*row['cell_only_envelope_radius_squared'])
        assert radius2 == Fraction(*row['Q_coefficient'])**2*Fraction(*row['Q_bound_squared'])
        assert Fraction(*row['envelope_radius_rational_upper'])**2 >= radius2
    large = {row['q']: row for row in envelopes}
    assert Fraction(*large[65537]['relative_radius_rational_upper']) < Fraction(34,10**13)
    assert Fraction(*large[1000033]['relative_radius_rational_upper']) < Fraction(12,10**16)

    six = read('scalar_fibers/degree_six_results.json')
    assert six['status'] == 'passed' and six['full_scalar_constant_codeword_pairs'] == 531441
    assert six['certificate']['complete_symbolic_node_count'] == 1457
    assert read('support_batches/verification_results.json')['status'] == 'passed'
    assert read('support_batches/iteration_results.json')['status'] == 'passed'

    formula_review = read('coefficient_structure_review/review_results.json')
    assert formula_review['status'] == 'passed' and formula_review['coefficient_formula_exactly_equal']
    bind_map('coefficient_structure', formula_review['source_sha256'])
    bind(HERE/'coefficient_structure_review/review_coefficient_structure.py', formula_review['review_sha256'])
    scaling_review = read('coefficient_structure_review/critical_scaling_review.json')
    assert scaling_review['status'] == 'passed'
    for name, expected in scaling_review['source_sha256'].items():
        bind(HERE/('coefficient_structure' if name=='symbolic_formula.json' else 'coefficient_structure_review')/name, expected)
    bind(HERE/'coefficient_structure_review/review_critical_scaling.py', scaling_review['review_sha256'])
    all_n = read('coefficient_structure/all_n_validation.json')
    assert all_n['total_exact_cardinality_checks'] == 1504
    for name,key in [('verify_all_n.py','source_sha256'),('closed_coefficient.py','closed_formula_sha256')]:
        bind(HERE/'coefficient_structure'/name, all_n[key])
    bind(HERE/'marked_moments/conditional_moments.py', all_n['compiler_sha256'])
    bind(HERE/'marked_moments/contraction_reduction.py', all_n['reduction_sha256'])
    for result, source in [('sensitivity.json','sensitivity.py'),('plot_metadata.json','plot_sensitivity.py'),
                           ('symbolic_formula.json','symbolic_formula.py'),('fixed_q_polynomials.json','coefficient_polynomial.py')]:
        bind(HERE/'coefficient_structure'/source, read('coefficient_structure/'+result)['source_sha256'])
    sys.path.insert(0,str(HERE/'coefficient_structure'))
    from closed_coefficient import contraction_coefficient_degree6
    for row in envelopes:
        assert contraction_coefficient_degree6(row['q'],row['n']) == Fraction(*row['Q_coefficient'])

    # Verify the prior mathematical snapshot. Its one documentation amendment
    # is explicit; it is not a claim of a new round3 mathematical replay.
    previous = json.loads((HERE.parent/'round3/manifest.json').read_text())
    for name, row in previous['artifacts'].items():
        bind(HERE.parent/'round3'/name, row['sha256'])

    broken = []
    for path in HERE.rglob('*.md'):
        prose = re.sub(r'```.*?```|\\\[.*?\\\]', '', path.read_text(), flags=re.S)
        prose = re.sub(r'`[^`\n]*`', '', prose)
        for target in re.findall(r'\]\(([^)]+)\)', prose):
            if target.startswith(('https://','http://','#','mailto:')):
                continue
            local_target = re.sub(r':\d+$', '', target.split('#')[0].strip('<>'))
            linked = path.parent/local_target
            if linked == HERE/'manifest.json':
                continue  # This successful invocation writes it below.
            if not linked.exists():
                broken.append((str(path.relative_to(HERE)),target))
    assert not broken, broken
    artifacts = {str(path.relative_to(HERE)):{'sha256':digest(path),'bytes':path.stat().st_size}
                 for path in sorted(HERE.rglob('*')) if path.is_file() and '__pycache__' not in path.parts
                 and path.name != 'verification.log' and path != HERE/'manifest.json'}
    manifest = {'created_utc':datetime.now(timezone.utc).isoformat(),
                'status':'saved independent exact evidence and current cross-lane bindings passed; no prize or historical-novelty certification',
                'verification_kind':'source/input binding and exact saved-evidence integration; full independent censuses were previously completed',
                'python_sources_compiled':len(sources), 'nested_snapshot_files_checked':nested_count,
                'source_and_input_bindings_checked':len(BINDINGS),
                'optimized_contraction_exports_cross_checked':count,
                'general_coefficient_independent_review_passed':True,
                'saved_every_cardinality_coefficient_checks':all_n['total_exact_cardinality_checks'],
                'markdown_local_links_passed':True,
                'artifacts':artifacts}
    (HERE/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(f'Round4 verification passed: {len(artifacts)} artifacts pinned, {len(BINDINGS)} source/input bindings, {count} optimized contraction exports.',flush=True)


if __name__ == '__main__':
    main()
