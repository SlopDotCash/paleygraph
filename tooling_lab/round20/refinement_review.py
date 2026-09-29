#!/usr/bin/env python3
"""Separate column-splitting and exact arithmetic review of adaptive suffixes."""
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def partition(rows):
    groups = [list(range(len(rows)))]; profile = []
    for col in range(len(rows[0])):
        refined = []
        for group in groups:
            for bit in '01':
                part = [i for i in group if rows[i][col] == bit]
                if part: refined.append(part)
        groups = refined
        if col+1 in (1, 2, 4, 8, 16, 32, 64, 128, 256): profile.append([col+1, len(groups)])
    return groups, profile


def membership(config, digits):
    p, g, f, N = config['p'], config['g'], config['relation'], config['N']; u = pow(g, -1, p)
    # Full dense convolution and direct modular powers; no producer imports.
    matrix = [[f[i-j] if i >= j else -f[N+i-j] for j in range(N)] for i in range(N)]
    def encode(a):
        F = [(a*pow(g, j, p)+p//2) % p-p//2 for j in range(N)]
        values = [sum(row[j]*F[j] for j in range(N)) for row in matrix]
        assert all(x % p == 0 for x in values)
        return [x//p for x in values]
    weights = [pow(u, j, p) for j in range(N)]
    beta = sum(x*y for x, y in zip(encode(1), weights)) % p; assert beta
    a = sum(x*y for x, y in zip(digits, weights))*pow(beta, -1, p) % p
    return encode(a) == digits


def verify_case(c, r):
    rows = r['cross_acceptance_rows']; groups, profile = partition(rows)
    reps = c['representative_prefix_indices']; assert len(groups) == len(reps) == c['all_supplied_suffixes_lower_bound']
    assert set(reps) == {min(group) for group in groups}
    cols = c['pair_witness_columns']; pos = 0
    for n, i in enumerate(reps):
        for j in reps[n+1:]:
            k = cols[pos]; assert 0 <= k < len(rows) and rows[i][k] != rows[j][k]; pos += 1
    assert pos == len(cols)
    a = c['adaptive_refinement']; arithmetic_words = 0
    if a['status'] == 'complete_for_supplied_prefixes':
        states = a['prefix_residual_state_ids']; observed = list(rows); got = [len(groups)]
        for step in a['additional_suffixes']:
            i, j = step['merged_prefix_indices']; assert observed[i] == observed[j] and states[i] != states[j]
            suffix = step['new_suffix_in_read_order']; assert len(suffix) == r['N']-r['cut']
            column = ''
            for original in r['words_in_original_coordinates']:
                word = original.copy()
                for k, value in zip(r['order'][r['cut']:], suffix): word[k] = value
                column += '1' if membership(r, word) else '0'; arithmetic_words += 1
            assert column == step['acceptance_column'] and column[i] != column[j]
            observed = [row+column[i] for i, row in enumerate(observed)]
            new_groups, _ = partition(observed); got.append(len(new_groups))
            assert got[-1] == step['observed_classes_after'] and got[-1] > got[-2]
        assert got == a['class_count_profile']
        final_groups, _ = partition(observed)
        assert all(len({states[i] for i in group}) == 1 for group in final_groups)
        assert len(final_groups) == len(set(states)) == a['distinct_true_residuals_of_supplied_prefixes']
    else: assert a['status'] == 'no_complete_residual_oracle' and not a['additional_suffixes']
    return {'name': r['name'], 'pair_witnesses_checked': len(cols), 'column_refinement_profile': profile,
            'adaptive_suffixes_checked': len(a['additional_suffixes']), 'adaptive_arithmetic_words': arithmetic_words,
            'final_supplied_prefix_classes': a.get('distinct_true_residuals_of_supplied_prefixes'),
            'certified_initial_lower_bound': len(groups)}


def main():
    source = HERE/'residual_refinement.json'; report = json.loads(source.read_text()); results = []
    bound = {source.name: sha256(source.read_bytes()).hexdigest()}; adaptive_records = []
    references = {'p65537_canonical': 'p65537_canonical.json', 'p65537_saved': 'p65537_saved.json', 'p6700417_canonical': 'symbolic_N32_k641.json'}
    for c in report['cases']:
        path = HERE/c['splice_artifact']; r = json.loads(path.read_text()); bound[path.name] = sha256(path.read_bytes()).hexdigest()
        if c['name'] in references:
            ref = HERE/references[c['name']]; data = json.loads(ref.read_text()); graph = data['graph']
            bound[ref.name] = sha256(ref.read_bytes()).hexdigest()
            assert c['exact_full_language_width'] == graph['widths'][r['cut']]
            true_states = []
            for word in r['words_in_original_coordinates']:
                s = graph['root']
                for d, k in enumerate(r['order'][:r['cut']]): s = dict(graph['layers'][d][s])[word[k]]
                true_states.append(s)
            assert true_states == c['adaptive_refinement']['prefix_residual_state_ids']
        else: assert c['exact_full_language_width'] is None
        B = r['p']//2*sum(map(abs, r['relation']))//r['p']
        upper = c['exact_full_language_width'] if c['exact_full_language_width'] is not None else min(r['p'], (2*B+1)**r['cut'])
        assert c['valid_width_upper_bound'] == upper and c['all_supplied_suffixes_lower_bound'] <= upper
        out = verify_case(c, r); results.append(out); print(json.dumps(out), flush=True)
        if c['adaptive_refinement']['additional_suffixes']: adaptive_records.append((c, r))
    assert adaptive_records
    bad = []
    c, r = adaptive_records[0]
    for label, target in [('wrong_distinguishing_column', 'witness'), ('wrong_adaptive_verdict', 'verdict')]:
        changed = deepcopy(c)
        if target == 'witness': changed['pair_witness_columns'][0] = len(r['scalars'])
        else:
            step = changed['adaptive_refinement']['additional_suffixes'][0]
            step['acceptance_column'] = ('1' if step['acceptance_column'][0] == '0' else '0')+step['acceptance_column'][1:]
        try: verify_case(changed, r)
        except (AssertionError, IndexError): bad.append(label)
        else: raise AssertionError('corruption accepted')
    out = {'status': 'passed', 'scope': 'All distinguishing columns checked by independent partition refinement. Every new adaptive suffix is replayed on every supplied prefix using direct scalar reconstruction and dense arithmetic. Complete diagrams certify only the supplied-prefix stopping criterion.',
           'cases': results, 'corrupt_artifacts_rejected': bad,
           'source_sha256': {'refinement_review.py': sha256(Path(__file__).read_bytes()).hexdigest()}, 'input_sha256': bound}
    (HERE/'refinement_review.json').write_text(json.dumps(out, indent=2)+'\n')


if __name__ == '__main__': main()
