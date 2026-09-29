#!/usr/bin/env python3
"""Refine residual classes using every supplied suffix, not just own suffixes.

Different row signatures are sufficient evidence; equal sampled signatures
never establish equality of the full residual languages.
"""
from hashlib import sha256
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def adaptive_suffixes(r, graph):
    layers = graph['layers']; N = r['N']; cut = r['cut']; order = r['order']
    def run(state, depth, word):
        for j, digit in enumerate(word, depth):
            state = dict(layers[j][state]).get(digit)
            if state is None: return None
        return state
    states = [run(0, 0, [w[j] for j in order[:cut]]) for w in r['words_in_original_coordinates']]
    assert all(s is not None for s in states)
    def any_suffix(depth, state):
        answer = []
        for j in range(depth, N):
            digit, state = layers[j][state][0]; answer.append(digit)
        return answer
    def witness(depth, a, b):
        assert a != b and depth < N
        left, right = dict(layers[depth][a]), dict(layers[depth][b])
        for digit in sorted(left.keys() | right.keys()):
            if digit not in left: return [digit]+any_suffix(depth+1, right[digit])
            if digit not in right: return [digit]+any_suffix(depth+1, left[digit])
            if left[digit] != right[digit]: return [digit]+witness(depth+1, left[digit], right[digit])
        raise AssertionError('different minimized nodes have identical residuals')
    observed = list(r['cross_acceptance_rows']); rounds = []; profile = [len(set(observed))]
    while True:
        groups = {}; pair = None
        for i, row in enumerate(observed):
            if row in groups and states[groups[row]] != states[i]: pair = (groups[row], i); break
            groups.setdefault(row, i)
        if pair is None: break
        i, j = pair; suffix = witness(cut, states[i], states[j])
        column = ''.join('1' if run(s, cut, suffix) is not None else '0' for s in states)
        assert column[i] != column[j]
        observed = [row+column[a] for a, row in enumerate(observed)]
        profile.append(len(set(observed)))
        rounds.append({'merged_prefix_indices': [i, j], 'new_suffix_in_read_order': suffix,
                       'acceptance_column': column, 'observed_classes_after': profile[-1]})
    assert profile[-1] == len(set(states))
    return {'status': 'complete_for_supplied_prefixes', 'prefix_residual_state_ids': states,
            'supplied_prefixes': len(states), 'distinct_true_residuals_of_supplied_prefixes': len(set(states)),
            'class_count_profile': profile, 'additional_suffixes': rounds,
            'scope': 'Equivalence oracle is the already certified complete diagram. Completeness here concerns supplied prefixes only.'}


def main():
    source = HERE/'splice_summary.json'; summary = json.loads(source.read_text())
    bound = {source.name: sha256(source.read_bytes()).hexdigest()}; results = []
    references = {'p65537_canonical': 'p65537_canonical.json', 'p65537_saved': 'p65537_saved.json',
                  'p6700417_canonical': 'symbolic_N32_k641.json'}
    for c in summary['cases']:
        path = HERE/c['artifact']; r = json.loads(path.read_text()); bound[path.name] = sha256(path.read_bytes()).hexdigest()
        rows = r['cross_acceptance_rows']; reps = {}; representatives = []
        for i, row in enumerate(rows):
            if row not in reps: reps[row] = i; representatives.append(i)
        witnesses = []
        for n, i in enumerate(representatives):
            for j in representatives[n+1:]:
                column = next(k for k in range(len(rows)) if rows[i][k] != rows[j][k])
                witnesses.append(column)
        exact_width = None; adaptive = {'status': 'no_complete_residual_oracle', 'additional_suffixes': []}
        if r['name'] in references:
            ref = HERE/references[r['name']]; graph = json.loads(ref.read_text())['graph']
            exact_width = graph['widths'][r['cut']]; bound[ref.name] = sha256(ref.read_bytes()).hexdigest()
            assert len(representatives) <= exact_width
            adaptive = adaptive_suffixes(r, graph)
        digit_bound = r['p']//2*sum(map(abs, r['relation']))//r['p']
        upper = exact_width if exact_width is not None else min(r['p'], (2*digit_bound+1)**r['cut'])
        assert r['certified_width_lower_bound'] <= len(representatives) <= upper
        result = {'name': r['name'], 'splice_artifact': path.name, 'cut': r['cut'],
                  'own_suffix_clique_lower_bound': r['certified_width_lower_bound'],
                  'all_supplied_suffixes_lower_bound': len(representatives),
                  'exact_full_language_width': exact_width, 'valid_width_upper_bound': upper,
                  'upper_bound_basis': 'complete minimal diagram' if exact_width is not None else 'at most p actual prefixes and at most (2B+1)^cut digit prefixes',
                  'representative_prefix_indices': representatives,
                  'pair_witness_columns': witnesses,
                  'adaptive_refinement': adaptive,
                  'witness_rule': 'Upper triangle in representative order; the indicated suffix column is accepted for exactly one of the two prefixes.'}
        results.append(result)
        print(json.dumps({k: result[k] for k in ['name', 'cut', 'own_suffix_clique_lower_bound', 'all_supplied_suffixes_lower_bound', 'exact_full_language_width', 'valid_width_upper_bound']}), flush=True)
    out = {'status': 'produced', 'scope': 'Finite observation-table refinement. Every differing row yields an explicit distinguishing suffix; no inference from equal observed rows to equal residual languages.',
           'cases': results, 'source_sha256': {'residual_refinement.py': sha256(Path(__file__).read_bytes()).hexdigest()},
           'input_sha256': bound}
    (HERE/'residual_refinement.json').write_text(json.dumps(out, separators=(',', ':'))+'\n')


if __name__ == '__main__': main()
