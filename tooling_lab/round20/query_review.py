#!/usr/bin/env python3
"""Exercise the usable query API against separate arithmetic expectations."""
from hashlib import sha256
import json
from pathlib import Path
from complete_query import query
from completion_review import binary_counter, direct_codebook

HERE = Path(__file__).resolve().parent


def check(path, fixed, limit):
    r = json.loads(path.read_text()); out = query(path, fixed, limit); c = r['config']; N = len(r['order'])
    if c.get('Q'):
        counter = binary_counter(N, c['k'], fixed)
        expected = counter(0, 0, 0, 0)+counter(0, 1, 1, 0)+int(all(v == 0 for v in fixed.values()))
        direct = None
    else:
        direct = direct_codebook(c, list(range(N))).tolist()
        expected = sum(all(w[j] == v for j, v in fixed.items()) for w in direct)
    assert out['count'] == expected and len(out['completions']) == min(expected, limit)
    assert out['truncated'] == (expected > limit)
    for item in out['completions']:
        a, word = item['scalar_residue'], item['digits']
        assert all(word[j] == v for j, v in fixed.items())
        if direct is not None: assert word == direct[a]
        else:
            # Direct centered orbit at the (possibly composite) scalar modulus.
            modulus = item['scalar_modulus']; f = c['relation']; m = modulus//2
            F = [(a*pow(2, j, modulus)+m) % modulus-m for j in range(N)]
            numerator = [2*F[j]-(F[j+1] if j+1 < N else -F[0]) for j in range(N)]
            assert all(v % modulus == 0 for v in numerator)
            assert word == [v//modulus for v in numerator]
    return {'artifact': path.name, 'fixed': fixed, 'limit': limit, 'count': expected,
            'listed_words_checked': len(out['completions'])}


def main():
    results = []; bound = {}
    for name in ['p65537_saved', 'p65537_saved_reordered', 'p41_general', 'symbolic_N32_k641', 'symbolic_N128_k1']:
        path = HERE/(name+'.json'); r = json.loads(path.read_text()); N = len(r['order'])
        bound[path.name] = sha256(path.read_bytes()).hexdigest()
        for fixed, limit in [({}, 0), ({0: 0, N//2: 1, N-1: -1}, 5), ({N-1: 1000}, 5)]:
            results.append(check(path, fixed, limit))
    example = json.loads((HERE/'query_example.json').read_text()); fixed = {int(k): v for k, v in example['fixed_original_coordinates'].items()}
    result = check(HERE/example['diagram'], fixed, 5); assert result['count'] == example['count'] == 209394
    live = query(HERE/example['diagram'], fixed, 5); live['fixed_original_coordinates'] = {str(k): v for k, v in fixed.items()}
    assert live == example; results.append(result)
    bound['query_example.json'] = sha256((HERE/'query_example.json').read_bytes()).hexdigest()
    invalid = []
    for name, fixed, limit in [('negative_coordinate', {-1: 0}, 1), ('out_of_range', {32: 0}, 1), ('boolean_digit', {0: True}, 1), ('negative_limit', {}, -1)]:
        try: query(HERE/'symbolic_N32_k641.json', fixed, limit)
        except ValueError: invalid.append(name)
        else: raise AssertionError('invalid query accepted')
    out = {'status': 'passed', 'cases': results, 'invalid_queries_rejected': invalid,
           'source_sha256': {n: sha256((HERE/n).read_bytes()).hexdigest() for n in ['query_review.py', 'complete_query.py', 'completion_review.py', 'decision_diagram.py']},
           'input_sha256': bound}
    (HERE/'query_review.json').write_text(json.dumps(out, indent=2)+'\n'); print(json.dumps(out), flush=True)


if __name__ == '__main__': main()
