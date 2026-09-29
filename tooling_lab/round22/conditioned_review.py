#!/usr/bin/env python3
"""Separate certificate review. No conditioned producer/decoder imports."""
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
from math import prod
from pathlib import Path
import json
import sys
import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'round21'))
from erasure_review import certify as certify_digit, rational
from pullback_review import certify_pullback


def certify(c):
    b = deepcopy(c)
    b['universal_candidate_box_cap'] = c['unconditioned_universal_cap']
    b['universal_unique_completion'] = c['unconditioned_unique_completion']
    if not c['erased']:
        # Avoid a shape ambiguity in the frozen empty SymPy matrix interface.
        assert c['centered_coordinate_forms'] == c['centered_coordinate_radii'] == []
        prepared = certify_digit(b)
    else:
        prepared = certify_pullback(b)
    N, f, U = c['N'], c['relation'], c['erased']; d = len(U)
    K = [i for i in range(N) if i not in U]
    assert c['conditioning']['known_coordinates'] == K
    assert c['conditioning']['rule'] == 'intersect_distinct_centered_integer_intervals'
    cuts = c['conditioning']['cuts']; assert len(cuts) == d
    M = sp.Matrix([[f[(i-t) % N]*(1 if i >= t else -1) for t in range(N)] for i in range(N)])
    T = sp.Matrix(d, d, lambda i,j: rational(c['inverse'][i][j]))
    Q = M[U, :].T*T if d else sp.zeros(N, 0)
    reduced = []
    for j, cut in enumerate(cuts):
        lam = [rational(v) for v in cut['lambda']]; assert len(lam) == len(K)
        residual = [rational(v) for v in cut['residual']]; assert len(residual) == N
        wanted = Q[:, j] - M[K, :].T*sp.Matrix(len(K), 1, lam)
        assert residual == [F(v) for v in wanted]
        norm = sum(map(abs, residual)); assert norm == rational(cut['l1'])
        radius = F(c['p']//2, c['p'])*norm; assert radius == rational(cut['radius'])
        reduced.append(min(prepared[1][j], radius))
        if cut['optimality_certified']:
            dual = [rational(v) for v in cut['dual']]; assert len(dual) == N
            assert all(abs(v) <= 1 for v in dual)
            assert M[K, :]*sp.Matrix(dual) == sp.zeros(len(K), 1)
            assert (Q[:, j].T*sp.Matrix(dual))[0] == norm
        else:
            assert cut['dual'] is None
    visible = c['visible_directions']
    assert c['universal_candidate_box_cap'] == prod(int(2*reduced[j])+1 for j in visible)
    assert c['universal_unique_completion'] == all(2*reduced[j] < 1 for j in visible)
    return prepared


def replay(c, known, budget, prepared):
    inv, radii, encoder, _ = prepared; U = c['erased']; p, k = c['p'], c['k']
    empty = {'status': 'complete', 'count': 0, 'completions': [], 'candidate_box_size': 0, 'scalar_candidates_checked': 0}
    if any(abs(v) > c['digit_bound'] for v in known.values()): return {**empty, 'reason': 'known_digit_height'}
    w = c['projection_weights']; dot = sum(w[i]*y for i,y in known.items()); rhs = -dot % k
    if rhs % c['kernel_gcd']: return {**empty, 'reason': 'syndrome_divisibility'}
    offset = [rhs//c['kernel_gcd']*v for v in c['bezout_lift']]
    boxes = []
    for j, cut in enumerate(c['conditioning']['cuts']):
        origin = sum(offset[i]*inv[i][j] for i in range(len(U)))
        shift = sum(rational(v)*known[i] for i, v in zip(c['conditioning']['known_coordinates'], cut['lambda']))
        r = rational(cut['radius'])
        low = max(-origin-radii[j], shift-origin-r)
        high = min(-origin+radii[j], shift-origin+r)
        lo = low.numerator//low.denominator
        if F(lo) < low: lo += 1
        boxes.append([lo, high.numerator//high.denominator])
    if any(lo > hi for lo, hi in boxes): return {**empty, 'reason': 'empty_integer_coordinate_interval', 'intervals': boxes}
    size = prod(boxes[j][1]-boxes[j][0]+1 for j in c['visible_directions'])
    if size > budget:
        return {'status': 'budget_exceeded', 'reason': 'candidate_box_too_large',
                'candidate_box_size': size, 'candidate_budget': budget, 'intervals': boxes}
    total = dot+sum(w[i]*v for i, v in zip(U, offset)); assert total % k == 0
    candidates = {total//k % p}
    for j in c['visible_directions']:
        candidates = {(a+t*c['scalar_steps'][j]) % p for a in candidates for t in range(boxes[j][0], boxes[j][1]+1)}
    scalars = sorted(candidates); words = encoder(scalars)
    found = [{'scalar': a, 'digits': word} for a, word in zip(scalars, words) if all(word[i] == y for i, y in known.items())]
    return {'status': 'complete', 'reason': 'exhausted_projected_candidate_box', 'count': len(found), 'completions': found,
            'candidate_box_size': size, 'scalar_candidates_checked': len(candidates), 'intervals': boxes}


def main():
    cases = []; inputs = {}; sample = None
    for summary_name in ['conditioned_summary.json', 'small_summary.json']:
        path = HERE/summary_name; summary = json.loads(path.read_text()); inputs[path.name] = sha256(path.read_bytes()).hexdigest()
        for row in summary['cases']:
            path = HERE/row['artifact']; data = json.loads(path.read_text()); inputs[path.name] = sha256(path.read_bytes()).hexdigest()
            c = data['certificate']; prepared = certify(c); scalar_count = 0; oracles = 0; changes = 0
            direct = prepared[2](list(range(c['p']))) if c['p'] < 100 else None
            for r in data['records']:
                known = dict(r['known']); expected = replay(c, known, data['candidate_budget'], prepared)
                assert expected == r['output']; scalar_count += expected.get('scalar_candidates_checked', 0)
                if direct is not None:
                    true = [a for a, word in enumerate(direct) if all(word[i] == y for i, y in known.items())]
                    assert expected['status'] == 'complete' and [v['scalar'] for v in expected['completions']] == true
                    oracles += 1
                changes += any(sum(rational(v)*known[i] for i, v in zip(c['conditioning']['known_coordinates'], cut['lambda'])) != 0 for cut in c['conditioning']['cuts'])
            out = {'name': data['name'], 'cuts_verified': len(c['erased']),
                   'optimality_certificates': sum(cut['optimality_certified'] for cut in c['conditioning']['cuts']),
                   'queries_replayed': len(data['records']), 'full_small_codebook_oracles': oracles,
                   'queries_with_nonzero_center_shift': changes, 'scalar_candidates_checked': scalar_count}
            cases.append(out); print(json.dumps(out), flush=True)
            if data['name'] == 'random32': sample = c
    rejected = []
    for name in ['false_lambda', 'false_residual', 'understated_radius', 'false_dual', 'wrong_known_order', 'false_cap']:
        c = deepcopy(sample); cut = c['conditioning']['cuts'][0]
        if name == 'false_lambda': cut['lambda'][0][0] += cut['lambda'][0][1]
        if name == 'false_residual': cut['residual'][0][0] += cut['residual'][0][1]
        if name == 'understated_radius': cut['radius'] = [0, 1]
        if name == 'false_dual': cut['dual'][0] = [2, 1]
        if name == 'wrong_known_order': c['conditioning']['known_coordinates'].reverse()
        if name == 'false_cap': c['universal_candidate_box_cap'] += 1
        try: certify(c)
        except AssertionError: rejected.append(name)
        else: raise AssertionError('corrupt certificate accepted: '+name)
    out = {'status': 'passed', 'cases': cases, 'corrupt_certificates_rejected': rejected,
           'input_sha256': inputs, 'source_sha256': {n: sha256((HERE/n).read_bytes()).hexdigest()
           for n in ['conditioned_review.py', '../round21/pullback_review.py', '../round21/erasure_review.py']}}
    (HERE/'conditioned_review.json').write_text(json.dumps(out, indent=2)+'\n')


if __name__ == '__main__': main()
