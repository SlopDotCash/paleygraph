#!/usr/bin/env python3
"""Compile joint difference cuts, then separate surviving integer directions."""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from math import floor
from pathlib import Path
import json
from conditioned_erasure import discover, pair

HERE = Path(__file__).resolve().parent


def main():
    path = HERE/'metric_consecutive32.json'; data = json.loads(path.read_text()); c = data['certificate']
    V = c['visible_directions']; N, f, U, p = c['N'], c['relation'], c['erased'], c['p']
    K = c['conditioning']['known_coordinates']
    rows = [[F(f[i-t] if i >= t else -f[N+i-t]) for t in range(N)] for i in K]
    Q = [[F(*v) for v in c['centered_coordinate_forms'][j]] for j in V]
    bounds = [floor(2*min(F(*c['coordinate_radii'][j]), F(*c['conditioning']['cuts'][j]['radius']))) for j in V]
    candidates = [list(t) for t in product(*(range(-b,b+1) for b in bounds))
                  if any(t) and next(v for v in t if v) > 0]
    initial = len(candidates); pairs = []
    for i,j in combinations(range(len(V)), 2):
        for sign in [1, -1]:
            direction = [0]*len(V); direction[i] = 1; direction[j] = sign
            cut = discover([Q[i][t]+sign*Q[j][t] for t in range(N)], rows)
            bound = floor(F(p-1,p)*F(*cut['l1']))
            before = len(candidates)
            candidates = [t for t in candidates if abs(t[i]+sign*t[j]) <= bound]
            pairs.append({'direction': direction, 'certificate': cut, 'integer_difference_bound': bound,
                          'candidates_removed': before-len(candidates)})
        if i % 4 == 0: print(json.dumps({'phase':'pair_cuts', 'pair':[i,j], 'remaining':len(candidates)}), flush=True)
    after_pairs = len(candidates); separators = []; feasible = []; unresolved = []
    # A bounded pilot; no unprocessed vector may be called impossible.
    budget = 128
    for index, t in enumerate(candidates):
        if index >= budget: unresolved.append(t); continue
        pivot = next(i for i, value in enumerate(t) if value)
        others = [j for j in range(len(V)) if j != pivot]
        q0 = [q/F(t[pivot]) for q in Q[pivot]]
        extra = [[Q[j][a]-F(t[j],t[pivot])*Q[pivot][a] for a in range(N)] for j in others]
        cut = discover(q0, rows+extra); multipliers = [F(*v) for v in cut['lambda']]
        mu = [F(0)]*len(V); mu[pivot] = F(1,t[pivot])
        for j, value in zip(others, multipliers[len(K):]):
            mu[j] = -value; mu[pivot] += F(t[j],t[pivot])*value
        assert sum(a*b for a,b in zip(mu,t)) == 1
        certificate = {'target':t, 'pivot':pivot, 'certificate':cut, 'direction':[pair(v) for v in mu]}
        norm = F(*cut['l1'])
        if F(p-1,p)*norm < 1:
            separators.append(certificate)
        elif cut['optimality_certified'] and norm:
            point = [F(*v)/norm for v in cut['dual']]
            assert all(abs(v) <= F(p-1,p) for v in point)
            assert all(sum(row[a]*point[a] for a in range(N)) == 0 for row in rows)
            assert all(sum(Q[j][a]*point[a] for a in range(N)) == t[j] for j in range(len(V)))
            certificate['normalized_centered_difference'] = [pair(v) for v in point]
            feasible.append(certificate)
        else:
            unresolved.append(t)
        print(json.dumps({'phase':'separation', 'index':index, 'separated':len(separators), 'continuous_feasible':len(feasible), 'unresolved':len(unresolved)}), flush=True)
    out = {'status':'produced', 'scope':'One representative of each nonzero signed visible integer difference; continuous witnesses are not actual codeword pairs.',
           'visible_directions':V, 'individual_difference_bounds':bounds, 'initial_sign_representatives':initial,
           'pair_cuts':pairs, 'remaining_after_pairs':after_pairs, 'target_budget':budget,
           'target_separators':separators, 'continuous_feasible_targets':feasible, 'unresolved_targets':unresolved,
           'universal_unique_completion':not feasible and not unresolved,
           'input_sha256':{path.name:sha256(path.read_bytes()).hexdigest()},
           'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['coupled_experiments.py','conditioned_erasure.py']}}
    (HERE/'coupled_summary.json').write_text(json.dumps(out, separators=(',', ':'))+'\n')
    print(json.dumps({k:out[k] for k in ['initial_sign_representatives','remaining_after_pairs','universal_unique_completion']}), flush=True)


if __name__ == '__main__': main()
