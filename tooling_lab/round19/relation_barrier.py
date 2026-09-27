#!/usr/bin/env python3
"""Certificates obstructing all bounded-L1 short-relation digit inequalities."""
from hashlib import sha256
import json
from pathlib import Path
import sys
import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'round17'))
from norm_carry import multiply, norm


def centered(a, p, N, g):
    return [(a*pow(g, j, p)+p//2) % p-p//2 for j in range(N)]


def compressed(f, F, p):
    numerator = multiply(f, F)
    assert all(v % p == 0 for v in numerator)
    return [v//p for v in numerator]


def main():
    source = HERE.parent/'round17/norm_compression.json'; old = json.loads(source.read_text())
    supplied = [{'name': 'tiny', 'p': 17, 'N': 4, 'n': 8, 'g': 2, 'relation': [2, 0, 0, 1]}]
    seen = set()
    for c in old['cases']:
        key = (c['p'], c['g'], tuple(c['relation']))
        if key not in seen: supplied.append(c); seen.add(key)
    cases = []
    for c in supplied:
        p, N, g, f = (c[key] for key in ('p', 'N', 'g', 'relation'))
        assert sp.isprime(p) and pow(g, N, p) == p-1
        m = (p-1)//2; a = m+1; F = centered(a, p, N, g); assert F[0] == -m
        escaped = F.copy(); escaped[0] += p
        delta = m-max(abs(v) for v in F[1:]); budget = delta+1
        assert 1 <= delta and budget < p and sum(map(abs, f)) <= budget
        D = compressed(f, F, p); bad = compressed(f, escaped, p)
        bound = m*sum(map(abs, f))//p
        assert max(map(abs, D)) <= bound and max(map(abs, bad)) <= bound
        nf = norm(f); assert nf % p == 0; k = nf//p
        nd, nb = norm(D), norm(bad); assert nd % k == nb % k == 0
        row = {'name': c['name'], 'p': p, 'N': N, 'n': 2*N, 'g': g, 'a': a,
               'centered_F': F, 'escaped_F': escaped, 'other_coordinate_slack': delta,
               'universal_relation_l1_budget': budget,
               'dual_multipliers': [delta+1, 1],
               'relation': f, 'relation_l1': sum(map(abs, f)), 'digit_bound': bound,
               'centered_digits': D, 'escaped_digits': bad,
               'relation_norm': nf, 'relation_cofactor': k,
               'centered_norm_defect': nd//k, 'escaped_norm_defect': nb//k,
               'coefficient_energy_difference': sum(v*v for v in escaped)-sum(v*v for v in F)}
        assert row['coefficient_energy_difference'] == p
        cases.append(row)
        print(json.dumps({key: row[key] for key in ('name', 'p', 'N', 'other_coordinate_slack',
            'universal_relation_l1_budget', 'relation_l1', 'digit_bound', 'centered_norm_defect', 'escaped_norm_defect')}), flush=True)

    # A finite implementation control: a16-element short-relation portfolio
    # barely changes128 single-coordinate aliases of the actual a1 vector.
    c = cases[-1]; p, N, g = c['p'], c['N'], c['g']; u = pow(g, -1, p); d = 16
    basis = [[p]+[0]*(d-1)]
    for j in range(1, d):
        r = [-pow(u, j, p)]+[0]*(d-1); r[j] = 1; basis.append(r)
    reduced = sp.Matrix(basis).lll()
    relations = sorted([[int(v) for v in reduced.row(j)]+[0]*(N-d) for j in range(d)],
                       key=lambda f: (sum(map(abs, f)), sum(v*v for v in f), f))
    original = centered(1, p, N, g); aliases = []
    for j in range(N):
        for sign in (-1, 1):
            F = original.copy(); F[j] += sign*p
            aliases.append({'coordinate': j, 'sign': sign, 'F': F})
    live = list(range(len(aliases))); stages = []
    for f in relations:
        assert sum(v*pow(u, j, p) for j, v in enumerate(f)) % p == 0 and any(f)
        bound = p//2*sum(map(abs, f))//p
        assert max(map(abs, compressed(f, original, p))) <= bound
        assert max(map(abs, compressed(f, c['escaped_F'], p))) <= bound
        live = [i for i in live if max(map(abs, compressed(f, aliases[i]['F'], p))) <= bound]
        stages.append({'relation': f, 'l1': sum(map(abs, f)), 'bound': bound, 'surviving_alias_indices': live.copy()})
    portfolio = {'p': p, 'N': N, 'g': g, 'actual_a': 1, 'actual_F': original,
                 'aliases': aliases, 'stages': stages, 'complete_scope': 'Only the128 specified aliases; no completeness claim for all nonrealizable words.'}
    out = {'status': 'produced', 'scope': 'An arithmetic near-boundary witness plus a two-variable nonnegative slack identity proves that every short-relation coefficient bound up to the certified L1 budget admits this nonrealizable vector. This is a limitation of the stated inequality interface, not of all arithmetic tools.',
           'cases': cases, 'finite_portfolio_control': portfolio,
           'source_sha256': {name: sha256((HERE/name).read_bytes()).hexdigest()
                             for name in ['relation_barrier.py', '../round17/norm_carry.py']},
           'input_sha256': {'../round17/norm_compression.json': sha256(source.read_bytes()).hexdigest()}}
    (HERE/'relation_barrier.json').write_text(json.dumps(out, separators=(',', ':'))+'\n')
    print(json.dumps({'portfolio_relations': len(stages), 'aliases': len(aliases),
                      'survivors': len(live), 'survivors_by_stage': [len(r['surviving_alias_indices']) for r in stages]}), flush=True)


if __name__ == '__main__': main()
