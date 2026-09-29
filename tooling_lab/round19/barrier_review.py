#!/usr/bin/env python3
"""Separate exact verification of the short-relation obstruction certificate."""
from hashlib import sha256
import json
from pathlib import Path
import sys
import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'round17'))
from carry_review import product, result_norm


def main():
    source = HERE/'relation_barrier.json'; data = json.loads(source.read_text()); cases = []
    t, s = sp.symbols('t s')
    for c in data['cases']:
        p, N, g, a = (c[key] for key in ('p', 'N', 'g', 'a')); m = p//2
        assert sp.isprime(p) and N >= 4 and N & (N-1) == 0 and pow(g, N, p) == p-1
        expected = []
        for j in range(N):
            v = a*pow(g, j, p) % p; expected.append(v if 2*v < p else v-p)
        assert a == m+1 and expected == c['centered_F'] and expected[0] == -m
        escaped = c['escaped_F']; assert escaped == [m+1]+expected[1:]
        assert all((v-a*pow(g, j, p)) % p == 0 for j, v in enumerate(escaped))
        delta = min(m-abs(v) for v in escaped[1:]); budget = c['universal_relation_l1_budget']
        assert delta == c['other_coordinate_slack'] >= 1 and budget == delta+1 < p
        lam, mu = c['dual_multipliers']; assert lam >= 0 and mu >= 0
        # For non-axis normals, s>=1, t+s<=budget. This exact polynomial
        # identity certifies delta*s-t>=0 for every such t,s.
        assert sp.expand(delta*s-t-lam*(s-1)-mu*(budget-t-s)) == 0
        # The only omitted case is an axis normal. Its modular kernel condition
        # forces its nonzero coefficient to be divisible by p, impossible below p.
        assert all(pow(g, j, p) != 0 for j in range(N))
        f, k = c['relation'], c['relation_cofactor']; u = pow(g, -1, p)
        assert any(f) and sum(v*pow(u, j, p) for j, v in enumerate(f)) % p == 0
        l1 = sum(map(abs, f)); assert l1 == c['relation_l1'] <= budget
        bound = m*l1//p; assert bound == c['digit_bound']
        nf = result_norm(f, N); assert nf == c['relation_norm'] == p*k
        for prefix, F in [('centered', expected), ('escaped', escaped)]:
            D = c[prefix+'_digits']; assert product(f, F, N) == [p*v for v in D]
            assert max(map(abs, D)) <= bound
            nd = result_norm(D, N); assert nd == k*c[prefix+'_norm_defect']
            assert result_norm(F, N) == p**(N-1)*c[prefix+'_norm_defect']
        assert sum(v*v for v in escaped)-sum(v*v for v in expected) == c['coefficient_energy_difference'] == p
        cases.append({'p': p, 'N': N, 'certified_all_relation_l1_budget': budget,
                      'other_coordinate_slack': delta, 'independent_resultants': 5,
                      'dual_identity_passed': True, 'axis_case_excluded': True})
        print(json.dumps(cases[-1]), flush=True)
    control = data['finite_portfolio_control']; p, N, g = control['p'], control['N'], control['g']
    F = control['actual_F']; assert all(2*abs(v) < p and (v-pow(g, j, p)) % p == 0 for j, v in enumerate(F))
    aliases = control['aliases']; assert len(aliases) == 2*N
    for i, alias in enumerate(aliases):
        j, sign = alias['coordinate'], alias['sign']; assert (j, sign) == (i//2, (-1, 1)[i % 2])
        expected = F.copy(); expected[j] += sign*p; assert alias['F'] == expected
        assert abs(expected[j])*2 > p
    live = list(range(len(aliases))); counts = []
    for stage in control['stages']:
        f = stage['relation']; assert sum(v*pow(pow(g, -1, p), j, p) for j, v in enumerate(f)) % p == 0
        assert sum(map(abs, f)) == stage['l1'] and stage['bound'] == p//2*stage['l1']//p
        retained = []
        for i in live:
            raw = product(f, aliases[i]['F'], N); assert all(v % p == 0 for v in raw)
            if all(abs(v//p) <= stage['bound'] for v in raw): retained.append(i)
        assert retained == stage['surviving_alias_indices']; live = retained; counts.append(len(live))
    out = {'status': 'passed', 'scope': 'Exact prime/order and coordinate checks, a symbolic nonnegative dual identity covering every allowed non-axis normal, modular exclusion of all axis normals, and separate polynomial products/resultants. Ordinary certificate verification, not Lean formalization.',
           'cases': cases, 'finite_portfolio_aliases_checked': len(aliases),
           'finite_portfolio_relations_checked': len(control['stages']), 'finite_portfolio_survivors': counts,
           'source_sha256': {name: sha256((HERE/name).read_bytes()).hexdigest()
                             for name in ['barrier_review.py', '../round17/carry_review.py']},
           'input_sha256': {source.name: sha256(source.read_bytes()).hexdigest()}}
    (HERE/'barrier_review.json').write_text(json.dumps(out, indent=2)+'\n')


if __name__ == '__main__': main()
