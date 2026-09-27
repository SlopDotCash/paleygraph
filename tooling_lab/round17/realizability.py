#!/usr/bin/env python3
"""Exact inverse certificates for small-digit centered coset encodings.

Discovery and small-cube census only. Separate matrix algebra in
realizability_review.py checks the inverse and every census decision.
"""
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import random
import sympy as sp
from norm_carry import multiply, norm

HERE = Path(__file__).resolve().parent


def norm_adjugate(a):
    """Quadratic tower descent; return n,A with a*A=n in Z[X]/(X^N+1)."""
    N = len(a)
    if N == 1:
        return a[0], [1]
    e, o = a[::2], a[1::2]
    ee, oo = multiply(e, e), multiply(o, o)
    yoo = [-oo[-1]] + oo[:-1]
    h = [x-y for x, y in zip(ee, yoo)]
    nh, ah = norm_adjugate(h)
    lifted = [0]*N
    lifted[::2] = ah
    conjugate = [v if j % 2 == 0 else -v for j, v in enumerate(a)]
    return nh, multiply(conjugate, lifted)


def prepare(p, g, f):
    N = len(f)
    if not (type(p) is int and sp.isprime(p) and p % 2 and
            N >= 4 and N & (N-1) == 0 and
            type(g) is int and 0 < g < p and pow(g, N, p) == p-1 and
            all(type(v) is int for v in f) and any(f)):
        raise ValueError('invalid arithmetic input')
    u = pow(g, -1, p)
    if sum(v*pow(u, j, p) for j, v in enumerate(f)) % p:
        raise ValueError('relation does not vanish at g inverse')
    nf, adj = norm_adjugate(f)
    assert nf > 0 and nf % p == 0
    assert multiply(f, adj) == [nf]+[0]*(N-1)
    return {'p': p, 'g': g, 'N': N, 'f': f, 'norm_f': nf,
            'k': nf//p, 'adj': adj,
            'bound': (p-1)*sum(map(abs, f))//(2*p)}


def classify(c, digits):
    N, p, k = c['N'], c['p'], c['k']
    if len(digits) != N or any(type(v) is not int for v in digits):
        return {'status': 'invalid_digits'}
    if any(abs(v) > c['bound'] for v in digits):
        return {'status': 'digit_height'}
    numerator = multiply(c['adj'], digits)
    if any(v % k for v in numerator):
        return {'status': 'nonintegral_inverse'}
    F = [v//k for v in numerator]
    a = F[0] % p
    scalar = all(v % p == a*pow(c['g'], j, p) % p for j, v in enumerate(F))
    centered = all(abs(v) <= p//2 for v in F)
    assert multiply(c['f'], F) == [p*v for v in digits]
    row = {'a': a, 'F': F, 'scalar_congruences': scalar, 'centered': centered}
    row['status'] = ('scalar_mismatch' if not scalar else 'outside_centered' if not centered
                     else 'zero' if a == 0 else 'nonzero')
    return row


def artifact_row(c, digits):
    return {'digits': list(digits), **classify(c, digits)}


def main():
    source = HERE/'norm_compression.json'
    compressed = json.loads(source.read_text())
    cases = []
    for case_index, case in enumerate(compressed['cases']):
        c = prepare(case['p'], case['g'], case['relation'])
        assert (c['k'], c['norm_f'], c['bound']) == tuple(case[k] for k in
                ('relation_cofactor', 'relation_norm', 'digit_height_bound'))
        positive = []
        for r in case['records']:
            out = artifact_row(c, r['digits'])
            assert out['status'] == 'nonzero' and out['a'] == r['a']
            positive.append(out)
        # Deterministic perturbations plus bounded random cube samples.
        probes = [[0]*c['N'], [c['bound']+1]+[0]*(c['N']-1)]
        for r in case['records'][:8]:
            for j in (0, c['N']//2, c['N']-1):
                d = list(r['digits']); d[j] = -d[j] if d[j] else 1
                probes.append(d)
        rng = random.Random(1729+case_index)
        probes += [[rng.randint(-c['bound'], c['bound']) for _ in range(c['N'])]
                   for _ in range(64)]
        negative = [artifact_row(c, d) for d in probes]
        cases.append({'name': case['name'], 'certificate': c, 'positives': positive,
                      'probes': negative,
                      'probe_status_counts': dict(sorted(Counter(r['status'] for r in negative).items()))})
        print(json.dumps({'case': case['name'], 'positive_vectors': len(positive),
                          'probe_status_counts': cases[-1]['probe_status_counts']}), flush=True)

    c = cases[0]['certificate']
    assert c['N'] == 8 and c['p'] == 257 and c['k'] == 1 and c['bound'] == 1
    census = []; by_norm = defaultdict(Counter); by_pair = defaultdict(Counter)
    witnesses = defaultdict(dict)
    for digits in product((-1, 0, 1), repeat=c['N']):
        row = artifact_row(c, digits)
        row['digit_norm'] = norm(digits)
        row['digit_energy'] = sum(v*v for v in digits)
        good = row['status'] in ('zero', 'nonzero')
        label = 'realized' if good else 'nonrealized'
        by_norm[row['digit_norm']][label] += 1
        key = (row['digit_norm'], row['digit_energy'])
        by_pair[key][label] += 1
        witnesses[key].setdefault(label, row)
        census.append(row)
    collisions = [v for k, v in sorted(witnesses.items()) if len(v) == 2]
    assert collisions
    realized_norms = {r['digit_norm'] for r in census if r['status'] in ('zero', 'nonzero')}
    realized_pairs = {(r['digit_norm'], r['digit_energy']) for r in census
                      if r['status'] in ('zero', 'nonzero')}
    realized = [r for r in census if r['status'] in ('zero', 'nonzero')]
    assert sorted(r['a'] for r in realized) == list(range(c['p']))
    stats = {'cube_words': len(census), 'actual_words_including_zero': len(realized),
             'status_counts': dict(sorted(Counter(r['status'] for r in census).items())),
             'actual_norm_values': sorted(realized_norms),
             'norm_divisible_by_k_accepts': sum(r['digit_norm'] % c['k'] == 0 for r in census),
             'oracle_actual_norm_set_accepts': sum(r['digit_norm'] in realized_norms for r in census),
             'oracle_actual_norm_energy_pair_set_accepts': sum(
                 (r['digit_norm'], r['digit_energy']) in realized_pairs for r in census),
             'by_norm': [{'norm': k, **dict(v)} for k, v in sorted(by_norm.items())],
             'by_norm_energy': [{'norm': k[0], 'energy': k[1], **dict(v)}
                                for k, v in sorted(by_pair.items())],
             'same_norm_and_energy_witness': collisions[0]}
    # f=p is valid but has all split roots zero; integrality alone does not
    # recover the required one-dimensional scalar section when p divides k.
    scaled = prepare(c['p'], c['g'], [c['p']]+[0]*(c['N']-1))
    generality = [artifact_row(scaled, [1]+[0]*(c['N']-1)),
                  artifact_row(scaled, [0]*c['N']),
                  artifact_row(scaled, cases[0]['positives'][0]['F'])]
    assert [r['status'] for r in generality] == ['scalar_mismatch', 'zero', 'nonzero']
    output = {'status': 'produced', 'scope': 'Exact inverse membership and full small digit-cube census. Norm-set relaxations use an oracle learned from the full small actual section, not a proposed large-scale algorithm.',
              'cases': cases, 'small_census': census, 'small_summary': stats,
              'p_divides_cofactor_control': {'certificate': scaled, 'records': generality},
              'source_sha256': {name: sha256((HERE/name).read_bytes()).hexdigest()
                                for name in ['realizability.py', 'norm_carry.py']},
              'input_sha256': {source.name: sha256(source.read_bytes()).hexdigest()}}
    (HERE/'realizability.json').write_text(json.dumps(output, separators=(',', ':'))+'\n')
    print(json.dumps({k: v for k, v in stats.items() if k not in ('by_norm', 'by_norm_energy')}), flush=True)


if __name__ == '__main__':
    main()
