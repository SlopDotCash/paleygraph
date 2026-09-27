#!/usr/bin/env python3
"""Verify an exact conditioned certificate, then list completions or report budget."""
import argparse
import json
from pathlib import Path
from conditioned_erasure import decode, encode
from conditioned_review import certify
from metric_review import check_transform


def query(c, digits, budget=50000, start=None):
    certify(c)
    if 'metric_change' in c: check_transform(c)
    N, U = c['N'], c['erased']; d = len(U)
    if not isinstance(digits, list) or len(digits) != N: raise ValueError('input must be an array of N digits')
    if any(v is not None and type(v) is not int for v in digits): raise ValueError('digits must be integers or null')
    if start is None:
        if [i for i,v in enumerate(digits) if v is None] != U: raise ValueError('input erasures do not match the certificate')
        known = {i:v for i,v in enumerate(digits) if v is not None}
    else:
        if U != list(range(d)) or type(start) is not int or start not in range(N): raise ValueError('invalid cyclic block')
        erased = {(start+j) % N for j in range(d)}
        if {i for i,v in enumerate(digits) if v is None} != erased: raise ValueError('input erasures do not match cyclic block')
        known = {j:digits[(j+start) % N]*(1 if j+start < N else -1) for j in range(d,N)}
    out = decode(c, known, budget)
    if start is not None and out['status'] == 'complete':
        completions = []
        for v in out['completions']:
            a = v['scalar']*pow(c['g'], -start, c['p']) % c['p']
            word = encode(c['p'], c['g'], c['relation'], a)[0]
            assert all(x is None or word[i] == x for i,x in enumerate(digits))
            completions.append({'scalar':a, 'digits':word})
        out['completions'] = sorted(completions, key=lambda x:x['scalar'])
    return {'certificate_verified':True, 'p':c['p'], 'N':N,
            'universal_candidate_box_cap':c['universal_candidate_box_cap'],
            'universal_unique_completion':c['universal_unique_completion'], **out}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--certificate', type=Path, required=True)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--candidate-budget', type=int, default=50000)
    parser.add_argument('--cyclic-start', type=int)
    args = parser.parse_args()
    try:
        data = json.loads(args.certificate.read_text()); c = data.get('certificate',data)
        output = query(c, json.loads(args.input.read_text()), args.candidate_budget, args.cyclic_start)
    except (ValueError, AssertionError, KeyError, TypeError, IndexError) as error:
        parser.exit(2, 'invalid input or certificate: '+str(error)+'\n')
    print(json.dumps(output, indent=2))


if __name__ == '__main__': main()
