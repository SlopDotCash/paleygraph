#!/usr/bin/env python3
"""Bracket the L1 cost of separating the near-boundary witness."""
from hashlib import sha256
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'round17'))
from carry_review import product


def main():
    source = HERE/'relation_barrier.json'; reviewed = HERE/'barrier_review.json'
    data = json.loads(source.read_text()); review = json.loads(reviewed.read_text())
    assert review['status'] == 'passed' and review['input_sha256'][source.name] == sha256(source.read_bytes()).hexdigest()
    cases = []
    for c in data['cases']:
        p, N, g = c['p'], c['N'], c['g']; m = p//2
        powers = [(pow(g, j, p), j) for j in range(2*N)]
        assert len({v for v, j in powers}) == 2*N
        odd, exponent = min((v, j) for v, j in powers if v != 1 and v % 2)
        assert odd == 2*c['other_coordinate_slack']+1
        # f=odd-X^(-exponent), expressed in the canonical negacyclic basis.
        f = [odd]+[0]*(N-1); at = (-exponent) % (2*N)
        f[at % N] += -1 if at < N else 1
        assert sum(v*pow(pow(g, -1, p), j, p) for j, v in enumerate(f)) % p == 0
        l1 = sum(map(abs, f)); bound = m*l1//p
        good = product(f, c['centered_F'], N); bad = product(f, c['escaped_F'], N)
        assert all(v % p == 0 for v in good+bad)
        good = [v//p for v in good]; bad = [v//p for v in bad]
        assert max(map(abs, good)) <= bound
        assert bad[0] == bound+1 == (odd+1)//2
        assert l1 == odd+1 == 2*c['universal_relation_l1_budget']
        cases.append({'p': p, 'N': N, 'g': g, 'least_nonidentity_odd_subgroup_residue': odd,
                      'exponent': exponent, 'separating_relation': f, 'separating_relation_l1': l1,
                      'digit_bound': bound, 'centered_encoding': good, 'escaped_encoding': bad,
                      'minimum_separating_l1_lower_bound': c['universal_relation_l1_budget']+1,
                      'minimum_separating_l1_upper_bound': l1})
        print(json.dumps({key: value for key, value in cases[-1].items()
                          if key not in ('separating_relation', 'centered_encoding', 'escaped_encoding')}), flush=True)
    out = {'status': 'passed', 'scope': 'The reviewed universal lower budget and an explicit binomial separator bracket the least coefficient L1 needed by this relation-bound interface. This is a factor-two bracket, not an exact optimum or a shortest-vector claim.',
           'cases': cases,
           'source_sha256': {name: sha256((HERE/name).read_bytes()).hexdigest()
                             for name in ['separation_bounds.py', '../round17/carry_review.py']},
           'input_sha256': {p.name: sha256(p.read_bytes()).hexdigest() for p in [source, reviewed]}}
    (HERE/'separation_bounds.json').write_text(json.dumps(out, indent=2)+'\n')


if __name__ == '__main__': main()
