#!/usr/bin/env python3
"""Finite, exact distinguishing-prefix witnesses from projection/re-encoding.

A greedy clique is a certified lower bound, not a maximum clique. Distinct
sampled prefixes alone are never counted as evidence of distinct states.
"""
from hashlib import sha256
import json
from pathlib import Path
import random
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'round19'))
from projection_codec import prepare, encode, decode


def clique(matrix):
    M = len(matrix)
    neighbors = [{j for j in range(M) if i != j and not (matrix[i][j] and matrix[j][i])} for i in range(M)]
    available = set(range(M)); selected = []
    while available:
        i = max(available, key=lambda j: (len(neighbors[j] & available), -j))
        selected.append(i); available &= neighbors[i]
    return selected


def main():
    source = HERE.parent/'round17/norm_compression.json'; saved = json.loads(source.read_text())['cases']
    small, mid, large = saved[1], saved[2], saved[-1]
    specs = []
    for label, c, g, f, order, cut in [
            ('p65537_canonical', small, 2, [2]+[0]*14+[1], list(range(16)), 8),
            ('p65537_saved', small, small['g'], small['relation'], list(range(16)), 8),
            ('p6700417_canonical', mid, 2, mid['relation'], list(range(32)), 16),
            *[(f'p2013265921_{label}_cut{cut}', large, large['g'], large['relation'], order, cut)
              for label, order, cut in [
                  ('natural', list(range(64)), 4), ('natural', list(range(64)), 8),
                  ('natural', list(range(64)), 32), ('reversed', list(range(63, -1, -1)), 32),
                  ('even_odd', list(range(0, 64, 2))+list(range(1, 64, 2)), 32)]]]:
        specs.append((label, c['p'], g, f, order, cut))
    reports = []
    for label, p, g, f, order, cut in specs:
        codec = prepare(p, g, f); N = len(f); seed = 20260906; rng = random.Random(seed)
        scalars = list(dict.fromkeys([0, 1, p//2, (p+1)//2, p//3, p//5]))
        while len(scalars) < 256:
            a = rng.randrange(p)
            if a not in scalars: scalars.append(a)
        words = [encode(p, g, f, a)[0] for a in scalars]; matrix = []
        for i, left in enumerate(words):
            row = []
            for right in words:
                mixed = right.copy()
                for j in order[:cut]: mixed[j] = left[j]
                row.append(decode(codec, mixed)['status'].startswith('accepted_'))
            assert row[i]; matrix.append(row)
        selected = clique(matrix); bits = []
        for n, i in enumerate(selected):
            for j in selected[n+1:]:
                assert not matrix[i][j] or not matrix[j][i]
                bits.append('0' if not matrix[i][j] else '1')
        r = {'name': label, 'p': p, 'g': g, 'relation': f, 'N': N, 'order': order, 'cut': cut,
             'seed': seed, 'sampling': 'six fixed scalar anchors followed by distinct seeded PRNG draws; no uniformity claim',
             'scalars': scalars, 'words_in_original_coordinates': words,
             'cross_acceptance_rows': [''.join('1' if v else '0' for v in row) for row in matrix],
             'clique_indices': selected, 'certified_width_lower_bound': len(selected),
             'pair_witness_selectors': ''.join(bits),
             'selector_rule': 'Pairs follow clique index order, upper triangle. 0: suffix of j accepts prefix j but rejects prefix i. 1: suffix of i accepts prefix i but rejects prefix j.',
             'distinct_sampled_prefixes': len({tuple(w[j] for j in order[:cut]) for w in words}),
             'cross_splices_checked': len(words)**2, 'accepted_cross_splices': sum(map(sum, matrix)),
             'certified_distinguished_pairs': len(bits)}
        path = HERE/(label+'_splice.json'); path.write_text(json.dumps(r, separators=(',', ':'))+'\n')
        row = {k: r[k] for k in ['name', 'N', 'cut', 'distinct_sampled_prefixes', 'certified_width_lower_bound', 'cross_splices_checked', 'accepted_cross_splices', 'certified_distinguished_pairs']}
        row['artifact'] = path.name; reports.append(row); print(json.dumps(row), flush=True)
    out = {'status': 'produced', 'scope': 'Exact membership matrix and explicit suffix selectors certify pairwise distinguishable prefixes at the stated cut and digit order. Finite lower bounds only; no all-orders or asymptotic lower bound.',
           'cases': reports,
           'source_sha256': {n: sha256((HERE/n).read_bytes()).hexdigest() for n in ['splice_certificates.py', '../round19/projection_codec.py']},
           'input_sha256': {'../round17/norm_compression.json': sha256(source.read_bytes()).hexdigest()}}
    (HERE/'splice_summary.json').write_text(json.dumps(out, indent=2)+'\n')


if __name__ == '__main__': main()
