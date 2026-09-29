#!/usr/bin/env python3
"""Small direct-product and permutation controls for the compiled six-set census."""
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json
import random
import subprocess
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'proximity'))
from extension_field_probe import Field


def kernel(S, support):
    out = 0
    for row in S:
        product = 1
        for x in support:
            product *= row[x]
        out += product
    return out


def census(S, T):
    q = len(S)
    masks = [[sum(1 << i for i, v in enumerate(row) if v < 0)
              for row in matrix] for matrix in [S, T]]
    payload = str(q) + '\n' + '\n'.join(' '.join(map(str, mask)) for mask in masks) + '\n'
    raw = subprocess.run([str(HERE / 'census_six')], input=payload,
                         text=True, capture_output=True, check=True).stdout
    hists = {'P': Counter(), 'Q': Counter(), 'J': Counter()}
    for line in raw.splitlines():
        kind, *rest = line.split()
        values = list(map(int, rest))
        if kind in ['P', 'Q']:
            hists[kind][values[0]] = values[1]
        elif kind == 'J':
            hists['J'][tuple(values[:2])] = values[2]
    direct = Counter((kernel(S, c), kernel(T, c)) for c in combinations(range(q), 6))
    assert direct == hists['J']
    assert Counter({x: sum(n for (a, _), n in direct.items() if a == x)
                    for x in hists['P']}) == hists['P']
    assert Counter({x: sum(n for (_, b), n in direct.items() if b == x)
                    for x in hists['Q']}) == hists['Q']
    return hists


def main():
    subprocess.run(['clang++', '-O3', '-std=c++17', str(HERE / 'census_six.cpp'),
                    '-o', str(HERE / 'census_six')], check=True)
    F = Field(3, [1, 0, 1])
    g = next(x for x in range(1, 9) if F.power(x, 4) != 1)
    classes = [{F.power(g, i) for i in range(8) if i % 4 == j} for j in range(4)]
    matrices = [[[0 if a == b else (1 if F.sub(a, b) in conn else -1)
                  for b in range(9)] for a in range(9)]
                for conn in [classes[0] | classes[2], classes[0] | classes[1]]]
    tiny = census(*matrices)
    assert tiny['P'] == tiny['Q']  # Order-nine Paley/Peisert isomorphic control.
    permutation = list(range(9))
    random.Random(84019).shuffle(permutation)
    permuted = [[matrices[0][permutation[i]][permutation[j]] for j in range(9)]
                for i in range(9)]
    null = census(matrices[0], permuted)
    assert null['P'] == null['Q']
    assert any(a != b for a, b in null['J'])  # Pointwise difference alone is insufficient.
    original = json.loads((HERE / 'results.json').read_text())
    large = original['sign_matrices']
    supports = list(combinations(range(9), 6))
    rng = random.Random(493061)
    supports += [sorted(rng.sample(range(49), 6)) for _ in range(1000)]
    large_checks = 0
    for S in large:
        masks = [sum(1 << i for i, v in enumerate(row) if v < 0) for row in S]
        for c in supports:
            parity = 0
            for x in c:
                parity ^= masks[x]
            live = ((1 << 49) - 1) ^ sum(1 << x for x in c)
            assert kernel(S, c) == 43 - 2 * (parity & live).bit_count()
            large_checks += 1
    witness = original['census']['largest_pointwise_difference']
    assert [kernel(S, witness['C']) for S in large] == [witness['paley_value'], witness['peisert_value']]
    output = {'status': 'passed', 'gf9_direct_complete_joint_censuses': 2,
              'sets_per_gf9_census': 84, 'gf9_histogram': dict(tiny['P']),
              'permutation_null_same_histogram': True,
              'permutation_null_has_pointwise_differences': True,
              'gf49_direct_product_vs_bit_parity_checks': large_checks,
              'largest_difference_witness_directly_rechecked': witness,
              'source_sha256': {name: sha256((HERE / name).read_bytes()).hexdigest()
                                for name in ['controls.py', 'census_six.cpp', 'results.json']}}
    (HERE / 'control_results.json').write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
