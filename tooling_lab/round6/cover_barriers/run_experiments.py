#!/usr/bin/env python3
"""Full support-mask census, alternate interpolation replay, and cap-family audit."""
from collections import defaultdict
from hashlib import sha256
from itertools import product
from math import comb
from pathlib import Path
import json
import random
import subprocess
import sys
import time
sys.dont_write_bytecode = True
from cover_barriers import *
from compressed_support import domain


def direct_book(F, dom, k):
    return [tuple(F.sum(F.mul(c, F.power(x, j)) for j, c in enumerate(cs)) for x in dom)
            for cs in product(range(F.q), repeat=k)]


def direct_maxima(book, words):
    supports = [[sum(1 << i for i, (x, y) in enumerate(zip(cw, w)) if x == y) for cw in book] for w in words]
    maxima = [max(x.bit_count() for x in masks) for masks in supports]
    joint = max((a & b).bit_count() for a in supports[0] for b in supports[1])
    return maxima + [joint]


def main():
    start = time.perf_counter()
    output = {'schema': 'cap-family-barriers/v1', 'tiny': [], 'hard_stacks': [], 'large_structured': [],
              'rate_quarter': [], 'scope': 'obstruction to the entire verified disjoint-piece root-count cap family, not to all decoding certificates or to the prize theorem'}
    partition_checks = 0
    for n in range(1, 41):
        for k in range(1, n + 1):
            for M in range(1, n + 1):
                dp = [0] + [10 ** 9] * n
                for size in range(1, n + 1):
                    dp[size] = min(dp[size - part] + min(part, k - 1) for part in range(1, min(M, size) + 1))
                result = minimum_cap(n, k, M)
                assert result['minimum'] == dp[n]
                assert sum(result['attaining_integer_partition']) == n
                partition_checks += 1
    output['partition_DP_checks'] = partition_checks
    print('partition DP', partition_checks, flush=True)

    cases = 0
    cpp_tiny_inputs, cpp_tiny_expected = [], []
    for F, n, k, count in [(PrimeField(3), 3, 1, 729), (PrimeField(5), 5, 2, 32),
                           (Field(3, [1, 0, 1]), 8, 2, 32)]:
        dom = list(range(n)); book = direct_book(F, dom, k); rng = random.Random(81531 + F.q)
        inputs = product(range(F.q), repeat=2 * n) if F.q == 3 else (
            [rng.randrange(F.q) for _ in range(2 * n)] for _ in range(count))
        max_hist = defaultdict(int)
        for values in inputs:
            words = [list(values[:n]), list(values[n:])]
            result = enumerate_bases(F, dom, k, words)
            expected = direct_maxima(book, words)
            assert result['maxima'] == expected
            if F.q == 5 and len(cpp_tiny_inputs) < 8:
                cpp_tiny_inputs.append(words)
                cpp_tiny_expected.append(expected)
            max_hist[tuple(expected)] += 1
            cases += 1
        output['tiny'].append({'field': field_record(F), 'n': n, 'k': k, 'cases': count,
                               'codewords_per_scalar_book': len(book), 'polynomial_pairs_compared_per_case': len(book) ** 2,
                               'maximum_histogram': [{'maxima': list(key), 'count': v} for key, v in sorted(max_hist.items())]})
    output['tiny_full_polynomial_pair_oracles'] = cases
    print('tiny maxima', cases, flush=True)
    tiny_input = HERE / 'tiny_prime.input.txt'
    tiny_ledger = HERE / 'tiny_prime.supports.bin'
    text = ['5 5 2 8', '0 1 2 3 4']
    for words in cpp_tiny_inputs:
        text.extend(' '.join(map(str, word)) for word in words)
    tiny_input.write_text('\n'.join(text) + '\n')
    for mode in ['generate', 'verify']:
        destination = HERE / f'tiny_prime.{mode}.json'
        subprocess.run([str(HERE / 'base_census'), mode, str(tiny_input), str(tiny_ledger), str(destination)], check=True)
        data = json.loads(destination.read_text())
        assert data['basis_count'] == 10
        assert [row['maxima_u0_u1_joint'] for row in data['samples']] == cpp_tiny_expected
    output['cpp_complete_tiny_polynomial_pair_comparisons'] = 8
    corrupt = HERE / 'tiny_corrupted.supports.bin'
    payload = bytearray(tiny_ledger.read_bytes()); payload[0] ^= 1; corrupt.write_bytes(payload)
    failed = subprocess.run([str(HERE / 'base_census'), 'verify', str(tiny_input), str(corrupt), str(HERE / 'invalid.json')], capture_output=True, text=True)
    assert failed.returncode != 0 and 'Newton support mismatch' in failed.stderr
    corrupt.unlink()
    output['corrupted_mask_ledger_rejected'] = True

    records = json.loads((LAB / 'proximity/stack_results.json').read_text())['records']
    grouped = defaultdict(list)
    for row in records:
        grouped[row['configuration']].append(row)
    total_bases = 0
    for name, rows in grouped.items():
        begin = time.perf_counter()
        p, n, k = rows[0]['p'], rows[0]['n'], rows[0]['k']; F = PrimeField(p); dom = domain(n, p)
        assert comb(n, k) <= 1000000
        source = HERE / f'{name}.input.txt'; ledger = HERE / f'{name}.supports.bin'
        values = [f'{p} {n} {k} {len(rows)}', ' '.join(map(str, dom))]
        for row in rows:
            values += [' '.join(map(str, row['u0'])), ' '.join(map(str, row['u1']))]
        source.write_text('\n'.join(values) + '\n')
        generated = HERE / f'{name}.enumerated.json'; checked = HERE / f'{name}.verified.json'
        if '--reuse-ledgers' not in sys.argv:
            subprocess.run([str(HERE / 'base_census'), 'generate', str(source), str(ledger), str(generated)], check=True)
        subprocess.run([str(HERE / 'base_census'), 'verify', str(source), str(ledger), str(checked)], check=True)
        a = json.loads(generated.read_text()); b = json.loads(checked.read_text())
        assert a['basis_count'] == b['basis_count'] == comb(n, k) and a['samples'] == b['samples']
        assert ledger.stat().st_size == comb(n, k) * len(rows) * 3 * 4
        for row, sample in zip(rows, a['samples']):
            witnesses = []
            for j in range(3):
                base = sample['witness_bases'][j]
                words = [row['u0'], row['u1']] if j == 2 else [row['u0' if j == 0 else 'u1']]
                cs = [interpolate(F, [dom[i] for i in base], [w[i] for i in base]) for w in words]
                supp = [i for i, x in enumerate(dom) if all(evaluate(F, h, x) == w[i] for h, w in zip(cs, words))]
                assert len(supp) == sample['maxima_u0_u1_joint'][j]
                assert sum(1 << i for i in supp) == sample['witness_support_masks'][j]
                witnesses.append({'mode': ['u0', 'u1', 'joint'][j], 'coefficients': cs, 'support': supp, 'base': base})
            barriers = [cap_barrier(n, k, row['s'], M, 'exact_maximum') for M in sample['maxima_u0_u1_joint']]
            attainments = [attain_from_maximum(F, dom, k, [row['u0'], row['u1']] if j == 2 else
                          [row['u0' if j == 0 else 'u1']], witness) for j, witness in enumerate(witnesses)]
            assert all(x is not None for x in attainments)
            output['hard_stacks'].append({'configuration': name, 'sample': row['sample'], 'field': field_record(F),
                'domain': dom, 'n': n, 'k': k, 's': row['s'], 'u0': row['u0'], 'u1': row['u1'],
                'exact_maxima_u0_u1_joint': sample['maxima_u0_u1_joint'], 'maximum_witnesses': witnesses,
                'actual_cap_minimum_attainment_u0_u1_joint': attainments,
                'barriers_u0_u1_joint': barriers, 'basis_count': comb(n, k),
                'ledger': {'path': ledger.name, 'sha256': sha256(ledger.read_bytes()).hexdigest(),
                           'sample_index': row['sample'], 'samples_per_basis': len(rows),
                           'encoding': 'lexicographic k-bases; for each sample u0,u1,joint full support masks, little-endian uint32'},
                'upper_bound_certificate': 'complete k-base enumeration, full support-mask ledger, independently re-evaluated by Newton interpolation',
                'candidate_only': False})
        total_bases += comb(n, k) * len(rows)
        print('hard', name, [s['maxima_u0_u1_joint'] for s in a['samples']], 'seconds', time.perf_counter() - begin, flush=True)
    output['hard_stack_bases_exhausted'] = total_bases

    frozen = json.loads((LAB / 'round6/pencil_tracks/results.json').read_text())
    for item in frozen['large']:
        cert = item['result']['certificate']; F = PrimeField(cert['field']['p'])
        bound = polynomial_partition_bound(F, cert['domain'], cert['k'], cert['u0'], cert['u1'], cert['tracks'])
        assert bound['proof_kind'] == 'exact_maximum'
        quarter_bound = polynomial_partition_bound(F, cert['domain'], 256, cert['u0'], cert['u1'], cert['tracks'])
        output['large_structured'].append({'name': item['name'], 'n': cert['n'], 'k': cert['k'], 's': cert['s'],
                'bound_certificate': bound, 'barrier': cap_barrier(cert['n'], cert['k'], cert['s'], bound['M'], bound['proof_kind']),
                'source': 'round6/pencil_tracks/results.json', 'field': cert['field'], 'domain': cert['domain'],
                'u0': cert['u0'], 'u1': cert['u1'], 'rate': '1/16, not the official 1/4 rate',
                'same_input_at_k256': {'bound_certificate': quarter_bound,
                    'barrier_at_s566': cap_barrier(1024, 256, 566, quarter_bound['M'], quarter_bound['proof_kind']),
                    'scope': 'simplified rate-quarter degree class, prime field; no official event certificate'}})

    for n in [16, 32, 64, 128, 256, 1024, 2 ** 30]:
        k = n // 4
        for s in sorted({n // 2, n // 2 + 1, (53 * n + 95) // 96}):
            lo, hi = 1, n
            while lo < hi:
                mid = (lo + hi) // 2
                if minimum_cap(n, k, mid)['minimum'] < s:
                    hi = mid
                else:
                    lo = mid + 1
            threshold = lo
            delta = s - n // 2
            if k >= 3 and 0 <= delta <= k - 3:
                assert threshold == n // 2 - (delta + 1) // 2
            output['rate_quarter'].append({'n': n, 'k': k, 's': s,
                    'minimum_M_not_blocked_by_abstract_cap': threshold,
                    'all_M_below_threshold_blocked': True,
                    'cap_at_threshold': minimum_cap(n, k, threshold)['minimum'],
                    'cap_one_below': minimum_cap(n, k, threshold - 1)['minimum'] if threshold > 1 else None,
                    'meaning': 'Necessary maximum joint agreement for this cap family; no claim an actual cover exists.'})
    # Changing k changes the polynomial class. The two-track M=512 bound was
    # explicitly re-proved at k256; the three-track M=342 comparison is hypothetical there.
    output['rate_comparison'] = [{'n': 1024, 'k': k, 's': s, 'hypothetical_M': M,
                                 'minimum_cap': minimum_cap(1024, k, M)['minimum'],
                                 'M_status': 'certified_exact' if k == 64 or M == 512 else 'hypothetical_only'}
                                for M in [342, 512] for k, s in [(64, 410), (256, 512), (256, 566)]]
    try:
        cap_barrier(16, 8, 11, 8, 'candidate_maximum')
    except ValueError:
        output['candidate_lower_bound_rejected_as_obstruction_evidence'] = True
    else:
        raise AssertionError('candidate-only maximum accepted')
    output['source_sha256'] = {name: sha256((HERE / name).read_bytes()).hexdigest()
                               for name in ['base_census.cpp', 'base_census', 'cover_barriers.py', 'run_experiments.py']}
    dependencies = ['proximity/stack_results.json', 'round6/pencil_tracks/results.json',
                    'round6/pencil_tracks/pencil_tracks.py', 'round5/piece_discovery/piece_discovery.py',
                    'round4/scalar_fibers/scalar_fibers.py', 'round2/proximity/compressed_support.py',
                    'proximity/extension_field_probe.py', 'proximity/deformation_microscope.py']
    output['dependency_sha256'] = {name: sha256((LAB / name).read_bytes()).hexdigest() for name in dependencies}
    output['artifact_sha256'] = {path.name: sha256(path.read_bytes()).hexdigest() for path in HERE.iterdir()
                                if path.name.endswith(('.input.txt', '.supports.bin', '.enumerated.json', '.verified.json'))
                                or path.name.startswith('tiny_prime.') and path.suffix == '.json'}
    output['elapsed_seconds'] = time.perf_counter() - start
    (HERE / 'results.json').write_text(json.dumps(output, indent=2) + '\n')
    print('complete', output['elapsed_seconds'], flush=True)


if __name__ == '__main__':
    main()
