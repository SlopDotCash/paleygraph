#!/usr/bin/env python3
"""Reuse frozen full-support ledgers, optimize two small covers, export results."""
import hashlib
import json
import struct
import sys
import time
from collections import Counter
from itertools import combinations
from math import comb, ceil
from pathlib import Path
sys.dont_write_bytecode = True
from overlap_cover import greedy_cover, track_from_base, compile_certificate, size_volume

HERE = Path(__file__).resolve().parent
LAB = HERE.parents[1]
BARRIERS = LAB / 'round6/cover_barriers'


def binding(path):
    return {'path': str(path.relative_to(LAB)), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def baseline_keys(nodes):
    return {(r['scalar'], tuple(r['codeword']), r['agreement_mask']) for r in nodes}


def main():
    started = time.perf_counter()
    baseline = LAB / 'proximity/stack_results.json'
    records = json.loads(baseline.read_text())['records']
    results = []
    for cfg, sample in [('coset_medium_characteristic', 0), ('coset_small_characteristic', 3)]:
        begin = time.perf_counter()
        record = next(r for r in records if r['configuration'] == cfg and r['sample'] == sample)
        p, n, k, s = [record[key] for key in ['p', 'n', 'k', 's']]
        input_path, ledger_path = [BARRIERS / (cfg + ext) for ext in ['.input.txt', '.supports.bin']]
        data = list(map(int, input_path.read_text().split()))
        assert data[:4] == [p, n, k, 4]
        domain = data[4:4 + n]
        assert data[4+n+sample*2*n:4+n+(sample+1)*2*n] == record['u0'] + record['u1']
        # Decode the independently barycentric-generated/Newton-replayed catalog.
        raw = ledger_path.read_bytes()
        assert len(raw) == comb(n, k) * 4 * 3 * 4
        words = struct.unpack('<' + str(len(raw) // 4) + 'I', raw)
        supports = words[3 * sample + 2::12]
        catalog = {}
        occurrences = Counter(supports)
        for base, mask in zip(combinations(range(n), k), supports):
            assert mask >> n == 0 and all(mask >> i & 1 for i in base)
            catalog.setdefault(mask, base)
        # Distinct polynomial pairs have joint-support intersection <k.
        # Hence each maximal support has exactly C(|S|,k) generating bases.
        assert all(occurrences[m] == comb(m.bit_count(), k) for m in catalog)
        chosen, stats = greedy_cover(n, k, s, catalog)
        tracks = [track_from_base(p, domain, k, record['u0'], record['u1'], catalog[m]) for m in chosen]
        assert [sum(1 << i for i in t['joint_support']) for t in tracks] == chosen
        cert = compile_certificate(p, domain, k, s, record['u0'], record['u1'], tracks)
        actual = {(r['scalar'], tuple(r['codeword']), sum(1 << i for i in r['agreement_support'])) for r in cert['nodes']}
        assert record['node_ledger_complete'] and actual == baseline_keys(record['nodes'])
        bindings = [binding(path) for path in [baseline, input_path, ledger_path, BARRIERS / 'base_census.cpp']]
        cert['input_provenance'] = bindings
        stem = cfg + '_sample' + str(sample)
        (HERE / (stem + '.certificate.json')).write_text(json.dumps(cert, indent=2) + '\n')
        cat_export = {'schema': 'complete_joint_support_catalog_v1', 'configuration': cfg,
                      'sample': sample, 'p': p, 'n': n, 'k': k, 'input_provenance': bindings,
                      'implicit_track_definition': 'interpolate u0 and u1 on the listed first basis',
                      'candidates': [{'support_mask': m, 'first_basis': list(catalog[m]),
                                      'generating_bases': occurrences[m]} for m in sorted(catalog)]}
        (HERE / (stem + '.catalog.json')).write_text(json.dumps(cat_export, separators=(',', ':')) + '\n')
        M = max(m.bit_count() for m in catalog)
        result = {'configuration': cfg, 'sample': sample, 'p': p, 'n': n, 'k': k, 's': s,
                  **stats, 'all_k_bases': comb(n, k), 'old_anchored_bases': comb(n - s + k, k),
                  'maximal_joint_support_size': M, 'one_largest_support_covers': size_volume(n, k, s, M),
                  'corrected_volume_lower_bound': ceil(comb(n, s) / size_volume(n, k, s, M)),
                  'k_block_counting_lower_bound': ceil(comb(n, k) / comb(s, k)),
                  'k_block_bound_applies_to_this_catalog': M == k,
                  'complete_output_equal_to_frozen_s_subset_oracle': True,
                  'node_count': cert['node_count'], 'qualifying_scalars': sorted(set(r['scalar'] for r in cert['nodes'])),
                  'whole_field_track_count': len(cert['whole_field_track_ids']),
                  'certificate': stem + '.certificate.json', 'catalog': stem + '.catalog.json',
                  'elapsed_seconds_including_export': time.perf_counter() - begin}
        print(json.dumps(result), flush=True)
        results.append(result)
    summary = {'schema': 'overlap_preflight_v1', 'results': results,
               'source_bindings': [binding(HERE / name) for name in ['overlap_cover.py', 'run_overlap.py']],
               'scope': 'Known cover mechanism; new local actual-support selection, finite exhaustive certificates, no optimum or scaling claim.',
               'elapsed_seconds': time.perf_counter() - started}
    (HERE / 'results.json').write_text(json.dumps(summary, indent=2) + '\n')


if __name__ == '__main__':
    main()
