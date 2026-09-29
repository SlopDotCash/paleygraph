#!/usr/bin/env python3
"""Source-bound readback of every saved full discovery record.

This verifier shares algebra helpers with the generator. The independent
review in round5/review is a separate source of evidence.
"""
from hashlib import sha256
from pathlib import Path
import json
import sys
sys.dont_write_bytecode = True
from piece_discovery import *


def main():
    data = json.loads((HERE / 'results.json').read_text())
    for name, expected in data['source_sha256'].items():
        assert sha256((HERE / name).read_bytes()).hexdigest() == expected, name
    counts = {'records': 0, 'linear_systems': 0, 'complete_static': 0, 'complete_scalar': 0,
              'minimal_piece_counts': 0, 'dual_inconsistency_witnesses': 0}

    def walk(obj):
        if isinstance(obj, dict):
            if obj.get('schema') == 'piece-discovery/v1':
                f = obj['field']
                F = PrimeField(f['p']) if f['e'] == 1 else Field(f['p'], f['modulus'])
                verify_export(F, obj)
                counts['records'] += 1
                counts['complete_static'] += obj['complete_static_list_certified']
                counts['complete_scalar'] += obj['complete_scalar_fibers_certified']
                counts['minimal_piece_counts'] += obj['minimum_polynomial_cover_size_certified'] is not None
                for a in obj['attempts']:
                    if 'linear_certificate' in a:
                        counts['linear_systems'] += 1
                        counts['dual_inconsistency_witnesses'] += a['linear_certificate']['status'] == 'inconsistent'
                return
            for v in obj.values():
                walk(v)
        elif isinstance(obj, list):
            for v in obj:
                walk(v)
    walk(data)
    output = {'status': 'passed', 'scope': 'source-bound readback sharing production algebra helpers',
              'counts': counts,
              'source_sha256': {name: sha256((HERE / name).read_bytes()).hexdigest()
                                for name in ['piece_discovery.py', 'run_experiments.py', 'verify_results.py']},
              'results_sha256': sha256((HERE / 'results.json').read_bytes()).hexdigest()}
    (HERE / 'verification.json').write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
