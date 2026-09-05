#!/usr/bin/env python3
"""Source-bound readback of saved discovery and symbolic track certificates."""
from hashlib import sha256
import json
from pathlib import Path
import sys
sys.dont_write_bytecode = True
from pencil_tracks import *


def main():
    data = json.loads((HERE / 'results.json').read_text())
    for name, expected in data['source_sha256'].items():
        assert sha256((HERE / name).read_bytes()).hexdigest() == expected
    counts = {'discoveries': 0, 'marginal_algebra_ledgers': 0, 'complete_discoveries': 0,
              'no_anchor_complete': 0, 'standalone_certificates': 0, 'symbolic_nodes': 0,
              'exceptional_scalar_records': 0}

    def F_for(obj):
        f = obj['field']
        return PrimeField(f['p']) if f['e'] == 1 else Field(f['p'], f['modulus'])

    def visit(obj):
        if isinstance(obj, dict):
            if obj.get('schema') == 'blind-pencil-discovery/v1':
                verify_discovery(F_for(obj), obj)
                counts['discoveries'] += 1
                counts['marginal_algebra_ledgers'] += 2
                if obj['complete']:
                    counts['complete_discoveries'] += 1
                    counts['no_anchor_complete'] += obj['code_anchor'] is None
                    counts['symbolic_nodes'] += obj['certificate']['complete_symbolic_node_count']
                    counts['exceptional_scalar_records'] += len(obj['certificate']['exceptional_scalars'])
                return
            if obj.get('schema') == 'polynomial-pencil-tracks/v1':
                verify_certificate(F_for(obj), obj)
                counts['standalone_certificates'] += 1
                return
            for x in obj.values():
                visit(x)
        elif isinstance(obj, list):
            for x in obj:
                visit(x)
    visit(data)
    own = ['pencil_tracks.py', 'run_experiments.py', 'verify_results.py']
    dependencies = ['round5/piece_discovery/piece_discovery.py', 'round4/scalar_fibers/scalar_fibers.py',
                    'round2/proximity/compressed_support.py', 'proximity/extension_field_probe.py',
                    'proximity/deformation_microscope.py']
    report = {'status': 'passed', 'scope': 'source-bound readback sharing production mathematical helpers',
              'counts': counts, 'source_sha256': {name: sha256((HERE / name).read_bytes()).hexdigest() for name in own},
              'dependency_sha256': {name: sha256((LAB / name).read_bytes()).hexdigest() for name in dependencies},
              'results_sha256': sha256((HERE / 'results.json').read_bytes()).hexdigest()}
    (HERE / 'verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
