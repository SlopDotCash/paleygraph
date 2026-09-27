#!/usr/bin/env python3
"""Recover every certified completion of a digit word with null erasures."""
import argparse
import json
from pathlib import Path
from lattice_erasure import decode, decode_cyclic
from erasure_review import certify
from pullback_review import certify_pullback


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--certificate', type=Path, required=True, help='saved case or certificate JSON')
    parser.add_argument('--input', type=Path, required=True, help='JSON array of integer digits or null at erased positions')
    parser.add_argument('--cyclic-start', type=int)
    parser.add_argument('--candidate-budget', type=int, default=50000)
    args = parser.parse_args()
    try:
        data = json.loads(args.certificate.read_text()); c = data.get('certificate', data)
        if 'radius_rule' in c: certify_pullback(c)
        else: certify(c)
        word = json.loads(args.input.read_text())
        if not isinstance(word, list) or len(word) != c['N'] or any(v is not None and type(v) is not int for v in word):
            raise ValueError('input must contain exactly N integer or null entries')
        known = {j: v for j, v in enumerate(word) if v is not None}
        if args.cyclic_start is None: result = decode(c, known, args.candidate_budget)
        else: result = decode_cyclic(c, args.cyclic_start, known, args.candidate_budget)
    except (AssertionError, ValueError, KeyError, IndexError) as exc:
        parser.error(str(exc) or 'certificate identity failed')
    print(json.dumps({'certificate_verified': True, 'prime': c['p'], 'dimension': c['N'],
                      'erased_coordinates': [j for j, v in enumerate(word) if v is None],
                      'universal_candidate_box_cap': c['universal_candidate_box_cap'],
                      'universal_unique_completion': c['universal_unique_completion'], **result}, indent=2))


if __name__ == '__main__': main()
