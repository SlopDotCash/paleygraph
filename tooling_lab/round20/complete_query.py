#!/usr/bin/env python3
"""Query complete diagrams with assignments in ORIGINAL digit coordinates."""
import argparse
import json
from pathlib import Path
import sys
from decision_diagram import complete

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'round19'))
from projection_codec import prepare, decode


def query(path, fixed, limit):
    r = json.loads(path.read_text()); order = r['order']; N = len(order)
    if any(type(i) is not int or i not in range(N) for i in fixed): raise ValueError('coordinate outside original word')
    assignments = {order.index(i): v for i, v in fixed.items()}
    result = complete(r['graph'], assignments, limit); c = r['config']; words = []
    codec = prepare(c['p'], c['g'], c['relation']) if c['p'] is not None else None
    for word in result.pop('words'):
        original = [0]*N
        for pos, j in enumerate(order): original[j] = word[pos]
        if codec is not None:
            decoded = decode(codec, original); assert decoded['status'].startswith('accepted_')
            scalar, modulus = decoded['extracted_scalar'], c['p']
        else:
            S = sum((1 << (N-1-j))*v for j, v in enumerate(original)); assert S % c['k'] == 0
            modulus = c['Q']//c['k']; scalar = S//c['k'] % modulus
        words.append({'scalar_residue': scalar, 'scalar_modulus': modulus, 'digits': original})
    return {'diagram': path.name, 'fixed_original_coordinates': fixed,
            'enumeration_order': 'lexicographic in the stored read order', **result, 'completions': words}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--diagram', type=Path, required=True)
    parser.add_argument('--fixed', default='', help='original coordinate=value pairs, comma-separated')
    parser.add_argument('--limit', type=int, default=12)
    args = parser.parse_args(); fixed = {}
    for item in filter(None, args.fixed.split(',')):
        i, v = map(int, item.split('='))
        if i in fixed: raise ValueError('repeated coordinate')
        fixed[i] = v
    print(json.dumps(query(args.diagram, fixed, args.limit), indent=2))


if __name__ == '__main__': main()
