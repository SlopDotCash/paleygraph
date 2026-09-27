#!/usr/bin/env python3
"""Scale a streaming-state lower bound without quadratic suffix testing.

Universal erasure uniqueness means different nonempty prefixes have disjoint
residual languages. Only actual encodings and prefix distinctness remain to
check; no pairwise cross-splice matrix is needed.
"""
from hashlib import sha256
import json
from pathlib import Path
from lattice_erasure import encode

HERE = Path(__file__).resolve().parent


def main():
    path = HERE/'pullback_erase20.json'; c = json.loads(path.read_text())['certificate']
    assert c['universal_unique_completion'] and c['erased'] == list(range(20))
    count = 100000; prefixes = {}; duplicates = []; stream = sha256()
    for a in range(count):
        word = encode(c['p'], c['g'], c['relation'], a)[0]; prefix = tuple(word[j] for j in c['erased'])
        if prefix in prefixes:
            if len(duplicates) < 12: duplicates.append([prefixes[prefix], a])
        else: prefixes[prefix] = a
        stream.update((json.dumps([a, *prefix], separators=(',', ':'))+'\n').encode())
    out = {'status': 'produced', 'scope': 'Finite state-width lower bound for the stated natural-order cut, obtained from universal erasure uniqueness and a complete specified scalar interval. No all-order or asymptotic lower bound.',
           'certificate_artifact': path.name, 'p': c['p'], 'N': c['N'], 'cut': len(c['erased']),
           'scalar_interval': [0, count-1], 'scalar_words_checked': count,
           'distinct_nonempty_prefixes': len(prefixes), 'certified_width_lower_bound': len(prefixes),
           'pairwise_distinctions_implied': len(prefixes)*(len(prefixes)-1)//2,
           'cross_splice_queries': 0, 'duplicate_prefix_examples': duplicates,
           'prefix_stream_sha256': stream.hexdigest(),
           'source_sha256': {n: sha256((HERE/n).read_bytes()).hexdigest() for n in ['prefix_certificate.py', '../round19/projection_codec.py']},
           'input_sha256': {path.name: sha256(path.read_bytes()).hexdigest()}}
    (HERE/'prefix_certificate.json').write_text(json.dumps(out, indent=2)+'\n'); print(json.dumps(out), flush=True)


if __name__ == '__main__': main()
