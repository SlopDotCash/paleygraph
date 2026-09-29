#!/usr/bin/env python3
"""Independent dense orbit census and uniqueness proof for the width bound."""
from hashlib import sha256
import json
from pathlib import Path
from pullback_review import certify_pullback

HERE = Path(__file__).resolve().parent


def main():
    source = HERE/'prefix_certificate.json'; r = json.loads(source.read_text()); path = HERE/r['certificate_artifact']
    c = json.loads(path.read_text())['certificate']; prepared = certify_pullback(c)
    assert c['universal_unique_completion'] and c['erased'] == list(range(r['cut']))
    assert (c['p'], c['N']) == (r['p'], r['N'])
    lo, hi = r['scalar_interval']; prefixes = {}; duplicates = []; stream = sha256()
    assert 0 <= lo <= hi < c['p']
    for start in range(lo, hi+1, 2048):
        scalars = list(range(start, min(hi+1, start+2048))); words = prepared[2](scalars)
        for a, word in zip(scalars, words):
            P = tuple(word[:r['cut']])
            if P in prefixes:
                if len(duplicates) < 12: duplicates.append([prefixes[P], a])
            else: prefixes[P] = a
            stream.update((json.dumps([a, *P], separators=(',', ':'))+'\n').encode())
    assert hi-lo+1 == r['scalar_words_checked']
    assert len(prefixes) == r['distinct_nonempty_prefixes'] == r['certified_width_lower_bound']
    assert len(prefixes)*(len(prefixes)-1)//2 == r['pairwise_distinctions_implied']
    assert duplicates == r['duplicate_prefix_examples'] and stream.hexdigest() == r['prefix_stream_sha256']
    assert r['cross_splice_queries'] == 0
    out = {'status': 'passed', 'scope': 'Universal exact erasure-uniqueness certificate plus direct-power/dense-matrix scalar census. Distinct nonempty prefixes are pairwise distinguished because their complete residual languages are disjoint.',
           'scalar_words_checked': hi-lo+1, 'certified_width_lower_bound': len(prefixes),
           'pairwise_distinctions_implied': r['pairwise_distinctions_implied'], 'prefix_stream_sha256': stream.hexdigest(),
           'source_sha256': {n: sha256((HERE/n).read_bytes()).hexdigest() for n in ['prefix_review.py', 'pullback_review.py', 'erasure_review.py']},
           'input_sha256': {source.name: sha256(source.read_bytes()).hexdigest(), path.name: sha256(path.read_bytes()).hexdigest()}}
    (HERE/'prefix_review.json').write_text(json.dumps(out, indent=2)+'\n'); print(json.dumps(out), flush=True)


if __name__ == '__main__': main()
