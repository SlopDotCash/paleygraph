#!/usr/bin/env python3
"""Complete small scalar families and symbolic full generator-2 families."""
from hashlib import sha256
import json
from pathlib import Path
import sys
from decision_diagram import from_words, generator_two, complete

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'round19'))
from projection_codec import encode


def queries(graph, anchors):
    N = graph['dimension']; specs = [{}, {0: 1000}]
    for w in anchors:
        for depth in (1, N//2, N): specs.append(dict(enumerate(w[:depth])))
        specs.append({i: w[i] for i in range(0, N, 2)})
        specs.append({0: w[0], N-1: w[-1]})
    result = []; seen = set()
    for spec in specs:
        key = tuple(sorted(spec.items()))
        if key in seen: continue
        seen.add(key)
        out = complete(graph, spec, 12)
        result.append({'fixed': [list(x) for x in key], **out})
    return result


def save_case(name, config, order, graph, qs, census):
    target = HERE/(name+'.json')
    target.write_text(json.dumps({'name': name, 'config': config, 'order': order,
                                 'graph': graph, 'queries': qs, 'complete_census': census},
                                separators=(',', ':'))+'\n')
    row = {'name': name, 'artifact': target.name, 'p': config.get('p'),
           'N': graph['dimension'], 'construction': graph['construction'],
           'word_count': graph['word_count'], 'widths': graph['widths'],
           'max_width': max(graph['widths']), 'nodes': graph['nodes'], 'edges': graph['edges'],
           'raw_states': graph['raw_states'], 'queries': len(qs)}
    print(json.dumps(row), flush=True); return row


def main():
    cases = []; saved_path = HERE.parent/'round17/norm_compression.json'
    saved = json.loads(saved_path.read_text())
    for c, e, h in zip(saved['cases'][:2], [11, 19], [5, 0]):
        p, N = c['p'], c['N']
        for label, g, f, order in [
                ('canonical', 2, [2]+[0]*(N-2)+[1], list(range(N))),
                ('saved', c['g'], c['relation'], list(range(N))),
                ('saved_reordered', c['g'], c['relation'], sorted(range(N), key=lambda j: (e*j-h) % N))]:
            words = []
            for a in range(p):
                digits, _ = encode(p, g, f, a)
                words.append(tuple(digits[i] for i in order))
            assert len(set(words)) == p
            graph = from_words(words)
            qs = queries(graph, [words[a] for a in [0, 1, (p+1)//2, p//3]])
            cases.append(save_case(f'p{p}_{label}', {'p': p, 'g': g, 'relation': f}, order, graph, qs, True))
    family_path = HERE.parent/'round19/codec_family_controls.json'
    family = json.loads(family_path.read_text())
    for c in family['cases'][1:]:
        p, N, g = c['p'], c['N'], c['g']
        f = min((r['relation'] for r in c['records']), key=lambda x: (sum(map(abs, x)), x))
        words = [encode(p, g, f, a)[0] for a in range(p)]
        graph = from_words(words); qs = queries(graph, [words[a] for a in [0, 1, (p+1)//2, p//3]])
        cases.append(save_case(f'p{p}_general', {'p': p, 'g': g, 'relation': f}, list(range(N)), graph, qs, True))
    # No scalar census and no enumeration of ternary ambient words here.
    for N, k in [(8, 1), (16, 1), (32, 641), (64, 1), (128, 1)]:
        graph = generator_two(N, k); Q = (1 << N)+1; p = Q//k
        config = {'p': p if (N, k) in [(8, 1), (16, 1), (32, 641)] else None,
                  'Q': Q, 'k': k, 'g': 2, 'relation': [2]+[0]*(N-2)+[1]}
        anchors = complete(graph, {}, 3)['words']+[ [0]*N ]
        qs = queries(graph, anchors)
        cases.append(save_case(f'symbolic_N{N}_k{k}', config, list(range(N)), graph, qs, False))
    out = {'status': 'produced', 'scope': 'Exact finite layered completion diagrams. Fixed read order, implicit rejecting state excluded. Full scalar census only where marked; symbolic Q families need not be prime.',
           'cases': cases,
           'source_sha256': {n: sha256((HERE/n).read_bytes()).hexdigest() for n in ['build_completions.py', 'decision_diagram.py', '../round19/projection_codec.py']},
           'input_sha256': {'../'+str(p.relative_to(HERE.parent)): sha256(p.read_bytes()).hexdigest() for p in [saved_path, family_path]}}
    (HERE/'completion_summary.json').write_text(json.dumps(out, indent=2)+'\n')


if __name__ == '__main__': main()
