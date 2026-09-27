#!/usr/bin/env python3
"""Command-line round trip, multiple completions, budget and corruption controls."""
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from pullback_review import certify_pullback

HERE = Path(__file__).resolve().parent


def run(certificate, word, directory, extra=()):
    path = directory/'input.json'; path.write_text(json.dumps(word)+'\n')
    result = subprocess.run([sys.executable, str(HERE/'erasure_query.py'), '--certificate', str(certificate),
                             '--input', str(path), *extra], text=True, capture_output=True, check=False)
    return result


def main():
    path = HERE/'pullback_erase20.json'; c = json.loads(path.read_text())['certificate']; prepared = certify_pullback(c)
    a, start = 1234567, 57; actual = prepared[2]([a])[0]
    erased = {(start+j) % c['N'] for j in range(20)}; input_word = [None if j in erased else v for j, v in enumerate(actual)]
    (HERE/'erasure_input.json').write_text(json.dumps(input_word, indent=2)+'\n')
    controls = []
    with TemporaryDirectory(prefix='round21-query-') as temp:
        directory = Path(temp)
        out = run(path, input_word, directory, ['--cyclic-start', str(start)])
        assert out.returncode == 0; example = json.loads(out.stdout)
        assert example['status'] == 'complete' and example['count'] == 1 and example['certificate_verified']
        assert example['completions'] == [{'scalar': a, 'digits': actual}]
        (HERE/'erasure_output.json').write_text(json.dumps(example, indent=2)+'\n')
        controls.append({'name': 'cyclic20_recovery', 'status': 'complete', 'count': 1})
        partial32 = [None]*32+actual[32:]
        for label, cert, args in [('half_recovery', HERE/'pullback_erase32.json', []),
                                  ('explicit_budget_failure', HERE/'p2013265921_erase32.json', ['--candidate-budget', '20000'])]:
            result = run(cert, partial32, directory, args); assert result.returncode == 0
            data = json.loads(result.stdout)
            if label == 'half_recovery': assert data['status'] == 'complete' and data['completions'] == [{'scalar': a, 'digits': actual}]
            else: assert data['status'] == 'budget_exceeded' and 'count' not in data and 'completions' not in data
            controls.append({'name': label, 'status': data['status'], 'candidate_box_size': data['candidate_box_size']})
        result = run(HERE/'p17_basic_erase4.json', [None]*4, directory)
        assert result.returncode == 0; data = json.loads(result.stdout)
        assert data['count'] == 17 and [r['scalar'] for r in data['completions']] == list(range(17))
        controls.append({'name': 'multiple_completions', 'status': 'complete', 'count': 17})
        bad_boolean = input_word.copy(); bad_boolean[next(j for j, v in enumerate(input_word) if v is not None)] = True
        wrong_mask = input_word.copy(); wrong_mask[next(j for j, v in enumerate(input_word) if v is None)] = 0
        for label, word, extra in [('wrong_dimension', input_word[:-1], ['--cyclic-start', str(start)]),
                                    ('boolean_digit', bad_boolean, ['--cyclic-start', str(start)]),
                                    ('wrong_erasure_mask', wrong_mask, ['--cyclic-start', str(start)]),
                                    ('invalid_budget', input_word, ['--cyclic-start', str(start), '--candidate-budget', '0'])]:
            result = run(path, word, directory, extra); assert result.returncode != 0
            controls.append({'name': label, 'status': 'rejected'})
        changed = deepcopy(c); changed['centered_coordinate_forms'][0][0][0] += 1
        corrupt = directory/'corrupt.json'; corrupt.write_text(json.dumps(changed))
        result = run(corrupt, input_word, directory, ['--cyclic-start', str(start)])
        assert result.returncode != 0; controls.append({'name': 'corrupt_certificate', 'status': 'rejected'})
    out = {'status': 'passed', 'controls': controls, 'expected_example_scalar': a,
           'source_sha256': {n: sha256((HERE/n).read_bytes()).hexdigest() for n in ['query_review.py', 'erasure_query.py', 'lattice_erasure.py', 'pullback_review.py', 'erasure_review.py']},
           'input_sha256': {n: sha256((HERE/n).read_bytes()).hexdigest() for n in ['pullback_erase20.json', 'pullback_erase32.json', 'p2013265921_erase32.json', 'p17_basic_erase4.json', 'erasure_input.json', 'erasure_output.json']}}
    (HERE/'query_review.json').write_text(json.dumps(out, indent=2)+'\n'); print(json.dumps(out), flush=True)


if __name__ == '__main__': main()
