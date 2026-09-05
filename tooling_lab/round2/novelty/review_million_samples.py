#!/usr/bin/env python3
"""Bounded independent direct-row checks for the million-prime NTT witness.

Adds diagnostic output only to a local copy of the inspected C++ backend.
The oracle computes each selected entry by ordinary integer summation.
This is not an independent full million-point convolution.
"""
import json
from pathlib import Path
import subprocess
from review_spectral_exact import HERE, SPEC, digest, load_witness


def main():
    saved = json.loads((SPEC / 'million_transport_check.json').read_text())
    path = SPEC / saved['file']
    assert digest(path) == saved['sha256']
    p, pairs, chi, z = load_witness(path)
    source = (SPEC / 'exact_ntt.cpp').read_text()
    source_hash = digest(SPEC / 'exact_ntt.cpp')
    row_hook = 'assert(abs(conv)<=absum);A+='
    corr_hook = 'assert(abs(corr)<=n);if(!t)'
    assert source.count(row_hook) == source.count(corr_hook) == 1
    source = source.replace(row_hook,
        'if(i==0 || i==m/7 || i==m/3 || i==m/2 || i==m-1)'
        'cerr<<"ROW "<<x<<" "<<conv<<"\\n"; ' + row_hook)
    source = source.replace(corr_hook,
        'if(t==0 || t==1 || t==2 || t==17 || t==p/7 || t==p/3 || t==p/2 || t==p-1)'
        'cerr<<"CORR "<<t<<" "<<corr<<"\\n"; ' + corr_hook)
    copied = HERE / 'exact_ntt_review_logged.cpp'
    binary = HERE / 'exact_ntt_review_logged'
    copied.write_text(source)
    subprocess.run(['clang++', '-std=c++17', '-O2', '-Wall', '-Wextra', str(copied), '-o', str(binary)], check=True)
    run = subprocess.run([str(binary), str(path)], check=True, capture_output=True, text=True)
    rec = json.loads(run.stdout)
    for key in ('norm_squared', 'sum_entries', 'sum_absolute_entries', 'signed_quadratic_form',
                'max_abs_nonzero_translation_numerator', 'max_translation_shift', 'count_translation_ge_3_5'):
        assert rec[key] == saved[key]
    rows, corrs = [], []
    for line in run.stderr.splitlines():
        kind, index, value = line.split()
        index, value = int(index), int(value)
        if kind == 'ROW':
            expected = sum(v * chi[(index - y) % p] for y, v in pairs)
            rows.append({'row': index, 'value': value})
        else:
            assert kind == 'CORR'
            expected = sum(v * z[(y + index) % p] for y, v in pairs)
            corrs.append({'shift': index, 'value': value})
        assert expected == value
    assert len(rows) == 5 and len(corrs) == 8
    n, sm, L = sum(v*v for v in z), sum(z), sum(abs(v) for v in z)
    assert (n, sm, L) == (rec['norm_squared'], rec['sum_entries'], rec['sum_absolute_entries'])
    a, b = saved['strict_directional_rayleigh_lower_bound']
    X, Y = b * rec['signed_quadratic_form'], a * p * n + b * sm * sm
    square_margin = X * X * p - Y * Y
    assert X > 0 and Y > 0 and square_margin > 0
    assert source_hash == digest(SPEC / 'exact_ntt.cpp')
    out = {'date': '2026-09-05', 'passed': True, 'file': path.name, 'sha256': digest(path),
           'cpp_source_sha256': source_hash, 'logged_cpp_sha256': digest(copied),
           'reviewer_sha256': digest(Path(__file__)), 'p': p, 'm': len(pairs),
           'direct_character_rows': rows, 'direct_autocorrelation_shifts': corrs,
           'directional_lower_bound': [a, b], 'strict_square_margin': square_margin,
           'scope': 'Full freshly compiled candidate replay plus independent sampled rows; not independent full million-point convolution'}
    (HERE / 'million_samples_review.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
