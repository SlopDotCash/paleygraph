#!/usr/bin/env python3
"""Independent integer-packing oracle for the exact spectral backend.

Expected convolutions use Python's arbitrary-precision integer multiplication,
not an FFT/NTT. The candidate C++ backend is freshly compiled in this folder.
"""
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import random
import struct
import subprocess
import sys
from decimal import Decimal, localcontext

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
SPEC = HERE.parent / 'spectral'
MOD = 998244353


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def packed_cyclic(a, b):
    """Exact cyclic convolution through nonnegative base-256^w packing.

    Offsets remove negative coefficients. The radix exceeds the largest
    possible linear coefficient, so integer multiplication has no carries
    between coefficients. Subtract the constant offset after cyclic folding.
    """
    p = len(a)
    assert len(b) == p
    oa, ob = max(0, -min(a)), max(0, -min(b))
    aa, bb = [x + oa for x in a], [x + ob for x in b]
    bound = p * max(aa) * max(bb)
    width = max(1, (bound.bit_length() + 7) // 8)
    assert 1 << (8 * width) > bound
    left = int.from_bytes(b''.join(x.to_bytes(width, 'little') for x in aa), 'little')
    right = int.from_bytes(b''.join(x.to_bytes(width, 'little') for x in bb), 'little')
    raw = (left * right).to_bytes((2 * p - 1) * width, 'little')
    coeff = [int.from_bytes(raw[i * width:(i + 1) * width], 'little')
             for i in range(2 * p - 1)]
    correction = ob * sum(a) + oa * sum(b) + p * oa * ob
    answer = [coeff[t] + (coeff[t + p] if t + p < len(coeff) else 0) - correction
              for t in range(p)]
    return answer, {'radix_bytes': width, 'linear_coefficient_bound': bound}


def load_witness(path):
    raw = path.read_bytes()
    p, m = struct.unpack('<II', raw[:8])
    assert len(raw) == 8 + 8 * m and m == (p - 5) // 4
    assert p % 4 == 1 and all(p % d for d in range(2, math.isqrt(p) + 1))
    pairs = list(struct.iter_unpack('<ii', raw[8:]))
    chi = [0 if x == 0 else (1 if pow(x, (p - 1) // 2, p) == 1 else -1)
           for x in range(p)]
    assert [x for x, _ in pairs] == [x for x in range(p) if chi[x] == chi[(x - 1) % p] == 1]
    assert max(abs(x) for _, x in pairs) <= 1024
    z = [0] * p
    for x, v in pairs:
        z[x] = v
    return p, pairs, chi, z


def sign_tests():
    spec = importlib.util.spec_from_file_location('review_candidate_above', SPEC / 'verify_exact.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    count = 0
    with localcontext() as ctx:
        ctx.prec = 90
        for p in (2, 5, 9, 101):
            for a in range(-3, 4):
                for b in range(1, 5):
                    for A in range(-3, 4):
                        for sm in range(-2, 3):
                            for n in (1, 2):
                                for sign in (-1, 1):
                                    X, Y = b * sign * A, a * p * n + b * sign * sm * sm
                                    # Handle equality algebraically, then a high-precision
                                    # independent real-arithmetic oracle on bounded integers.
                                    equal = X * X * p == Y * Y and (X > 0) - (X < 0) == (Y > 0) - (Y < 0)
                                    expected = False if equal else Decimal(X) * Decimal(p).sqrt() > Decimal(Y)
                                    assert module.above(a, b, p, A, sm, n, sign) == expected
                                    count += 1
    return count


def replay(path, binary):
    p, pairs, chi, z = load_witness(path)
    rec = json.loads(subprocess.check_output([str(binary), str(path)], text=True))
    n, sm, L = sum(v * v for v in z), sum(z), sum(abs(v) for v in z)
    assert (rec['norm_squared'], rec['sum_entries'], rec['sum_absolute_entries']) == (n, sm, L)
    conv, packing = packed_cyclic(chi, z)
    assert all(abs(x) <= L for x in conv)
    A = sum(v * conv[x] for x, v in pairs)
    assert A == rec['signed_quadratic_form']
    if path.name.endswith('negative_512.bin'):
        stored = json.loads((SPEC / 'extra_transport_check.json').read_text())
        assert digest(path) == stored['sha256']
        a, b = stored['strict_directional_rayleigh_lower_bound']
        sign = -1
    else:
        stored = json.loads((SPEC / f'results_{p}.json').read_text())
        side = next(side for side in stored['sides'] if side['witness']['file'] == path.name)
        assert digest(path) == side['witness']['sha256']
        sign = 1 if side['side'] == 'positive' else -1
        b = 10 ** 6
        a = math.floor(sign * side['witness']['rounded_rayleigh_H_float'] * b) - 1
    X, Y = b * sign * A, a * p * n + b * sign * sm * sm
    # Every stored witness lies in this sign case; no floating arithmetic
    # participates in the following acceptance check.
    assert X > 0 and Y > 0
    square_margin = X * X * p - Y * Y
    assert square_margin > 0
    rng = random.Random(20260905)
    indices = list(range(len(pairs))) if p <= 4001 else sorted({0, len(pairs) - 1, *rng.sample(range(len(pairs)), 30)})
    for i in indices:
        x = pairs[i][0]
        assert conv[x] == sum(chi[(x - y) % p] * v for y, v in pairs)
    autocorr = None
    if 2 * n < MOD:
        corr, cp = packed_cyclic(z, [z[-t % p] for t in range(p)])
        assert corr[0] == n and all(abs(v) <= n for v in corr)
        maximum = max(abs(v) for v in corr[1:])
        shift = next(t for t in range(1, p) if abs(corr[t]) == maximum)
        count = sum(5 * abs(v) >= 3 * n for v in corr)
        assert (maximum, shift, count) == (rec['max_abs_nonzero_translation_numerator'],
                                          rec['max_translation_shift'], rec['count_translation_ge_3_5'])
        assert rec.get('autocorrelation_computed', True)
        sampled_shifts = sorted({0, 1, p - 1, shift, p - shift, *rng.sample(range(p), min(10, p))})
        for t in sampled_shifts:
            assert corr[t] == sum(v * z[(x + t) % p] for x, v in pairs)
        autocorr = {'all_shifts': p, 'max_abs_nonzero_numerator': maximum, 'max_shift': shift,
                    'count_ge_3_5_including_zero': count, 'direct_sampled_shifts': len(sampled_shifts),
                    'packing': cp}
    else:
        assert not rec.get('autocorrelation_computed', False)
        assert rec['max_abs_nonzero_translation_numerator'] in (None, -1)
    return {'file': path.name, 'sha256': digest(path), 'p': p, 'm': len(pairs),
            'norm_squared': n, 'sum_entries': sm, 'signed_quadratic_form': A,
            'directional_lower_bound': [a, b], 'strict_square_margin': square_margin,
            'direct_rows_checked': len(indices), 'packing': packing, 'autocorrelation': autocorr}


def main():
    binary = HERE / 'exact_ntt_review'
    cpp_hash = digest(SPEC / 'exact_ntt.cpp')
    subprocess.run(['clang++', '-std=c++17', '-O2', '-Wall', '-Wextra', str(SPEC / 'exact_ntt.cpp'),
                    '-o', str(binary)], check=True)
    tiny = []
    for p in (101, 401, 1009, 4001):
        for side in ('positive', 'negative'):
            tiny.append(replay(SPEC / f'witness_{p}_{side}.bin', binary))
    large = [replay(SPEC / f'witness_65537_{suffix}.bin', binary)
             for suffix in ('positive', 'negative', 'negative_512')]
    # Worst-case integer bounds for the explicitly reviewed input envelope.
    maxp, maxcoef = 1000033, 1024
    maxm = (maxp - 5) // 4
    maxL, maxn = maxm * maxcoef, maxm * maxcoef ** 2
    bounds = {'p_at_most': maxp, 'coefficient_at_most': maxcoef, 'm_at_most': maxm,
              'sum_absolute_at_most': maxL, 'norm_squared_at_most': maxn,
              'quadratic_form_abs_at_most': maxL ** 2,
              'ntt_residue_product_at_most': (MOD - 1) ** 2, 'int64_max': 2 ** 63 - 1}
    assert maxL ** 2 < 2 ** 63 and (MOD - 1) ** 2 < 2 ** 63 and 5 * maxn < 2 ** 63
    assert 2 * maxL < MOD and 2 * maxp - 1 < 2 ** 21
    assert cpp_hash == digest(SPEC / 'exact_ntt.cpp'), 'Source changed during review; replay required'
    out = {'date': '2026-09-05', 'passed': True, 'method': 'independent Python integer packing plus direct row sums',
           'cpp_source_sha256': digest(SPEC / 'exact_ntt.cpp'), 'python_source_sha256': digest(SPEC / 'verify_exact.py'),
           'reviewer_sha256': digest(Path(__file__)), 'sign_cases_checked': sign_tests(),
           'small_witnesses': tiny, 'large_witnesses': large, 'int64_bounds': bounds,
           'scope': 'Exact finite witnesses; million-prime numerical solver not invoked; no spectral upper bound'}
    (HERE / 'spectral_exact_review.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({'passed': True, 'sign_cases': out['sign_cases_checked'],
                      'witnesses': len(tiny) + len(large), 'large_witnesses': large, 'int64_bounds': bounds}, indent=2))


if __name__ == '__main__':
    main()
