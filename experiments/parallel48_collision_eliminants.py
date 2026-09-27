#!/usr/bin/env python3
"""Complete fixed-order triple-fiber certificates, not a uniform Paley proof.

Requires python-flint==0.9.0. An isolated installation can be supplied
with --flint-path; no changes to any Lean environment are needed.
"""
import argparse
from collections import Counter
from hashlib import sha256
from math import gcd as integer_gcd, isqrt, prod
from pathlib import Path
import json
import sys
import time

from parallel2_subgroup_2026_09_04 import (
    primitive_polynomial, derivative, gcd as poly_gcd, from_roots,
    generator, primitive_labels, primes_through,
)
from cyclotomic_norm_audit import determinant

ROOT = Path(__file__).resolve().parents[1]


def resultant_matrix(f, b):
    """Multiplication by b in the monic quotient Z[Y]/f."""
    d = len(f) - 1
    matrix = [[0] * d for _ in range(d)]
    for j in range(d):
        column = [0] * j + list(b)
        for k in range(len(column) - 1, d - 1, -1):
            leading = column[k]
            for i, coefficient in enumerate(f):
                column[k - d + i] -= leading * coefficient
        for i in range(d):
            matrix[i][j] = column[i] if i < len(column) else 0
    return matrix


def prime_certificate(p, trial_primes):
    assert p >= 2
    last = 1
    for q in trial_primes:
        if q * q > p:
            return dict(method='trial division through sqrt(p)', last_trial_prime=last)
        assert p % q != 0
        last = q
    assert trial_primes[-1] >= isqrt(p)
    return dict(method='trial division through sqrt(p)', last_trial_prime=last)


def field_check(p, n, f, exponent):
    assert p % n == 1
    g = generator(p, n)
    roots = primitive_labels(p, n, g)
    mult = Counter(roots)
    assert from_roots(roots, p) == [c % p for c in f]
    fp = derivative(f)
    common = poly_gcd(poly_gcd(f, fp, p), derivative(fp), p)
    defect = sum(max(c - 2, 0) for c in mult.values())
    assert len(common) - 1 == defect
    Y = sum(c * (c - 2) for c in mult.values() if c >= 3)
    assert Y <= exponent
    singles = sum(c == 1 for c in mult.values())
    k = n // 2
    B = n * sum(c * c for c in mult.values())
    assert B == 2 * k * k + n * (Y - singles)
    assert B <= 2 * k * k + n * exponent

    K = [pow(g, 2 * j, p) for j in range(k)]
    L = [g * a % p for a in K]
    pairs = Counter((a + b) % p for a in K for b in L)
    assert sum(c * c for c in pairs.values()) == B
    assert max(pairs.values()) == max(mult.values())
    expected = Counter()
    for c in mult.values():
        expected[c] += n
    assert Counter(pairs.values()) == expected
    H = K + L
    total = Counter((a + b) % p for a in H for b in H)
    child = Counter((a + b) % p for a in K for b in K)
    E = sum(c * c for c in total.values())
    E_child = sum(c * c for c in child.values())
    T = sum(c * pairs.get(x, 0) for x, c in child.items())
    assert E == 2 * E_child + 6 * B + 8 * T
    assert T * T <= E_child * B
    assert E <= 4 * E_child + 14 * B and E <= 3 * E_child + 22 * B
    out = dict(p=p, n=n, generator=g, maximum_fiber=max(mult.values()),
        multiplicity_histogram=dict(sorted(Counter(mult.values()).items())),
        triple_gcd_degree=defect, triple_mass_Y=Y, eliminant_valuation=exponent,
        balanced_energy=B, child_energy=E_child, parent_energy=E, T=T,
        quartic_window=n**4 <= 4*p <= 4*n**4,
        literal_pair_checks=True, literal_pair_counts=k*k+n*n,
        roots_sha256=sha256(json.dumps(roots).encode()).hexdigest())
    if out['quartic_window'] and defect:
        x = min(x for x, c in pairs.items() if c == max(pairs.values()))
        Lset = set(L)
        representations = [(a, (x-a) % p) for a in K if (x-a) % p in Lset]
        assert len(representations) == out['maximum_fiber']
        out['largest_fiber_witness'] = dict(x=x, representations=representations)
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--flint-path')
    args = parser.parse_args()
    if args.flint_path:
        sys.path.insert(0, args.flint_path)
    import flint
    assert flint.__version__ == '0.9.0'
    candidates = json.loads((ROOT/'results/parallel48_factorization_candidates_2026_09_06.json').read_text())
    factorizations = candidates['factorizations']
    maximum_prime = max(p for factors in factorizations.values() for p, _ in factors)
    trial_primes = primes_through(isqrt(maximum_prime) + 100)
    prime_certificates = {
        str(p): prime_certificate(p, trial_primes)
        for p in sorted({p for factors in factorizations.values() for p, _ in factors})
    }
    records = []
    for n in [8, 16, 32, 64, 128]:
        started = time.monotonic()
        f = primitive_polynomial(n)
        F = flint.fmpz_poly(f)
        r1 = int(F.resultant(F.derivative()))
        r2 = int(F.resultant(F.derivative().derivative()))
        assert r1 and r2
        matrices = [resultant_matrix(f, b) for b in [derivative(f), derivative(derivative(f))]]
        # Different exact algorithm: matrix determinants versus polynomial resultants.
        matrix_values = [int(flint.fmpz_mat(M).det()) for M in matrices]
        assert matrix_values == [r1, r2]
        bareiss_checked = n <= 64
        if bareiss_checked:
            assert [determinant(M) for M in matrices] == [r1, r2]
        G = integer_gcd(r1, r2)
        factors = factorizations[str(n)]
        assert prod(p**e for p, e in factors) == G
        # No factor is accepted from the factoring library without a separate
        # deterministic primality certificate and this exact product equality.
        fields = [field_check(p, n, f, e) for p, e in factors if p != 2]
        assert all(row['triple_gcd_degree'] > 0 for row in fields)
        quartic = [row for row in fields if row['quartic_window']]
        assert all(row['balanced_energy'] <= 2 * (n//2)**2 for row in quartic)
        records.append(dict(n=n, degree=n//4, coefficients_ascending=f,
            resultant_first_hex=hex(r1), resultant_second_hex=hex(r2),
            gcd_hex=hex(G), complete_factorization=factors,
            matrix_determinants_match=True, stdlib_bareiss_match=bareiss_checked,
            candidate_primes_are_exact_triple_primes=True, fields=fields,
            quartic_triple_primes=[row['p'] for row in quartic],
            seconds=time.monotonic()-started))
        print(json.dumps({'n':n,'certified_triple_primes':len(fields),
            'quartic_triple_primes':records[-1]['quartic_triple_primes']}),flush=True)
    assert records[-1]['quartic_triple_primes'] == [77796353,118593281,181312129]
    assert max(row['p'] for rec in records[:-1] for row in rec['fields']) == 697601
    package_meta = Path(flint.__file__).resolve().parent.parent/'python_flint-0.9.0.dist-info'/'METADATA'
    out = dict(scope='Complete triple-fiber classification at fixed dyadic orders through128; three quartic counterexamples to a maximum-fiber-two hypothesis. No improved uniform energy or Paley theorem.',
        python_version=sys.version.split()[0], python_flint_version=flint.__version__,
        package_metadata_sha256=sha256(package_meta.read_bytes()).hexdigest(),
        prime_certificates=prime_certificates, orders=records,
        total_split_prime_checks=sum(len(rec['fields']) for rec in records),
        literal_pair_counts=sum(row['literal_pair_counts'] for rec in records for row in rec['fields']),
        no_lean_process_started=True, all_passed=True)
    p=ROOT/'results/parallel48_collision_eliminants_2026_09_06.json'
    p.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'output':str(p)}))


if __name__ == '__main__':
    main()
