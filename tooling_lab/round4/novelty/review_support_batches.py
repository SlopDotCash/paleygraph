#!/usr/bin/env python3
"""Independent finite-field polynomial/scalar and selected-base batch oracle."""
import sys
sys.dont_write_bytecode = True
from collections import defaultdict
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, product
from math import comb
from pathlib import Path
import json
import random

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'support_batches' / 'symbolic_batches.py'
sys.path.insert(0, str(SOURCE.parent))


def load(path):
    spec = spec_from_file_location(path.stem, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class OracleField:
    """Prime fields or F3[X]/(X^2+1), with inverses by exhaustive search."""
    def __init__(self, q):
        self.q = q
        self.p = 3 if q == 9 else q
    def add(self, a, b):
        if self.q == 9:
            return (a%3+b%3)%3 + 3*((a//3+b//3)%3)
        return (a+b)%self.q
    def neg(self, a):
        if self.q == 9:
            return (-(a%3))%3 + 3*(-(a//3)%3)
        return -a % self.q
    def sub(self, a, b):
        return self.add(a, self.neg(b))
    def mul(self, a, b):
        if self.q == 9:
            return ((a%3)*(b%3)-(a//3)*(b//3))%3 + 3*(((a%3)*(b//3)+(a//3)*(b%3))%3)
        return a*b % self.q
    def inv(self, a):
        return next(x for x in range(self.q) if self.mul(a, x) == 1)
    def eval(self, coeffs, x):
        out = 0
        for c in reversed(coeffs):
            out = self.add(self.mul(out, x), c)
        return out
    def lagrange(self, xs, ys, domain):
        result = []
        for z in domain:
            answer = 0
            for i, x in enumerate(xs):
                numerator = denominator = 1
                for j, y in enumerate(xs):
                    if i != j:
                        numerator = self.mul(numerator, self.sub(z, y))
                        denominator = self.mul(denominator, self.sub(x, y))
                answer = self.add(answer, self.mul(ys[i], self.mul(numerator, self.inv(denominator))))
            result.append(answer)
        return tuple(result)


def complete_nodes(F, domain, k, s, u0, u1):
    out = set()
    for coeffs in product(range(F.q), repeat=k):
        word = tuple(F.eval(coeffs, x) for x in domain)
        for scalar in range(F.q):
            support = tuple(i for i, a in enumerate(word) if a == F.add(u0[i], F.mul(scalar, u1[i])))
            if len(support) >= s:
                out.add((scalar, word, support))
    return out


def expand_record(F, r):
    out = {(x['scalar'], tuple(x['codeword']), tuple(x['agreement_support'])) for x in r['nodes']}
    for track in r['all_field_track_certificates']:
        for scalar in range(F.q):
            word = tuple(F.add(a, F.mul(scalar, b)) for a, b in zip(track['intercept'], track['slope']))
            support = tuple(i for i, a in enumerate(word) if a == F.add(r['u0'][i], F.mul(scalar, r['u1'][i])))
            assert len(support) >= r['s']
            out.add((scalar, word, support))
    return out


def check_batches(F, r):
    plan = r['cover_certificate']
    assert (plan['n'], plan['k'], plan['s']) == (r['n'], r['k'], r['s'])
    selected = {B for block in plan['blocks'] for B in combinations(block, r['k'])}
    assert len(selected) == plan['basis_count']
    for support in combinations(range(r['n']), r['s']):
        assert any(set(B) <= set(support) for B in selected)
    actual = defaultdict(set)
    for base in selected:
        xs = [r['domain'][i] for i in base]
        A = F.lagrange(xs, [r['u0'][i] for i in base], r['domain'])
        B = F.lagrange(xs, [r['u1'][i] for i in base], r['domain'])
        actual[A, B].add(base)
    assert len(actual) == len(r['batch_certificates'])
    seen = set()
    for batch in r['batch_certificates']:
        track = (tuple(batch['intercept']), tuple(batch['slope']))
        assert track in actual and track not in seen
        seen.add(track)
        common = {i for i in range(r['n']) if track[0][i] == r['u0'][i] and track[1][i] == r['u1'][i]}
        bases = {B for B in selected if set(B) <= common}
        assert bases == actual[track]
        assert len(bases) == batch['base_multiplicity']
        assert tuple(batch['first_base']) in bases
        assert sorted(common) == batch['common_coordinates']
    assert sum(map(len, actual.values())) == len(selected)
    return len(selected)


def family_checks(candidate):
    rng = random.Random(47702)
    count = 0
    for n in range(1, 9):
        for k in range(n+1):
            zdd = candidate.FamilyZDD()
            root = zdd.choose(range(n), k)
            exact = set(combinations(range(n), k))
            for _ in range(50):
                allowed = {i for i in range(n) if rng.randrange(2)}
                root = zdd.subtract_subsets(root, sum(1 << i for i in allowed))
                exact = {B for B in exact if not set(B) <= allowed}
                assert zdd.materialize(root) == exact and zdd.counts[root] == len(exact)
                if exact:
                    assert zdd.first(root) in exact
                count += 1
    return count


def main():
    candidate = load(SOURCE)
    verifier = load(SOURCE.parent / 'verify_exports.py')
    hashes = {p.name: sha256(p.read_bytes()).hexdigest() for p in (SOURCE, SOURCE.parent/'verify_exports.py')}
    family_count = family_checks(candidate)
    rng = random.Random(89182)
    checked = bases = whole = symbolic = 0
    examples = []
    for q, n, max_k in ((3, 3, 3), (5, 5, 3), (9, 6, 2)):
        F = OracleField(q)
        actual_field = candidate.Field(3, [1, 0, 1]) if q == 9 else candidate.PrimeField(q)
        domain = list(range(n))
        pencils = [([0]*n, [0]*n),
                   ([rng.randrange(q) for _ in domain], [rng.randrange(q) for _ in domain]),
                   ([i % q for i in domain], [(2*i+1) % q for i in domain])]
        for k in range(1, max_k+1):
            for s in range(k, n+1):
                for pencil_id, (u0, u1) in enumerate(pencils):
                    expected = complete_nodes(F, domain, k, s, u0, u1)
                    for mode in ('partition', 'anchor', 'full_k'):
                        limit = 0 if pencil_id == 0 else 1000
                        r = candidate.run_census(actual_field, domain, k, s, u0, u1, mode, limit)
                        verifier.verify_export(actual_field, r)
                        assert expand_record(F, r) == expected, (q, n, k, s, pencil_id, mode)
                        scalar_set = sorted({x[0] for x in expected})
                        assert r['finite_bad_scalar_count'] == len(scalar_set)
                        assert r['bad_scalars'] == (None if r['whole_field_correlated'] and limit == 0 else scalar_set)
                        bases += check_batches(F, r)
                        whole += bool(r['whole_field_correlated'])
                        symbolic += not r['node_ledger_materialized']
                        checked += 1
        examples.append({'q': q, 'n': n, 'maximum_k': max_k})
    for filename, digest in hashes.items():
        assert sha256((SOURCE.parent/filename).read_bytes()).hexdigest() == digest
    out = {'status': 'passed', 'date': '2026-09-05', 'configurations': examples,
           'complete_polynomial_scalar_comparisons': checked,
           'independent_selected_base_interpolations': bases,
           'whole_field_cases': whole, 'symbolic_node_ledger_cases': symbolic,
           'sequential_family_subtraction_checks': family_count,
           'candidate_sha256': hashes, 'reviewer_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
           'scope': 'Independent prime/GF9 field arithmetic, complete coefficient/scalar enumeration, explicit selected bases; does not rerun large structured cases.'}
    (HERE/'support_batches_review.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
