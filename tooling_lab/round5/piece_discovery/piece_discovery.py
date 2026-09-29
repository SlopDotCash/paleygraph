#!/usr/bin/env python3
"""Bounded monic polynomial reconstruction feeding frozen piece certificates.

Interpolation, finite-field splitting and Hensel lifting are known methods.
Discovery is deliberately incomplete; every successful certificate is exact.
"""
from pathlib import Path
import sys
import random
import time
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
LAB = HERE.parents[1]
sys.path.insert(0, str(LAB / 'round4' / 'scalar_fibers'))
from scalar_fibers import (PrimeField, Field, evaluate, canonical, field_record,
    find_code_anchor, certify_static_word, certify_pencil, verify_certificate,
    UnavailableCertificate)


def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a or [0]


def add(F, a, b):
    return trim([F.add(a[i] if i < len(a) else 0, b[i] if i < len(b) else 0)
                 for i in range(max(len(a), len(b)))])


def scale(F, a, c):
    return trim([F.mul(x, c) for x in a])


def mul(F, a, b, limit=None):
    size = len(a) + len(b) - 1
    if limit is not None:
        size = min(size, limit)
    c = [0] * max(1, size)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b[:max(0, size - i)]):
                if y:
                    c[i + j] = F.add(c[i + j], F.mul(x, y))
    return trim(c)


def divmod_poly(F, a, b):
    a, b = trim(a), trim(b)
    assert b != [0]
    q = [0] * max(1, len(a) - len(b) + 1)
    inverse = F.inv(b[-1])
    while a != [0] and len(a) >= len(b):
        i = len(a) - len(b)
        c = F.mul(a[-1], inverse)
        q[i] = c
        for j, v in enumerate(b):
            a[i + j] = F.sub(a[i + j], F.mul(c, v))
        a = trim(a)
    return trim(q), a


def gcd(F, a, b):
    a, b = trim(a), trim(b)
    while trim(b) != [0]:
        a, b = b, divmod_poly(F, a, b)[1]
    return scale(F, a, F.inv(a[-1])) if trim(a) != [0] else [0]


def powmod_poly(F, a, exponent, modulus):
    out = [1]
    while exponent:
        if exponent & 1:
            out = divmod_poly(F, mul(F, out, a), modulus)[1]
        a = divmod_poly(F, mul(F, a, a), modulus)[1]
        exponent >>= 1
    return out


def simple_roots(F, poly, seed=0, budget=64):
    """Bounded Cantor-Zassenhaus splitting, not an enumeration of field values."""
    poly = trim(poly)
    assert poly[-1] == 1
    derivative = [F.mul(j % F.p, poly[j]) for j in range(1, len(poly))] or [0]
    if len(gcd(F, poly, derivative)) > 1:
        return None, {'status': 'not_squarefree', 'splitting_attempts': 0}
    if divmod_poly(F, add(F, powmod_poly(F, [0, 1], F.q, poly), [0, F.sub(0, 1)]), poly)[1] != [0]:
        return None, {'status': 'not_completely_split', 'splitting_attempts': 0}
    rng = random.Random(seed)
    attempts = 0

    def split(f):
        nonlocal attempts
        if len(f) == 2:
            return [F.mul(F.sub(0, f[0]), F.inv(f[1]))]
        if len(f) <= 1:
            return []
        for _ in range(budget):
            attempts += 1
            a = [rng.randrange(F.q) for _ in range(len(f) - 1)]
            d = gcd(F, f, a)
            if not 1 < len(d) < len(f):
                # Odd characteristic is an explicit supported domain.
                b = powmod_poly(F, a, (F.q - 1) // 2, f)
                d = gcd(F, f, add(F, b, [F.sub(0, 1)]))
            if 1 < len(d) < len(f):
                left = split(d)
                right = split(divmod_poly(F, f, d)[0])
                return None if left is None or right is None else left + right
        return None

    roots = split(poly)
    if roots is None:
        return None, {'status': 'split_budget_exhausted', 'splitting_attempts': attempts}
    roots.sort()
    check = [1]
    for root in roots:
        check = mul(F, check, [F.sub(0, root), 1])
    assert check == poly
    return roots, {'status': 'all_simple_roots_verified', 'splitting_attempts': attempts}


def shift(F, a, center, limit=None):
    """Return a(center+T), including in small characteristic."""
    out = [0]
    for c in reversed(a):
        out = add(F, mul(F, out, [center, 1], limit), [c])
    return out[:limit] if limit else out


def substitute_y(F, relation, h, limit=None):
    out = [0]
    for a in reversed(relation):
        out = add(F, mul(F, out, h, limit), a[:limit] if limit else a)
        if limit:
            out = out[:limit]
    return trim(out)


def recover_branches(F, relation, k, center_budget=24, split_budget=64):
    """Lift simple roots, then verify global polynomial factors exactly."""
    if F.p == 2:
        return {'status': 'unsupported_even_characteristic', 'branches': [], 'centers': []}
    centers = list(range(min(F.q, center_budget, 8)))
    rng = random.Random(70399 + F.q + k)
    # The search budget stays constant as q grows.
    sampled = rng.sample(range(F.q), min(F.q, center_budget))
    centers = list(dict.fromkeys(centers + sampled))[:center_budget]
    records = []
    for center in centers:
        at_center = trim([evaluate(F, a, center) for a in relation])
        roots, record = simple_roots(F, at_center, seed=9317 + center, budget=split_budget)
        records.append({'center': center, **record})
        if roots is None:
            continue
        shifted = [shift(F, a, center, k) for a in relation]
        branches = []
        for root in roots:
            derivative = F.sum(F.mul(j % F.p, F.mul(at_center[j], F.power(root, j - 1)))
                               for j in range(1, len(at_center)))
            assert derivative
            inverse = F.inv(derivative)
            h = [root]
            for ell in range(1, k):
                value = substitute_y(F, shifted, h, ell + 1)
                coefficient = value[ell] if ell < len(value) else 0
                h.append(F.mul(F.sub(0, coefficient), inverse))
            global_h = shift(F, h, F.sub(0, center))
            if substitute_y(F, relation, global_h) == [0]:
                branches.append(global_h + [0] * (k - len(global_h)))
        # Each branch was checked as a complete identity, not merely a jet.
        factor_product = [[1]]
        for h in branches:
            new = [[0] for _ in range(len(factor_product) + 1)]
            for j, a in enumerate(factor_product):
                new[j] = add(F, new[j], scale(F, mul(F, a, h), F.sub(0, 1)))
                new[j + 1] = add(F, new[j + 1], a)
            factor_product = new
        complete = [trim(a) for a in factor_product] == [trim(a) for a in relation]
        return {'status': 'complete_polynomial_factorization' if complete else 'some_simple_branches_nonpolynomial',
                'branches': branches, 'centers': records, 'chosen_center': center,
                'factor_product_verified': complete,
                'root_method': 'bounded Cantor-Zassenhaus, then Hensel coefficients and exact identities'}
    return {'status': 'no_simple_split_center_within_budget', 'branches': [], 'centers': records}


def system(F, dom, word, k, r):
    monomials = [(i, j) for j in range(r) for i in range((r - j) * (k - 1) + 1)]
    matrix = []
    rhs = []
    for x, y in zip(dom, word):
        xp, yp = [1], [1]
        for _ in range(r * (k - 1)):
            xp.append(F.mul(xp[-1], x))
        for _ in range(r):
            yp.append(F.mul(yp[-1], y))
        matrix.append([F.mul(xp[i], yp[j]) for i, j in monomials])
        rhs.append(F.sub(0, yp[r]))
    return monomials, matrix, rhs


def solve_exact(F, matrix, rhs):
    """Echelon solve with particular/nullspace and independent-row witnesses."""
    rows, cols = len(matrix), len(matrix[0])
    use_numpy = F.e == 1 and F.p <= 2147483713
    if use_numpy:
        import numpy as np
        a = np.array([row + [b] for row, b in zip(matrix, rhs)], dtype=np.int64)
    else:
        a = [list(row) + [b] for row, b in zip(matrix, rhs)]
    row_ids = list(range(rows))
    pivots = []
    rank = 0
    for col in range(cols):
        pivot = next((i for i in range(rank, rows) if a[i][col]), None)
        if pivot is None:
            continue
        if use_numpy:
            a[[rank, pivot]] = a[[pivot, rank]]
        else:
            a[rank], a[pivot] = a[pivot], a[rank]
        row_ids[rank], row_ids[pivot] = row_ids[pivot], row_ids[rank]
        inverse = F.inv(int(a[rank][col]))
        if use_numpy:
            a[rank, col:] = a[rank, col:] * inverse % F.p
            if rank + 1 < rows:
                factors = a[rank + 1:, col].copy()
                a[rank + 1:, col:] = (a[rank + 1:, col:] - factors[:, None] * a[rank, col:]) % F.p
        else:
            a[rank][col:] = [F.mul(x, inverse) for x in a[rank][col:]]
            for i in range(rank + 1, rows):
                factor = a[i][col]
                if factor:
                    for j in range(col, cols + 1):
                        a[i][j] = F.sub(a[i][j], F.mul(factor, a[rank][j]))
        pivots.append(col)
        rank += 1
        if rank == rows:
            break
    if use_numpy:
        a = a.tolist()
    free = [j for j in range(cols) if j not in pivots]

    def backsolve(free_index=None, homogeneous=False):
        x = [0] * cols
        if free_index is not None:
            x[free_index] = 1
        for i in range(rank - 1, -1, -1):
            col = pivots[i]
            value = 0 if homogeneous else a[i][-1]
            for j in range(col + 1, cols):
                if x[j] and a[i][j]:
                    value = F.sub(value, F.mul(a[i][j], x[j]))
            x[col] = value
        return x

    nullspace = [backsolve(j, True) for j in free]
    inconsistent = next((i for i in range(rank, rows) if a[i][-1]), None)
    record = {'rank': rank, 'unknowns': cols, 'equations': rows,
              'pivot_columns': pivots, 'independent_rows': row_ids[:rank],
              'free_columns': free, 'kernel_basis': nullspace,
              'status': 'inconsistent' if inconsistent is not None else 'consistent',
              'backend': 'prime NumPy int64 echelon' if use_numpy else 'generic exact field echelon'}
    if inconsistent is None:
        record['particular'] = backsolve()
    else:
        # A selected nonsingular minor expresses the conflicting original row
        # in the certified independent rows. This gives a short dual witness.
        original_i = row_ids[inconsistent]
        pivot_rows = row_ids[:rank]
        if rank:
            transposed = [[matrix[i][j] for i in pivot_rows] for j in pivots]
            target = [matrix[original_i][j] for j in pivots]
            auxiliary = solve_exact(F, transposed, target)
            assert auxiliary['status'] == 'consistent' and auxiliary['rank'] == rank
            lam = auxiliary['particular']
        else:
            lam = []
        record['inconsistency_witness'] = {
            'row': original_i, 'basis_rows': pivot_rows, 'coefficients': lam,
            'rhs_residual': F.sub(rhs[original_i], F.sum(F.mul(c, rhs[i]) for c, i in zip(lam, pivot_rows)))}
    return record


def verify_linear_certificate(F, matrix, rhs, record, verify_rank=True):
    rows, cols = len(matrix), len(matrix[0])
    assert record['equations'] == rows and record['unknowns'] == cols
    rank = record['rank']
    pivots, selected = record['pivot_columns'], record['independent_rows']
    assert len(pivots) == len(selected) == rank
    assert len(set(pivots)) == len(set(selected)) == rank
    assert all(0 <= j < cols for j in pivots) and all(0 <= i < rows for i in selected)
    free = [j for j in range(cols) if j not in pivots]
    assert record['free_columns'] == free and len(record['kernel_basis']) == len(free)

    def product_vector(x):
        assert len(x) == cols
        canonical(F, x)
        # Bounded int64 dot only when the whole accumulated sum is safe.
        if F.e == 1 and cols * (F.p - 1) ** 2 < 2 ** 63:
            import numpy as np
            return (np.asarray(matrix, dtype=np.int64) @ np.asarray(x, dtype=np.int64) % F.p).tolist()
        return [F.sum(F.mul(a, b) for a, b in zip(row, x)) for row in matrix]

    for j, vector in zip(free, record['kernel_basis']):
        assert all(vector[z] == int(z == j) for z in free)
        assert product_vector(vector) == [0] * rows
    if verify_rank and rank:
        minor = [[matrix[i][j] for j in pivots] for i in selected]
        # A fresh elimination of the exported original-row minor establishes
        # the lower bound; independent null vectors establish the upper bound.
        assert solve_exact(F, minor, [0] * rank)['rank'] == rank
    if record['status'] == 'consistent':
        assert product_vector(record['particular']) == list(rhs)
    else:
        w = record['inconsistency_witness']
        i, basis, lam = w['row'], w['basis_rows'], w['coefficients']
        assert basis == selected and len(lam) == rank
        canonical(F, lam)
        for j in range(cols):
            assert matrix[i][j] == F.sum(F.mul(c, matrix[z][j]) for c, z in zip(lam, basis))
        residual = F.sub(rhs[i], F.sum(F.mul(c, rhs[z]) for c, z in zip(lam, basis)))
        assert residual and residual == w['rhs_residual']
    return True


def relation_from_solution(F, monomials, solution, r):
    relation = [[0] for _ in range(r)] + [[1]]
    for (i, j), c in zip(monomials, solution):
        if len(relation[j]) <= i:
            relation[j].extend([0] * (i + 1 - len(relation[j])))
        relation[j][i] = c
    return [trim(a) for a in relation]


def discover(F, dom, k, s, word, u0=None, max_pieces=3, center_budget=24,
             split_budget=64, unknown_budget=512, verify_rank=True):
    """No planted labels or piece coefficients enter this interface."""
    start = time.perf_counter()
    n = len(dom)
    assert 1 <= k <= s <= n and len(set(dom)) == n and len(word) == n
    assert max_pieces >= 1 and center_budget >= 1 and split_budget >= 1 and unknown_budget >= 1
    canonical(F, dom); canonical(F, word)
    if u0 is not None:
        assert len(u0) == n
        canonical(F, u0)
    out = {'schema': 'piece-discovery/v1', 'field': field_record(F), 'domain': list(dom),
           'n': n, 'k': k, 's': s, 'word': list(word), 'u0': u0,
           'discovery_policy': 'one deterministic monic solution per degree; bounded simple-center search; no kernel-combination enumeration',
           'limits': {'max_pieces': max_pieces, 'center_budget': center_budget,
                      'split_budget_per_factor': split_budget, 'unknown_budget': unknown_budget},
           'attempts': [], 'anchor': find_code_anchor(F, dom, k, u0, word) if u0 is not None else None,
           'status': 'discovery_incomplete', 'complete_static_list_certified': False,
           'complete_scalar_fibers_certified': False,
           'polynomial_cover_size_lower_bound': 1,
           'minimum_polynomial_cover_size_certified': None}
    for r in range(1, max_pieces + 1):
        unknowns = (k - 1) * r * (r + 1) // 2 + r
        if unknowns > unknown_budget:
            out['attempts'].append({'piece_bound': r, 'status': 'linear_system_budget_exceeded', 'unknowns': unknowns})
            break
        begin = time.perf_counter()
        monomials, matrix, rhs = system(F, dom, word, k, r)
        linear = solve_exact(F, matrix, rhs)
        verify_linear_certificate(F, matrix, rhs, linear, verify_rank)
        attempt = {'piece_bound': r, 'weighted_degree': r * (k - 1),
                   'linear_certificate': linear, 'linear_seconds': time.perf_counter() - begin}
        out['attempts'].append(attempt)
        if linear['status'] == 'inconsistent':
            attempt['status'] = 'no_monic_relation_in_this_ansatz'
            if out['polynomial_cover_size_lower_bound'] == r:
                out['polynomial_cover_size_lower_bound'] = r + 1
            continue
        relation = relation_from_solution(F, monomials, linear['particular'], r)
        attempt['relation_y_coefficients'] = relation
        begin = time.perf_counter()
        factors = recover_branches(F, relation, k, center_budget, split_budget)
        attempt['factor_recovery'] = factors
        attempt['factor_seconds'] = time.perf_counter() - begin
        candidates = sorted(set(tuple(h) for h in factors['branches']))
        supports = [[i for i, x in enumerate(dom) if evaluate(F, h, x) == word[i]] for h in candidates]
        attempt['recovered_maximal_supports'] = supports
        assignments = [[] for _ in candidates]
        uncovered = []
        for i in range(n):
            j = next((j for j, supp in enumerate(supports) if i in supp), None)
            if j is None:
                uncovered.append(i)
            else:
                assignments[j].append(i)
        attempt['uncovered_coordinates'] = uncovered
        if uncovered:
            attempt['status'] = 'recovered_branches_do_not_cover_input'
            continue
        pieces = [{'coefficients': list(h), 'coordinates': indices}
                  for h, indices in zip(candidates, assignments) if indices]
        attempt['discovered_pieces'] = pieces
        attempt['identifiability'] = {
            'each_maximal_support_exceeds_weighted_degree': all(len(supp) > r * (k - 1) for supp in supports),
            'unique_monic_relation_by_rank': linear['rank'] == unknowns,
            'reason': 'If each of r distinct covering polynomials has more than r(k-1) agreement points, substitution forces every linear Y factor; monicity fixes the product.'}
        out['status'] = 'piece_cover_discovered'
        out['discovered_pieces'] = pieces
        out['minimum_polynomial_cover_size_certified'] = None
        if out['polynomial_cover_size_lower_bound'] == len(pieces):
            out['minimum_polynomial_cover_size_certified'] = len(pieces)
        try:
            static = certify_static_word(F, dom, k, s, word, pieces)
        except UnavailableCertificate as exc:
            attempt['status'] = 'cover_verified_but_root_count_list_certificate_unavailable'
            attempt['list_failure'] = str(exc)
            continue
        out['static_certificate'] = static
        out['complete_static_list_certified'] = True
        attempt['status'] = 'complete_static_list_certified'
        if u0 is not None and out['anchor'] is not None:
            certificate = certify_pencil(F, dom, k, s, u0, word, pieces)
            verify_certificate(F, certificate)
            out['scalar_certificate'] = certificate
            out['complete_scalar_fibers_certified'] = True
            out['status'] = 'complete_scalar_fibers_certified'
        else:
            out['status'] = 'complete_static_list_certified'
        break
    out['seconds'] = time.perf_counter() - start
    return out


def verify_export(F, record):
    """Read back algebra witnesses and frozen certificates without rediscovery."""
    assert record['schema'] == 'piece-discovery/v1' and record['field'] == field_record(F)
    dom, word = record['domain'], record['word']
    n, k, s = record['n'], record['k'], record['s']
    assert n == len(dom) == len(word) and 1 <= k <= s <= n and len(set(dom)) == n
    canonical(F, dom); canonical(F, word)
    if record['u0'] is not None:
        assert record['anchor'] == find_code_anchor(F, dom, k, record['u0'], word)
    else:
        assert record['anchor'] is None
    certified_pieces = []
    lower_bound = 1
    for attempt in record['attempts']:
        r = attempt['piece_bound']
        if attempt['status'] == 'linear_system_budget_exceeded':
            assert (k - 1) * r * (r + 1) // 2 + r > record['limits']['unknown_budget']
            continue
        assert attempt['weighted_degree'] == r * (k - 1)
        mons, matrix, rhs = system(F, dom, word, k, r)
        linear = attempt['linear_certificate']
        verify_linear_certificate(F, matrix, rhs, linear)
        if linear['status'] == 'inconsistent':
            assert attempt['status'] == 'no_monic_relation_in_this_ansatz'
            if lower_bound == r:
                lower_bound = r + 1
            continue
        relation = attempt['relation_y_coefficients']
        assert relation == relation_from_solution(F, mons, linear['particular'], r)
        factors = attempt['factor_recovery']
        branches = factors['branches']
        assert len({tuple(h) for h in branches}) == len(branches)
        factor_product = [[1]]
        for h in branches:
            assert len(h) == k
            canonical(F, h)
            assert substitute_y(F, relation, h) == [0]
            new = [[0] for _ in range(len(factor_product) + 1)]
            for j, a in enumerate(factor_product):
                new[j] = add(F, new[j], scale(F, mul(F, a, h), F.sub(0, 1)))
                new[j + 1] = add(F, new[j + 1], a)
            factor_product = new
        if factors.get('factor_product_verified'):
            assert [trim(a) for a in factor_product] == [trim(a) for a in relation]
            assert factors['status'] == 'complete_polynomial_factorization'
        candidates = sorted(set(tuple(h) for h in branches))
        supports = [[i for i, x in enumerate(dom) if evaluate(F, h, x) == word[i]] for h in candidates]
        assert attempt['recovered_maximal_supports'] == supports
        uncovered = [i for i in range(n) if not any(i in supp for supp in supports)]
        assert attempt['uncovered_coordinates'] == uncovered
        if not uncovered:
            assignments = [[] for _ in candidates]
            for i in range(n):
                assignments[next(j for j, supp in enumerate(supports) if i in supp)].append(i)
            pieces = [{'coefficients': list(h), 'coordinates': indices}
                      for h, indices in zip(candidates, assignments) if indices]
            assert attempt['discovered_pieces'] == pieces
            certified_pieces.append(pieces)
            ident = attempt['identifiability']
            assert ident['each_maximal_support_exceeds_weighted_degree'] == all(len(x) > r * (k - 1) for x in supports)
            assert ident['unique_monic_relation_by_rank'] == (linear['rank'] == len(mons))
    assert record['polynomial_cover_size_lower_bound'] == lower_bound
    if record['minimum_polynomial_cover_size_certified'] is not None:
        pieces = record['discovered_pieces']
        assert pieces in certified_pieces
        assert record['minimum_polynomial_cover_size_certified'] == len(pieces) == lower_bound
    if record['complete_static_list_certified']:
        pieces = record['discovered_pieces']
        assert pieces in certified_pieces
        assert record['static_certificate'] == certify_static_word(F, dom, k, s, word, pieces)
    if record['complete_scalar_fibers_certified']:
        assert record['complete_static_list_certified'] and record['anchor'] is not None
        cert = record['scalar_certificate']
        assert cert['domain'] == dom and cert['u0'] == record['u0'] and cert['u1'] == word
        assert (cert['n'], cert['k'], cert['s']) == (n, k, s)
        verify_certificate(F, cert)
        assert cert['static_direction_certificate'] == record['static_certificate']
    return True
