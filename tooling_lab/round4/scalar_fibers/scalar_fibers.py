#!/usr/bin/env python3
"""Exact through-code scalar transport with a verified piecewise direction list.

Known linear-code symmetry and polynomial root counting; no basis-cover search.
The piece polynomials are supplied witnesses, not discovered by this tool.
"""
from pathlib import Path
import sys
import time

sys.dont_write_bytecode = True
LAB = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(LAB / 'round2' / 'proximity'))
from compressed_support import PrimeField, Field


class UnavailableCertificate(ValueError):
    pass


def field_record(F):
    return {'p': F.p, 'e': F.e, 'q': F.q, 'modulus': F.f if F.e > 1 else None}


def canonical(F, values):
    assert all(isinstance(x, int) and 0 <= x < F.q for x in values)


def coefficients(F, cs, k):
    assert len(cs) <= k
    canonical(F, cs)
    return tuple(cs) + (0,) * (k - len(cs))


def evaluate(F, cs, x):
    value = 0
    for c in reversed(cs):
        value = F.add(F.mul(value, x), c)
    return value


def interpolate(F, xs, ys):
    """Newton divided differences, then conversion to monomial coefficients."""
    assert len(xs) == len(ys) >= 1 and len(set(xs)) == len(xs)
    dd = list(ys)
    for j in range(1, len(xs)):
        for i in range(len(xs) - 1, j - 1, -1):
            dd[i] = F.mul(F.sub(dd[i], dd[i - 1]), F.inv(F.sub(xs[i], xs[i - j])))
    out = [dd[-1]]
    for i in range(len(xs) - 2, -1, -1):
        next_out = [0] * (len(out) + 1)
        for j, c in enumerate(out):
            next_out[j] = F.sub(next_out[j], F.mul(xs[i], c))
            next_out[j + 1] = F.add(next_out[j + 1], c)
        next_out[0] = F.add(next_out[0], dd[i])
        out = next_out
    return out


def find_code_anchor(F, dom, k, u0, u1):
    """Solve the two residual vectors' common affine zero, without scalar scan."""
    n = len(dom)
    assert 1 <= k <= n and len(set(dom)) == n and len(u0) == len(u1) == n
    for xs in [dom, u0, u1]:
        canonical(F, xs)
    c0 = interpolate(F, dom[:k], u0[:k])
    c1 = interpolate(F, dom[:k], u1[:k])
    r0 = [F.sub(y, evaluate(F, c0, x)) for x, y in zip(dom, u0)]
    r1 = [F.sub(y, evaluate(F, c1, x)) for x, y in zip(dom, u1)]
    pivot = next((i for i, value in enumerate(r1) if value), None)
    if pivot is None:
        if any(r0):
            return None
        scalar, mode = 0, 'entire_pencil_in_code'
    else:
        scalar = F.mul(F.sub(0, r0[pivot]), F.inv(r1[pivot]))
        if any(F.add(a, F.mul(scalar, b)) for a, b in zip(r0, r1)):
            return None
        mode = 'unique_code_intersection'
    cs = [F.add(a, F.mul(scalar, b)) for a, b in zip(c0, c1)]
    return {'scalar': scalar, 'coefficients': cs, 'discovery': mode,
            'method': 'two interpolants and a common affine residual zero'}


def certify_static_word(F, dom, k, s, word, pieces):
    """Verify a supplied piece partition and certify the entire static list."""
    n = len(dom)
    assert 1 <= k <= s <= n and len(set(dom)) == n and len(word) == n
    canonical(F, dom); canonical(F, word)
    flat, grouped = [], {}
    for piece in pieces:
        cs = coefficients(F, piece['coefficients'], k)
        indices = piece['coordinates']
        assert all(isinstance(i, int) and 0 <= i < n for i in indices)
        flat.extend(indices)
        for i in indices:
            assert evaluate(F, cs, dom[i]) == word[i]
        if indices:
            grouped.setdefault(cs, []).extend(indices)
    assert len(flat) == n and sorted(flat) == list(range(n))
    normalized = [{'coefficients': list(cs), 'coordinates': sorted(indices)}
                  for cs, indices in sorted(grouped.items())]
    # Any other degree<k polynomial has at most k-1 roots of its difference
    # with each distinct piece polynomial, restricted to that piece's region.
    cap = sum(min(len(piece['coordinates']), k - 1) for piece in normalized)
    if s <= cap:
        raise UnavailableCertificate(f'piece root-count certificate needs s>{cap}, got s={s}')
    candidates = []
    for piece in normalized:
        cs = piece['coefficients']
        cw = [evaluate(F, cs, x) for x in dom]
        support = [i for i in range(n) if cw[i] == word[i]]
        candidates.append({'coefficients': cs, 'codeword': cw, 'agreement_support': support,
                           'qualifies': len(support) >= s})
    qualified = [x for x in candidates if x['qualifies']]
    return {'complete': True, 'kind': 'disjoint_piece_root_count',
            'piece_discovery': 'supplied polynomial/partition witness; only verification is performed',
            'pieces': normalized, 'distinct_piece_polynomials': len(normalized),
            'outside_piece_agreement_cap': cap, 'threshold': s,
            'piece_candidates': candidates, 'complete_static_list': qualified,
            'complete_static_list_size': len(qualified)}


def certify_pencil(F, dom, k, s, u0, u1, pieces, anchor=None):
    start = time.perf_counter()
    if s < k:
        raise UnavailableCertificate('s>=k is required; the anchor can otherwise have additional codewords')
    n = len(dom)
    assert 1 <= k <= s <= n and len(set(dom)) == n and len(u0) == len(u1) == n
    for xs in [dom, u0, u1]:
        canonical(F, xs)
    if anchor is None:
        anchor = find_code_anchor(F, dom, k, u0, u1)
    if anchor is None:
        raise UnavailableCertificate('the pencil has no complete codeword intersection')
    zstar = anchor['scalar']; canonical(F, [zstar])
    fcs = coefficients(F, anchor['coefficients'], k)
    f = [evaluate(F, fcs, x) for x in dom]
    assert f == [F.add(a, F.mul(zstar, b)) for a, b in zip(u0, u1)]
    static = certify_static_word(F, dom, k, s, u1, pieces)
    fibers = []
    for h in static['complete_static_list']:
        hcs, cw = h['coefficients'], h['codeword']
        fibers.append({'scalar_domain': {'kind': 'all_except', 'excluded': [zstar]},
            'direction_coefficients': hcs, 'direction_codeword': cw,
            'intercept_coefficients': [F.sub(a, F.mul(zstar, b)) for a, b in zip(fcs, hcs)],
            'intercept_codeword': [F.sub(a, F.mul(zstar, b)) for a, b in zip(f, cw)],
            'slope_coefficients': hcs, 'slope_codeword': cw,
            'agreement_support_for_every_scalar_in_domain': h['agreement_support']})
    L = len(fibers)
    return {'schema': 'exact-scalar-fiber-certificate/v1', 'status': 'certified_complete',
            'event_definition': 'exists a codeword with at least s matches; MCA non-jointness is not imposed',
            'field': field_record(F), 'n': n, 'k': k, 's': s, 'domain': dom, 'u0': u0, 'u1': u1,
            'anchor': {'scalar': zstar, 'coefficients': list(fcs), 'codeword': f,
                       'agreement_support': list(range(n)),
                       'discovery': anchor.get('discovery', 'supplied witness')},
            'static_direction_certificate': static, 'nonanchor_fibers': fibers,
            'transport': 'g(z)=f+(z-zstar)*h; inverse h=(g-f)/(z-zstar)',
            'list_size_at_anchor': 1, 'list_size_at_every_nonanchor_scalar': L,
            'bad_scalar_count': F.q if L else 1,
            'bad_scalar_set': {'kind': 'whole_field'} if L else {'kind': 'singleton', 'scalar': zstar},
            'complete_symbolic_node_count': 1 + (F.q - 1) * L,
            'elapsed_seconds': round(time.perf_counter() - start, 6)}


def materialize_scalar(F, record, scalar):
    canonical(F, [scalar])
    anchor = record['anchor']
    if scalar == anchor['scalar']:
        return [{'scalar': scalar, 'codeword': anchor['codeword'],
                 'agreement_support': anchor['agreement_support']}]
    out = []
    for fiber in record['nonanchor_fibers']:
        cw = [F.add(a, F.mul(scalar, b)) for a, b in zip(fiber['intercept_codeword'], fiber['slope_codeword'])]
        out.append({'scalar': scalar, 'codeword': cw,
                    'agreement_support': fiber['agreement_support_for_every_scalar_in_domain']})
    return out


def verify_certificate(F, record):
    """Rebuild from explicit polynomial/partition evidence, never enumerate bases."""
    assert record['field'] == field_record(F)
    rebuilt = certify_pencil(F, record['domain'], record['k'], record['s'],
        record['u0'], record['u1'], record['static_direction_certificate']['pieces'], record['anchor'])
    assert record['n'] == len(record['domain'])
    for key in rebuilt:
        if key != 'elapsed_seconds':
            assert rebuilt[key] == record[key], key
    return {'all_symbolic_fibers_verified': True,
            'complete_static_list_size': record['list_size_at_every_nonanchor_scalar'],
            'symbolic_node_count': record['complete_symbolic_node_count']}
