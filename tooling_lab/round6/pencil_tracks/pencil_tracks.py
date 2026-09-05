#!/usr/bin/env python3
"""Certified affine polynomial-track lists without a codeword anchor.

Known interpolation/root-counting and affine-incidence foundations. Discovery
is bounded and incomplete; verified track partitions give exact whole-field lists.
"""
from pathlib import Path
import sys
import time
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
LAB = HERE.parents[1]
sys.path.insert(0, str(LAB / 'round5/piece_discovery'))
import piece_discovery as pd
from scalar_fibers import (PrimeField, Field, evaluate, canonical, coefficients,
                          field_record, find_code_anchor, UnavailableCertificate)


def affine_zero(F, intercept, slope):
    """Exact zero set of a vector-valued affine function."""
    assert len(intercept) == len(slope)
    pivot = next((i for i, b in enumerate(slope) if b), None)
    if pivot is None:
        return {'kind': 'all'} if not any(intercept) else {'kind': 'empty'}
    scalar = F.mul(F.sub(0, intercept[pivot]), F.inv(slope[pivot]))
    if all(F.add(a, F.mul(scalar, b)) == 0 for a, b in zip(intercept, slope)):
        return {'kind': 'singleton', 'scalar': scalar}
    return {'kind': 'empty'}


def normalize_tracks(F, dom, k, u0, u1, tracks):
    n = len(dom)
    assert 1 <= k <= n and len(u0) == len(u1) == n and len(set(dom)) == n
    for xs in [dom, u0, u1]:
        canonical(F, xs)
    flat = []
    grouped = {}
    for track in tracks:
        a = coefficients(F, track['intercept_coefficients'], k)
        b = coefficients(F, track['slope_coefficients'], k)
        indices = track['coordinates']
        assert all(isinstance(i, int) and 0 <= i < n for i in indices)
        flat.extend(indices)
        for i in indices:
            assert evaluate(F, a, dom[i]) == u0[i]
            assert evaluate(F, b, dom[i]) == u1[i]
        if indices:
            grouped.setdefault((a, b), []).extend(indices)
    assert sorted(flat) == list(range(n))
    return [{'intercept_coefficients': list(a), 'slope_coefficients': list(b),
             'coordinates': sorted(indices)} for (a, b), indices in sorted(grouped.items())]


def support_at(track, scalar):
    return sorted(track['always_coordinates'] + track['scalar_buckets'].get(str(scalar), []))


def track_coefficients(F, track, scalar):
    return [F.add(a, F.mul(scalar, b)) for a, b in
            zip(track['intercept_coefficients'], track['slope_coefficients'])]


def compile_certificate(F, dom, k, s, u0, u1, tracks):
    """Completeness requires a supplied, fully verified polynomial-track partition."""
    start = time.perf_counter()
    n = len(dom)
    assert 1 <= k <= s <= n
    normalized = normalize_tracks(F, dom, k, u0, u1, tracks)
    cap = sum(min(len(t['coordinates']), k - 1) for t in normalized)
    if s <= cap:
        raise UnavailableCertificate(f'joint track root-count certificate needs s>{cap}, got s={s}')
    events = set()
    compiled = []
    for index, track in enumerate(normalized):
        a, b = track['intercept_coefficients'], track['slope_coefficients']
        av = [evaluate(F, a, x) for x in dom]
        bv = [evaluate(F, b, x) for x in dom]
        always, never, buckets = [], [], {}
        for i in range(n):
            zero = affine_zero(F, [F.sub(av[i], u0[i])], [F.sub(bv[i], u1[i])])
            if zero['kind'] == 'all':
                always.append(i)
            elif zero['kind'] == 'empty':
                never.append(i)
            else:
                z = zero['scalar']
                buckets.setdefault(str(z), []).append(i)
                events.add(z)
        compiled.append({'id': index, **track, 'intercept_codeword': av, 'slope_codeword': bv,
                         'always_coordinates': always, 'never_coordinates': never,
                         'scalar_buckets': {z: buckets[z] for z in sorted(buckets, key=int)}})
    collisions = []
    for i, left in enumerate(compiled):
        for j in range(i + 1, len(compiled)):
            right = compiled[j]
            da = [F.sub(a, b) for a, b in zip(left['intercept_coefficients'], right['intercept_coefficients'])]
            db = [F.sub(a, b) for a, b in zip(left['slope_coefficients'], right['slope_coefficients'])]
            zero = affine_zero(F, da, db)
            assert zero['kind'] != 'all'  # Coincident affine tracks were merged.
            if zero['kind'] == 'singleton':
                z = zero['scalar']; events.add(z)
                collisions.append({'tracks': [i, j], 'scalar': z})
    generic = [t['id'] for t in compiled if len(t['always_coordinates']) >= s]
    finite = []
    for z in sorted(events):
        groups = {}
        for t in compiled:
            if len(support_at(t, z)) >= s:
                key = tuple(track_coefficients(F, t, z))
                groups.setdefault(key, []).append(t['id'])
        nodes = []
        for cs, ids in sorted(groups.items()):
            support = support_at(compiled[ids[0]], z)
            assert all(support_at(compiled[j], z) == support for j in ids)
            nodes.append({'representative_track': ids[0], 'equal_track_ids': ids})
        finite.append({'scalar': z, 'nodes': nodes, 'list_size': len(nodes)})
    count = (F.q - len(events)) * len(generic) + sum(x['list_size'] for x in finite)
    bad = (F.q - len(events) if generic else 0) + sum(bool(x['nodes']) for x in finite)
    return {'schema': 'polynomial-pencil-tracks/v1', 'status': 'certified_complete',
            'scope': 'all ordinary >=s agreement codewords and maximal supports; MCA nonjointness is not imposed',
            'field': field_record(F), 'n': n, 'k': k, 's': s, 'domain': list(dom),
            'u0': list(u0), 'u1': list(u1), 'tracks': compiled,
            'root_count_certificate': {'cap': cap, 'threshold': s, 'strict': True,
                'reason': 'At each scalar a polynomial distinct from every track polynomial has at most k-1 matches on each verified region.'},
            'track_collisions': collisions,
            'generic_scalar_domain': {'kind': 'all_except', 'excluded': sorted(events), 'cardinality': F.q - len(events)},
            'generic_qualifying_track_ids': generic, 'generic_list_size': len(generic),
            'exceptional_scalars': finite, 'complete_symbolic_node_count': count,
            'bad_scalar_count': bad, 'every_scalar_has_a_qualifying_codeword': bad == F.q,
            'node_encoding': 'At z, representative j means coefficients a_j+z*b_j, codeword on the saved domain, and support always_j union bucket_j[z]. equal_track_ids are precisely qualifying tracks giving that same polynomial.',
            'code_anchor': find_code_anchor(F, dom, k, u0, u1),
            'seconds': time.perf_counter() - start}


def materialize_scalar(F, certificate, z):
    canonical(F, [z])
    event = next((item for item in certificate['exceptional_scalars'] if item['scalar'] == z), None)
    references = event['nodes'] if event is not None else [
        {'representative_track': j, 'equal_track_ids': [j]} for j in certificate['generic_qualifying_track_ids']]
    nodes = []
    for ref in references:
        t = certificate['tracks'][ref['representative_track']]
        cs = track_coefficients(F, t, z)
        nodes.append({'scalar': z, 'coefficients': cs,
                      'codeword': [F.add(a, F.mul(z, b)) for a, b in zip(t['intercept_codeword'], t['slope_codeword'])],
                      'agreement_support': support_at(t, z), 'equal_track_ids': ref['equal_track_ids']})
    return nodes


def verify_certificate(F, certificate):
    assert certificate['field'] == field_record(F)
    assert certificate['n'] == len(certificate['domain'])
    rebuilt = compile_certificate(F, certificate['domain'], certificate['k'], certificate['s'],
                                  certificate['u0'], certificate['u1'], certificate['tracks'])
    for key in rebuilt:
        if key != 'seconds':
            assert rebuilt[key] == certificate[key], key
    return True


def intersect_covers(F, dom, k, u0, u1, left, right, slice_scalars=None):
    """Intersect two certified point partitions, then convert to affine tracks."""
    tracks = []
    if slice_scalars is not None:
        z0, z1 = slice_scalars
        canonical(F, [z0, z1]); assert z0 != z1
        inverse = F.inv(F.sub(z1, z0))
    for a in left:
        for b in right:
            indices = sorted(set(a['coordinates']).intersection(b['coordinates']))
            if not indices:
                continue
            ca = coefficients(F, a['coefficients'], k)
            cb = coefficients(F, b['coefficients'], k)
            if slice_scalars is None:
                intercept, slope = list(ca), list(cb)
            else:
                slope = [F.mul(F.sub(v, u), inverse) for u, v in zip(ca, cb)]
                intercept = [F.sub(u, F.mul(z0, v)) for u, v in zip(ca, slope)]
            tracks.append({'intercept_coefficients': intercept, 'slope_coefficients': slope, 'coordinates': indices})
    return normalize_tracks(F, dom, k, u0, u1, tracks)


def greedy_reassign(F, dom, k, u0, u1, tracks):
    """Known greedy set cover on verified joint supports; no minimality claim."""
    tracks = normalize_tracks(F, dom, k, u0, u1, tracks)
    full = []
    for t in tracks:
        full.append({i for i, x in enumerate(dom)
                     if evaluate(F, t['intercept_coefficients'], x) == u0[i]
                     and evaluate(F, t['slope_coefficients'], x) == u1[i]})
    remaining = set(range(len(dom)))
    out = []
    while remaining:
        j = min(range(len(tracks)), key=lambda j: (-len(full[j] & remaining), j))
        indices = sorted(full[j] & remaining)
        assert indices
        out.append({**tracks[j], 'coordinates': indices})
        remaining.difference_update(indices)
    return normalize_tracks(F, dom, k, u0, u1, out)


def discover_pencil(F, dom, k, s, u0, u1, max_pieces=3, slice_scalars=None, reassign=True):
    start = time.perf_counter()
    if slice_scalars is None:
        words = [u0, u1]
    else:
        z0, z1 = slice_scalars
        canonical(F, [z0, z1]); assert z0 != z1
        words = [[F.add(a, F.mul(z, b)) for a, b in zip(u0, u1)] for z in [z0, z1]]
    discoveries = [pd.discover(F, dom, k, s, w, max_pieces=max_pieces) for w in words]
    result = {'schema': 'blind-pencil-discovery/v1', 'field': field_record(F), 'n': len(dom),
              'k': k, 's': s, 'domain': list(dom), 'u0': list(u0), 'u1': list(u1),
              'mode': 'intercept_and_direction' if slice_scalars is None else 'two_scalar_slices',
              'slice_scalars': list(slice_scalars) if slice_scalars is not None else None,
              'marginal_discoveries': discoveries, 'status': 'discovery_incomplete',
              'complete': False, 'code_anchor': find_code_anchor(F, dom, k, u0, u1),
              'reassignment_policy': 'greedy verified joint-support cover' if reassign else 'raw partition intersections'}
    if any('discovered_pieces' not in d for d in discoveries):
        result['failure'] = 'at least one marginal polynomial cover was not discovered'
        result['seconds'] = time.perf_counter() - start
        return result
    tracks = intersect_covers(F, dom, k, u0, u1, discoveries[0]['discovered_pieces'],
                              discoveries[1]['discovered_pieces'], slice_scalars)
    result['raw_intersection_tracks'] = tracks
    if reassign:
        tracks = greedy_reassign(F, dom, k, u0, u1, tracks)
    result['verified_tracks'] = tracks
    result['joint_root_count_cap'] = sum(min(len(t['coordinates']), k - 1) for t in tracks)
    try:
        certificate = compile_certificate(F, dom, k, s, u0, u1, tracks)
    except UnavailableCertificate as exc:
        result['status'] = 'joint_cover_found_certificate_unavailable'
        result['failure'] = str(exc)
    else:
        verify_certificate(F, certificate)
        result['complete'] = True
        result['status'] = 'certified_complete'
        result['certificate'] = certificate
    result['seconds'] = time.perf_counter() - start
    return result


def verify_discovery(F, result, verify_marginals=True):
    assert result['field'] == field_record(F)
    assert result['reassignment_policy'] in ['greedy verified joint-support cover', 'raw partition intersections']
    dom, k, s, u0, u1 = (result[x] for x in ['domain', 'k', 's', 'u0', 'u1'])
    assert result['n'] == len(dom) and 1 <= k <= s <= len(dom)
    if result['slice_scalars'] is None:
        assert result['mode'] == 'intercept_and_direction'
        words = [u0, u1]
    else:
        assert result['mode'] == 'two_scalar_slices'
        zs = result['slice_scalars']
        assert len(zs) == 2 and zs[0] != zs[1]; canonical(F, zs)
        words = [[F.add(a, F.mul(z, b)) for a, b in zip(u0, u1)] for z in zs]
    discoveries = result['marginal_discoveries']
    assert len(discoveries) == 2
    for d, w in zip(discoveries, words):
        assert d['field'] == field_record(F) and d['domain'] == dom and d['word'] == w and d['u0'] is None
        assert (d['n'], d['k'], d['s']) == (len(dom), k, s)
        if verify_marginals:
            pd.verify_export(F, d)
    assert result['code_anchor'] == find_code_anchor(F, dom, k, u0, u1)
    if 'raw_intersection_tracks' in result:
        raw = intersect_covers(F, dom, k, u0, u1, discoveries[0]['discovered_pieces'],
                               discoveries[1]['discovered_pieces'], result['slice_scalars'])
        assert result['raw_intersection_tracks'] == raw
        tracks = greedy_reassign(F, dom, k, u0, u1, raw) if result['reassignment_policy'] == 'greedy verified joint-support cover' else raw
        assert result['verified_tracks'] == tracks
        assert result['joint_root_count_cap'] == sum(min(len(t['coordinates']), k - 1) for t in tracks)
    if result['complete']:
        assert result['status'] == 'certified_complete'
        cert = result['certificate']
        assert (cert['domain'], cert['k'], cert['s'], cert['u0'], cert['u1']) == (dom, k, s, u0, u1)
        assert normalize_tracks(F, dom, k, u0, u1, cert['tracks']) == result['verified_tracks']
        verify_certificate(F, cert)
    elif result['status'] == 'joint_cover_found_certificate_unavailable':
        assert result['joint_root_count_cap'] >= s
    else:
        assert any('discovered_pieces' not in d for d in discoveries)
    return True
