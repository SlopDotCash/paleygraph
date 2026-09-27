#!/usr/bin/env python3
"""Ablate the new mean matrix into old target incidence and induced boundaries."""
from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from decomposition import components,project_edges
from review import inward
from query import query

HERE=Path(__file__).resolve().parent


def rat(v):v=F(v);return [v.numerator,v.denominator]


def analyze(decomp):
    n=len(decomp['selected']);N=decomp['ordered_insertion_pairs_per_deletion']
    fields={}
    for name in ('target_sum','high_degree_part','lower_degree_part','induced_boundary_part','pair_target_high','pair_target_low','pair_target_lower'):
        fields[name],den=project_edges(n,{tuple(r['indices']):r[name] for r in decomp['pairs']})
    H,L,J=(decomp[k] for k in ('high_pair_coefficient','low_pair_coefficient','lower_pair_coefficient'))
    for edge in fields['target_sum']:
        assert fields['high_degree_part'][edge]==H*fields['pair_target_high'][edge]
        assert fields['lower_degree_part'][edge]==L*fields['pair_target_low'][edge]+J*fields['pair_target_lower'][edge]
        assert fields['target_sum'][edge]==sum(fields[name][edge] for name in ('high_degree_part','lower_degree_part','induced_boundary_part'))
    norm=lambda name:sum(v*v for v in fields[name].values())
    total=norm('target_sum');high=norm('high_degree_part');low=norm('lower_degree_part');boundary=norm('induced_boundary_part')
    error=sum((fields['target_sum'][e]-fields['high_degree_part'][e])**2 for e in fields['target_sum'])
    return {'projected_mean_norm_squared':rat(F(total,(N*den)**2)),
            'high_target_norm_squared':rat(F(high,(N*den)**2)),
            'lower_degree_norm_squared':rat(F(low,(N*den)**2)),
            'boundary_norm_squared':rat(F(boundary,(N*den)**2)),
            'remainder_relative_to_high_norm_squared':rat(F(error,high)) if high else None,
            'high_target_projection_identically_zero':not bool(high),
            'all_projected_means_explained_by_lower_degrees_and_induced_graph':not bool(high),
            'projected_T2_pair_targets_identically_zero':not bool(norm('pair_target_lower'))}


def toy():
    p=HERE.parent/'round7/local_edits/orbit_results.json';old=json.loads(p.read_text());cases=[]
    for case in old['cases']:
        n=case['n'];decks=defaultdict(list);norms=defaultdict(list);rows=[]
        for c in case['canonical_orbit_representatives']:
            t=inward(17,c,6);decomp=components(17,c,6,t['induced'],t);r=query(17,c,6)
            assert [x['target_sum'] for x in r['pairs']]==[x['target_sum'] for x in decomp['pairs']]
            a=analyze(decomp);deck=tuple(sorted(F(*x['mean']) for x in r['pairs']))
            decks[deck].append(c);norms[tuple(a['projected_mean_norm_squared'])].append(c)
            rows.append({'selected':c,'analysis':a})
        def canonical(c):return min(tuple(sorted(a*(x-pivot)%17 for x in c)) for a in range(1,17) for pivot in c)
        full_decks=[sorted({canonical(c) for c in group}) for group in decks.values()]
        full_norms=[sorted({canonical(c) for c in group}) for group in norms.values()]
        cases.append({'q':17,'n':n,'square_affine_representatives':len(rows),'mean_deck_fibres':len(decks),
                      'mean_deck_fibres_joining_distinct_square_affine_orbits':sum(len(v)>1 for v in decks.values()),
                      'mean_deck_max_orbits_per_fibre':max(map(len,decks.values())),
                      'full_affine_orbits':len({canonical(r['selected']) for r in rows}),
                      'mean_deck_fibres_joining_distinct_full_affine_orbits':sum(len(v)>1 for v in full_decks),
                      'mean_deck_full_affine_orbit_collisions':[v for v in full_decks if len(v)>1],
                      'projected_norm_fibres_joining_distinct_full_affine_orbits':sum(len(v)>1 for v in full_norms),
                      'projected_norm_fibres':len(norms),'all_projected_means_explained_by_lower_degrees_and_induced_graph':all(r['analysis']['all_projected_means_explained_by_lower_degrees_and_induced_graph'] for r in rows),
                      'rows':rows})
        print(json.dumps({k:v for k,v in cases[-1].items() if k!='rows'}),flush=True)
    result={'status':'passed','scope':'Complete prior square-affine representative lists for Paley17. The n6/n7 projected separation is exactly lower-degree incidence plus induced boundary, not new arithmetic.',
            'cases':cases,'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ('projection_analysis.py','decomposition.py','review.py','query.py','inward_targets.cpp','inward_targets','pair_means_backend.cpp','pair_means_backend')},
            'input_sha256':{'../round7/local_edits/orbit_results.json':sha256(p.read_bytes()).hexdigest()}}
    (HERE/'projection_toy.json').write_text(json.dumps(result,indent=2)+'\n')


def scale():
    p=HERE/'scale_verification.json';old=json.loads(p.read_text());rows=[]
    for r in old['cases']:
        row={'q':r['q'],'n':r['n'],'family':r['family'],**analyze(r['decomposition'])};rows.append(row)
        print(json.dumps(row),flush=True)
    result={'status':'passed','cases':rows,'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ('projection_analysis.py','decomposition.py')},
            'input_sha256':{'scale_verification.json':sha256(p.read_bytes()).hexdigest()}}
    (HERE/'projection_scale.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    import sys
    if sys.argv[1:]==['scale']:scale()
    else:toy()
