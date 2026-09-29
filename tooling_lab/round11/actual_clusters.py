#!/usr/bin/env python3
"""Apply the rank-two compiler to actual saved polynomial candidates."""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from rank_two_pinning import compile_pinning

HERE=Path(__file__).resolve().parent


def evaluate(f,x,p):return sum(v*pow(x,j,p) for j,v in enumerate(f))%p


def difference(a,b,p):
    return [((a[j] if j<len(a) else 0)-(b[j] if j<len(b) else 0))%p for j in range(max(len(a),len(b)))]


def matrix(f,g,domain,p):return [(evaluate(f,x,p),evaluate(g,x,p)) for x in domain]


def main():
    smallpath=HERE.parent/'round6/overlap_preflight/coset_small_characteristic_sample3.certificate.json'
    source=json.loads(smallpath.read_text());p,n,k,s=(source[x] for x in ('p','n','k','s'));nodes=source['nodes'];candidates=[];degenerate=0
    for node in nodes:
        word=[evaluate(node['coefficients'],x,p) for x in source['domain']]
        assert word==node['codeword'] and len(node['coefficients'])<=k
        assert [i for i,(a,b,v) in enumerate(zip(source['u0'],source['u1'],word)) if (a+node['scalar']*b)%p==v]==node['agreement_support']
    for indices in combinations(range(len(nodes)),3):
        a,b,c=(nodes[i]['coefficients'] for i in indices)
        f=difference(b,a,p);g=difference(c,a,p);columns=matrix(f,g,source['domain'],p)
        try:r=compile_pinning(p,columns,s)
        except ValueError:degenerate+=1;continue
        candidates.append((r['minimum_injective_pairs'],indices,f,g,r))
    minimum,indices,f,g,result=min(candidates,key=lambda a:(a[0],a[1]))
    columns=matrix(f,g,source['domain'],p)
    edges=[(i,j) for i,j in combinations(range(n),2) if (columns[i][0]*columns[j][1]-columns[j][0]*columns[i][1])%p]
    brute_min=min(sum(i in A and j in A for i,j in edges) for A in map(set,combinations(range(n),s)))
    assert brute_min==minimum
    small={'p':p,'n':n,'k':k,'s':s,'candidate_triples':len(candidates)+degenerate,'rank_two_triples':len(candidates),
           'dependent_triples':degenerate,'worst_triple_node_indices':indices,'directions':[f,g],
           'result':result,'complete_agreement_sets_checked':4368}
    largepath=HERE.parent/'round6/pencil_tracks/results.json';old=json.loads(largepath.read_text())
    case=next(c['result'] for c in old['large'] if c['name']=='skew-two')
    p=case['field']['p'];n,k,s=(case[x] for x in ('n','k','s'));tracks=case['verified_tracks']
    origin=tracks[0]['intercept_coefficients'];f=tracks[0]['slope_coefficients'];g=difference(tracks[1]['intercept_coefficients'],origin,p)
    assert all(len(poly)<=k for poly in (origin,f,g))
    agreement_counts=[]
    # Three actual near-codeword points determine the declared affine plane.
    summed=difference(origin,[-x%p for x in f],p)
    polynomials=[origin,summed,tracks[1]['intercept_coefficients']]
    for scalar,poly in zip((0,1,0),polynomials):
        count=sum(evaluate(poly,x,p)==(u+scalar*v)%p for x,u,v in zip(case['domain'],case['u0'],case['u1']))
        assert count>=s;agreement_counts.append(count)
    result=compile_pinning(p,matrix(f,g,case['domain'],p),s)
    A=result['worst_agreement_set'];cols=matrix(f,g,case['domain'],p)
    actual=sum((cols[i][0]*cols[j][1]-cols[j][0]*cols[i][1])%p!=0 for i,j in combinations(A,2))
    assert actual==result['minimum_injective_pairs']
    large={'p':p,'n':n,'k':k,'s':s,'origin':origin,'directions':[f,g],'source_node_scalars':[0,1,0],
           'source_node_agreement_counts':agreement_counts,'result':result,'attaining_witness_pairs_checked':len(A)*(len(A)-1)//2,
           'universal_RS_rank_two_lower_bound':(s-k+2)*(s-k+1)//2}
    out={'status':'initial_actual_inputs_passed','scope':'Declared affine clusters from saved actual polynomial candidates; no discovery of all required clusters. Large minimum uses the proved rank-two occupancy formula and a fully checked attaining witness, not enumeration of all agreement sets.',
         'small':small,'large':large,'source_sha256':{name:sha256((HERE/name).read_bytes()).hexdigest() for name in ('actual_clusters.py','rank_two_pinning.py')},
         'input_sha256':{'../round6/overlap_preflight/coset_small_characteristic_sample3.certificate.json':sha256(smallpath.read_bytes()).hexdigest(),
                         '../round6/pencil_tracks/results.json':sha256(largepath.read_bytes()).hexdigest()}}
    (HERE/'actual_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'small_rank_two_triples':small['rank_two_triples'],'small_minimum_pairs':minimum,
                      'large_minimum_pairs':result['minimum_injective_pairs'],'large_success_probability':result['uniform_pair_success_probability'],
                      'large_universal_bound':large['universal_RS_rank_two_lower_bound']}),flush=True)


if __name__=='__main__':main()
