#!/usr/bin/env python3
"""Audit the simple-space premise on actual saved near-codeword candidates."""
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from rank_three_pinning import certify
from rank_three_verifier import rank,validate
from run_rank_three import literal

HERE=Path(__file__).resolve().parent


def difference(a,b,p):
    return [((a[j] if j<len(a) else 0)-(b[j] if j<len(b) else 0))%p for j in range(max(len(a),len(b)))]


def evaluate(b,x,p):return sum(a*pow(x,j,p) for j,a in enumerate(b))%p


def simple_profile(columns,p):
    zeros=[];directions={}
    for i,v in enumerate(columns):
        first=next((a for a in v if a),None)
        if first is None:zeros.append(i);continue
        inv=pow(first,-1,p);key=tuple(a*inv%p for a in v);directions.setdefault(key,[]).append(i)
    return {'rank':rank(columns,p),'zero_columns':len(zeros),'projective_classes':len(directions),
            'parallel_class_sizes':sorted(map(len,directions.values())),
            'simple':not zeros and len(directions)==len(columns)}


def main():
    sourcepath=HERE.parent/'round6/overlap_preflight/coset_small_characteristic_sample3.certificate.json'
    source=json.loads(sourcepath.read_text());nodes=source['nodes'];p=source['p'];rows=[];counts={'dependent':0,'simple_rank_three':0,'nonsimple_rank_three':0}
    for node in nodes:
        assert len(node['coefficients'])<=source['k']
        word=[evaluate(node['coefficients'],x,p) for x in source['domain']];assert word==node['codeword']
        support=[j for j,(a,b,v) in enumerate(zip(source['u0'],source['u1'],word)) if (a+node['scalar']*b)%p==v]
        assert support==node['agreement_support'] and len(support)>=source['s']
    for ids in combinations(range(len(nodes)),4):
        words=[nodes[i]['codeword'] for i in ids]
        columns=[tuple((word[j]-words[0][j])%p for word in words[1:]) for j in range(source['n'])]
        profile=simple_profile(columns,p)
        if profile['rank']<3:counts['dependent']+=1;continue
        if not profile['simple']:counts['nonsimple_rank_three']+=1;continue
        counts['simple_rank_three']+=1
        basis=[difference(nodes[i]['coefficients'],nodes[ids[0]]['coefficients'],p) for i in ids[1:]]
        data={'p':p,'k':source['k'],'s':source['s'],'domain':source['domain'],'basis':basis}
        cert=certify(data,200_000);review=validate(cert);brute=literal(data)
        assert review['interval']==[brute['minimum']]*2
        rows.append({'node_indices':list(ids),'certificate':cert,'independent_review':review,'literal_review':brute})
    assert sum(counts.values())==3060 and counts['simple_rank_three']==55
    largepath=HERE.parent/'round6/pencil_tracks/results.json'
    large=next(x['result'] for x in json.loads(largepath.read_text())['large'] if x['name']=='skew-two')
    p=large['field']['p'];t0,t1=large['verified_tracks'];a=t0['intercept_coefficients'];b=t0['slope_coefficients'];c=t1['intercept_coefficients'];d=t1['slope_coefficients']
    basis=[b,difference(c,a,p),d];columns=[tuple(evaluate(poly,x,p) for poly in basis) for x in large['domain']]
    profile=simple_profile(columns,p);node_checks=[]
    for scalar,poly in [(0,a),(1,difference(a,[-v%p for v in b],p)),(0,c),(1,difference(c,[-v%p for v in d],p))]:
        agreement=sum(evaluate(poly,x,p)==(u+scalar*v)%p for x,u,v in zip(large['domain'],large['u0'],large['u1']))
        assert len(poly)<=large['k'] and agreement>=large['s'];node_checks.append({'scalar':scalar,'agreement':agreement})
    out={'status':'passed','scope':'All 3060 quadruples of the saved18 nodes classified; only the55 simple rank-three spaces enter this version of the compiler. These nodes come from possibly different received-word scalars; no single-scalar cluster coverage is asserted.',
         'small_classification':counts,'small_cases':rows,'small_agreement_sets_enumerated':sum(r['literal_review']['agreement_sets'] for r in rows),
         'small_minimum_over_simple_clusters':min(r['independent_review']['interval'][0] for r in rows),
         'large_declared_space':{'p':p,'n':large['n'],'k':large['k'],'s':large['s'],'domain':large['domain'],'origin':a,'basis':basis,
                                 'source_node_checks':node_checks,'evaluation_profile':profile,
                                 'optimizer_applied':False,'next_requirement':'Full-length support and scalable verification; current public compiler is capped at256 coordinates.'},
         'source_sha256':{name:sha256((HERE/name).read_bytes()).hexdigest() for name in ('rank_three_actual.py','rank_three_pinning.py','rank_three_verifier.py','run_rank_three.py')},
         'input_sha256':{'../round6/overlap_preflight/coset_small_characteristic_sample3.certificate.json':sha256(sourcepath.read_bytes()).hexdigest(),
                         '../round6/pencil_tracks/results.json':sha256(largepath.read_bytes()).hexdigest()}}
    (HERE/'rank_three_actual.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'passed','small_classification':counts,'simple_cluster_minimum':out['small_minimum_over_simple_clusters'],
                      'large_profile':{k:v for k,v in profile.items() if k!='parallel_class_sizes'},'large_source_nodes':node_checks}),flush=True)


if __name__=='__main__':main()
