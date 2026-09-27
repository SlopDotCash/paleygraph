#!/usr/bin/env python3
"""Received-word-dependent split on one projective equation class.

If m_b coordinates in the chosen class ask for theta.w=b, separate the values
with m_b>h. Each such fibre has dimension one less and needs s-m_b agreements
outside the class. All other candidates need at least s-h outside agreements.
The covers and parameter lifts certify completeness within the supplied space.
This prototype currently requires original dimension three and residual
thresholds at least the residual dimension; unsupported inputs are explicit.
"""
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'round15'))
from capacity_search import find_capacity_cover,projective_groups,inputs
sys.path.insert(0,str(HERE.parent/'round14'))
from arc_profile import profile


def add_scaled(origin,basis,scalar,k,p):
    return [((origin[i] if i<len(origin) else 0)+scalar*(basis[i] if i<len(basis) else 0))%p for i in range(k)]


def reduce_polynomial(poly,domain,p):
    """Canonical equivalent representative of degree below punctured length."""
    if len(poly)<=len(domain):return poly[:]
    modulus=[1]
    for x in domain:
        new=[0]*(len(modulus)+1)
        for j,a in enumerate(modulus):new[j]=(new[j]-x*a)%p;new[j+1]=(new[j+1]+a)%p
        modulus=new
    result=poly[:]
    while len(result)>len(domain):
        top=result[-1];offset=len(result)-len(modulus)
        for j,a in enumerate(modulus):result[offset+j]=(result[offset+j]-top*a)%p
        assert result[-1]==0;result.pop()
    return result


def make_punctured(data,outside,s,origin,basis):
    domain=[data['domain'][i] for i in outside];p=data['p']
    return {'p':p,'n':len(outside),'k':min(data['k'],len(outside)),'s':s,'dimension':len(basis),
            'domain':domain,'origin':reduce_polynomial(origin,domain,p),
            'basis':[reduce_polynomial(b,domain,p) for b in basis]}


def plan(data,received):
    columns=inputs(data);p,n,s,d,k=(data[x] for x in ('p','n','s','dimension','k'))
    if d!=3:raise ValueError('this conditional prototype requires dimension three')
    if len(received)!=n or any(type(v) is not int or not 0<=v<p for v in received):raise ValueError('canonical received word required')
    groups,_=projective_groups(columns,p);chosen=max(groups,key=lambda g:(len(g),-g[0]));chosen_set=set(chosen)
    outside=[i for i in range(n) if i not in chosen_set]
    pivot=next(j for j,v in enumerate(columns[chosen[0]]) if v)
    w=[v*pow(columns[chosen[0]][pivot],-1,p)%p for v in columns[chosen[0]]]
    origin_values=[sum(a*pow(x,j,p) for j,a in enumerate(data['origin']))%p for x in data['domain']]
    labels={i:(received[i]-origin_values[i])*pow(columns[i][pivot],-1,p)%p for i in chosen}
    counts=Counter(labels.values());rows=[]
    for h in range(max(counts.values())+1):
        heavy=sorted(b for b,m in counts.items() if m>h);light_threshold=s-h
        if light_threshold>len(outside):light_cost=0
        elif light_threshold<3:continue
        else:
            shape=profile(len(outside),light_threshold,3)
            if shape['status']!='profile_found':continue
            light_cost=shape['queries']
        fibre_cost=0;valid=True
        for b in heavy:
            threshold=s-counts[b]
            if threshold>len(outside):continue
            if threshold<2:valid=False;break
            shape=profile(len(outside),threshold,2)
            if shape['status']!='profile_found':valid=False;break
            fibre_cost+=shape['queries']
        if valid:rows.append({'cutoff':h,'heavy_values':heavy,'light_agreement_threshold':light_threshold,
                              'light_queries_lower_bound':light_cost,'fibre_queries_lower_bound':fibre_cost,
                              'total_queries_lower_bound':light_cost+fibre_cost})
    if not rows:raise ValueError('no supported threshold split for this prototype')
    best=min(rows,key=lambda r:(r['total_queries_lower_bound'],len(r['heavy_values']),r['cutoff']))
    return {'input':data,'received':received,'class_coordinates':chosen,'outside_coordinates':outside,
            'projective_direction':w,'pivot':pivot,'normalized_received_labels':[[i,labels[i]] for i in chosen],
            'label_counts':sorted(counts.items()),'candidate_plans':rows,'best':best,
            'scope':'Cost-optimal cutoff for this chosen projective class under balanced-profile estimates. Arithmetic feasibility and different class choices are separate.'}


def compile_plan(planned,seed=160906,max_attempts=64,max_queries=250000):
    data=planned['input'];p,k=(data[x] for x in ('p','k'));outside=planned['outside_coordinates'];best=planned['best']
    if best['total_queries_lower_bound']>max_queries:return {'status':'query_budget_exceeded','plan':planned}
    parts=[];searches=[];pivot=planned['pivot'];w=planned['projective_direction'];counts=dict(planned['label_counts'])
    threshold=best['light_agreement_threshold']
    if threshold<=len(outside):
        restricted=make_punctured(data,outside,threshold,data['origin'],data['basis'])
        search=find_capacity_cover(restricted,seed=seed,max_attempts=max_attempts,max_queries=max_queries)
        searches.append({'part':'light','search':search})
        if search['status']!='found':return {'status':'arithmetic_search_failed','plan':planned,'searches':searches}
        parts.append({'kind':'light','excluded_values':best['heavy_values'],'certificate':search['certificate']})
    free=[j for j in range(3) if j!=pivot]
    fibre_basis=[add_scaled(data['basis'][j],data['basis'][pivot],-w[j],k,p) for j in free]
    for b in best['heavy_values']:
        threshold=data['s']-counts[b]
        if threshold>len(outside):continue
        origin=add_scaled(data['origin'],data['basis'][pivot],b,k,p)
        restricted=make_punctured(data,outside,threshold,origin,fibre_basis)
        search=find_capacity_cover(restricted,seed=seed+b,max_attempts=max_attempts,max_queries=max_queries)
        searches.append({'part':f'fibre_{b}','search':search})
        if search['status']!='found':return {'status':'arithmetic_search_failed','plan':planned,'searches':searches}
        parts.append({'kind':'fibre','value':b,'class_agreements':counts[b],'free_parameter_indices':free,'certificate':search['certificate']})
    total=sum(part['certificate']['query_count'] for part in parts)
    if total>max_queries:return {'status':'query_budget_exceeded','plan':planned,'searches':searches}
    return {'status':'found','schema':'single_projective_class_conditional_cover_v1','plan':planned,'parts':parts,
            'searches':searches,'query_count':total,
            'scope':'Every qualifying parameter is covered by its heavy fibre or the light residual cover. This is a received-word-dependent certificate in the supplied space, not an all-s-coordinate-set cover or cluster discovery.'}


def main():
    prepath=HERE.parent/'round15/root_class_results.json';data=json.loads(prepath.read_text())['input']
    records=[];bindings={'../round15/root_class_results.json':sha256(prepath.read_bytes()).hexdigest()}
    for name in ['boundary','late_near_miss','true_guard_mismatch']:
        path=HERE.parent/'round15'/f'{name}.run.json'
        if not path.exists():raise SystemExit(f'required completed producer artifact missing: {path}')
        received=json.loads(path.read_text())['transcript']['received'];bindings[f'../round15/{path.name}']=sha256(path.read_bytes()).hexdigest()
        result=compile_plan(plan(data,received));assert result['status']=='found'
        output=HERE/f'{name}.conditional.json';output.write_text(json.dumps(result,separators=(',',':'))+'\n')
        records.append({'case':name,'file':output.name,'query_count':result['query_count'],'best':result['plan']['best']})
        print(json.dumps(records[-1]),flush=True)
    out={'status':'conditional_covers_compiled','scope':'All residual query covers certified by Gaussian and integer determinant implementations. Routing proof and candidate decoding have not yet been independently reviewed.',
         'cases':records,'source_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['conditional_cover.py','../round15/capacity_search.py','../round15/projective_profile.py','../round14/arc_profile.py','../round14/arc_certificate.py','../round14/dimension_preflight.py','../round14/arc_verifier.py']},
         'input_sha256':bindings,'artifact_sha256':{r['file']:sha256((HERE/r['file']).read_bytes()).hexdigest() for r in records}}
    (HERE/'conditional_results.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
