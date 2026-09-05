#!/usr/bin/env python3
"""Explicit Paley/Peisert49 counterfactual with exact pair types and six-set census."""
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
from math import comb
from pathlib import Path
import json
import subprocess
import sys
import time

sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
LAB=HERE.parents[1]
sys.path.insert(0,str(LAB/'proximity'))
sys.path.insert(0,str(LAB/'round2'/'observables'))
from extension_field_probe import Field
from exchange_grouped import grouped_coefficient,containment_falling
from review_exchange import eberlein,solve,encode,C


def matrices():
    F=Field(7,[1,0,1])
    generator=next(g for g in range(1,49) if all(F.power(g,k)!=1 for k in [16,24]))
    powers=[F.power(generator,i) for i in range(48)]
    assert len(set(powers))==48
    cosets=[sorted(powers[i] for i in range(48) if i%4==j) for j in range(4)]
    paley=set(cosets[0]+cosets[2]);peisert=set(cosets[0]+cosets[1])
    signs=[]
    for conn in [paley,peisert]:
        signs.append([[0 if x==y else (1 if F.sub(x,y) in conn else -1) for y in range(49)] for x in range(49)])
    removed=[list(e) for e in combinations(range(49),2) if signs[0][e[0]][e[1]]==1 and signs[1][e[0]][e[1]]==-1]
    added=[list(e) for e in combinations(range(49),2) if signs[0][e[0]][e[1]]==-1 and signs[1][e[0]][e[1]]==1]
    return signs,{'field_p':7,'field_degree':2,'modulus':[1,0,1],'encoding':'a+7b represents a+bX',
        'primitive_element':generator,'primitive_element_coefficients':F.decode(generator),
        'cyclotomic_classes':cosets,'paley_connection_classes':[0,2],'peisert_connection_classes':[0,1],
        'trade':'remove differences in class2, add differences in class1; explicit cyclotomic trade, not a Godsil-McKay switching sequence',
        'removed_edges':removed,'added_edges':added,'nonzero_inverses_checked':48}


def validate(S):
    q=len(S);pairs=[];pattern_counts=Counter()
    for i in range(q):
        assert S[i][i]==0 and sum(S[i])==0 and Counter(S[i])=={0:1,1:24,-1:24}
        for j in range(q):
            assert S[i][j]==S[j][i]
            inner=sum(a*b for a,b in zip(S[i],S[j]))
            assert inner==(q-1 if i==j else -1)
            hist=Counter(zip(S[i],S[j]))
            if i==j:
                expected={(0,0):1,(1,1):24,(-1,-1):24}
            else:
                eps=S[i][j]
                expected={(0,eps):1,(eps,0):1,(1,1):(q-3-2*eps)//4,
                          (-1,-1):(q-3+2*eps)//4,(1,-1):(q-1)//4,(-1,1):(q-1)//4}
                common=sum(a==b==1 for a,b in zip(S[i],S[j]))
                assert common==(11 if eps==1 else 12)
            assert dict(hist)==expected
            pattern=tuple(sorted(hist.items()));pattern_counts[pattern]+=1
            pairs.append({'rows':[i,j],'inner_product':inner,'types':[[a,b,count] for (a,b),count in sorted(hist.items())]})
    return pairs,pattern_counts


def spectrum_from_patterns(pattern_counts,q,n,d=6):
    K=[sum(multiplicity*grouped_coefficient([(a,b,count) for (a,b),count in pattern],r,d)
           for pattern,multiplicity in pattern_counts.items()) for r in range(d+1)]
    corr=[sum(Q(K[r])*containment_falling(q,n,j,r,d) for r in range(d+1)) for j in range(d+1)]
    theta=[[eberlein(q,n,z,j) for z in range(d+1)] for j in range(d+1)]
    energies=solve(theta,corr)
    assert all(v>=0 for v in energies) and sum(energies)==corr[0]
    return {'n':n,'degree':d,'kernel_overlap_sums':K,'distance_correlations':corr,'energies':energies}


def direct_kernel(S,c):
    out=0
    for row in S:
        v=1
        for i in c:v*=row[i]
        out+=v
    return out


def compiled_census(signs):
    compiler='clang++';binary=HERE/'census_six'
    subprocess.run([compiler,'-O3','-std=c++17',str(HERE/'census_six.cpp'),'-o',str(binary)],check=True)
    masks=[[sum(1<<i for i in range(49) if row[i]<0) for row in S] for S in signs]
    # Symmetry turns row masks into column masks needed by the XOR formula.
    payload='49\n'+'\n'.join(' '.join(map(str,m)) for m in masks)+'\n'
    started=time.perf_counter()
    raw=subprocess.run([str(binary)],input=payload,text=True,capture_output=True,check=True).stdout
    result={'paley':{},'peisert':{},'joint':[],'witnesses':{'paley':{},'peisert':{}}}
    for line in raw.splitlines():
        pieces=line.split();kind=pieces[0];xs=list(map(int,pieces[1:]))
        if kind=='COUNT':result['six_sets']=xs[0]
        elif kind in ['P','Q']:
            name='paley' if kind=='P' else 'peisert'
            value,count,*w=xs;result[name][value]=count;result['witnesses'][name][value]=w
            assert direct_kernel(signs[0 if kind=='P' else 1],w)==value
        elif kind=='J':result['joint'].append({'paley_value':xs[0],'peisert_value':xs[1],'count':xs[2]})
        elif kind=='D':result['largest_pointwise_difference']={'paley_value':xs[0],'peisert_value':xs[1],'C':xs[2:]}
    assert result['six_sets']==comb(49,6)
    assert sum(result['paley'].values())==sum(result['peisert'].values())==result['six_sets']
    assert sum(x['count'] for x in result['joint'])==result['six_sets']
    result['elapsed_seconds']=round(time.perf_counter()-started,6)
    return result


def main():
    started=time.perf_counter();signs,construction=matrices()
    # Explicit inverse check for the field certificate.
    F=Field(7,[1,0,1]);assert all(F.mul(x,F.inv(x))==1 for x in range(1,49))
    validations=[validate(S) for S in signs]
    assert validations[0][1]==validations[1][1]
    trade=[row[:] for row in signs[0]]
    for a,b in construction['removed_edges']:trade[a][b]=trade[b][a]=-1
    for a,b in construction['added_edges']:trade[a][b]=trade[b][a]=1
    assert trade==signs[1]
    spectra=[spectrum_from_patterns(validations[0][1],49,n) for n in [6,7,8]]
    result=compiled_census(signs)
    moments={}
    for name in ['paley','peisert']:
        moments[name]={str(k):Q(sum(count*value**k for value,count in result[name].items()),result['six_sets']) for k in [1,2,3,4,6]}
    assert moments['paley']['1']==moments['peisert']['1']
    assert moments['paley']['2']==moments['peisert']['2']==sum(spectra[0]['energies'])
    result['moments']=moments
    result['histograms_differ']=result['paley']!=result['peisert']
    result['isomorphism_certificate']='different complete six-kernel histograms' if result['histograms_differ'] else 'not established by six-kernel census'
    metadata={'status':'exact finite counterfactual audit','construction':construction,'sign_matrices':signs,
              'pair_validation':{'each_matrix_ordered_pairs_checked':49**2,'all_pair_types_equal_as_multiset':True,
                                 'row_sums':[0]*49,'srg_parameters':[49,24,11,12],
                                 'patterns':[[[[a,b,count] for (a,b),count in pattern],multiplicity] for pattern,multiplicity in validations[0][1].items()]},
              'all_input_L2_spectra_identical':spectra,'census':result,
              'source_sha256':{name:sha256((HERE/name).read_bytes()).hexdigest() for name in ['build_and_census.py','census_six.cpp']},
              'elapsed_seconds':round(time.perf_counter()-started,3)}
    (HERE/'results.json').write_text(json.dumps(encode(metadata),indent=2)+'\n')
    (HERE/'pair_type_certificates.json').write_text(json.dumps({'paley':validations[0][0],'peisert':validations[1][0]},indent=2)+'\n')
    print('primitive',construction['primitive_element'],'trade',len(construction['removed_edges']),len(construction['added_edges']))
    print('count',result['six_sets'],'hist different',result['histograms_differ'],'census seconds',result['elapsed_seconds'])
    print('Paley histogram',result['paley']);print('Peisert histogram',result['peisert'])
    print('pointwise witness',result['largest_pointwise_difference'])
    print('moments',encode(moments))


if __name__=='__main__':main()
