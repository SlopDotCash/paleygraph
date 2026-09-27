#!/usr/bin/env python3
"""Independent Gaussian query replay and complete-oracle review of pruning runs."""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from certificate_verifier import validate

HERE=Path(__file__).resolve().parent


def solve(rows,rhs,p):
    matrix=[list(row)+[v] for row,v in zip(rows,rhs)]
    for j in range(3):
        pivot=next((i for i in range(j,3) if matrix[i][j]%p),None)
        if pivot is None:return None
        matrix[j],matrix[pivot]=matrix[pivot],matrix[j]
        inv=pow(matrix[j][j]%p,-1,p);matrix[j]=[v*inv%p for v in matrix[j]]
        for i in range(3):
            if i==j:continue
            a=matrix[i][j];matrix[i]=[(v-a*w)%p for v,w in zip(matrix[i],matrix[j])]
    return tuple(row[3] for row in matrix)


def evaluate(poly,x,p):return sum(v*pow(x,j,p) for j,v in enumerate(poly))%p


def replay(c,verified,run):
    data=c['input'];p,n,k,s=(data[x] for x in ('p','n','k','s'));domain=data['domain'];basis=data['basis'];origin=data['origin'];received=run['received']
    assert run['schema']=='declared_space_pruning_run_v1' and len(received)==n and all(type(v) is int and 0<=v<p for v in received)
    budget=run['sampling_budget'];num,den=verified['success_probability_lower'];T=budget['trials'];bits=budget['failure_bits']
    assert type(T) is int and 1<=T<=100_000 and type(bits) is int and 1<=bits<=128
    assert budget['candidate_count_upper']==p**3 and budget['per_candidate_success_lower']==[num,den]
    assert budget['failure_probability_at_most']==[1,1<<bits]
    assert p**3*pow(den-num,T)*(1<<bits)<=pow(den,T)
    assert T==1 or p**3*pow(den-num,T-1)*(1<<bits)>pow(den,T-1)
    queries=run['queries'];assert len(queries)==T
    assert all(B==sorted(set(B)) and len(B)==3 and all(type(i) is int and 0<=i<n for i in B) for B in queries)
    columns=[tuple(evaluate(poly,x,p) for poly in basis) for x in domain];origin_values=[evaluate(origin,x,p) for x in domain]
    seen=set();accepted={};singular=duplicates=0
    for index,B in enumerate(queries):
        parameters=solve([columns[i] for i in B],[(received[i]-origin_values[i])%p for i in B],p)
        if parameters is None:singular+=1;continue
        if parameters in seen:duplicates+=1;continue
        seen.add(parameters)
        word=[(a+parameters[0]*v[0]+parameters[1]*v[1]+parameters[2]*v[2])%p for a,v in zip(origin_values,columns)]
        support=[i for i in range(n) if word[i]==received[i]]
        if len(support)>=s:accepted[parameters]=(index,support,word)
    assert len(seen)==run['distinct_candidates_evaluated'] and singular==run['singular_queries'] and duplicates==run['duplicate_candidates']
    assert [tuple(row['parameters']) for row in run['outputs']]==sorted(accepted)
    for row in run['outputs']:
        params=tuple(row['parameters']);first,support,word=accepted[params]
        assert row['first_query_index']==first and row['agreement_support']==support and row['codeword']==word
        expected=[((origin[j] if j<len(origin) else 0)+sum(params[i]*(basis[i][j] if j<len(basis[i]) else 0) for i in range(3)))%p for j in range(k)]
        assert row['coefficients']==expected and [evaluate(expected,x,p) for x in domain]==word
    return {'p':p,'n':n,'s':s,'queries_replayed':T,'singular_queries':singular,'distinct_candidates_recomputed':len(seen),
            'outputs_checked':len(accepted),'failure_probability_at_most_under_uniform_sampling':[1,1<<bits]}


def main():
    path=HERE/'pruning_results.json';source=json.loads(path.read_text());assert source['status']=='passed'
    smallcpath=HERE/'hard_rank3.certificate.json';largecpath=HERE/'full_length.certificate.json';smallc=json.loads(smallcpath.read_text());largec=json.loads(largecpath.read_text())
    smallverified=validate(smallc);largeverified=validate(largec)
    smallsourcepath=HERE.parent/'round6/overlap_preflight/coset_small_characteristic_sample3.certificate.json';smallsource=json.loads(smallsourcepath.read_text())
    largesourcepath=HERE.parent/'round6/pencil_tracks/results.json';largesource=next(x['result'] for x in json.loads(largesourcepath.read_text())['large'] if x['name']=='skew-two')
    d=smallc['input'];p=d['p'];origin=[evaluate(d['origin'],x,p) for x in d['domain']];columns=[tuple(evaluate(b,x,p) for b in d['basis']) for x in d['domain']]
    all_words=[(parameters,[(a+parameters[0]*v[0]+parameters[1]*v[1]+parameters[2]*v[2])%p for a,v in zip(origin,columns)]) for parameters in product(range(p),repeat=3)]
    assert len(all_words)==p**3==source['small_full_space_members_enumerated'] and len({tuple(word) for parameters,word in all_words})==p**3
    rows=[];prob=F(0)
    assert sorted(r['scalar'] for r in source['small_runs'])==list(range(p))
    for record in source['small_runs']:
        run=record['run'];scalar=record['scalar'];assert run['received']==[(u+scalar*v)%p for u,v in zip(smallsource['u0'],smallsource['u1'])]
        review=replay(smallc,smallverified,run)
        truth=[list(parameters) for parameters,word in all_words if sum(a==b for a,b in zip(word,run['received']))>=d['s']]
        assert record['oracle_parameters']==truth==[row['parameters'] for row in run['outputs']]
        rows.append({'case':'small','scalar':scalar,**review});prob+=F(*review['failure_probability_at_most_under_uniform_sampling'])
    assert sorted(r['scalar'] for r in source['large_runs'])==[0,1]
    for record in source['large_runs']:
        run=record['run'];scalar=record['scalar'];d=largec['input'];p=d['p'];assert run['received']==[(u+scalar*v)%p for u,v in zip(largesource['u0'],largesource['u1'])]
        review=replay(largec,largeverified,run);pieces=record['oracle_polynomials'];assert len(pieces)==2 and pieces[0]!=pieces[1]
        assert all(len(poly)<=d['k'] and all(type(v) is int and 0<=v<p for v in poly) for poly in pieces)
        supports=[[i for i,x in enumerate(d['domain']) if evaluate(poly,x,p)==run['received'][i]] for poly in pieces]
        assert supports==record['oracle_agreement_supports'] and all(len(S)>=d['s'] for S in supports)
        assert set(supports[0])|set(supports[1])==set(range(d['n']))
        assert record['other_degree_bounded_polynomial_agreement_cap']==2*(d['k']-1)<d['s']
        assert sorted(pieces)==sorted(row['coefficients'] for row in run['outputs'])
        rows.append({'case':'large','scalar':scalar,**review});prob+=F(*review['failure_probability_at_most_under_uniform_sampling'])
    assert prob<=F(*source['bundle_failure_probability_at_most_under_uniform_sampling'])==F(1,1<<40)
    out={'status':'passed','scope':'Every saved query replayed using Gaussian elimination, all emitted candidates checked, complete small-space and two-piece large oracles verified. Probability statement remains conditional on independent uniform draws; the saved outputs themselves are completely checked on these inputs.',
         'cases':rows,'queries_replayed':sum(r['queries_replayed'] for r in rows),'candidate_evaluations_recomputed':sum(r['distinct_candidates_recomputed'] for r in rows),
         'small_space_members_enumerated':len(all_words),'large_global_degree_root_bound_checked':True,
         'sum_of_run_failure_bounds':[prob.numerator,prob.denominator],
         'source_sha256':{'pruning_verifier.py':sha256(Path(__file__).read_bytes()).hexdigest(),'certificate_verifier.py':sha256((HERE/'certificate_verifier.py').read_bytes()).hexdigest()},
         'input_sha256':{'pruning_results.json':sha256(path.read_bytes()).hexdigest(),'hard_rank3.certificate.json':sha256(smallcpath.read_bytes()).hexdigest(),'full_length.certificate.json':sha256(largecpath.read_bytes()).hexdigest(),
                         '../round6/overlap_preflight/coset_small_characteristic_sample3.certificate.json':sha256(smallsourcepath.read_bytes()).hexdigest(),
                         '../round6/pencil_tracks/results.json':sha256(largesourcepath.read_bytes()).hexdigest()}}
    (HERE/'pruning_verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'passed','queries_replayed':out['queries_replayed'],'candidate_evaluations_recomputed':out['candidate_evaluations_recomputed'],'cases':len(rows)}),flush=True)


if __name__=='__main__':main()
