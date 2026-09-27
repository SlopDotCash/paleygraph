#!/usr/bin/env python3
"""Independent Gaussian execution review with complete finite/root-bound oracles."""
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from cover_verifier import validate,evaluate

HERE=Path(__file__).resolve().parent


def solve(rows,rhs,p):
    matrix=[list(row)+[v] for row,v in zip(rows,rhs)]
    for j in range(3):
        pivot=next((i for i in range(j,3) if matrix[i][j]%p),None)
        assert pivot is not None
        matrix[j],matrix[pivot]=matrix[pivot],matrix[j]
        inv=pow(matrix[j][j]%p,-1,p);matrix[j]=[v*inv%p for v in matrix[j]]
        for i in range(3):
            if i==j:continue
            a=matrix[i][j];matrix[i]=[(x-a*y)%p for x,y in zip(matrix[i],matrix[j])]
    return tuple(row[3] for row in matrix)


def replay(c,record):
    d=c['input'];p,n,k,s=(d[x] for x in ('p','n','k','s'));received=record['received']
    assert record['schema']=='deterministic_declared_space_pruning_v1'
    assert len(received)==n and all(type(v) is int and 0<=v<p for v in received)
    columns=[tuple(evaluate(b,x,p) for b in d['basis']) for x in d['domain']]
    origin=[evaluate(d['origin'],x,p) for x in d['domain']]
    queries=[q for row in c['blocks'] for q in row['queries']]
    assert record['queries_executed']==len(queries)
    seen=set();accepted={};duplicates=0
    for index,Q in enumerate(queries):
        params=solve([columns[i] for i in Q],[(received[i]-origin[i])%p for i in Q],p)
        if params in seen:duplicates+=1;continue
        seen.add(params)
        word=[(a+params[0]*v[0]+params[1]*v[1]+params[2]*v[2])%p for a,v in zip(origin,columns)]
        support=[i for i in range(n) if word[i]==received[i]]
        if len(support)>=s:accepted[params]=(index,support,word)
    assert record['distinct_candidates_evaluated']==len(seen) and record['duplicate_candidates']==duplicates
    assert [tuple(row['parameters']) for row in record['outputs']]==sorted(accepted)
    for row in record['outputs']:
        params=tuple(row['parameters']);index,support,word=accepted[params]
        assert row['first_query_index']==index and row['agreement_support']==support and row['codeword']==word
        coeffs=[((d['origin'][j] if j<len(d['origin']) else 0)+sum(params[i]*(d['basis'][i][j] if j<len(d['basis'][i]) else 0) for i in range(3)))%p for j in range(k)]
        assert row['coefficients']==coeffs and [evaluate(coeffs,x,p) for x in d['domain']]==word
    return {'queries_replayed':len(queries),'distinct_candidate_evaluations':len(seen),'outputs_checked':len(accepted)}


def main():
    path=HERE/'actual_results.json';actual=json.loads(path.read_text());assert actual['status']=='passed'
    smallpath=HERE.parent/'round6/overlap_preflight/coset_small_characteristic_sample3.certificate.json';small=json.loads(smallpath.read_text())
    largepath=HERE.parent/'round6/pencil_tracks/results.json';large=next(x['result'] for x in json.loads(largepath.read_text())['large'] if x['name']=='skew-two')
    rows=[]
    assert [case['name'] for case in actual['cases']]==['small','full_length']
    for case in actual['cases']:
        c=json.loads((HERE/case['certificate']).read_text());review=validate(c);assert review==case['coverage_review']
        d=c['input'];p=d['p'];n=d['n'];s=d['s'];k=d['k']
        origin=[evaluate(d['origin'],x,p) for x in d['domain']]
        columns=[tuple(evaluate(b,x,p) for b in d['basis']) for x in d['domain']]
        if case['name']=='small':
            source=small
            all_words=[(params,[(a+params[0]*v[0]+params[1]*v[1]+params[2]*v[2])%p for a,v in zip(origin,columns)]) for params in product(range(p),repeat=3)]
            assert len(all_words)==len({tuple(w) for _,w in all_words})==4913
            assert [r['scalar'] for r in case['runs']]==list(range(17))
        else:
            source=large;assert [r['scalar'] for r in case['runs']]==[0,1]
        for row in case['runs']:
            record=row['run'];scalar=row['scalar'];received=record['received']
            assert received==[(a+scalar*b)%p for a,b in zip(source['u0'],source['u1'])]
            result=replay(c,record)
            if case['name']=='small':
                truth=[list(params) for params,word in all_words if sum(a==b for a,b in zip(word,received))>=s]
                assert truth==[o['parameters'] for o in record['outputs']]
                oracle={'full_space_members_enumerated':len(all_words)}
            else:
                pieces=[]
                for track in source['verified_tracks']:
                    a=track['intercept_coefficients'];b=track['slope_coefficients']
                    assert len(a)<=k and len(b)<=k
                    pieces.append([((a[i] if i<len(a) else 0)+scalar*(b[i] if i<len(b) else 0))%p for i in range(k)])
                assert len(pieces)==2 and pieces[0]!=pieces[1]
                supports=[[i for i,x in enumerate(d['domain']) if evaluate(poly,x,p)==received[i]] for poly in pieces]
                assert set(supports[0])|set(supports[1])==set(range(n)) and all(len(S)>=s for S in supports)
                assert 2*(k-1)<s and sorted(pieces)==sorted(o['coefficients'] for o in record['outputs'])
                oracle={'two_piece_agreement_sizes':[len(S) for S in supports],'other_polynomial_agreement_cap':2*(k-1)}
            rows.append({'case':case['name'],'scalar':scalar,**result,'complete_oracle':oracle})
    out={'status':'passed','scope':'All saved queries independently replayed; actual output completeness checked by exhaustive small-space enumeration and a full ambient two-piece root bound on the two special large words. Universal supplied-space completeness follows separately from checked query-cover certificates.',
         'cases':rows,'queries_replayed':sum(row['queries_replayed'] for row in rows),
         'candidate_evaluations_recomputed':sum(row['distinct_candidate_evaluations'] for row in rows),
         'source_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['pruning_review.py','cover_verifier.py']},
         'input_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['actual_results.json','small.certificate.json','full_length.certificate.json','../round6/overlap_preflight/coset_small_characteristic_sample3.certificate.json','../round6/pencil_tracks/results.json']}}
    (HERE/'pruning_review.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'passed','cases':len(rows),'queries_replayed':out['queries_replayed'],'candidate_evaluations_recomputed':out['candidate_evaluations_recomputed']}),flush=True)


if __name__=='__main__':main()
