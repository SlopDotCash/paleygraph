#!/usr/bin/env python3
"""Independent query replay and full-word review of every candidate decision.

Queries use integer Cramer determinants, not the producer's modular Gaussian
solver. All candidate words are independently evaluated in bounded int64 batches
after proving that the unreduced dot products cannot overflow.
"""
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations,product
import json
from math import comb
from pathlib import Path
from time import perf_counter
import numpy as np
from arc_verifier import validate,determinant,evaluate

HERE=Path(__file__).resolve().parent


def solve(rows,rhs,p):
    det=determinant(rows)%p;assert det;inv=pow(det,-1,p);answer=[]
    for j in range(len(rows)):
        changed=[list(v) for v in rows]
        for i in range(len(rows)):changed[i][j]=rhs[i]
        answer.append(determinant(changed)*inv%p)
    return tuple(answer)


def query_replay(c,transcript,columns,origin):
    data=c['input'];p,n,d=(data[x] for x in ('p','n','dimension'));received=transcript['received']
    assert transcript['schema']=='complete_arc_query_transcript_v1'
    assert len(received)==n and all(type(v) is int and 0<=v<p for v in received)
    reconstructed={};index=0
    for j,B in enumerate(c['blocks']):
        labels=[];local_queries=list(combinations(range(len(B)),d))
        for local in local_queries:
            Q=[B[t] for t in local];params=solve([columns[i] for i in Q],[(received[i]-origin[i])%p for i in Q],p);labels.append(params)
            if params not in reconstructed:reconstructed[params]={'first_query_index':index,'positive_blocks':[]}
            index+=1
        counts=Counter(labels);unions=defaultdict(set)
        for label,Q in zip(labels,local_queries):unions[label].update(Q)
        for params,count in counts.items():
            S=unions[params];assert count==comb(len(S),d)
            reconstructed[params]['positive_blocks'].append({'block':j,'count':count,'support_mask':sum(2**i for i in S)})
    expected=[{'parameters':list(params),**reconstructed[params]} for params in sorted(reconstructed)]
    assert expected==transcript['candidates'] and index==transcript['queries_executed']==c['query_count']
    return index


def review_run(c,result):
    data=c['input'];p,n,k,d,s=(data[x] for x in ('p','n','k','dimension','s'));blocks=c['blocks']
    transcript=result['transcript'];filtered=result['filter'];received=transcript['received'];candidates=transcript['candidates']
    assert filtered['schema']=='arc_transcript_pruning_v1'
    columns=[tuple(evaluate(b,x,p) for b in data['basis']) for x in data['domain']];origin=[evaluate(data['origin'],x,p) for x in data['domain']]
    queries=query_replay(c,transcript,columns,origin)
    assert all(len(row['parameters'])==d and all(type(v) is int and 0<=v<p for v in row['parameters']) for row in candidates)
    assert [tuple(row['parameters']) for row in candidates]==sorted({tuple(row['parameters']) for row in candidates})
    decisions=filtered['decisions'];assert len(decisions)==len(candidates)
    # All terms are nonnegative canonical residues, and the origin addition is
    # bounded too. This is an exact integer oracle, not floating-point numerics.
    assert d*(p-1)**2+(p-1)<2**63
    matrix=np.array(columns,dtype=np.int64).T;origin_array=np.array(origin,dtype=np.int64);received_array=np.array(received,dtype=np.int64)
    owner={i:j for j,B in enumerate(blocks) for i in B};base=[min(len(B),d-1) for B in blocks]
    accepted=[];evaluations=0;reject_hist=Counter();output_index=0
    for start in range(0,len(candidates),256):
        chunk=candidates[start:start+256];params=np.array([r['parameters'] for r in chunk],dtype=np.int64)
        words=(params@matrix+origin_array)%p;equal=words==received_array
        for local,row in enumerate(chunk):
            index=start+local;decision=decisions[index];assert decision['candidate_index']==index
            positives={v['block']:v for v in row['positive_blocks']};caps=base[:]
            for j,B in enumerate(blocks):
                actual=[t for t,i in enumerate(B) if bool(equal[local,i])]
                if j in positives:
                    r=positives[j];assert r['support_mask']==sum(2**t for t in actual) and r['count']==comb(len(actual),d)>0
                    caps[j]=len(actual)
                else:assert len(actual)<=base[j]
            upper=sum(caps);assert upper==decision['initial_agreement_upper']
            unknown=sorted((j for j in range(len(blocks)) if j not in positives),key=lambda j:(len(blocks[j]),j))
            order=[i for j in unknown for i in blocks[j]];tested=decision['tested_coordinates']
            assert [record[0] for record in tested]==order[:len(tested)] and len(tested)<=len(order)
            done=Counter();matches=Counter()
            for i,agrees in tested:
                assert upper>=s and type(i) is int and type(agrees) is int and agrees in (0,1)
                assert bool(agrees)==bool(equal[local,i]);j=owner[i];done[j]+=1;matches[j]+=agrees
                revised=min(base[j],matches[j]+len(blocks[j])-done[j]);upper+=revised-caps[j];caps[j]=revised
            evaluations+=len(tested);assert upper==decision['final_agreement_upper']
            actual_count=int(equal[local].sum());assert actual_count<=upper
            assert type(decision['accepted']) is bool and decision['accepted']==(actual_count>=s)
            if decision['accepted']:
                assert len(tested)==len(order) and actual_count==upper
                out=filtered['outputs'][output_index];output_index+=1
                assert out['parameters']==row['parameters'] and out['first_query_index']==row['first_query_index']
                word=[int(v) for v in words[local]];support=[i for i in range(n) if bool(equal[local,i])]
                assert out['codeword']==word and out['agreement_support']==support
                params=row['parameters'];coeffs=[((data['origin'][j] if j<len(data['origin']) else 0)+sum(params[i]*(data['basis'][i][j] if j<len(data['basis'][i]) else 0) for i in range(d)))%p for j in range(k)]
                assert coeffs==out['coefficients'] and [evaluate(coeffs,x,p) for x in data['domain']]==word
                accepted.append(params)
            else:
                assert upper<s;reject_hist[len(tested)]+=1
    assert output_index==len(filtered['outputs'])
    assert filtered['candidate_coordinate_evaluations']==evaluations
    assert filtered['output_word_coordinate_evaluations']==len(accepted)*n
    assert filtered['total_explicit_coordinate_evaluations']==evaluations+len(accepted)*n
    assert filtered['full_scan_coordinate_evaluations']==len(candidates)*n
    return {'queries_independently_replayed':queries,'candidate_words_fully_evaluated':len(candidates),
            'full_word_coordinate_values_checked':len(candidates)*n,
            'adaptive_coordinate_evaluations_checked':evaluations,'output_materialization_evaluations':len(accepted)*n,
            'outputs_checked':len(accepted),'rejection_evaluation_histogram':sorted(reject_hist.items())}


def main():
    actualpath=HERE/'actual_arc_results.json';actual=json.loads(actualpath.read_text());assert actual['status']=='passed'
    smallpath=HERE.parent/'round6/overlap_preflight/coset_small_characteristic_sample3.certificate.json';small=json.loads(smallpath.read_text())
    largepath=HERE.parent/'round6/pencil_tracks/results.json';large=next(r['result'] for r in json.loads(largepath.read_text())['large'] if r['name']=='skew-two')
    cases=[];bindings={'actual_arc_results.json':sha256(actualpath.read_bytes()).hexdigest()}
    for case in actual['cases']:
        cpath=HERE/case['certificate'];c=json.loads(cpath.read_text());assert validate(c)==case['review'];bindings[cpath.name]=sha256(cpath.read_bytes()).hexdigest()
        data=c['input'];p,n,k,d,s=(data[x] for x in ('p','n','k','dimension','s'));raw=small if case['name']=='small' else large
        if case['name']=='small':
            columns=[tuple(evaluate(b,x,p) for b in data['basis']) for x in data['domain']];origin=[evaluate(data['origin'],x,p) for x in data['domain']]
            all_words=[(params,[(a+sum(t*v for t,v in zip(params,col)))%p for a,col in zip(origin,columns)]) for params in product(range(p),repeat=d)]
            assert len(all_words)==4913
        for record in case['runs']:
            path=HERE/record['run_file'];result=json.loads(path.read_text());bindings[path.name]=sha256(path.read_bytes()).hexdigest();scalar=record['scalar'];received=result['transcript']['received']
            assert received==[(a+scalar*b)%p for a,b in zip(raw['u0'],raw['u1'])]
            started=perf_counter();checked=review_run(c,result)
            if case['name']=='small':
                truth=[list(params) for params,word in all_words if sum(a==b for a,b in zip(word,received))>=s]
                assert truth==[o['parameters'] for o in result['filter']['outputs']]
                oracle={'kind':'complete_space','members':4913}
            else:
                pieces=[]
                for track in large['verified_tracks']:
                    a=track['intercept_coefficients'];b=track['slope_coefficients'];assert len(a)<=k and len(b)<=k
                    pieces.append([((a[j] if j<len(a) else 0)+scalar*(b[j] if j<len(b) else 0))%p for j in range(k)])
                assert len(pieces)==2 and pieces[0]!=pieces[1]
                supports=[[i for i,x in enumerate(data['domain']) if evaluate(poly,x,p)==received[i]] for poly in pieces]
                assert set(supports[0])|set(supports[1])==set(range(n)) and all(len(A)>=s for A in supports)
                assert 2*(k-1)<s and sorted(pieces)==sorted(o['coefficients'] for o in result['filter']['outputs'])
                oracle={'kind':'two_piece_complete_ambient','agreement_sizes':[len(A) for A in supports],'other_polynomial_agreement_cap':2*(k-1)}
            row={'case':case['name'],'scalar':scalar,'review_seconds':perf_counter()-started,'review':checked,'complete_oracle':oracle};cases.append(row)
            print(json.dumps(row),flush=True)
    out={'status':'passed','scope':'Every query independently solved with integer Cramer determinants; every candidate word fully evaluated with overflow-bounded integer arrays. All local support inferences, adaptive rejection bounds, outputs and special complete oracles checked. No timing speedup or global prize theorem is inferred.',
         'cases':cases,'queries_independently_replayed':sum(r['review']['queries_independently_replayed'] for r in cases),
         'candidate_words_fully_evaluated':sum(r['review']['candidate_words_fully_evaluated'] for r in cases),
         'full_word_coordinate_values_checked':sum(r['review']['full_word_coordinate_values_checked'] for r in cases),
         'source_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['arc_review.py','arc_verifier.py']},
         'input_sha256':{**bindings,'../round6/pencil_tracks/results.json':sha256(largepath.read_bytes()).hexdigest(),'../round6/overlap_preflight/coset_small_characteristic_sample3.certificate.json':sha256(smallpath.read_bytes()).hexdigest()}}
    (HERE/'arc_review.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('cases','source_sha256','input_sha256')}),flush=True)


if __name__=='__main__':main()
