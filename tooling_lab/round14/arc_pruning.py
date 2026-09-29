#!/usr/bin/env python3
"""Retain query occurrence and absence as exact local agreement information.

Every d-subset of an arc block is queried. A parameter vector agreeing on h
coordinates appears C(h,d) times. If it appears, the union of its query supports
is its entire local agreement set. Otherwise h<=d-1. These facts provide global
upper bounds for candidate filtering before full-word materialization.
"""
from itertools import combinations
from math import comb
from time import perf_counter
from arc_verifier import validate


def evaluate(poly,x,p):
    value=0
    for a in poly[::-1]:value=(value*x+a)%p
    return value


def solve(rows,rhs,p):
    a=[list(row)+[v] for row,v in zip(rows,rhs)];d=len(a)
    for j in range(d):
        pivot=next((i for i in range(j,d) if a[i][j]%p),None)
        if pivot is None:raise ValueError('certified query lost rank')
        a[j],a[pivot]=a[pivot],a[j]
        inv=pow(a[j][j]%p,-1,p)
        for k in range(j,d+1):a[j][k]=a[j][k]*inv%p
        for i in range(j+1,d):
            v=a[i][j]
            for k in range(j,d+1):a[i][k]=(a[i][k]-v*a[j][k])%p
    answer=[0]*d
    for j in range(d-1,-1,-1):answer[j]=(a[j][d]-sum(a[j][k]*answer[k] for k in range(j+1,d)))%p
    return tuple(answer)


def prepare(c):
    review=validate(c);data=c['input'];p=data['p']
    return {'data':data,'blocks':c['blocks'],'review':review,
            'columns':[tuple(evaluate(b,x,p) for b in data['basis']) for x in data['domain']],
            'origin':[evaluate(data['origin'],x,p) for x in data['domain']]}


def collect(prepared,received):
    data=prepared['data'];p,n,d=(data[x] for x in ('p','n','dimension'));columns=prepared['columns'];origin=prepared['origin']
    if len(received)!=n or any(type(v) is not int or not 0<=v<p for v in received):raise ValueError('canonical field-valued received word of certified length required')
    found={};index=0
    for j,B in enumerate(prepared['blocks']):
        for local in combinations(range(len(B)),d):
            Q=[B[t] for t in local];params=solve([columns[i] for i in Q],[(received[i]-origin[i])%p for i in Q],p)
            if params not in found:found[params]={'first_query_index':index,'positive_blocks':{}}
            row=found[params]['positive_blocks'].setdefault(j,[0,0]);row[0]+=1;row[1]|=sum(1<<t for t in local);index+=1
    candidates=[]
    for params in sorted(found):
        old=found[params];positives=[]
        for j,(count,mask) in sorted(old['positive_blocks'].items()):
            assert count==comb(mask.bit_count(),d)
            positives.append({'block':j,'count':count,'support_mask':mask})
        candidates.append({'parameters':list(params),'first_query_index':old['first_query_index'],'positive_blocks':positives})
    assert index==prepared['review']['queries_rank_checked']
    return {'schema':'complete_arc_query_transcript_v1','received':list(received),'queries_executed':index,'candidates':candidates}


def filter_candidates(prepared,transcript):
    data=prepared['data'];p,n,k,d,s=(data[x] for x in ('p','n','k','dimension','s'));blocks=prepared['blocks'];received=transcript['received'];columns=prepared['columns'];origin=prepared['origin']
    base_caps=[min(len(B),d-1) for B in blocks];base=sum(base_caps)
    outputs=[];decisions=[];tests=materialized=0
    for index,row in enumerate(transcript['candidates']):
        params=row['parameters'];positive={v['block']:v for v in row['positive_blocks']}
        support=[];upper=base
        for j,record in positive.items():
            mask=record['support_mask'];size=mask.bit_count();assert record['count']==comb(size,d)
            upper+=size-base_caps[j]
            support.extend(i for t,i in enumerate(blocks[j]) if mask>>t&1)
        initial=upper;tested=[]
        # Unqueried singletons are often decisive. This fixed ordering is a
        # simple policy, not a claim of optimal adaptive testing.
        unknown=sorted((j for j in range(len(blocks)) if j not in positive),key=lambda j:(len(blocks[j]),j))
        for j in unknown:
            if upper<s:break
            B=blocks[j];cap=base_caps[j];matches=0
            for done,i in enumerate(B,1):
                agrees=(origin[i]+sum(a*b for a,b in zip(params,columns[i])))%p==received[i]
                tested.append([i,int(agrees)]);tests+=1
                if agrees:matches+=1;support.append(i)
                next_cap=min(base_caps[j],matches+len(B)-done);upper+=next_cap-cap;cap=next_cap
                if upper<s:break
            assert matches<=base_caps[j]
        accepted=upper>=s
        if accepted:
            assert upper==len(support)>=s
            word=[(a+sum(c*v for c,v in zip(params,column)))%p for a,column in zip(origin,columns)];materialized+=n
            support.sort();assert [i for i in range(n) if word[i]==received[i]]==support
            coefficients=[((data['origin'][j] if j<len(data['origin']) else 0)+sum(params[i]*(data['basis'][i][j] if j<len(data['basis'][i]) else 0) for i in range(d)))%p for j in range(k)]
            outputs.append({'parameters':params,'coefficients':coefficients,'codeword':word,'agreement_support':support,'first_query_index':row['first_query_index']})
        decisions.append({'candidate_index':index,'initial_agreement_upper':initial,'final_agreement_upper':upper,'tested_coordinates':tested,'accepted':accepted})
    return {'schema':'arc_transcript_pruning_v1','outputs':outputs,'decisions':decisions,
            'candidate_coordinate_evaluations':tests,'output_word_coordinate_evaluations':materialized,
            'total_explicit_coordinate_evaluations':tests+materialized,
            'full_scan_coordinate_evaluations':len(transcript['candidates'])*n,
            'scope':'All query labels, including absence from each block, give deterministic agreement bounds. Filtering is exact inside the supplied space, conditional on a checked cover and complete query transcript.'}


def run(prepared,received):
    start=perf_counter();transcript=collect(prepared,received);query_seconds=perf_counter()-start
    start=perf_counter();filtered=filter_candidates(prepared,transcript);filter_seconds=perf_counter()-start
    return {'transcript':transcript,'filter':filtered,'query_seconds':query_seconds,'filter_seconds':filter_seconds}


if __name__=='__main__':
    import json,sys
    from pathlib import Path
    if len(sys.argv)!=3:raise SystemExit('Usage: arc_pruning.py certificate.json received.json')
    print(json.dumps(run(prepare(json.loads(Path(sys.argv[1]).read_text())),json.loads(Path(sys.argv[2]).read_text())),separators=(',',':')))
