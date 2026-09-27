#!/usr/bin/env python3
"""Prototype: all conditional two-insertion means in one weighted-Gram pass.

The first round10 artifact only. Prime-field streaming and scale review remain.
"""
from fractions import Fraction as F
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent


def compile_means(matrix,selected,degree=6):
    original=np.asarray(matrix);S=np.asarray(matrix,dtype=np.int64);q=len(S);c=sorted(selected);n=len(c)
    if not np.array_equal(original,S) or S.shape!=(q,q) or q>257:
        raise ValueError('exact square matrix of order at most257 required')
    if not (2<=n<=min(64,q-2) and len(set(c))==n and 0<=min(c)<=max(c)<q and 0<=degree<=min(6,n)):
        raise ValueError('invalid query')
    if not (np.array_equal(S,S.T) and np.all(np.diag(S)==0) and np.all(np.isin(S+np.eye(q,dtype=np.int64),[-1,1]))
            and np.all(S.sum(axis=0)==0) and np.array_equal(S@S,q*np.eye(q,dtype=np.int64)-1)):
        raise ValueError('balanced conference identities required')
    m=q-n;N=m*(m-1);total0=0;linear=[0]*n;gram=[[0]*n for _ in c];correction=[[0]*n for _ in c]
    indices={a:i for i,a in enumerate(c)}
    for x in range(q):
        signs=[int(S[x,a]) for a in c];e=[1]+[0]*degree
        for sign in signs:
            for j in range(degree,0,-1):e[j]+=sign*e[j-1]
        R=sum(signs);inside=int(x in indices)
        def phi(u,v):
            first=[e[0]]
            for j in range(1,degree+1):first.append(e[j]-u*first[-1])
            second=[first[0]]
            for j in range(1,degree+1):second.append(first[j]-v*second[-1])
            get=lambda d:second[d] if d>=0 else 0
            return N*get(degree)-2*(m-1)*R*get(degree-1)+(R*R-(m-1)-inside)*get(degree-2)
        pp,mm,pm=phi(1,1),phi(-1,-1),phi(1,-1)
        w0=pp+mm+2*pm;w1=pp-mm;w2=pp+mm-2*pm;total0+=w0
        for i in range(n):
            linear[i]+=w1*signs[i]
            for j in range(i+1,n):gram[i][j]+=w2*signs[i]*signs[j]
        if inside:
            i=indices[x]
            for j in range(n):
                if i!=j:
                    s=signs[j];correction[min(i,j)][max(i,j)]+=4*phi(0,s)-w0-w1*s
    means=[]
    for i in range(n):
        for j in range(i+1,n):
            four=total0+linear[i]+linear[j]+gram[i][j]+correction[i][j]
            assert four%4==0
            numerator=four//4;mean=F(numerator,N)
            means.append({'deleted':[c[i],c[j]],'target_sum':numerator,'mean':[mean.numerator,mean.denominator]})
    return {'q':q,'selected':c,'degree':degree,'ordered_insertion_pairs_per_deletion':N,'pairs':means}


def main():
    old=HERE.parent/'round7/local_edits/review_results.py'
    spec=importlib.util.spec_from_file_location('literal',old);literal=importlib.util.module_from_spec(spec);spec.loader.exec_module(literal)
    path=HERE.parent/'round9/preflight_results.json';data=json.loads(path.read_text());cache={};checked=0
    def get(q,c,d):
        key=(q,tuple(c),d)
        if key not in cache:cache[key]=compile_means(literal.literal_signs(q),c,d)
        return {tuple(r['deleted']):r for r in cache[key]['pairs']}
    for case in data['boundary_cases']:
        r=get(case['q'],case['selected'],case['degree'])[tuple(case['deleted'])]
        assert r['mean']==case['mean'];checked+=1
    pairs=[]
    for pair in data['pairs']:
        members=[]
        for member in pair['members']:
            values=get(17,member['selected'],6)
            for fibre in member['fibres']:
                r=values[tuple(fibre['deleted'])]
                assert r['mean']==fibre['result']['mean'] and r['target_sum']==fibre['result']['target_sum'];checked+=1
            members.append({'selected':member['selected'],'mean_deck':sorted([r['mean'] for r in values.values()],key=lambda a:F(*a))})
        assert members[0]['mean_deck']!=members[1]['mean_deck'];pairs.append(members)
    out={'status':'toy_passed_scale_pending','scope':'One weighted-Gram pass computes every deletion-pair mean; no convolution. This is an initial small-matrix prototype, not a completed round.',
         'distinct_compiler_queries':len(cache),'independent_mean_comparisons':checked,'pairs':pairs,
         'source_sha256':{'deletion_pair_means.py':sha256(Path(__file__).read_bytes()).hexdigest()},
         'input_sha256':{'../round9/preflight_results.json':sha256(path.read_bytes()).hexdigest(),
                         '../round7/local_edits/review_results.py':sha256(old.read_bytes()).hexdigest()}}
    (HERE/'preflight_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'compiler_queries':len(cache),'exact_mean_comparisons':checked}),flush=True)


if __name__=='__main__':main()
