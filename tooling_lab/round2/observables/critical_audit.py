#!/usr/bin/env python3
"""Critical-scale audit of inexpensive quartic-incidence contractions.

Exact arithmetic features, shuffled-incidence negative controls, disjoint seeded
cohorts, and simple explicitly empirical linear prediction. No uniform theorem.
"""
from collections import Counter
from hashlib import sha256
from itertools import combinations,islice
import json
from math import comb,isqrt,sqrt
from pathlib import Path
import random
import sys,time
import numpy as np

HERE=Path(__file__).resolve().parent
LAB=HERE.parents[1]
sys.path.insert(0,str(LAB/'observables'))
from fiber_audit import direct_record


def chi_table(p):
    assert p>=3 and all(p%d for d in range(2,isqrt(p)+1))
    c=np.full(p,-1,dtype=np.int64);c[0]=0
    c[np.array([i*i%p for i in range(1,p)])]=1
    return c


def normalized(c,p):
    c=sorted(c);a=pow(c[1]-c[0],-1,p)
    return tuple(sorted((x-c[0])*a%p for x in c))


def cohort(p,n,count,seed,planted=False):
    rng=random.Random(seed);seen=set();chi=chi_table(p)
    assert count<=comb(p-2,n-2)
    while len(seen)<count:
        if planted:
            anchors=rng.sample(range(p),3)
            pool=[x for x in range(p) if all(chi[(x-a)%p]==1 for a in anchors)]
            if len(pool)<n:continue
            c=normalized(rng.sample(pool,n),p)
        else:c=(0,1)+tuple(sorted(rng.sample(range(2,p),n-2)))
        if c not in seen:seen.add(c);yield c


def contractions(graph):
    n=len(graph);w=np.array(graph,dtype=object)
    t=int(np.trace(w@w@w));v=w.sum(axis=1)
    s1=int(v.sum());s2=int(v@v);s3=int(v@w@v)
    centered=n**3*t-3*n*n*s3+3*n*s1*s2-s1**3
    # Exact check of n^3 tr((PWP)^3) without rational arithmetic.
    return t,centered


def features(sets,p,chi,seed):
    c=np.array(sets,dtype=np.int64);b,n=c.shape
    assert p*n**6*1000<np.iinfo(np.int64).max
    rows=chi[(np.arange(p)[None,None,:]-c[:,:,None])%p]
    f=rows.sum(axis=1);f2=f*f;m2=f2.sum(axis=1);m4=(f2*f2).sum(axis=1);m6=(f2*f2*f2).sum(axis=1)
    bd=np.take_along_axis(f,c,axis=1);b2=(bd*bd).sum(axis=1);b4=(bd**4).sum(axis=1)
    assert np.all(m2==n*(p-n))
    # Compute all elementary coefficients through six in O(6np), avoiding C(n,6).
    elementary=np.zeros((7,b,p),dtype=np.int64);elementary[0]=1
    for i in range(n):
        for j in range(min(6,i+1),0,-1):elementary[j]+=rows[:,i,:]*elementary[j-1]
    t6=elementary[6].sum(axis=1)
    qs=list(combinations(range(n),4))
    ks=np.array([(rows[:,i,:]*rows[:,j,:]*rows[:,k,:]*rows[:,l,:]).sum(axis=1)
                 for i,j,k,l in qs],dtype=np.int64).T
    t4=ks.sum(axis=1)
    a0=15*n**3-30*n*n+16*n;a2=90*n*n-300*n+272;a4=360*n-960
    assert np.array_equal(m6,p*a0-a2*comb(n,2)+a4*t4+720*t6-15*b4-15*b2-n)
    pairs=((c[:,:,None]+c[:,None,:])%p).reshape(b,-1)
    paircounts=np.zeros((b,p),dtype=np.int64);np.add.at(paircounts,(np.arange(b)[:,None],pairs),1)
    energies=(paircounts*paircounts).sum(axis=1)
    ix=[[k for k,q in enumerate(qs) if i not in q and j not in q] for i,j in combinations(range(n),2)]
    rng=np.random.default_rng(seed)
    out=[]
    for z in range(b):
        variants=[ks[z]]+[rng.permutation(ks[z]) for _ in range(3)]
        vals=[]
        for deck in variants:
            g=[[0]*n for _ in range(n)]
            for (i,j),inds in zip(combinations(range(n),2),ix):g[i][j]=g[j][i]=sum(int(deck[t]) for t in inds)
            vals.append(contractions(g))
        out.append({'C':sets[z],'m2':int(m2[z]),'m4':int(m4[z]),'m6':int(m6[z]),
                    'b2':int(b2[z]),'b4':int(b4[z]),'energy':int(energies[z]),'t4':int(t4[z]),'t6':int(t6[z]),
                    'quartic_deck':ks[z].tolist(),'triangle':vals[0][0],'centered_triangle_numerator':vals[0][1],
                    'shuffled_triangle':[x[0] for x in vals[1:]],'shuffled_centered':[x[1] for x in vals[1:]]})
    return out


def generate(p,n,count,seed,planted=False):
    chi=chi_table(p);source=iter(cohort(p,n,count,seed,planted));records=[];batchno=0
    while batch:=list(islice(source,64)):
        records.extend(features(batch,p,chi,seed+batchno));batchno+=1
    for r in records[:2]:
        check=direct_record(r['C'],p)
        for name in ['m2','m4','m6','b2','b4','energy']:assert r[name]==check[name]
        assert r['t6']==check['t6_excluded']
    return records


def design(records,p,n,mode):
    x=[];y=[]
    scale=(p*comb(n-2,4))**1.5*n**3
    for r in records:
        row=[1,r['m4']/(p*n*n),r['b2']/n**3,r['b4']/n**5,r['energy']/n**2]
        if mode=='actual':row += [r['triangle']/scale,r['centered_triangle_numerator']/(n**3*scale)]
        if mode=='shuffled':row += [r['shuffled_triangle'][0]/scale,r['shuffled_centered'][0]/(n**3*scale)]
        x.append(row);y.append(r['t6']/sqrt(p*comb(n,6)))
    return np.array(x),np.array(y)


def metrics(actual,pred):
    err=np.mean((actual-pred)**2);base=np.mean((actual-np.mean(actual))**2)
    return {'mse':float(err),'r2_against_cohort_mean':float(1-err/base),
            'max_abs_error':float(max(abs(actual-pred)))}


def correlations(records):
    y=np.array([r['t6'] for r in records]);a=np.array([r['triangle'] for r in records]);b=np.array([r['centered_triangle_numerator'] for r in records])
    shuffle=[float(np.corrcoef(y,[r['shuffled_triangle'][j] for r in records])[0,1]) for j in range(3)]
    return {'target_raw_triangle':float(np.corrcoef(y,a)[0,1]),
            'target_centered_triangle':float(np.corrcoef(y,b)[0,1]),'target_shuffled_triangles':shuffle}


def compact_fibers(records):
    groups={}
    for r in records:groups.setdefault(tuple(r[k] for k in ['m4','b2','b4','energy']),[]).append(r)
    ambiguous=[rs for rs in groups.values() if len({r['m6'] for r in rs})>1]
    w=max(ambiguous,key=lambda rs:max(r['m6'] for r in rs)-min(r['m6'] for r in rs)) if ambiguous else None
    return {'fibers':len(groups),'singleton_records':sum(len(rs)==1 for rs in groups.values()),
            'ambiguous_fibers':len(ambiguous),'witness':None if w is None else [min(w,key=lambda r:r['m6']),max(w,key=lambda r:r['m6'])]}


def main():
    started=time.monotonic();configs=[(61,6,2048),(1297,6,8192),(2437,7,4096),(4129,8,2048)]
    allcases=[];models=None
    for p,n,size in configs:
        rs=generate(p,n,size,200000+p);train,test=rs[:size//2],rs[size//2:]
        overlap=set(r['C'] for r in train)&set(r['C'] for r in test);assert not overlap
        planted=generate(p,n,128,300000+p,True)
        local={};transport={}
        for mode in ['baseline','actual','shuffled']:
            xt,yt=design(train,p,n,mode);beta=np.linalg.lstsq(xt,yt,rcond=None)[0]
            local[mode]={}
            for name,coh in [('random_holdout',test),('three_anchor_holdout',planted)]:
                x,y=design(coh,p,n,mode);local[mode][name]=metrics(y,x@beta)
            if models is not None:
                x,y=design(test,p,n,mode);transport[mode]=metrics(y,x@models[mode])
            if models is None:pass
        if models is None:
            models={mode:np.linalg.lstsq(*design(train,p,n,mode),rcond=None)[0] for mode in ['baseline','actual','shuffled']}
        result={'p':p,'n':n,'random_count':size,'planted_count':128,'critical_size':n**4<=p,
                'low_cost_fibers':compact_fibers(rs),'correlations':correlations(rs),
                'local_fit':local,'toy_model_transferred':transport,
                'target_range':[min(r['t6'] for r in rs),max(r['t6'] for r in rs)],
                'maximum_m6_over_gaussian':max(r['m6'] for r in rs)/(15*p*n**3)}
        (HERE/f'data_{p}.json').write_text(json.dumps({'p':p,'n':n,'random':rs,'planted':planted},separators=(',',':'))+'\n')
        allcases.append(result)
        print(json.dumps({'p':p,'count':size,'correlations':result['correlations'],'random_holdout_r2':{k:v['random_holdout']['r2_against_cohort_mean'] for k,v in local.items()}}),flush=True)
    report={'status':'finite evidence about candidate informativeness; empirical regression is not a certificate or uniform bound',
            'cases':allcases,'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
            'elapsed_seconds':round(time.monotonic()-started,3)}
    (HERE/'results.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
