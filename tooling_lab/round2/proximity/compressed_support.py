#!/usr/bin/env python3
"""Complete agreement-support census via a certified interpolation-base cover.

Known ingredients: Turan covers, Lagrange interpolation, affine pencil counting.
The new local instrument retains exact support/track provenance and cover proofs.
"""
from collections import Counter
from itertools import combinations
from math import comb
import json
from pathlib import Path
import sys
import time

sys.dont_write_bytecode=True  # Keep imported round1 sources read-only, including caches.
ROUND1=Path(__file__).resolve().parents[2]/'proximity'
sys.path.insert(0,str(ROUND1))
from extension_field_probe import Field
from deformation_microscope import domain, is_prime


class PrimeField:
    def __init__(self,p):
        assert is_prime(p)
        self.p=self.q=p
        self.e=1
    def add(self,a,b): return (a+b)%self.p
    def sub(self,a,b): return (a-b)%self.p
    def mul(self,a,b): return a*b%self.p
    def inv(self,a): return pow(a,-1,self.p)
    def sum(self,xs): return sum(xs)%self.p
    def power(self,a,n): return pow(a,n,self.p)


def cover_plan(n,k,s,mode='partition'):
    assert 1<=k<=s<=n
    if mode=='full_k':
        return {'mode':mode,'n':n,'k':k,'s':s,'blocks':[list(range(n))],
                'basis_count':comb(n,k),'ignored_count':0,'outside_plus_block_capacity':k-1}
    best=None
    limit=(s-1)//(k-1) if k>1 else 1
    for m in range(1,limit+1):
        if mode=='anchor' and m>1:
            break
        t=n-s+m*(k-1)+1
        if t>n or t<m*k:
            continue
        sizes=[t//m+int(i<t%m) for i in range(m)]
        cost=sum(comb(z,k) for z in sizes)
        if best is None or cost<best[0]:
            best=cost,sizes,t
    assert best is not None
    cost,sizes,t=best
    blocks=[]
    offset=0
    for size in sizes:
        blocks.append(list(range(offset,offset+size)))
        offset+=size
    return {'mode':mode,'n':n,'k':k,'s':s,'blocks':blocks,'basis_count':cost,
            'ignored_count':n-t,'outside_plus_block_capacity':n-t+len(blocks)*(k-1)}


def verify_cover(plan):
    """Check a sufficient universal proof certificate, not sampled coverage."""
    n,k,s=plan['n'],plan['k'],plan['s']
    flat=sum(plan['blocks'],[])
    assert len(set(flat))==len(flat) and all(0<=i<n for i in flat)
    assert plan['ignored_count']==n-len(flat)
    capacity=n-len(flat)+sum(min(k-1,len(block)) for block in plan['blocks'])
    assert capacity<s
    assert plan['outside_plus_block_capacity']==capacity
    assert plan['basis_count']==sum(comb(len(block),k) for block in plan['blocks'])
    return True


def interpolation_matrix_reference(F,xs,dom):
    rows=[]
    den=[]
    for i,x in enumerate(xs):
        p=1
        for j,y in enumerate(xs):
            if i!=j:
                p=F.mul(p,F.sub(x,y))
        den.append(F.inv(p))
    for x in dom:
        if x in xs:
            j=xs.index(x)
            rows.append([int(i==j) for i in range(len(xs))])
            continue
        row=[]
        for i,w in enumerate(den):
            v=w
            for j,y in enumerate(xs):
                if i!=j:
                    v=F.mul(v,F.sub(x,y))
            row.append(v)
        rows.append(row)
    return rows


def interpolation_matrix(F,xs,dom):
    """Same Lagrange matrix, O(n*k+k^2) products via prefix/suffix exclusion.

    Added after the first measured run showed fewer bases could still lose
    to the old subset oracle because naive global interpolation cost O(n*k^2).
    """
    den=[]
    for i,x in enumerate(xs):
        v=1
        for j,y in enumerate(xs):
            if i!=j:v=F.mul(v,F.sub(x,y))
        den.append(F.inv(v))
    known={x:i for i,x in enumerate(xs)}
    rows=[]
    for x in dom:
        if x in known:
            rows.append([int(i==known[x]) for i in range(len(xs))])
            continue
        diffs=[F.sub(x,y) for y in xs]
        pre=[1]
        for d in diffs:pre.append(F.mul(pre[-1],d))
        suf=[1]*(len(xs)+1)
        for i in range(len(xs)-1,-1,-1):suf[i]=F.mul(suf[i+1],diffs[i])
        rows.append([F.mul(F.mul(pre[i],suf[i+1]),den[i]) for i in range(len(xs))])
    return rows


def decode(F,matrix,ys):
    return tuple(F.sum(F.mul(x,y) for x,y in zip(row,ys)) for row in matrix)


def run_census(F,dom,k,s,u0,u1,mode='partition',materialize_limit=1000):
    start=time.perf_counter()
    n=len(dom)
    assert len(set(dom))==n and len(u0)==len(u1)==n
    plan=cover_plan(n,k,s,mode)
    verify_cover(plan)
    tracks={}
    stats=Counter()
    nodes={}
    all_field_tracks=[]
    finite_provenance=[]
    for block in plan['blocks']:
        for base in combinations(block,k):
            stats['bases_processed']+=1
            matrix=interpolation_matrix(F,[dom[i] for i in base],dom)
            A=decode(F,matrix,[u0[i] for i in base])
            B=decode(F,matrix,[u1[i] for i in base])
            track=(A,B)
            if track in tracks:
                tracks[track]['base_multiplicity']+=1
                continue
            common=[]
            buckets={}
            for i in range(n):
                d0=F.sub(u0[i],A[i])
                d1=F.sub(u1[i],B[i])
                if d1==0:
                    if d0==0:
                        common.append(i)
                    continue
                z=F.mul(F.sub(0,d0),F.inv(d1))
                buckets.setdefault(z,[]).append(i)
            record={'first_base':list(base),'base_multiplicity':1,'intercept':list(A),
                    'slope':list(B),'common_coordinates':common}
            tracks[track]=record
            if len(common)>=s:
                record['all_field']=True
                all_field_tracks.append(record)
                if F.q<=materialize_limit:
                    chosen=range(F.q)
                else:
                    chosen=[]
            else:
                chosen=[z for z,indices in buckets.items() if len(indices)+len(common)>=s]
                if chosen:
                    record['all_field']=False
                    record['qualifying_scalar_buckets']={str(z):buckets[z] for z in sorted(chosen)}
                    finite_provenance.append(record)
            for z in chosen:
                cw=tuple(F.add(a,F.mul(z,b)) for a,b in zip(A,B))
                support=sorted(common+buckets.get(z,[]))
                assert len(support)>=s
                nodes[(z,cw)]={'scalar':z,'codeword':list(cw),'agreement_support':support}
    stats['distinct_tracks_examined']=len(tracks)
    stats['duplicate_tracks_skipped']=stats['bases_processed']-len(tracks)
    assert stats['bases_processed']==plan['basis_count']
    node_list=[node for _,node in sorted(nodes.items())]
    whole=bool(all_field_tracks)
    return {'field':{'p':F.p,'e':F.e,'q':F.q},'n':n,'k':k,'s':s,'domain':dom,
            'u0':u0,'u1':u1,'cover_certificate':plan,'cover_verified':True,
            'completeness':'all finite affine scalars and all qualifying codewords; all-field tracks may be symbolic',
            'finite_bad_scalar_count':F.q if whole else len({node['scalar'] for node in node_list}),
            'whole_field_correlated':whole,
            'bad_scalars':list(range(F.q)) if whole and F.q<=materialize_limit else
                (None if whole else sorted({node['scalar'] for node in node_list})),
            'nodes':node_list,'node_ledger_materialized':not whole or F.q<=materialize_limit,
            'finite_track_certificates':finite_provenance,'all_field_track_certificates':all_field_tracks,
            'stats':dict(stats),'elapsed_seconds':round(time.perf_counter()-start,6)}


def node_set(record):
    return {(node['scalar'],tuple(node['codeword']),tuple(node['agreement_support'] if 'agreement_support' in node else
            [i for i in range(record['n']) if node['agreement_mask']&(1<<i)]))
            for node in record['nodes']}


if __name__=='__main__':
    print(json.dumps([cover_plan(*cfg) for cfg in [(16,8,11),(20,10,14),(32,8,20),(64,4,40)]],indent=2))
