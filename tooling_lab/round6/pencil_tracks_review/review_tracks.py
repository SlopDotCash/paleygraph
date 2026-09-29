#!/usr/bin/env python3
"""Independent exact arithmetic review of whole-field polynomial-track lists."""
import sys
sys.dont_write_bytecode=True
from collections import Counter
from copy import deepcopy
from hashlib import sha256
from itertools import product
from pathlib import Path
import importlib.util
import json
import random
import time

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'pencil_tracks'/'pencil_tracks.py'
ORACLE=HERE.parents[1]/'round4'/'novelty'/'review_support_batches.py'


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m


def zero_set(F,a,b):
    """Intersect scalar constraints coordinate by coordinate, unlike production's pivot solve."""
    root=None
    for x,y in zip(a,b):
        if not y:
            if x:return ('empty',None)
        else:
            z=F.mul(F.neg(x),F.inv(y))
            if root is not None and root!=z:return ('empty',None)
            root=z
    return ('all',None) if root is None else ('singleton',root)


def direct_nodes(F,dom,k,s,u0,u1):
    nodes=set();anchors=[]
    for cs in product(range(F.q),repeat=k):
        cw=tuple(F.eval(cs,x) for x in dom)
        for z in range(F.q):
            support=tuple(i for i,v in enumerate(cw) if v==F.add(u0[i],F.mul(z,u1[i])))
            if len(support)>=s:nodes.add((z,tuple(cs),cw,support))
            if len(support)==len(dom):anchors.append((z,tuple(cs)))
    return nodes,anchors


def audit(F,cert,full=False):
    dom,u0,u1=cert['domain'],cert['u0'],cert['u1'];n=len(dom);k,s=cert['k'],cert['s']
    assert n==cert['n'] and len(set(dom))==n and 1<=k<=s<=n
    assert all(0<=x<F.q for x in dom+u0+u1)
    tracks=cert['tracks'];r=len(tracks)
    assert sorted(i for t in tracks for i in t['coordinates'])==list(range(n))
    assert len({(tuple(t['intercept_coefficients']),tuple(t['slope_coefficients'])) for t in tracks})==r
    cap=sum(min(len(t['coordinates']),k-1) for t in tracks)
    assert cap<s and cert['root_count_certificate']['cap']==cap
    events=set();ledger=[]
    for index,t in enumerate(tracks):
        assert t['id']==index
        a,b=t['intercept_coefficients'],t['slope_coefficients']
        assert len(a)==len(b)==k and all(0<=x<F.q for x in a+b)
        av=[F.eval(a,x) for x in dom];bv=[F.eval(b,x) for x in dom]
        assert av==t['intercept_codeword'] and bv==t['slope_codeword']
        assert all(av[i]==u0[i] and bv[i]==u1[i] for i in t['coordinates'])
        always=[];never=[];buckets={}
        for i in range(n):
            kind,z=zero_set(F,[F.sub(av[i],u0[i])],[F.sub(bv[i],u1[i])])
            if kind=='all':always.append(i)
            elif kind=='empty':never.append(i)
            else:buckets.setdefault(z,[]).append(i);events.add(z)
        assert t['always_coordinates']==always and t['never_coordinates']==never
        assert t['scalar_buckets']=={str(z):buckets[z] for z in sorted(buckets)}
        ledger.append((always,buckets))
    collisions=[];pair_sets=[]
    for i in range(r):
        for j in range(i+1,r):
            da=[F.sub(a,b) for a,b in zip(tracks[i]['intercept_coefficients'],tracks[j]['intercept_coefficients'])]
            db=[F.sub(a,b) for a,b in zip(tracks[i]['slope_coefficients'],tracks[j]['slope_coefficients'])]
            kind,z=zero_set(F,da,db);assert kind!='all'
            pair_sets.append((i,j,kind,z))
            if kind=='singleton':collisions.append({'tracks':[i,j],'scalar':z});events.add(z)
    assert cert['track_collisions']==collisions
    assert cert['generic_scalar_domain']=={'kind':'all_except','excluded':sorted(events),'cardinality':F.q-len(events)}
    generic=[j for j,(a,b) in enumerate(ledger) if len(a)>=s]
    assert cert['generic_qualifying_track_ids']==generic and cert['generic_list_size']==len(generic)
    finite=[]
    for z in sorted(events):
        bypoly={}
        for j,t in enumerate(tracks):
            support=sorted(ledger[j][0]+ledger[j][1].get(z,[]))
            if len(support)>=s:
                cs=tuple(F.add(a,F.mul(z,b)) for a,b in zip(t['intercept_coefficients'],t['slope_coefficients']))
                bypoly.setdefault(cs,[]).append(j)
        nodes=[{'representative_track':ids[0],'equal_track_ids':ids} for cs,ids in sorted(bypoly.items())]
        for ids in bypoly.values():
            assert all(sorted(ledger[j][0]+ledger[j][1].get(z,[]))==sorted(ledger[ids[0]][0]+ledger[ids[0]][1].get(z,[])) for j in ids)
        finite.append({'scalar':z,'nodes':nodes,'list_size':len(nodes)})
    assert cert['exceptional_scalars']==finite
    node_count=(F.q-len(events))*len(generic)+sum(x['list_size'] for x in finite)
    bad=(F.q-len(events) if generic else 0)+sum(bool(x['nodes']) for x in finite)
    assert cert['complete_symbolic_node_count']==node_count and cert['bad_scalar_count']==bad
    assert cert['every_scalar_has_a_qualifying_codeword']==(bad==F.q)
    no_anchor_proof=None
    if all(len(t['coordinates'])>=k for t in tracks) and r>=2:
        # An exact codeword must equal every track polynomial: >=k roots on
        # each assigned region. Distinct scalar collision constraints may
        # then certify nonexistence without scanning the field.
        possibilities=None
        for i,j,kind,z in pair_sets:
            allowed=set() if kind=='empty' else {z}
            possibilities=allowed if possibilities is None else possibilities&allowed
        if not possibilities:
            assert cert['code_anchor'] is None
            no_anchor_proof={'all_regions_have_at_least_k_coordinates':True,
                             'pair_collision_constraints':[[i,j,kind,z] for i,j,kind,z in pair_sets],
                             'common_scalar_solution_set':'empty'}
    if full:
        expected,anchors=direct_nodes(F,dom,k,s,u0,u1);actual=set()
        for z in range(F.q):
            bypoly={}
            for j,t in enumerate(tracks):
                cs=tuple(F.add(a,F.mul(z,b)) for a,b in zip(t['intercept_coefficients'],t['slope_coefficients']))
                cw=tuple(F.eval(cs,x) for x in dom)
                support=tuple(i for i,v in enumerate(cw) if v==F.add(u0[i],F.mul(z,u1[i])))
                if len(support)>=s:bypoly[cs]=(z,cs,cw,support)
            actual.update(bypoly.values())
        assert actual==expected and len(expected)==node_count
        if cert['code_anchor'] is None:assert not anchors
        else:
            a=cert['code_anchor'];assert (a['scalar'],tuple(a['coefficients'])) in anchors
    return {'q':F.q,'n':n,'k':k,'s':s,'tracks':r,'cap':cap,'exceptional_scalars':len(events),
            'collisions':len(collisions),'generic_list_size':len(generic),'nodes':node_count,
            'qualifying_scalars':bad,'complete_codebook_and_scalar_oracle':full,'no_anchor_certificate':no_anchor_proof}


def main():
    started=time.perf_counter();digest=sha256(SOURCE.read_bytes()).hexdigest()
    candidate=load('reviewed_pencil_tracks',SOURCE);oracle=load('independent_field',ORACLE)
    counts=Counter();rng=random.Random(61049);sample=None
    for k,s in ((1,2),(2,3)):
        F=oracle.OracleField(3);CF=candidate.PrimeField(3);dom=[0,1,2]
        for u0 in product(range(3),repeat=3):
            for u1 in product(range(3),repeat=3):
                if k==1:
                    tracks=[{'intercept_coefficients':[u0[i]],'slope_coefficients':[u1[i]],'coordinates':[i]} for i in range(3)]
                else:
                    tracks=[{'intercept_coefficients':[u0[0],F.sub(u0[1],u0[0])],
                             'slope_coefficients':[u1[0],F.sub(u1[1],u1[0])],'coordinates':[0,1]},
                            {'intercept_coefficients':[u0[2]],'slope_coefficients':[u1[2]],'coordinates':[2]}]
                cert=candidate.compile_certificate(CF,dom,k,s,list(u0),list(u1),tracks)
                row=audit(F,cert,True);counts['complete_tiny_oracles']+=1;counts['nodes_compared']+=row['nodes']
                counts['no_anchor_tiny_cases']+=cert['code_anchor'] is None
                counts['whole_field_tiny_cases']+=cert['every_scalar_has_a_qualifying_codeword']
                if cert['track_collisions'] and cert['exceptional_scalars']:sample=cert
    for q,n,r,number in ((5,5,2,64),(9,7,3,64)):
        F=oracle.OracleField(q);CF=candidate.Field(3,[1,0,1]) if q==9 else candidate.PrimeField(q);dom=list(range(n));k=2;s=r+1
        for trial in range(number):
            labels=[rng.randrange(r) for _ in dom]
            polys=[([rng.randrange(q) for _ in range(k)],[rng.randrange(q) for _ in range(k)]) for _ in range(r)]
            u0=[F.eval(polys[j][0],x) for x,j in zip(dom,labels)];u1=[F.eval(polys[j][1],x) for x,j in zip(dom,labels)]
            tracks=[{'intercept_coefficients':a,'slope_coefficients':b,'coordinates':[i for i,j in enumerate(labels) if j==z]} for z,(a,b) in enumerate(polys)]
            cert=candidate.compile_certificate(CF,dom,k,s,u0,u1,tracks);row=audit(F,cert,True)
            counts['complete_tiny_oracles']+=1;counts['nodes_compared']+=row['nodes'];counts[f'field_{q}_oracles']+=1
            counts['no_anchor_tiny_cases']+=cert['code_anchor'] is None
            counts['whole_field_tiny_cases']+=cert['every_scalar_has_a_qualifying_codeword']
    saved_path=SOURCE.parent/'results.json';saved=json.loads(saved_path.read_text());readbacks=[]
    for entry in saved['medium']+saved['large']:
        cert=entry['result']['certificate'];F=oracle.OracleField(cert['field']['q'])
        F.inv=lambda a,p=F.q:pow(a,-1,p)
        row=audit(F,cert);assert row['no_anchor_certificate'] is not None
        readbacks.append({'name':entry['name'],**row})
    extension=saved['saved_extension_certificate'];readbacks.append({'name':'saved GF9 certificate',**audit(oracle.OracleField(9),extension,True)})
    for entry in saved['iterations']:
        initial,final=entry['initial'],entry['iteration'];cert=final['certificate'];F=oracle.OracleField(cert['field']['q'])
        F.inv=lambda a,p=F.q:pow(a,-1,p)
        readbacks.append({'name':entry['name'],**audit(F,cert,cert['k']<=2)})
        assert (initial['domain'],initial['u0'],initial['u1'])==(final['domain'],final['u0'],final['u1'])
        if 'joint_root_count_cap' in initial:
            assert initial['joint_root_count_cap']>=initial['s'] and not initial['complete']
        else:assert not initial['complete'] and initial['status']=='discovery_incomplete'
    rejected=[]
    assert sample is not None
    mutations=[]
    bad=deepcopy(sample);bad['complete_symbolic_node_count']+=1;mutations.append(('false_total',bad))
    bad=deepcopy(sample);bad['track_collisions'].pop();mutations.append(('missing_collision',bad))
    bad=deepcopy(sample);bad['generic_scalar_domain']['excluded']=[];mutations.append(('false_generic_domain',bad))
    bad=deepcopy(sample);bad['tracks'][0]['always_coordinates']=[];mutations.append(('false_support_ledger',bad))
    bad=deepcopy(sample);bad['code_anchor']={'scalar':0,'coefficients':[0]};mutations.append(('false_anchor',bad))
    for name,bad in mutations:
        try:candidate.verify_certificate(candidate.PrimeField(3),bad)
        except AssertionError:rejected.append(name)
        else:raise AssertionError(('corrupt certificate accepted',name))
    assert sha256(SOURCE.read_bytes()).hexdigest()==digest
    out={'status':'whole-field track compiler independently verified','counts':dict(counts),'saved_readbacks':readbacks,
         'malformed_certificates_rejected':rejected,'candidate_sha256':digest,'results_sha256':sha256(saved_path.read_bytes()).hexdigest(),
         'reviewer_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'field_oracle_sha256':sha256(ORACLE.read_bytes()).hexdigest(),
         'elapsed_seconds':time.perf_counter()-started,
         'scope':'Tiny lists enumerate every polynomial and scalar independently. Larger cases independently reconstruct the complete event ledger and prove the generic complement symbolically; no large codebook or full-field scan.'}
    (HERE/'review_tracks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))


if __name__=='__main__':main()
