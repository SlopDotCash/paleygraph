#!/usr/bin/env python3
"""Independent symmetry and exact-oracle checks of the folded adapters."""
import sys
sys.dont_write_bytecode=True
from array import array
from collections import Counter
from copy import deepcopy
from fractions import Fraction
from hashlib import sha256
from itertools import permutations,product
from pathlib import Path
import importlib.util
import json
import random
import time

HERE=Path(__file__).resolve().parent
FOLD=HERE.parent/'trace_invariants'
sys.path.insert(0,str(FOLD))


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m


def primitive_orbit(p,t):
    """BFS from two generators, independently of the six closed formulas."""
    orbit={t};frontier=[t]
    while frontier:
        a=frontier.pop()
        for b in ((1-a)%p,pow(a,-1,p)):
            if b not in orbit:orbit.add(b);frontier.append(b)
    return orbit


def main():
    start=time.perf_counter()
    source_paths=[FOLD/'folded_moment.py',FOLD/'residue_moment.py',FOLD/'verify_symmetry.py']
    hashes={p.name:sha256(p.read_bytes()).hexdigest() for p in source_paths}
    folded=load('reviewed_fold',source_paths[0]);residue=load('reviewed_residue',source_paths[1])
    oracle=load('fold_literal_oracle',HERE/'direct_third_oracle.py')
    literal=load('fold_literal_column_oracle',HERE/'review_union_backend.py')
    derivation=load('fold_difference_derivation',HERE/'fold_review_derivation.py')
    counts=Counter()
    # All eight mutual sign patterns, not only the normalized four.
    for q in (13,17,49):
        for edges in product((-1,1),repeat=3):
            edge_matrix=[[0,edges[0],edges[1]],[edges[0],0,edges[2]],[edges[1],edges[2],0]]
            for tau in range(-(q-3),q-2):
                try:hist=oracle.triple_histogram(q,edges,tau)
                except AssertionError:continue
                theta=edges[0]*edges[1]*edges[2]*tau;mono=len(set(edges))==1
                assert theta%8==((15-q)%8 if mono else (q-7)%8)
                assert folded.theta_family(q,theta)==mono
                for order in permutations(range(3)):
                    for sign in (-1,1):
                        e=tuple(sign*edge_matrix[order[a]][order[b]] for a,b in ((0,1),(0,2),(1,2)))
                        changed=Counter({tuple(sign*v[i] for i in order):count for v,count in hist.items()})
                        assert changed==oracle.triple_histogram(q,e,sign*tau)
                        assert e[0]*e[1]*e[2]*sign*tau==theta and (len(set(e))==1)==mono
                        counts['literal_histogram_symmetry_checks']+=1
                counts['admissible_arbitrary_edge_histograms']+=1
    # Literal per-column multiplication also checks row/global-sign symmetry
    # without any conference-type or row-balance assumption.
    rng=random.Random(56913)
    for trial in range(2):
        columns=[tuple(rng.choice((-1,0,1)) for _ in range(3)) for _ in range(8)]
        expected=literal.literal_columns(columns,6,8)
        for order in permutations(range(3)):
            for sign in (-1,1):
                transformed=[tuple(sign*row[i] for i in order) for row in columns]
                assert literal.literal_columns(transformed,6,8)==expected
                counts['arbitrary_literal_column_coefficient_symmetries']+=1
    base=json.loads((HERE/'direct_third_oracle.json').read_text())
    inventories={r['graph']:r for r in base['normalized_joint_records']}
    cases=[r for r in base['records'] if r['q']>=13]
    for n in (6,7,8):cases+=json.loads((HERE/f'twins_n{n}_third_oracle.json').read_text())['records']
    for d in (0,2,4):
        direct=oracle.direct(oracle.prime(13),7,[d])[0]
        cases.append({'graph':'Paley13','q':13,**direct})
    reports=[]
    for row in cases:
        inventory=inventories[row['graph']]
        summary=folded.fold_inventory(inventory)
        assert folded.unfold_inventory(summary)==sorted(inventory['records'],key=lambda x:(x['edges'],x['tau']))
        counts['lossless_inventory_roundtrips']+=1
        powers_only=deepcopy(summary);del powers_only['records']
        for supplied in (summary,powers_only):
            out=folded.folded_global_third_moment(supplied,row['n'],row['degree'])
            assert out['third_moment']==row['moments'][2] and out['coefficient_evaluations']<=19
            counts['folded_complete_subset_comparisons']+=1
        if row['degree']==6:
            statistics=residue.seven_statistics(summary)
            for supplied in (summary,statistics):
                out=residue.seven_stat_global_third_moment(supplied,row['n'])
                assert out['third_moment']==row['moments'][2] and out['coefficient_evaluations']<=15
                counts['seven_stat_complete_subset_comparisons']+=1
        reports.append({key:row[key] for key in ('graph','q','n','degree')}|{'third_moment':row['moments'][2]})
    # Independently form K(z)^6/720 and compare every raw union coefficient.
    universal=[Fraction(1)]
    for _ in range(6):
        out=[Fraction() for _ in range(len(universal)+3)]
        for i,a in enumerate(universal):
            for j,b in ((1,1),(2,-3),(3,2)):out[i+j]+=a*b
        universal=out
    universal=[v/720 for v in universal]
    assert residue.universal_sixth_union()==universal
    hs=[];plans=[]
    for mono in (False,True):
        edges=(1,1,1) if mono else (1,1,-1);sign=1 if mono else -1
        chosen=[]
        for tau in range(-98,99):
            try:hist=oracle.triple_histogram(101,edges,tau)
            except AssertionError:continue
            chosen.append((sign*tau,hist))
        chosen=chosen[:8];plans.append((len(hs),[t for t,h in chosen]));hs.extend(h for t,h in chosen)
    raw,_=folded.compiler.union_coefficients_batch(hs,6,max_union=18)
    for index,nodes in plans:
        for k in range(19):
            values=[r['coefficients'][k] for r in raw[index:index+8]]
            poly=derivation.interpolation(nodes[:7],values[:7])
            assert poly[6]==universal[k] and poly[5]==0
            assert derivation.value(poly,nodes[-1])==values[-1]
            counts['universal_high_coefficients_checked']+=1
    # Independent parameter-orbit handling: complete tiny fields, then bounded
    # large-prime samples plus every saved exceptional orbit.
    saved_sym=json.loads((FOLD/'symmetry_validation.json').read_text())
    orbit_reports=[]
    for report in saved_sym['prime_orbits']:
        p=report['p'];path=HERE.parent/'trace_backend'/f'tau_p{p}.bin'
        data=array('i');data.frombytes(path.read_bytes())
        if sys.byteorder!='little':data.byteswap()
        assert data[0]==p and len(data)==p+1;taus=data[1:]
        sign=lambda x:0 if x%p==0 else 1 if pow(x%p,(p-1)//2,p)==1 else -1
        samples=list(range(2,p)) if p<=1297 else sorted({rng.randrange(2,p) for _ in range(128)})
        samples+= [t for o in report['special_orbits'] for t in o['parameters']]
        seen=set();census=Counter()
        for t in samples:
            orbit=primitive_orbit(p,t);assert len(orbit) in (2,3,6)
            theta=sign(t)*sign(t-1)*taus[t]
            for u in orbit:
                assert sign(u)*sign(u-1)*taus[u]==theta
                assert taus[(1-u)%p]==taus[u] and taus[pow(u,-1,p)]==sign(u)*taus[u]
                counts['prime_orbit_parameter_checks']+=1
            if p<=1297 and t not in seen:census[len(orbit)]+=1;seen.update(orbit)
        harmonic=primitive_orbit(p,p-1)
        assert harmonic=={p-1,2,pow(2,-1,p)} and len(harmonic)==3
        for o in report['special_orbits']:
            orbit=primitive_orbit(p,o['parameters'][0]);assert orbit==set(o['parameters'])
            assert o['stabilizer_size']*len(orbit)==6
            if len(orbit)==2:assert all((t*t-t+1)%p==0 for t in orbit)
        eq=1 if p%3==1 else 0
        expected={3:1,6:(p-5-2*eq)//6}
        if eq:expected[2]=1
        assert {int(k):v for k,v in report['orbit_census'].items()}==expected
        if p<=1297:assert census==Counter(expected) and len(seen)==p-2
        orbit_reports.append({'p':p,'complete_census_independently_replayed':p<=1297,'sampled_starting_parameters':len(samples),
                              'exact_orbit_count_formula':expected,'tau_sha256':sha256(path.read_bytes()).hexdigest()})
    # Corruptions cannot be passed off as a checked full folded inventory.
    good=folded.fold_inventory(inventories['Paley13']);bad_records=[]
    bad=deepcopy(good);bad['family_power_sums'][0]['powers'][3]+=1;bad_records.append(('false_attached_cube',bad))
    bad=deepcopy(good);bad['family_power_sums'][-1]=bad['family_power_sums'][0];bad_records.append(('duplicate_family',bad))
    bad=deepcopy(good);bad['family_power_sums'][0]['monochromatic']=0;bad_records.append(('nonboolean_family',bad))
    bad=deepcopy(good);bad['records'][0]['theta']+=800;bad_records.append(('inadmissible_same_residue_theta',bad))
    rejections=[]
    for name,bad in bad_records:
        for label,consumer in (('folded',lambda x:folded.folded_global_third_moment(x,7)),('seven_converter',residue.seven_statistics)):
            try:consumer(bad)
            except AssertionError:rejections.append(label+':'+name)
            else:raise AssertionError(('corrupt folded record accepted',label,name))
    for p in source_paths:assert sha256(p.read_bytes()).hexdigest()==hashes[p.name],'reviewed source changed'
    out={'status':'folded and seven-stat adapters pass independent exact checks','counts':dict(counts),'records':reports,
         'prime_orbit_reviews':orbit_reports,'malformed_records_rejected':rejections,
         'source_sha256':hashes|{Path(__file__).name:sha256(Path(__file__).read_bytes()).hexdigest()},
         'elapsed_seconds':time.perf_counter()-start,
         'scope':'Full theta histogram retains family through residue. Powers-only and seven-stat objects require certified provenance and are not graph-realizability certificates.'}
    (HERE/'fold_review.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('records','prime_orbit_reviews')},indent=2))


if __name__=='__main__':main()
