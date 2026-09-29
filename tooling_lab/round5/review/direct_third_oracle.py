#!/usr/bin/env python3
"""Independent complete-set T_d third moments and normalized triple records."""
import sys
sys.dont_write_bytecode=True
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations,product
from math import comb
from pathlib import Path
import argparse
import json
import subprocess
import time

HERE=Path(__file__).resolve().parent


def choose(n,k):return comb(n,k) if 0<=k<=n else 0


def prime(p):
    sign=[0 if x==0 else (1 if pow(x,(p-1)//2,p)==1 else -1) for x in range(p)]
    return [[sign[(x-y)%p] for y in range(p)] for x in range(p)]


def mul(x,y):
    return ((x%7)*(y%7)-(x//7)*(y//7))%7+7*(((x%7)*(y//7)+(x//7)*(y%7))%7)


def sub(x,y):return (x%7-y%7)%7+7*((x//7-y//7)%7)


def power(x,k):
    ans=1
    for _ in range(k):ans=mul(ans,x)
    return ans


def twins():
    powers=[power(9,j) for j in range(48)]
    assert len(set(powers))==48 and power(9,48)==1
    log={x:j for j,x in enumerate(powers)}
    result={}
    for name,classes in (('Paley49',(0,2)),('Peisert49',(0,1))):
        sign=[0]+[1 if log[x]%4 in classes else -1 for x in range(1,49)]
        result[name]=[[sign[sub(x,y)] for y in range(49)] for x in range(49)]
    return result


def direct(S,n,degrees):
    started=time.perf_counter();q=len(S)
    tables={d:{(a,b):sum((-1)**j*choose(b,j)*choose(a,d-j) for j in range(d+1))
               for a in range(n+1) for b in range(n+1-a)} for d in degrees}
    totals={d:[0,0,0] for d in degrees};count=0
    for C in combinations(range(q),n):
        counts=[]
        for row in S:
            signs=[row[x] for x in C]
            counts.append((signs.count(1),signs.count(-1)))
        for d in degrees:
            value=sum(tables[d][a,b] for a,b in counts)
            totals[d][0]+=value;totals[d][1]+=value**2;totals[d][2]+=value**3
        count+=1
    assert count==comb(q,n)
    return [{'n':n,'degree':d,'subset_count':count,
             'power_sums':sums,'moments':[[Fraction(s,count).numerator,Fraction(s,count).denominator] for s in sums],
             'method':'all subsets directly enumerated','elapsed_seconds_shared':time.perf_counter()-started}
            for d,sums in totals.items()]


def triple_histogram(q,edges,tau):
    a,b,c=edges
    out=Counter(((0,a,b),(a,0,c),(b,c,0)))
    moments=(q-3,-a-b,-a-c,-b-c,-1-b*c,-1-a*c,-1-a*b,tau)
    for x,y,z in product((-1,1),repeat=3):
        basis=(1,x,y,z,x*y,x*z,y*z,x*y*z)
        total=sum(v*w for v,w in zip(moments,basis))
        assert total>=0 and total%8==0
        if total:out[x,y,z]=total//8
    return out


def normalized_record(name,S):
    q=len(S);hist=Counter();checked=0
    for t in range(2,q):
        literal=Counter((S[0][x],S[1][x],S[t][x]) for x in range(q))
        tau=sum(a*b*c*count for (a,b,c),count in literal.items())
        edges=S[0][1],S[0][t],S[1][t]
        assert literal==triple_histogram(q,edges,tau)
        hist[edges,tau]+=1;checked+=1
    # Every unordered row triple gets a literal validation, not only normalized ones.
    for a,b,c in combinations(range(q),3):
        literal=Counter((S[a][x],S[b][x],S[c][x]) for x in range(q))
        tau=sum(x*y*z*count for (x,y,z),count in literal.items())
        assert literal==triple_histogram(q,(S[a][b],S[a][c],S[b][c]),tau)
        checked+=1
    return {'graph':name,'q':q,'normalization':'first ordered pair maps to 0,1, up to common sign; semilinear maps for Peisert49 verified in frozen preflight',
            'records':[{'edges':list(edges),'tau':tau,'count':count} for (edges,tau),count in sorted(hist.items())],
            'literal_triple_type_checks':checked,'normalized_parameter_count':sum(hist.values())}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--twins-n',type=int,choices=[6,7,8]);args=parser.parse_args()
    started=time.perf_counter();T=twins()
    if args.twins_n is not None:
        n=args.twins_n;binary=HERE/'origin_third_census'
        subprocess.run(['clang++','-std=c++17','-O3',str(HERE/'origin_third_census.cpp'),'-o',str(binary)],check=True)
        payload='\n'.join(' '.join(map(str,row)) for S in T.values() for row in S)+'\n'
        raw=subprocess.check_output([str(binary),str(n)],input=payload,text=True)
        lines=raw.splitlines();count=int(lines[0].split()[1]);assert count==comb(48,n-1)
        records=[]
        for name,line in zip(T,lines[1:3]):
            graph,first,second,third,minimum,maximum=map(int,line.split())
            sums=[first,second,third]
            assert all(s*49%n==0 for s in sums)
            global_sums=[s*49//n for s in sums]
            total=comb(49,n)
            records.append({'graph':name,'q':49,'n':n,'degree':6,'subset_count':total,
                            'origin_subset_count':count,'origin_power_sums':sums,'power_sums':global_sums,
                            'moments':[[Fraction(s,total).numerator,Fraction(s,total).denominator] for s in global_sums],
                            'minimum':minimum,'maximum':maximum,
                            'method':'all origin-containing subsets directly evaluated, global translation double count49/n'})
        out={'status':'complete direct finite subset oracle','records':records,
             'cpp_elapsed_seconds':float(lines[3].split()[1]),'elapsed_seconds':time.perf_counter()-started,
             'source_sha256':{name:sha256((HERE/name).read_bytes()).hexdigest() for name in ('direct_third_oracle.py','origin_third_census.cpp')}}
        (HERE/f'twins_n{n}_third_oracle.json').write_text(json.dumps(out,indent=2)+'\n')
        print(json.dumps(out,indent=2));return
    records=[];inventories=[]
    for p in (13,17):
        S=prime(p);inventories.append(normalized_record(f'Paley{p}',S))
        for n in (6,7,8):
            for row in direct(S,n,[6]):records.append({'graph':f'Paley{p}','q':p,**row})
    for n in range(10):
        # Paley9 is generated independently in the quadratic extension F3[i].
        def subtract(x,y):return (x%3-y%3)%3+3*((x//3-y//3)%3)
        def multiply(x,y):return ((x%3)*(y%3)-(x//3)*(y//3))%3+3*(((x%3)*(y//3)+(x//3)*(y%3))%3)
        squares={multiply(x,x) for x in range(1,9)}
        S=[[0 if x==y else (1 if subtract(x,y) in squares else -1) for y in range(9)] for x in range(9)]
        for row in direct(S,n,list(range(n+1))):records.append({'graph':'Paley9','q':9,**row})
    for name,S in T.items():inventories.append(normalized_record(name,S))
    preflight=HERE.parent/'preflight'/'third_moment_preflight.json'
    out={'status':'complete independent third-moment oracles','date':'2026-09-05','records':records,'normalized_joint_records':inventories,
         'preflight_normalization_evidence_sha256':sha256(preflight.read_bytes()).hexdigest(),
         'source_sha256':{Path(__file__).name:sha256(Path(__file__).read_bytes()).hexdigest()},
         'elapsed_seconds':time.perf_counter()-started}
    (HERE/'direct_third_oracle.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'records':len(records),'normalized_type_checks':sum(r['literal_triple_type_checks'] for r in inventories),'seconds':out['elapsed_seconds']},indent=2))


if __name__=='__main__':main()
