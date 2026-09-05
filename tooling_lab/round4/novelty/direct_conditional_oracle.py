#!/usr/bin/env python3
"""Direct conditional T6 moments; no imports of any moment or cell compiler."""
from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import comb
from pathlib import Path
import json
import subprocess

HERE=Path(__file__).resolve().parent


def choose(n,k):return comb(n,k) if 0<=k<=n else 0


def six_table(n):
    return {(a,b):sum((-1)**j*choose(b,j)*choose(a,6-j) for j in range(7))
            for a in range(n+1) for b in range(n+1-a)}


def target(S,C,table):
    out=0
    for row in S:
        a=b=0
        for x in C:
            if row[x]==1:a+=1
            elif row[x]==-1:b+=1
        out+=table[a,b]
    return out


def record(marks,count,total,squares):
    mean,second=Fraction(total,count),Fraction(squares,count)
    return {'marks':list(marks),'count':count,'sum_T6':total,'sum_T6_squared':squares,
            'mean':[mean.numerator,mean.denominator],
            'second_moment':[second.numerator,second.denominator],
            'variance':[(second-mean*mean).numerator,(second-mean*mean).denominator]}


def prime_case(p,n):
    h=[0 if x==0 else (1 if pow(x,(p-1)//2,p)==1 else -1) for x in range(p)]
    S=[[h[(y-x)%p] for x in range(p)] for y in range(p)]
    table=six_table(n);sums=defaultdict(lambda:[0,0,0]);seen=0
    for C in combinations(range(p),n):
        value=target(S,C,table);seen+=1
        for m in range(4):
            for M in combinations(C,m):
                a=sums[M];a[0]+=1;a[1]+=value;a[2]+=value*value
    rows=[]
    for M,(count,total,squares) in sorted(sums.items(),key=lambda t:(len(t[0]),t[0])):
        assert count==comb(p-len(M),n-len(M))
        rows.append(record(M,count,total,squares))
    baseline=rows[0]
    for r in rows:
        if len(r['marks'])<=2:
            assert r['mean']==baseline['mean'] and r['second_moment']==baseline['second_moment']
    boundary=list(range(n));value=target(S,boundary,table)
    return {'graph':f'Paley{p}','q':p,'n':n,'enumerated_subsets':seen,'records':rows,
            'full_mark_boundary':record(boundary,1,value,value*value)}


def mul(x,y):
    a,b,c,d=x%7,x//7,y%7,y//7
    return (a*c-b*d)%7+7*((a*d+b*c)%7)


def sub(x,y):return (x%7-y%7)%7+7*((x//7-y//7)%7)


def twins():
    powers=[1]
    for _ in range(47):powers.append(mul(powers[-1],9))
    assert len(set(powers))==48 and mul(powers[-1],9)==1
    connections=[{powers[k] for k in range(48) if k%4 in kinds} for kinds in ((0,2),(0,1))]
    matrices=[[[0 if x==y else (1 if sub(x,y) in conn else -1) for y in range(49)]
               for x in range(49)] for conn in connections]
    marks=[[],[0],[0,1],[0,2],[0,7],[0,1,2],[0,1,7],[0,1,8],
           [0,7,14],[0,1,41],[0,1,47],[0,7,41]]
    payload='\n'.join(' '.join(map(str,row)) for S in matrices for row in S)+'\n'
    payload+=str(len(marks))+'\n'+' '.join(str(sum(1<<x for x in M)) for M in marks)+'\n'
    binary=HERE/'origin_marked_oracle'
    subprocess.run(['clang++','-O3','-std=c++17',str(HERE/'origin_marked_oracle.cpp'),'-o',str(binary)],check=True)
    output=subprocess.check_output([str(binary)],input=payload,text=True)
    rows=[[],[]]
    for line in output.splitlines():
        i,count,sp,pp,sq,qq=map(int,line.split());M=marks[i]
        if not M:
            assert all(x*49%6==0 for x in (count,sp,pp,sq,qq))
            count,sp,pp,sq,qq=[x*49//6 for x in (count,sp,pp,sq,qq)]
        assert count==comb(49-len(M),6-len(M))
        rows[0].append(record(M,count,sp,pp));rows[1].append(record(M,count,sq,qq))
    out=[]
    for name,recs in zip(('Paley49','Peisert49'),rows):
        for r in recs:
            if len(r['marks'])<=2:assert r['mean']==recs[0]['mean'] and r['second_moment']==recs[0]['second_moment']
        out.append({'graph':name,'q':49,'n':6,'enumerated_origin_subsets':comb(48,5),
                    'records':recs,'baseline_note':'Empty mark moments recovered by translation double counting; all other marks contain0.'})
    selected=[[0,1,2],[0,1,7],[0,7,14]]
    payload='\n'.join(' '.join(map(str,row)) for S in matrices for row in S)+'\n'
    payload+=str(len(selected))+'\n'+'\n'.join('7 3 '+' '.join(map(str,M)) for M in selected)+'\n'
    binary=HERE/'selected_marked_oracle'
    subprocess.run(['clang++','-O3','-std=c++17',str(HERE/'selected_marked_oracle.cpp'),'-o',str(binary)],check=True)
    raw=subprocess.check_output([str(binary)],input=payload,text=True)
    selected_rows=[[],[]]
    for line in raw.splitlines():
        i,count,sp,pp,sq,qq=map(int,line.split());assert count==comb(46,4)
        selected_rows[0].append(record(selected[i],count,sp,pp))
        selected_rows[1].append(record(selected[i],count,sq,qq))
    for name,recs in zip(('Paley49','Peisert49'),selected_rows):
        out.append({'graph':name,'q':49,'n':7,'records':recs,
                    'direct_completions_per_mark':comb(46,4),'scope':'generic n>d case, all completions of each selected mark'})
    return out


def main():
    cases=[prime_case(p,n) for p,n in ((13,6),(13,7),(17,6),(17,8))]
    cases.extend(twins())
    out={'date':'2026-09-05','status':'direct exact conditional oracle, independent of compiler',
         'cases':cases,'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest()
                                      for n in ('direct_conditional_oracle.py','origin_marked_oracle.cpp','selected_marked_oracle.cpp')}}
    (HERE/'direct_conditional_oracle.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps([{'graph':c['graph'],'n':c['n'],'conditional_records':len(c['records']),
                       'distinct_three_mark_mean_second_pairs':len({(tuple(r['mean']),tuple(r['second_moment']))
                                                                  for r in c['records'] if len(r['marks'])==3})}
                      for c in cases],indent=2))


if __name__=='__main__':main()
