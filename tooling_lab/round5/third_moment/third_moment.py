#!/usr/bin/env python3
"""Exact global third moments from normalized conference row-triple inventory.

The normalized API requires even degree and a separately verified transitive
normalization by automorphisms/complementation. The integer polynomial backend
itself accepts arbitrary nonnegative sign-type histograms and degrees0..6.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from math import comb
from pathlib import Path
import subprocess
import time

HERE=Path(__file__).resolve().parent
SIGNS=tuple(product((-1,1),repeat=3))


def encode(value):
    return [value.numerator,value.denominator]


def compile_backend():
    source,binary=HERE/'union_coefficients.cpp',HERE/'union_coefficients'
    if not binary.exists() or binary.stat().st_mtime_ns<source.stat().st_mtime_ns:
        subprocess.run(['clang++','-O3','-std=c++17','-I/opt/homebrew/include',str(source),'-o',str(binary)],check=True)
    return binary


def union_coefficients_batch(histograms,degree=6,max_union=None):
    assert 0<=degree<=6 and histograms
    binary=compile_backend()
    text=[str(len(histograms))]
    for hist in histograms:
        assert all(len(pattern)==3 and set(pattern)<={-1,0,1} and isinstance(count,int) and count>=0
                   for pattern,count in hist.items())
        q=sum(hist.values())
        limit=min(3*degree,q) if max_union is None else min(max_union,3*degree,q)
        entries=[(*pattern,count) for pattern,count in sorted(hist.items()) if count]
        text.append(f'{degree} {limit} {len(entries)}')
        text.extend(' '.join(map(str,row)) for row in entries)
    started=time.perf_counter()
    output=subprocess.check_output([str(binary)],input='\n'.join(text)+'\n',text=True)
    elapsed=time.perf_counter()-started
    rows=[json.loads(line) for line in output.splitlines()]
    assert len(rows)==len(histograms)
    for i,row in enumerate(rows):
        assert row['case_id']==i and row['degree']==degree and len(row['coefficients'])==row['max_union']+1
    return rows,elapsed


def inclusion_moment(coefficients,q,n):
    assert 0<=n<=q
    probability=Fraction(1)
    value=Fraction()
    for k,coefficient in enumerate(coefficients):
        if k>n:break
        if k:probability*=Fraction(n-k+1,q-k+1)
        value+=coefficient*probability
    return value


def distinct_histogram(q,edges,tau):
    assert isinstance(q,int) and q>=5 and q%4==1 and len(edges)==3 and set(edges)<={-1,1}
    e01,e02,e12=edges
    out=Counter([(0,e01,e02),(e01,0,e12),(e02,e12,0)])
    linear=(-e01-e02,-e01-e12,-e02-e12)
    paired=(-1-e02*e12,-1-e01*e12,-1-e01*e02)
    for a,b,c in SIGNS:
        numerator=q-3+a*linear[0]+b*linear[1]+c*linear[2]+a*b*paired[0]+a*c*paired[1]+b*c*paired[2]+a*b*c*tau
        assert numerator>=0 and numerator%8==0
        if numerator:out[a,b,c]+=numerator//8
    assert sum(out.values())==q
    return out


def pair_histogram(q,sign):
    assert sign in (-1,1)
    return {(0,sign):1,(sign,0):1,(1,1):(q-3-2*sign)//4,
            (-1,-1):(q-3+2*sign)//4,(1,-1):(q-1)//4,(-1,1):(q-1)//4}


def repeated_histograms(q):
    same={(0,0,0):1,(1,1,1):(q-1)//2,(-1,-1,-1):(q-1)//2}
    twice=[{(a,a,b):count for (a,b),count in pair_histogram(q,s).items() if count} for s in (-1,1)]
    return [same,*twice]


def admissible_tau_nodes(q,edges):
    """Integral nonnegative synthetic inventories; no graph-realizability claim."""
    e01,e02,e12=edges
    linear=(-e01-e02,-e01-e12,-e02-e12)
    paired=(-1-e02*e12,-1-e01*e12,-1-e01*e02)
    lower,upper=-(q-3),q-3
    residue=None
    for a,b,c in SIGNS:
        constant=q-3+a*linear[0]+b*linear[1]+c*linear[2]+a*b*paired[0]+a*c*paired[1]+b*c*paired[2]
        slope=a*b*c
        if slope==1:lower=max(lower,-constant)
        else:upper=min(upper,constant)
        this=(-constant*slope)%8
        assert residue is None or residue==this
        residue=this
    first=lower+(residue-lower)%8
    if first>upper:return []
    count=(upper-first)//8+1
    if count<=8:
        return list(range(first,upper+1,8))
    # Eight central admissible values: seven interpolation nodes and one
    # additional implementation check. Degree is established algebraically.
    midpoint=min(max(((-first)//8)-3,0),count-8)
    return [first+8*(midpoint+j) for j in range(8)]


def lagrange_weights(nodes,x):
    result=[]
    for i,node in enumerate(nodes):
        value=Fraction(1)
        for j,other in enumerate(nodes):
            if i!=j:value*=Fraction(x-other,node-other)
        result.append(value)
    return result


def validate_inventory(record):
    q=record.get('q',record.get('p'))
    assert isinstance(q,int) and q>=5 and q%4==1
    records=record.get('records',record.get('joint_edges_tau_histogram'))
    assert isinstance(records,list) and records
    seen=set()
    for row in records:
        edges,tau,count=tuple(row['edges']),row['tau'],row['count']
        assert edges[0]==1 and set(edges)<={-1,1} and len(edges)==3
        assert isinstance(tau,int) and isinstance(count,int) and count>0 and (edges,tau) not in seen
        seen.add((edges,tau))
        distinct_histogram(q,edges,tau)
    assert sum(x['count'] for x in records)==q-2
    class_counts=Counter()
    for row in records:class_counts[tuple(row['edges'])]+=row['count']
    for a,b in product((-1,1),repeat=2):
        assert class_counts[1,a,b]==((q-5)//4 if a==b==1 else (q-1)//4)
    # These identities follow from S1=0 and S²=qI-J for the actual
    # normalized pair, with t=0,1 removed; they are inexpensive diagnostics.
    assert sum(row['count']*row['tau'] for row in records)==2
    assert sum(row['count']*row['tau']**2 for row in records)==(q-3)*(q+1)
    for i in (1,2):
        assert sum(row['count']*row['edges'][i]*row['tau'] for row in records)==2
    if 'edge_power_sums' in record:
        power_classes=set()
        for row in record['edge_power_sums']:
            edges=tuple(row['edges'])
            assert edges not in power_classes
            power_classes.add(edges)
            assert row['powers']==[sum(v['count']*v['tau']**j for v in records if tuple(v['edges'])==edges)
                                    for j in range(7)]
        assert power_classes=={(1,a,b) for a,b in product((-1,1),repeat=2)}
    if 'ordered_distinct_row_weight' in record:
        assert record['ordered_distinct_row_weight']==q*(q-1)
    return q,records


def global_third_moment(record,n,degree=6,method='interpolate'):
    started=time.perf_counter()
    q,records=validate_inventory(record)
    assert degree in (0,2,4,6) and degree<=n<=q
    assert method in ('direct_inventory','interpolate')
    classes=defaultdict(list)
    for row in records:classes[tuple(row['edges'])].append(row)
    histograms=repeated_histograms(q)
    plans=[]
    for edges,rows in sorted(classes.items()):
        nodes=admissible_tau_nodes(q,edges)
        interpolate=method=='interpolate' and len(nodes)>=degree+1
        # The tau-degree bound is degree, by consuming >=1 of each target
        # variable per tau-dependent logarithm factor.
        if interpolate:
            nodes=nodes[:degree+2]
        else:
            nodes=sorted(row['tau'] for row in rows)
        index=len(histograms)
        histograms.extend(distinct_histogram(q,edges,tau) for tau in nodes)
        plans.append({'edges':edges,'rows':rows,'nodes':nodes,'index':index,'interpolate':interpolate})
    raw,backend_elapsed=union_coefficients_batch(histograms,degree,max_union=min(n,3*degree))
    moments=[inclusion_moment(row['coefficients'],q,n) for row in raw]
    assert moments[1]==moments[2]  # even degree, global sign-complement symmetry
    repeated=q*moments[0]+3*q*(q-1)*moments[1]
    distinct=Fraction()
    plan_records=[]
    for plan in plans:
        nodes=plan['nodes']
        values=moments[plan['index']:plan['index']+len(nodes)]
        if plan['interpolate']:
            basis_nodes=nodes[:degree+1]
            basis_values=values[:degree+1]
            checked_extra=False
            if len(nodes)>degree+1:
                assert sum((w*v for w,v in zip(lagrange_weights(basis_nodes,nodes[-1]),basis_values)),Fraction())==values[-1]
                checked_extra=True
            # Aggregate interpolation weights over the actual inventory.
            # This is equivalent to its first degree+1 power sums, and uses
            # actual counts, not the synthetic interpolation-node frequencies.
            aggregate_weights=[Fraction() for _ in basis_nodes]
            for row in plan['rows']:
                for i,w in enumerate(lagrange_weights(basis_nodes,row['tau'])):
                    aggregate_weights[i]+=row['count']*w
            contribution=sum((w*v for w,v in zip(aggregate_weights,basis_values)),Fraction())
        else:
            actual_values=dict(zip(nodes,values))
            contribution=sum((row['count']*actual_values[row['tau']] for row in plan['rows']),Fraction())
            checked_extra=False
            aggregate_weights=[]
        distinct+=q*(q-1)*contribution
        plan_records.append({'edges':list(plan['edges']),'method':'proved-degree interpolation' if plan['interpolate'] else 'direct actual inventory',
                             'synthetic_nodes_claimed_realizable':False,'nodes':nodes,
                             'row_triple_product_moments':[encode(v) for v in values],
                             'aggregate_interpolation_weights':[encode(v) for v in aggregate_weights],
                             'extra_node_checked':checked_extra,
                             'normalized_class_sum':encode(contribution)})
    total=repeated+distinct
    return {'q':q,'n':n,'degree':degree,'third_moment':encode(total),'third_moment_float':float(total),
            'all_equal_rows_contribution':encode(q*moments[0]),
            'exactly_two_equal_rows_contribution':encode(3*q*(q-1)*moments[1]),
            'distinct_rows_contribution':encode(distinct),
            'normalized_parameter_count':q-2,'normalization_contract':'supplied inventory requires verified automorphism/complement transitivity; current intended instances are prime Paley and checked GF49 twins',
            'class_plans':plan_records,'coefficient_evaluations':len(histograms),
            'backend_integer_products':sum(row['integer_products'] for row in raw),
            'maximum_intermediate_bits':max(row['maximum_intermediate_bits'] for row in raw),
            'backend_elapsed_seconds':backend_elapsed,'elapsed_seconds':time.perf_counter()-started,
            'source_sha256':{name:sha256((HERE/name).read_bytes()).hexdigest()
                             for name in ('third_moment.py','union_coefficients.cpp','union_coefficients')}}


def main():
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inventory',type=Path)
    parser.add_argument('--n',type=int,required=True)
    parser.add_argument('--degree',type=int,default=6)
    parser.add_argument('--method',choices=('direct_inventory','interpolate'),default='interpolate')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    result=global_third_moment(json.loads(args.inventory.read_text()),args.n,args.degree,args.method)
    text=json.dumps(result,indent=2)+'\n'
    if args.output:args.output.write_text(text)
    else:print(text,end='')


if __name__=='__main__':main()
