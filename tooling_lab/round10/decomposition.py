#!/usr/bin/env python3
"""Exact decomposition into inward targets and the induced-graph boundary."""
from math import comb


def elementary_count(size,positive,zero,degree):
    if degree<0:return 0
    negative=size-positive-zero
    return sum((-1)**j*comb(negative,j)*comb(positive,degree-j)
               for j in range(max(0,degree-positive),min(degree,negative)+1))


def components(q,selected,degree,induced,targets):
    n=len(selected);m=q-n;r=m-degree;v=n-degree+1
    c0=-(2*r+1)*(n-degree)+1-r
    H=(r+1)*(r+2);L=c0-4*v;J=v*(v+1)
    base_counts=[(sum(s==1 for s in row),sum(s==0 for s in row)) for row in induced]
    rows=[]
    for record in targets['pairs']:
        i,j=record['indices'];A=set(range(n))-{i,j};b=0
        for x,row in enumerate(induced):
            p,z=base_counts[x];p-=int(row[i]==1)+int(row[j]==1);z-=int(row[i]==0)+int(row[j]==0)
            coeff=lambda d:elementary_count(n-2,p,z,d)
            if x in A:b+=2*r*coeff(degree-2)-2*(row[i]+row[j])*coeff(degree-3)-2*v*coeff(degree-4)
            else:b-=2*coeff(degree-2)
        high=H*record['values'][0]-2*(r+1)*(targets['singles'][i][0]+targets['singles'][j][0])+2*targets['full'][0]
        low=L*record['values'][1]+2*v*(targets['singles'][i][1]+targets['singles'][j][1])+J*record['values'][2]
        rows.append({'indices':[i,j],'deleted':[selected[i],selected[j]],'high_degree_part':high,
                     'lower_degree_part':low,'induced_boundary_part':b,'target_sum':high+low+b,
                     'pair_target_high':record['values'][0],'pair_target_low':record['values'][1],
                     'pair_target_lower':record['values'][2]})
    return {'q':q,'selected':selected,'degree':degree,'high_pair_coefficient':H,'low_pair_coefficient':L,
            'lower_pair_coefficient':J,'ordered_insertion_pairs_per_deletion':m*(m-1),'pairs':rows}


def project_edges(n,edges):
    """Remove all vertex-additive effects; return a common integer denominator."""
    if n<3:raise ValueError('at least3 vertices required')
    degree=[0]*n;total=0
    for (i,j),value in edges.items():degree[i]+=value;degree[j]+=value;total+=value
    den=(n-1)*(n-2)
    residual={(i,j):den*value-(n-1)*(degree[i]+degree[j])+2*total for (i,j),value in edges.items()}
    assert all(sum(v for (i,j),v in residual.items() if a in (i,j))==0 for a in range(n))
    return residual,den
