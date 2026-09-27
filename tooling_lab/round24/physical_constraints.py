#!/usr/bin/env python3
"""Compile exact inverse-coordinate inequalities and eliminate disjoint carries."""
from fractions import Fraction as F
from itertools import combinations
from math import gcd,lcm,floor


def pair(v):
    v=F(v);return [v.numerator,v.denominator]


def integerize(direction,radius):
    scale=lcm(*(v.denominator for v in direction)) if direction else 1
    raw=[int(v*scale) for v in direction];divisor=gcd(*raw) if raw else 0
    divisor=divisor or 1
    return [v//divisor for v in raw],floor(radius*scale/divisor)


def compile_constraints(c):
    N,p,k,U,adj=c['N'],c['p'],c['k'],c['erased'],c['adj']
    V,I=c['visible_directions'],c['invisible_directions'];P=[]
    for row in c['basis']:
        P.append([F(sum((adj[i-u] if i>=u else -adj[N+i-u])*x for u,x in zip(U,row)),k) for i in range(N)])
    supports={j:[i for i,x in enumerate(P[j]) if x] for j in I}
    owners=[[j for j in I if P[j][i]] for i in range(N)]
    disjoint=all(len(js)<=1 for js in owners);cuts=[]
    def append(kind,direction,radius,**metadata):
        ints,bound=integerize(direction,radius)
        cuts.append({'kind':kind,**metadata,'direction':[pair(x) for x in direction],'radius':pair(radius),
                     'integer_direction':ints,'integer_radius':bound})
    for i,js in enumerate(owners):
        if not js:append('fixed_coordinate',[P[j][i] for j in V],F(p-1),coordinate=i)
    if disjoint:
        for j in I:
            for i,h in combinations(supports[j],2):
                direction=[P[v][i]/P[j][i]-P[v][h]/P[j][h] for v in V]
                radius=F(p-1)*(1/abs(P[j][i])+1/abs(P[j][h]))
                append('single_carry_elimination',direction,radius,invisible_direction=j,coordinates=[i,h])
    return {'image_rows':[[pair(x) for x in row] for row in P],
            'visible_directions':V,'invisible_directions':I,
            'carry_supports':[{'basis_index':j,'coordinates':supports[j]} for j in I],
            'projection_mode':'disjoint_invisible_supports' if disjoint else 'direct_rows_only',
            'continuous_projection_complete':disjoint,
            'integer_rounding':'clear_direction_denominators_then_divide_gcd_and_floor_radius',
            'cuts':cuts}
