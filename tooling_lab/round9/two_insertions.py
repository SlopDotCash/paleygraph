#!/usr/bin/env python3
"""Exact two-insertion moments. All final contractions use unbounded integers."""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json
import subprocess
import time
import numpy as np

HERE = Path(__file__).resolve().parent


def rat(x):
    x = F(x)
    return [x.numerator, x.denominator]


def validate(q, selected, deleted, degree):
    c = sorted(selected)
    a = sorted(deleted)
    if not (5 <= q <= 10_000_000 and 2 <= len(c) <= min(64, q-2)
            and len(set(c)) == len(c) and all(type(x) is int and 0 <= x < q for x in c)
            and len(a) == 2 and len(set(a)) == 2 and set(a) <= set(c)
            and type(degree) is int and 0 <= degree <= min(6, len(c))):
        raise ValueError('invalid set, deletion pair, degree or field size')
    return c, a


def moments(s):
    q, c, degree = s['q'], s['selected'], s['degree']
    n = len(c); m = q-n; N = m*(m-1)
    g, h, f, v, w = (s[k] for k in ('g_selected','h_selected','f_selected','v_selected','w_selected'))
    K = s['K_selected']; G0, G2, H0, H2, Q, t0 = (s[k] for k in ('G0','G2','H0','H2','Q','t0'))
    outside_f = -sum(f)
    outside_f2 = q*G2-G0*G0-sum(x*x for x in f)
    outside_h = H0-sum(h)
    outside_h2 = H2-sum(x*x for x in h)
    outside_fh = Q-sum(x*y for x,y in zip(f,h))
    ksum = sum(map(sum,K))
    knorm = (q*q-2*q)*H2+H0*H0
    rownorm = [q*(H2-x*x)-y*y for x,y in zip(h,v)]
    outside_k2 = knorm-2*sum(rownorm)+sum(x*x for row in K for x in row)
    L = -q*sum(w)+G0*sum(v)+sum(f[i]*sum(K[i]) for i in range(n))
    z1 = 2*(m-1)*outside_f+ksum-m*H0+outside_h
    z2 = (2*(m-2)*outside_f2+2*outside_f**2+outside_k2+4*L
          -m*H0**2-outside_h2-4*H0*outside_f+4*outside_fh+2*H0*outside_h)
    total = N*t0+z1
    squares = N*t0*t0+2*t0*z1+z2
    mean = F(total,N); second = F(squares,N); variance = second-mean*mean
    qb = F((q*G2-G0*G0)*(q*H2-H0*H0),q)
    assert outside_f2 >= 0 and outside_k2 >= 0 and all(x >= 0 for x in rownorm)
    assert Q*Q <= qb and variance >= 0
    return {'ordered_pair_count':N,'distinct_final_sets':N//2,'target_sum':total,
            'target_square_sum':squares,'mean':rat(mean),'second_moment':rat(second),
            'variance':rat(variance),'Q_variance_coefficient':rat(F(4,N)),
            'variance_without_Q':rat(variance-F(4*Q,N)),
            'Q_variance_correction':rat(F(4*Q,N)),
            'Q_cauchy_bound_squared':rat(qb),
            'Q_free_variance_radius_squared':rat(F(16,N*N)*qb)}


def coefficient(signs, degree):
    if degree < 0: return 0
    p = sum(x == 1 for x in signs); m = sum(x == -1 for x in signs)
    return sum((-1)**j*comb(m,j)*comb(p,degree-j)
               for j in range(max(0,degree-p),min(degree,m)+1))


def from_matrix(matrix, selected, deleted, degree=6):
    original = np.asarray(matrix)
    S = np.asarray(matrix, dtype=np.int64)
    if not np.array_equal(original, S):
        raise ValueError('matrix entries must be exact integers')
    q = len(S); c, a = validate(q, selected, deleted, degree)
    if q > 257 or S.shape != (q,q) or not np.array_equal(S,S.T):
        raise ValueError('matrix adapter requires a symmetric matrix of size at most257')
    if (not np.all(np.diag(S)==0) or not np.all(np.isin(S+np.eye(q,dtype=np.int64),[-1,1]))
            or not np.all(S.sum(axis=0)==0) or not np.array_equal(S@S,q*np.eye(q,dtype=np.int64)-1)):
        raise ValueError('matrix fails balanced conference identities')
    remaining = sorted(set(c)-set(a))
    g = np.array([coefficient(row[remaining],degree-1) for row in S],dtype=np.int64)
    h = np.array([coefficient(row[remaining],degree-2) for row in S],dtype=np.int64)
    f = S@g; v = S@h; w = S@(g*h)
    K = (S[c,:]*h)@S[:,c]
    s = {'q':q,'selected':c,'deleted':a,'degree':degree,
         't0':sum(coefficient(row[remaining],degree) for row in S),
         'G0':sum(map(int,g)),'G2':sum(int(x)**2 for x in g),
         'H0':sum(map(int,h)),'H2':sum(int(x)**2 for x in h),
         'Q':sum(int(x)*int(y) for x,y in zip(g,v)),
         'K_selected':K.tolist()}
    for key, field in [('g',g),('h',h),('f',f),('v',v),('w',w)]:s[key+'_selected']=field[c].tolist()
    return {'statistics':s,**moments(s)}


def query(q, selected, deleted, degree=6):
    c, a = validate(q,selected,deleted,degree)
    payload = ' '.join(map(str,[q,len(c),degree,*a,*c]))+'\n'
    start = time.monotonic()
    proc = subprocess.run([str(HERE/'two_insertions_backend')],input=payload,text=True,capture_output=True)
    if proc.returncode:raise ValueError(proc.stderr.strip())
    s = json.loads(proc.stdout)
    return {'statistics':s,**moments(s),'seconds':time.monotonic()-start}


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('q',type=int);parser.add_argument('--selected',type=int,nargs='+',required=True)
    parser.add_argument('--deleted',type=int,nargs=2,required=True);parser.add_argument('--degree',type=int,default=6)
    args=parser.parse_args()
    print(json.dumps(query(args.q,args.selected,args.deleted,args.degree),indent=2))
