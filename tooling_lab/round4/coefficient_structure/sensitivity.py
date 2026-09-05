#!/usr/bin/env python3
"""Symbolically verified scaling limits and a finite contribution audit."""
from fractions import Fraction
import hashlib,json,math,sys
from pathlib import Path
import sympy as sp
from closed_coefficient import quartic,contraction_coefficient_degree6
HERE=Path(__file__).resolve().parent

def encode(x):return [x.numerator,x.denominator]

def main():
    q,n,r,alpha,lam,s,beta=sp.symbols('q n r alpha lambda s beta',positive=True)
    c=-sp.ff(r,4)*sp.ff(q-r-3,4)*quartic(q,r)/(2*sp.ff(q-3,12))
    fixed=sp.factor(sp.limit((q**3*c).subs(r,n-3),q,sp.oo))
    bulk=sp.factor(sp.limit((c/q**2).subs(r,alpha*q-3),q,sp.oo))
    transition=sp.factor(sp.limit((q*c).subs({q:s*s,r:lam*s},simultaneous=True),s,sp.oo))
    middle=sp.factor(sp.limit((c/q).subs({q:s*s,r:s*s/2+beta*s-3},simultaneous=True),s,sp.oo))
    critical=sp.factor(sp.limit((q*q*c).subs({q:s**4,r:alpha*s-3},simultaneous=True),s,sp.oo))
    assert sp.expand(fixed+sp.ff(n-3,4))==0
    assert sp.expand(bulk-alpha**6*(1-alpha)**4*(1-2*alpha)**2/2)==0
    assert sp.expand(transition-lam**4*(lam**2-2)/2)==0
    assert sp.expand(middle-(beta**2-sp.Rational(1,4))/512)==0
    assert sp.expand(critical+alpha**4)==0
    low={}
    for card in [6,7,8]:low[card]=str(sp.factor(c.subs(r,card-3)))
    source=HERE.parent/'marked_counts_ablation/contraction_p1000033_marks_0_1_2.json';data=json.loads(source.read_text());p=data['p'];card=31
    coefficient=contraction_coefficient_degree6(p,card);term=coefficient*data['Q'];second=Fraction(*next(row for row in data['reduced_moments'] if row['n']==31)['second_moment'])
    norms=[0,0,0,0]
    for cell in data['cells']:
        a,b,d=cell['pattern'];size=cell['size']
        if 0 in [a,b,d]:continue
        h2=a*b+a*d+b*d;h3=a*b*d
        norms[0]+=size*h2;norms[1]+=size*h2*h2;norms[2]+=size*h3;norms[3]+=size*h3*h3
    sum2,norm2,sum3,norm3=norms
    bound_Q_squared=(Fraction(norm2)-Fraction(sum2*sum2,p))*(p*norm3-sum3*sum3)
    assert data['Q']*data['Q']<=bound_Q_squared
    out={'status':'all five scaling limits verified symbolically from the exact formula',
         'low_cardinality_formulas':low,
         'fixed_n_q3_times_c_limit':str(fixed),
         'n_equal_alpha_q_c_over_q2_limit':str(bulk),
         'r_equal_lambda_sqrtq_q_times_c_limit':str(transition),
         'n_equal_q_half_plus_beta_sqrtq_c_over_q_limit':str(middle),
         'n_equal_alpha_q_quarter_q2_times_c_limit':str(critical),
         'uniform_critical_scalar_correction_bound':'For n~alpha*q^(1/4), |cQ| <= (sqrt(3)*alpha^4+o(1))*q^(-1/2), using ||h2||²<=3q, ||h3||²<q and ||S||=sqrt(q). Only the conditional-average scalar correction is bounded.',
         'scope':'limits follow from the exact rational function; integer set sizes may use rounding without changing leading limits',
         'million_prime_n31_example':{'p':p,'n':card,'Q':data['Q'],'coefficient':encode(coefficient),
             'coefficient_float':float(coefficient),'scalar_contribution_cQ':encode(term),'scalar_contribution_float':float(term),
             'full_second_moment_float':float(second),'scalar_contribution_fraction_of_second':encode(term/second),
             'scalar_contribution_fraction_float':float(term/second),
             'conference_Cauchy_Q_bound_squared':encode(bound_Q_squared),
             'conference_Cauchy_cQ_absolute_bound_display':math.sqrt(float(bound_Q_squared))*abs(float(coefficient)),
             'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest()},
         'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'sympy_version':sp.__version__,'python':sys.executable}
    (HERE/'sensitivity.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2),flush=True)
if __name__=='__main__':main()
