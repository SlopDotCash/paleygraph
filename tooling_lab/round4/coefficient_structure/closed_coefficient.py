#!/usr/bin/env python3
"""Exact closed coefficient, ONLY for degree6 and exactly three marks.
Derived in symbolic_formula.py by a polynomial recurrence, not fitted samples.
"""
from fractions import Fraction
from math import prod

def quartic(q,r):
    a4=-4*(q-21)*(q-19)
    a3=4*(q**3-52*q*q+715*q-1792)
    a2=-q**4+87*q**3-1689*q*q+8537*q-13846
    a1=-15*q**4+445*q**3-3683*q*q+13299*q-19838
    a0=2*q**5-74*q**4+942*q**3-6158*q*q+21760*q-32984
    return ((((a4*r+a3)*r+a2)*r+a1)*r+a0)

def contraction_coefficient_degree6(q,n):
    assert isinstance(q,int) and q>=17 and q%4==1
    assert isinstance(n,int) and 6<=n<=q
    r=n-3
    return Fraction(-prod(r-i for i in range(4))*prod(q-n-i for i in range(4))*quartic(q,r),
                    2*prod(q-3-i for i in range(12)))
