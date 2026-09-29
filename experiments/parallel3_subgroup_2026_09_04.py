"""Exact order-index and reduction certificates; no prime scan or factorization."""

from collections import Counter, defaultdict
from hashlib import sha256
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def multiply(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def reduce_poly(a, modulus, p=None):
    a = a[:]
    d = len(modulus) - 1
    assert modulus[-1] == 1
    if p is not None:
        a = [x % p for x in a]
    for j in range(len(a) - 1, d - 1, -1):
        c = a[j]
        for i, x in enumerate(modulus):
            a[j - d + i] -= c * x
            if p is not None:
                a[j - d + i] %= p
    return (a + [0] * d)[:d]


def power_mod(a, exponent, modulus, p=None):
    out = [1]
    while exponent:
        if exponent & 1:
            out = reduce_poly(multiply(out, a), modulus, p)
        a = reduce_poly(multiply(a, a), modulus, p)
        exponent //= 2
    return reduce_poly(out, modulus, p)


def chebyshev(d):
    out, degree = [0, 1], 1
    while degree < d:
        out = multiply(out, out)
        out[0] -= 2
        degree *= 2
    assert degree == d
    return out


def multiplication_matrix(a, modulus):
    d = len(modulus) - 1
    columns = [reduce_poly([0]*j + a, modulus) for j in range(d)]
    return [[columns[j][i] for j in range(d)] for i in range(d)]


def power_basis_matrix(a, modulus, p=None):
    d = len(modulus)-1
    columns, power = [], [1]
    for _ in range(d):
        columns.append(reduce_poly(power, modulus, p))
        power = reduce_poly(multiply(power, a), modulus, p)
    return [[columns[j][i] for j in range(d)] for i in range(d)]


def determinant(a):
    a = [row[:] for row in a]
    d, previous, sign = len(a), 1, 1
    for j in range(d-1):
        if not a[j][j]:
            pivot = next((i for i in range(j+1,d) if a[i][j]), None)
            if pivot is None:
                return 0
            a[j], a[pivot] = a[pivot], a[j]
            sign *= -1
        pivot = a[j][j]
        for i in range(j+1,d):
            for h in range(j+1,d):
                value = pivot*a[i][h]-a[i][j]*a[j][h]
                assert value % previous == 0
                a[i][h] = value//previous
        for i in range(j+1,d):
            a[i][j] = 0
        previous = pivot
    return sign*a[-1][-1]


def derivative(a):
    return [i*a[i] for i in range(1,len(a))]


def discriminant(a):
    d = len(a)-1
    return (-1)**(d*(d-1)//2)*determinant(multiplication_matrix(derivative(a),a))


def rank_mod(a, p):
    a = [[x % p for x in row] for row in a]
    rank = 0
    for j in range(len(a[0])):
        pivot = next((i for i in range(rank,len(a)) if a[i][j]),None)
        if pivot is None:
            continue
        a[rank],a[pivot]=a[pivot],a[rank]
        inverse=pow(a[rank][j],-1,p)
        a[rank]=[x*inverse%p for x in a[rank]]
        for i in range(rank+1,len(a)):
            c=a[i][j]
            if c:
                a[i]=[(x-c*y)%p for x,y in zip(a[i],a[rank])]
        rank+=1
    return rank


def evaluate(a,x,p):
    out=0
    for c in reversed(a):
        out=(out*x+c)%p
    return out


def prime(p):
    return p>=2 and all(p%d for d in range(2,math.isqrt(p)+1))


def index_certificate(n,R):
    d,k=n//4,n//2
    C=chebyshev(d)
    alpha=[-x for x in power_mod([2,1],k,C)]
    M=power_basis_matrix(alpha,C)
    index=abs(determinant(M))
    assert index>0
    dc=discriminant(C)
    dr=discriminant(R)
    assert dc==2**(d-1)*d**d
    assert dr==index*index*dc
    A_derivative=[-k*x for x in power_mod([2,1],k-1,C)]
    derivative_norm=abs(determinant(multiplication_matrix(A_derivative,C)))
    assert derivative_norm==k**d*2**(k-1)
    return {"parent_order":n,"degree":d,"ambient_polynomial":C,"alpha_remainder":alpha,
            "ambient_discriminant_hex":hex(dc),"primitive_discriminant_hex":hex(dr),
            "order_index_hex":hex(index),"order_index_bits":index.bit_length(),
            "unreduced_derivative_norm_hex":hex(derivative_norm)}


def modular_certificate(p,n,g):
    assert prime(p) and (p-1)%n==0 and pow(g,n,p)==1 and pow(g,n//2,p)==p-1
    d,k=n//4,n//2
    C=chebyshev(d)
    alpha=[-x%p for x in power_mod([2,1],k,C,p)]
    M=power_basis_matrix(alpha,C,p)
    rank=rank_mod(M,p)
    fibers=defaultdict(list)
    t_values=[]
    for j in range(1,n//2,2):
        h=pow(g,j,p)
        t=(h+pow(h,-1,p))%p
        assert evaluate(C,t,p)==0
        assert -k*pow(2+t,k-1,p)%p!=0
        label=pow(1+h,n,p)
        assert label==evaluate(alpha,t,p)==-pow(2+t,k,p)%p
        fibers[label].append(j)
        t_values.append(t)
    assert len(set(t_values))==d
    assert rank==len(fibers)
    first_pairs=sum(len(js)*(len(js)-1)//2 for js in fibers.values())
    lifted_g=pow(g,p,p*p)
    assert lifted_g%p==g and pow(lifted_g,n,p*p)==1
    lifted=Counter(pow(1+pow(lifted_g,j,p*p),n,p*p) for j in range(1,n//2,2))
    assert max(lifted.values())==1, "This selected witness needs further precision"
    return {"p":p,"parent_order":n,"generator":g,"rank":rank,"nullity":d-rank,
            "ambient_discriminant_mod_p":2**(d-1)*pow(d,d,p)%p,
            "all_unreduced_derivatives_nonzero":True,
            "fiber_histogram":dict(sorted(Counter(map(len,fibers.values())).items())),
            "colliding_exponents":{str(label):js for label,js in fibers.items() if len(js)>1},
            "lifted_mod_p_squared_squarefree":True,"exact_index_valuation":first_pairs,
            "balanced_energy_from_fibers":2*k*sum(len(js)**2 for js in fibers.values())}


def main():
    prior=ROOT/'results/parallel2_subgroup_2026_09_04.json'
    r=json.loads(prior.read_text())
    polys={int(n):list(map(int,a)) for n,a in r['integer_primitive_polynomials'].items()}
    integers=[]
    for n in (4,8,16,32,64):
        item=index_certificate(n,polys[n])
        integers.append(item)
        print(json.dumps({"parent_order":n,"index_bits":item['order_index_bits'],"discriminant_identity":"passed"}),flush=True)
    modular=[modular_certificate(17,16,3),modular_certificate(67403009,128,64701253),
             modular_certificate(17189277697,512,730170861)]
    p,n,g=67403009,128,64701253
    before=[pow(1+pow(g,j,p),n,p) for j in (13,25)]
    after=[pow(1+pow(g,3*j,p),n,p) for j in (13,25)]
    assert before[0]==before[1] and after[0]!=after[1]
    # The source's exponential hypothesis is outside every dyadic quartic window k>=16.
    for k in (16,32,64,128,256,512):
        assert 14**(k//2)>16*k**4
    script=ROOT/'experiments/parallel3_subgroup_2026_09_04.py'
    note=ROOT/'research/parallel3-subgroup-2026-09-04.md'
    result={"status":"Uniform algebraic obstruction and sufficient index estimate; uniform energy bound unproved.",
            "integer_index_certificates":integers,"modular_certificates":modular,
            "galois_partition_counterexample":{"p":p,"n":n,"g":g,"exponents":[13,25],"multiplier":3,"before":before,"after":after},
            "source_sha256":{str(path.relative_to(ROOT)):sha256(path.read_bytes()).hexdigest() for path in (script,prior,note)}}
    (ROOT/'results/parallel3_subgroup_2026_09_04.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({"modular_certificates":len(modular),"status":"passed"}),flush=True)


if __name__=='__main__':
    main()
