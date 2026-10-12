"""Exact full-parameter determinant of the protected Cartesian third-jet symbol.

This verifies only the linear symbol, not the missing conic exclusion.
The determinant rows use d2 tr(C), d2 tr(C^2)/2, and d2 det(F).
"""
from sage.all import *
from pathlib import Path
import signal, time, json, resource

OUT = Path(__file__).resolve().parent
start = time.time()
signal.alarm(110)
resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
P = PolynomialRing(QQ, names=('s', 't', 'p', 'q', 'c'))
s,t,p,q,c = P.gens()
K = P.fraction_field()
F = matrix(P, [[s,0,0], [0,t,0], [p,q,c]])
C = F.transpose()*F
J = F.det()
adjF = F.adjugate()
R = PolynomialRing(P, names=tuple('z'+str(j) for j in range(12)))
z = R.gens()
third = {}
for ix, triple in enumerate([(0,0,0),(0,0,1),(0,1,1),(1,1,1)]):
    third[triple] = vector(R,z[3*ix:3*ix+3])

def T(a,b,k):
    return third.get(tuple(sorted((a,b,k))), vector(R,[0,0,0]))

ddF = {}
ddC = {}
for a in range(3):
    for b in range(a,3):
        ddf = matrix(R,3,3,lambda i,j:T(j,a,b)[i])
        ddF[a,b] = ddf
        ddC[a,b] = F.transpose()*ddf + ddf.transpose()*F

def dc(a,b):
    return ddC[tuple(sorted((a,b)))]

equations = []
for a,b in [(0,0),(0,1),(1,1)]:
    equations += [dc(a,b).trace(), (C*dc(a,b)).trace(),
                  (adjF*ddF[a,b]).trace()]
for i,j in [(0,1),(0,2),(1,2)]:
    equations.append(sum(dc(i,a)[a,j]-dc(j,a)[a,i] for a in range(3)))
M = matrix(P, [[f.monomial_coefficient(v) for v in z] for f in equations])
save({'ring':P, 'matrix':M, 'equations':equations},
     str(OUT/'rr_cartesian_twelve_symbol_source.sobj'))
print('symbol saved',time.time()-start,flush=True)
det = M.det()
factors = det.factor()
gap = s*s-t*t
left = t*t*(gap+p*p)-s*s*q*q-gap*c*c
right = s*s*(gap-q*q)+t*t*p*p-gap*c*c
positive = (p*p*left*left + q*q*right*right
            +p*p*q*q*((s*s+t*t)**2*(p*p+q*q)
                       +2*gap**2*(s*s+t*t)))
expected = 128*p*q*s*t*gap*c**5*positive
assert det == expected
real_part = p*(left-(s*s+t*t)*q*q)
imag_part = q*(right+(s*s+t*t)*p*p)
assert positive == real_part**2+imag_part**2
assert ((right+(s*s+t*t)*p*p)-(left-(s*s+t*t)*q*q)
        == gap**2+(s*s+t*t)*(p*p+q*q))
report = {'status':'FULL_PARAMETER_SYMBOL_ONLY_NOT_RR_EXCLUSION',
          'shape':list(M.dimensions()),'determinant':str(factors),
          'exact_positive_factor_identity':True,
          'exact_complex_norm_factor_identity':True,
          'positive_factor':'p^2 A^2+q^2 B^2+p^2 q^2[(s^2+t^2)^2(p^2+q^2)+2(s^2-t^2)^2(s^2+t^2)]',
          'A':str(left),'B':str(right),
          'seconds':time.time()-start,
          'domain':'s,t,c>0 and p*q*(s^2-t^2)!=0',
          'note':'The missing projective no-root theorem is not proved here.'}
(OUT/'rr_cartesian_twelve_symbol_reader.json').write_text(
    json.dumps(report,indent=2),encoding='utf-8')
save({'matrix':M,'determinant':det,'factorization':factors,
      'positive_factor':positive,'A':left,'B':right,
      'real_part':real_part,'imag_part':imag_part},
     str(OUT/'rr_cartesian_twelve_symbol_certificate.sobj'))
print(json.dumps(report,indent=2),flush=True)
