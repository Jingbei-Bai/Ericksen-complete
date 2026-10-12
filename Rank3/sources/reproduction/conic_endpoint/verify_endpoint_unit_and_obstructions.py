"""Exact symbolic endpoint Jacobian and obstruction proof inputs.

All material parameters k,r stay indeterminate. This file establishes the
unit and Taylor-coefficient inputs of the all-orders relative-ideal proof;
it does NOT turn a finite jet test into a no-lift theorem.
No historical normal matrix, tangent matrix, or quadratic certificate is read.
"""
from sage.all import *
from pathlib import Path
import argparse, hashlib, json, time

ap=argparse.ArgumentParser()
ap.add_argument('--source-directory',type=Path,required=True)
ap.add_argument('--output-directory',type=Path,default=Path(__file__).absolute().parent)
args=ap.parse_args(); srcdir=args.source_directory; out=args.output_directory
out.mkdir(parents=True,exist_ok=True); start=time.monotonic()
names=['independent_subcritical_single_RR_full_regular_chart_H1.sobj',
       'independent_subcritical_single_RR_full_complete_E_chart_H1.sobj']
d,g=(load(str(srcdir/name)) for name in names)
P=d['ring']; PK=P.base_ring()
K=FractionField(PolynomialRing(QQ,names=('k','r'))); k,r=K.gens(); m=1+k; q=1+r
UE=PolynomialRing(K,'xb'); E=UE.quotient(UE.gen()**2-2*m,'xb'); xb=E.gen()
ph=PK.ring().hom([k,r],K)
def scalar(c):
    c=PK(c); return E(ph(c.numerator())/ph(c.denominator()))
base=(xb,-xb,E.zero(),E.zero(),E.zero(),2*xb,E.zero(),E.zero(),E.zero())
def ev(f):
    return sum((scalar(c)*prod(z**int(n) for z,n in zip(base,e))
                for e,c in f.dict().items() if not any(n and not z for n,z in zip(e,base))),E.zero())
def norm(c):
    c=E(c); return K((c*E([c[0],-c[1]]))[0])
checks=[]
def ck(name,value):
    assert bool(value),name; checks.append(name)
sources=list(d['sources'])+[it['source'] for it in g['items']]
ck('37 complete sources',len(sources)==37)
ck('sixth pivot is full own-p21 differentiated Q6',g['items'][5]['gauge']==2 and g['items'][5]['Q_index']==6)
ck('all 37 values vanish at both symbolic roots',all(ev(f)==0 for f in sources))
M=matrix(E,37,9,lambda i,j:ev(sources[i].derivative(P.gen(j))))
# Range coordinates: x at fixed (s=x+z,v=e2-(2x+10s)), then e0,e1,delta,t,kappa.
R=matrix(E,9,6); R[0,0]=1; R[1,0]=-1; R[5,0]=2
for j,i in enumerate((3,4,6,7,8),1): R[i,j]=1
# Before implicit compensation: free coordinates v,s,w.
N=matrix(E,9,3); N[5,0]=1; N[1,1]=1; N[5,1]=10; N[2,2]=1; N[3,2]=-2*xb
rows=(2,3,6,25,27,33); A=(M*R).matrix_from_rows(rows)
D=A.det(); expectedD=QQ(36)/7*k*k*q*m*m*xb
ck('six-source unit determinant',D==expectedD)
ck('two-root norm determinant',norm(D)==-QQ(2592)/49*k**4*q*q*m**5)
minor3=M.matrix_from_rows_and_columns((2,3,6),(3,4,6)).det()
ck('three-source unit determinant',minor3==m*q/3584)
T=N-R*A.inverse()*M.matrix_from_rows(rows)*N
ck('compensated tangent kernel',M*T==0)
targetT=matrix(E,[[-QQ(1)/16,-QQ(1)/2,0],[QQ(1)/16,QQ(3)/2,0],[0,0,1],
                  [0,0,-2*xb],[0,0,0],[QQ(7)/8,9,0],[0,0,0],[0,0,0],[0,0,0]])
ck('explicit tangent columns',T==targetT)
B=vector(E,37); B[6]=1; B[13]=QQ(27)/49/(m*q); B[21]=-QQ(12)/49/(r*m); B[23]=QQ(39)/98/m*xb
F=vector(E,37); F[25]=1; F[26]=-1/k
ck('B annihilates every ambient linear derivative',B*M==0)
ck('F annihilates every ambient linear derivative',F*M==0)
pairs=((0,0),(0,1),(1,1),(0,2),(1,2),(2,2))
obstructions={}
for name,V,expected in [('B',B,[0,0,0,-QQ(3)/224/m*xb,0,0]),
                       ('F',F,[QQ(1)/8,QQ(7)/2,32,0,0,0])]:
    H=matrix(E,9,9)
    for i in V.nonzero_positions():
        H+=V[i]*matrix(E,9,9,lambda j,l:ev(sources[i].derivative(P.gen(j)).derivative(P.gen(l))))
    reduced=T.transpose()*H*T
    obs=vector(E,[reduced[a,b]/(2 if a==b else 1) for a,b in pairs])
    ck(name+' true implicit-graph second coefficients',obs==vector(E,expected))
    obstructions[name]=obs
U=-QQ(3)/224/m*xb
ck('B quotient coefficient norm',norm(U)==-QQ(9)/25088/m)
# Every coefficient denominator of every input polynomial uses protected
# material factors only. Each symbolic factor is tested, not sampled.
denominators=set()
for f in sources:
    for c in f.coefficients():
        cc=scalar(c)[0]
        for factor,power in K(cc).denominator().factor():
            if factor.is_constant(): continue
            monic=factor/factor.leading_coefficient()
            denominators.add(str(monic))
            ck('protected source denominator '+str(monic),any(monic==z for z in [K.ring()(k),K.ring()(r),K.ring()(m),K.ring()(q)]))
payload={'status':'PASS_EXACT_SYMBOLIC_INPUTS','scope':__doc__,
         'source_hashes':{name:hashlib.sha256((srcdir/name).read_bytes()).hexdigest() for name in names},
         'k_r_indeterminate':True,'endpoint_root_equation':'xb^2=2*(1+k)',
         'six_source_rows_zero_based':list(rows),'six_source_determinant':str(D),
         'six_source_determinant_norm':str(norm(D)),'three_source_determinant':str(minor3),
         'tangent_columns':[[str(c) for c in row] for row in T.rows()],
         'B_weights':{str(i):str(B[i]) for i in B.nonzero_positions()},
         'F_weights':{str(i):str(F[i]) for i in F.nonzero_positions()},
         'second_monomials':['v^2','v*s','s^2','v*w','s*w','w^2'],
         'obstructions':{key:[str(c) for c in obs] for key,obs in obstructions.items()},
         'quotient_U_constant':str(U),'quotient_U_norm':str(norm(U)),
         'source_denominator_factors':sorted(denominators),
         'checks_count':len(checks),'seconds':time.monotonic()-start}
(out/'endpoint_unit_and_obstructions.json').write_text(json.dumps(payload,indent=2)+'\n')
(out/'endpoint_normal_matrix.txt').write_text('\n'.join(', '.join(str(c) for c in row) for row in M.rows())+'\n')
save({'ring':E,'matrix':M,'range':R,'free':N,'tangents':T,'B':B,'F':F,'obstructions':obstructions,
      'proof_inputs':payload},str(out/'endpoint_unit_and_obstructions.sobj'),compress=False)
print(json.dumps(payload,indent=2),flush=True)
