"""Small exact scalar-monodromy graph in one auxiliary amplitude lambda.

No lambda or entry pivot is inverted. Equivalence to scalar monodromy is
conditional only on the two explicitly recorded bilinear block determinants.
"""
from pathlib import Path
import json,time
from sage.all import load,save,QQ,PolynomialRing,matrix,vector

out=Path(__file__).resolve().parent;start=time.monotonic()
source=load(str(out/'general_curve_dimensionless_leading.sobj'))
C=PolynomialRing(QQ,names=('H','Z','S','N'));H,Z,S,N=C.gens()
def cv(poly):
    ans=C.zero()
    for mm,c in poly.dict().items():
        assert mm[0]==0 and mm[1]==0
        ans+=c*H**mm[2]*Z**mm[3]*S**mm[4]*N**mm[5]
    return ans
B={ij:matrix(C,2,2,[cv(c) for c in M.list()]) for ij,M in source['blocks'].items()}
P=B[0,1];Q=B[0,2];R=B[1,2];J=matrix(C,[[0,-1],[1,0]])
checks={}
for ij,M in B.items():
    checks['only_22_N_'+str(ij)]=all(M[i,j].degree(N)<=int(i==1 and j==1) for i in range(2) for j in range(2))
D=P.transpose()*J*Q
checks['D_at_most_linear_N']=all(c.degree(N)<=1 for c in D.list())
checks['D00_independent_N']=D[0,0].degree(N)<=0

# Four rows multiplying (1,lambda,N,lambda*N).
L=matrix(C,4,4)
for row,(i,j) in enumerate([(0,0),(0,1),(1,0),(1,1)]):
    r0=R[i,j].subs({N:0});r1=R[i,j].derivative(N)
    d0=D[i,j].subs({N:0});d1=D[i,j].derivative(N)
    L[row]=vector(C,[r0,-d0,r1,-d1])
checks['lift_coefficients_independent_N']=all(c.degree(N)<=0 for c in L.list())
checks['first_three_no_bare_N']=all(L[i,2]==0 for i in range(3))
U=PolynomialRing(C,'lambda');lam=U.gen()
eq=[U(R[i,j])-lam*U(D[i,j]) for i,j in [(0,0),(0,1),(1,0),(1,1)]]
checks['four_lift_rows_replay']=all(sum((U(L[i,j])*v for j,v in enumerate([U.one(),lam,U(N),lam*U(N)])),U.zero())==eq[i] for i in range(4))

# Source monodromy and an exact division-free implication of the four rows.
# If R=lambda*P^T J Q then M=-lambda det(P)det(Q) I.
M=J*Q*J*R.transpose()*J*P.transpose()
residual=matrix(U,2,2,[U(R[i,j])-lam*U(D[i,j]) for i in range(2) for j in range(2)])
Mu=matrix(U,M)
checks['monodromy_source_identity']=Mu+lam*U(P.det()*Q.det())*matrix.identity(U,2)==matrix(U,J*Q*J)*residual.transpose()*matrix(U,J*P.transpose())

# Converse encoded without actually inverting P,Q: the two outer matrices
# have determinants det(Q),det(P), hence are invertible precisely on the
# stated chart. Algebraically check the proposed lambda from trace M.
K=C.fraction_field();PK=matrix(K,P);QK=matrix(K,Q);RK=matrix(K,R)
CK=PK.inverse().transpose()*RK*QK.inverse()
checks['skew_to_scalar_linear_map']=matrix(K,M)==-K(P.det()*Q.det())*CK.transpose()*matrix(K,J)

assert all(checks.values()),[key for key,val in checks.items() if not val]
factor_data={'D00':D[0,0].factor(),'B12_00':R[0,0].factor(),
             'det_B01':P.det().factor(),'det_B02':Q.det().factor(),
             'B12_N_coefficient':R[1,1].derivative(N).factor()}
payload={'ring':C,'variables':C.gens(),'blocks':B,'J':J,'D':D,'lift_matrix':L,
 'lift_monomials':('1','lambda','N','lambda*N'),'lambda_ring':U,'equations':eq,
 'source_monodromy':M,'factor_data':factor_data,'checks':checks,
 'source_file':'general_curve_dimensionless_leading.sobj',
 'scope':'Exact source equations with auxiliary lambda; no lambda/D00 division. Equivalent to scalar monodromy only on det(B01)*det(B02)!=0. The scalar condition at singular block charts is retained outside this equivalence and requires separate coverage.'}
save(payload,str(out/'independent_scalar_lambda_lift.sobj'))
report={'status':'EXACT_LAMBDA_LIFT_PASS','checks':checks,
 'D_terms':[len(c.dict()) for c in D.list()],
 'lift_terms':[[len(c.dict()) for c in row] for row in L.rows()],
 'factors':{key:str(value) for key,value in factor_data.items()},
 'elapsed_seconds':time.monotonic()-start,'scope':payload['scope']}
(out/'independent_scalar_lambda_lift.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2),flush=True)
