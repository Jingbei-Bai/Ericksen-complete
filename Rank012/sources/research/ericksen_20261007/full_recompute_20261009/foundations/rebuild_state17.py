"""Independent reconstruction from Section 2; no historical coefficient input.
Run with SageMath's Python. Every spectral calculation is over QQ(a,b).
"""
from sage.all import *
from pathlib import Path
import json, pickle, hashlib, time
OUT=Path(__file__).absolute().parent
START=time.monotonic(); checks=[]
def check(name,passed,**extra):
    checks.append(dict(name=name,passed=bool(passed),**extra))
    print(name, bool(passed), '%.2fs'%(time.monotonic()-START),flush=True)
    (OUT/'state17_checks.json').write_text(json.dumps(dict(status='RUNNING',checks=checks),indent=2))
    assert passed,name
AB=PolynomialRing(QQ,names=('a','b')); a,b=AB.gens(); K=AB.fraction_field()
names=tuple('p%d%d'%(i+1,j+1) for i in range(3) for j in range(3))+tuple('t%d'%i for i in range(6))+('u0','u1')
extras=tuple('d%d'%i for i in range(27))+tuple('T%d'%i for i in range(18))+tuple('U%d'%i for i in range(6))
R=PolynomialRing(K,names=names+extras,order='lex'); gens=R.gens()
p=list(gens[:9]); t=list(gens[9:15]); u=list(gens[15:17]); d=list(gens[17:44]); T=list(gens[44:62]); U=list(gens[62:68])
state=p+t+u
def eps(i,j,k):
    if len({i,j,k})<3:return 0
    return 1 if (i,j,k) in ((0,1,2),(1,2,0),(2,0,1)) else -1
def C(k,i,j):return -sum(eps(i,j,l)*p[3*k+l] for l in range(3))
def derivative(f,variables,values):return sum((f.derivative(v)*z for v,z in zip(variables,values)),R.zero())
def deriv0(i,f):return derivative(f,p,d[9*i:9*i+9])
def gamma(w,k,i,j):return (w[k]*C(k,i,j)-w[i]*C(i,j,k)+w[j]*C(j,k,i))/(2*w[k])
pairs=((1,2),(2,0),(0,1))
def curvatures(w):
    g={(k,i,j):gamma(w,k,i,j) for k in range(3) for i in range(3) for j in range(3)}
    return [w[l]*(deriv0(i,g[l,j,k])-deriv0(j,g[l,i,k])+sum(g[m,j,k]*g[l,i,m]-g[m,i,k]*g[l,j,m]-C(m,i,j)*g[l,m,k] for m in range(3))) for l,k in pairs for i,j in pairs]
def curls(w):
    aa=[(w[i]-w[(i+1)%3])*p[3*((i+1)%3)+(i+2)%3]-(w[i]-w[(i+2)%3])*p[3*((i+2)%3)+(i+1)%3] for i in range(3)]
    return [deriv0((i+1)%3,aa[(i+2)%3])-deriv0((i+2)%3,aa[(i+1)%3])+sum(aa[l]*p[3*l+i] for l in range(3)) for i in range(3)]
def linear_system(expressions,variables):
    M=matrix(K,[[K(f.derivative(v)) for v in variables] for f in expressions])
    rhs=vector(R,[-f+sum(M[i,j]*v for j,v in enumerate(variables)) for i,f in enumerate(expressions)])
    return M,rhs
def selected_solve(M,rhs,rows,cols,free,free_values):
    print('solving selected block',len(rows),flush=True)
    B=M.matrix_from_rows_and_columns(rows,cols)
    rhs0=vector(R,[rhs[i] for i in rows])-M.matrix_from_rows_and_columns(rows,free)*vector(R,free_values)
    # Sparse Gauss--Jordan avoids a dense characteristic-polynomial inverse.
    n=B.nrows(); augmented=[list(B.row(i))+[K(int(i==j)) for j in range(n)] for i in range(n)]
    determinant=K.one()
    for col in range(n):
        pivot=next(i for i in range(col,n) if augmented[i][col])
        if pivot!=col: augmented[col],augmented[pivot]=augmented[pivot],augmented[col]; determinant=-determinant
        value=augmented[col][col]; determinant*=value; determinant.reduce()
        augmented[col]=[v/value for v in augmented[col]]
        for i in range(n):
            if i==col or not augmented[i][col]:continue
            fac=augmented[i][col]
            augmented[i]=[v-fac*w for v,w in zip(augmented[i],augmented[col])]
            for v in augmented[i]:v.reduce()
    inverse=matrix(K,[row[n:] for row in augmented])
    print('inverse complete',n,flush=True)
    solution=inverse*rhs0
    answer=[R.zero()]*M.ncols()
    for j,v in zip(cols,solution):answer[j]=v
    for j,v in zip(free,free_values):answer[j]=v
    return determinant,answer
E=curvatures([K.one()]*3)+curvatures([K.one(),1/K(a),1/K(b)])+curls([K.one(),K(a),K(b)])+curls([K.one(),1/K(a),1/K(b)])
print('24 curvature/pressure sources constructed',flush=True)
M0,rhs0=linear_system(E,d)
print('first coefficient matrix constructed',flush=True)
I0=(0,1,2,3,4,5,6,7,8,9,10,11,13,14,17,18,19,20,21,22,23)
J0=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,15,16,17,18,19,21,26); K0=(14,20,22,23,24,25)
det0,Q=selected_solve(M0,rhs0,I0,J0,K0,t)
check('fnd:firstminor determinant',det0==-(a-1)**4*(a-b)**5*(b-1)**6/(32*K(a)**6*b**8),value=str(det0.factor()))
check('all 24 original equations after first solve',M0*vector(R,Q)==rhs0,count=24)
check('three diagonal contact derivatives',all(Q[9*i+3*i+i]+2*p[3*i+i]*(p[3*((i+1)%3)+(i+2)%3]-p[3*((i+2)%3)+(i+1)%3])==0 for i in range(3)),count=3)
def deriv1(i,f):return derivative(f,p+t,Q[9*i:9*i+9]+T[6*i:6*i+6])
E1=[deriv1(i,Q[9*j+k])-deriv1(j,Q[9*i+k])-sum(C(m,i,j)*Q[9*m+k] for m in range(3)) for i,j in pairs for k in range(9)]
M1,rhs1=linear_system(E1,T)
I1=(0,1,2,3,4,5,6,7,9,10,11,14,15,16,20,23)
J1=(0,1,2,3,4,5,6,7,8,9,10,11,13,15,16,17); K1=(12,14)
det1,Tsol=selected_solve(M1,rhs1,I1,J1,K1,u)
check('fnd:secondminor determinant',det1==-2*(a-b)**2/(K(a)**2*(b-1)**2),value=str(det1.factor()))
raw3=list(M1*vector(R,Tsol)-rhs1); C3=[]
for f in raw3:
    if f:
        ff=f/f.leading_coefficient()
        if ff not in C3:C3.append(R(ff))
check('exactly three distinct monic cubic residuals',len(C3)==3,count=len(C3))
def deriv2(i,f):return derivative(f,state,Q[9*i:9*i+9]+Tsol[6*i:6*i+6]+U[2*i:2*i+2])
E2=[deriv2(i,Tsol[6*j+k])-deriv2(j,Tsol[6*i+k])-sum(C(m,i,j)*Tsol[6*m+k] for m in range(3)) for i,j in pairs for k in range(6)]
M2,rhs2=linear_system(E2,U)
I2=(0,2,3,5,6,8)
det2,Usol=selected_solve(M2,rhs2,I2,tuple(range(6)),(),())
check('six by six prolongation determinant',det2==b*(a-1)/(2*K(a)*(b-1)),value=str(det2.factor()))
F=[Q[9*i:9*i+9]+Tsol[6*i:6*i+6]+Usol[2*i:2*i+2] for i in range(3)]
def D(i,f):return derivative(f,state,F[i])
Q4=[f for f in M2*vector(R,Usol)-rhs2 if f]+[D(i,f) for f in C3 for i in range(3)]
Q4=[f for f in Q4 if f]
check('21 original quartic residuals',len(Q4)==21,count=len(Q4))
weights=[1]*9+[2]*6+[3]*2
check('all 51 field weights',all(all(sum(e*w for e,w in zip(mon,weights))==weights[j]+1 for mon in f.dict()) for row in F for j,f in enumerate(row)),count=51)
guard=AB(a*b*(a-1)*(b-1)*(a-b))
denominators={K(c).denominator() for f in sum(F,[])+C3+Q4 for c in f.coefficients()}
bad=[str(q) for q in denominators if any(guard%factor for factor,power in q.factor())]
check('all field and cubic/quartic spectral denominators protected',not bad,denominator_count=len(denominators),unprotected=bad)
# Audit the combined two-minor u reconstruction without assuming either minor nonzero.
Mu=matrix(R,[[f.derivative(v) for v in u] for f in Q4])
B0=Mu.matrix_from_rows((0,3)); B1=Mu.matrix_from_rows((12,16))
A=K(a*(b-1)); CC=K(b*(a-1)); DD=K(a-b)
X0=CC*p[0]*p[4];Y0=A*p[0]*p[8];Z0=DD*p[4]*p[8]
check('fnd:uleftinverse first determinant',B0.det()==2*(6*X0+Y0+Z0)/A)
check('fnd:uleftinverse second determinant',B1.det()==DD*(X0-Y0-Z0)/(2*A**2))
check('division free combined left inverse',A/2*B0.adjugate()*B0+2*A**2/DD*B1.adjugate()*B1==7*CC*p[0]*p[4]*identity_matrix(R,2))
P=PolynomialRing(K,names=names,order='lex'); morph=R.hom(list(P.gens())+[P.zero()]*len(extras),P)
Fp=[[morph(f) for f in row] for row in F]; Cp=list(map(morph,C3)); Qp=list(map(morph,Q4))
def sparse(f):return {tuple(m):K(c) for m,c in f.dict().items()}
# The indices are part of the published finite coefficient-source convention.
selectedQ=[0,1,3,5,6,10,11,12,13,15,18,19]
selectedD=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,21,22,25,27,30,31,34]
Dsel=[morph(D(i,Q4[j])) for j in selectedQ for i in range(3)]
G5=[Dsel[j] for j in selectedD]
defects=[[morph(D(i,F[j][k])-D(j,F[i][k])-sum(C(m,i,j)*F[m][k] for m in range(3))) for k in range(17)] for i,j in pairs]
prepared={'Fraw':[[sparse(f) for f in row] for row in Fp],'C3':list(map(sparse,Cp)),'Q4':list(map(sparse,Qp)),
          'basis':list(map(sparse,G5)),'Dsel':list(map(sparse,Dsel)),'top_defects':[sparse(row[k]) for row in defects for k in (15,16)],
          'source':'Fresh Section 2 Koszul/curvature/curl reconstruction. basis has ONLY the first 27 selected degree-five sources.',
          'selected_quartic_indices':selectedQ,'selected_derivative_indices':selectedD}
(OUT/'generic_degree5_sparse_prepared.sobj').write_bytes(pickle.dumps(prepared))
save({'ring':P,'F':Fp,'C3':Cp,'Q4':Qp,'G5':G5,'defects':defects,'Q':list(map(morph,Q)),'T':list(map(morph,Tsol)),'U':list(map(morph,Usol))},str(OUT/'state17_fresh.sobj'))
(OUT/'state17_fields.txt').write_text('\n'.join('D%d %s = %s'%(i+1,names[j],f) for i,row in enumerate(Fp) for j,f in enumerate(row)))
report=dict(status='PASS',checks=checks,elapsed_seconds=time.monotonic()-START,field_count=51,
            archived_coefficient_inputs=0,symbolic_parameter_field='QQ(a,b)',
            fresh_source_sha256=hashlib.sha256((OUT/'generic_degree5_sparse_prepared.sobj').read_bytes()).hexdigest())
(OUT/'state17_checks.json').write_text(json.dumps(report,indent=2))
print('COMPLETE',report['elapsed_seconds'],flush=True)
