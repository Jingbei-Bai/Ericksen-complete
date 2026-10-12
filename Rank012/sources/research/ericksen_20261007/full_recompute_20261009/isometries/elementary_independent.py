"""Independent exact reconstruction of the small all-parameter Section 7 identities."""
from pathlib import Path
import json,time
import sympy as S
OUT=Path(__file__).absolute().parent;start=time.time();checks={};data={}
def ck(name,expr):
    ans=S.factor(expr);checks[name]=ans==0
    if ans!=0:data[name+'_residual']=str(ans)
s,t,p,q,c,v,B,j,I1,I2,lam1,lam2,lam3=S.symbols('s t p q c v B j I1 I2 lam1 lam2 lam3')
g11,g12,g22,v1,v2,kap=S.symbols('g11 g12 g22 v1 v2 kap')
GG=S.Matrix([[g11,g12],[g12,g22]]);vv=S.Matrix([v1,v2]);GI=GG.inv()
CTT=(GG+vv*vv.T).row_join(kap*vv).col_join((kap*vv.T).row_join(S.Matrix([[kap**2]])))
ITT=GI.row_join(-GI*vv/kap).col_join((-vv.T*GI/kap).row_join(S.Matrix([[(1+(vv.T*GI*vv)[0])/kap**2]])))
checks['TT_full_block_inverse']=all(S.factor(x)==0 for x in CTT*ITT-S.eye(3))
ck('TT_full_block_determinant',CTT.det()-kap**2*GG.det())
RR,uq,uz=S.symbols('RR uq uz')
FT=S.Matrix([[1/RR,0,0],[RR*uq,RR*uz,RR],[p*uq,1+p*uz,p]])
ck('TRTS_coordinate_Jacobian',FT.det()+1)
aa=S.Matrix([uq,uz+p/(RR**2+p**2)])
KS=S.diag(1/RR**2,RR**2/(RR**2+p**2));vt=RR**2+p**2
cov=(KS+vt*aa*aa.T).row_join(vt*aa).col_join((vt*aa.T).row_join(S.Matrix([[vt]])))
checks['TRTS_complete_covariant_metric']=all(S.factor(x)==0 for x in FT.T*FT-cov)
ck('TRTS_covariant_determinant',cov.det()-1)
F=S.Matrix([[s,0,0],[0,t,0],[p,q,c]]);C=F.T*F;ei=S.Matrix([0,0,1]);CI=C.inv()
ck('RR_physical_chart_direction_determinant',S.Matrix.hstack(ei,CI*ei,CI*CI*ei).det()-p*q*(s*s-t*t)/(s**4*t**4*c*c))
P=v**3-I1*v**2+I2*v-j;A=j/v-P/(v*B);R=I2-A-j*(1+B)/v
theta=v*A*R-j*v*v
ck('TRTS_theta_full_expansion',theta-(-j*j*B/v+2*j*I2-j*I1*v-3*j*j/v+P*(3*j/v-I2)/B-P*P/(v*B*B)))
subst={I1:lam1+lam2+lam3,I2:lam1*lam2+lam1*lam3+lam2*lam3,j:lam1*lam2*lam3}
ck('TRTS_spectral_moment_product',theta.subs(subst)+S.prod((v-u)*(v-w)+u*w*B for u,w in [(lam2,lam3),(lam1,lam3),(lam1,lam2)])/(v*B*B))
Q=v-p*p
ck('TRTS_Lzero_pressure_norm',4*p*p+(v-2*p*p)**2/Q-v*v/Q)
ck('TRTS_Lzero_constant_norm',(v-2*p*p)**2+4*p*p*Q-v*v)
ck('TRTS_wall_pressure_norm',p*p*(3*v-4*p*p)**2/Q**2+(v-4*p*p)**2/Q-v**3/Q**2)
ck('TRTS_wall_constant_norm',(v-4*p*p)**2+p*p*(3*v-4*p*p)**2/Q-v**3/Q)
L=S.symbols('L');D=L*L+v
bE=-p*(p*p*L*L+4*v*p*p-3*v*v);bO=-Q*(3*p*p*L*L+4*v*p*p-v*v)
ck('TRTS_finite_angle_KC_norm',Q*D*bE*bE+v*bO*bO-Q*(p*p*L*L+v*v)**3)
M=p*p*L*L+2*v*p*p-v*v
ck('TRTS_finite_angle_KCzero_norm',M*M+4*p*p*Q*v*D-(p*p*L*L+v*v)**2)
SS,TT,PP,QQ,CC=S.symbols('SS TT PP QQ CC');de=SS-TT
f20=de*(TT*de*PP+TT*PP**2+5*SS*de*QQ+5*(SS+TT)*PP*QQ-3*SS*QQ**2-de*(PP+5*QQ)*CC)
f13=de*(-5*TT*de*PP-3*TT*PP**2-SS*de*QQ+5*(SS+TT)*PP*QQ+SS*QQ**2+de*(5*PP+QQ)*CC)
psi=TT*PP**3+(15*SS+8*TT)*PP**2*QQ+(8*SS+15*TT)*PP*QQ**2+SS*QQ**3+12*(SS-TT)**2*PP*QQ
ck('RR_bad_fibre_positive_combination',(PP+5*QQ)*f13+(5*PP+QQ)*f20-2*de*psi)
E,k,w=S.symbols('E k w');Delta=E*E-E*(1+w*w+k*k)+k*k
Z=3*(E-1)*(k*k-1)-(E+1)*w*w
T=(E-1)*(k*k-E)*(k*k-1)+w*w*((5*E+1)*k*k-E*(E+5))
Hstar=(E-1)*(k*k-1)+(E+1)*w*w
# Identity modulo Z, without dividing on E=1, E=k^2, or Z=0.
numer=S.expand((E-1)*T+Delta*Hstar)
quot,rem=S.div(numer,Z,w)
ck('RR_conic_nonzero_on_Zwall',rem)
data['RR_Zwall_exact_division_quotient']=str(S.factor(quot))
alpha,beta,u,vs,Kk=S.symbols('alpha beta u vs Kk')
Fp=S.Matrix([[alpha,0],[alpha*beta,alpha*k]])
Cyl=Fp.T*Fp
subc={alpha**2:u*vs/Kk,k**2:Kk**2,beta**2:(Kk-u/vs)*(vs/u-Kk)}
ck('cylinder_planar_trace',S.expand(Cyl.trace()).subs(subc)-(u*u+vs*vs))
ck('cylinder_planar_determinant',S.expand(Cyl.det()).subs({alpha**4:u*u*vs*vs/(Kk*Kk),k**2:Kk*Kk})-u*u*vs*vs)
r=S.symbols('r',positive=True);a,b0,cc,z=S.symbols('a b0 cc z')
stress_div=S.Matrix([(a-cc)/r,2*b0/r,0])
phi=(a-cc)*S.log(r)+2*b0*S.Symbol('theta')
ck('cylinder_local_pressure_r',S.diff(phi,r)-stress_div[0])
ck('cylinder_local_pressure_theta',S.diff(phi,S.Symbol('theta'))/r-stress_div[1])
x1,x2,x3,aa,tau=S.symbols('x1 x2 x3 aa tau');xx=S.Matrix([x1,x2,x3]);co=S.symbols('q0:6')
qform=co[0]*x1*x1+co[1]*x1*x2+co[2]*x1*x3+co[3]*x2*x2+co[4]*x2*x3+co[5]*x3*x3
xt=S.exp(-tau)*xx;at=S.exp(-2*tau)*(aa+tau*qform)
ck('SS_resonant_stable_flow',S.diff(at,tau)+2*at-qform.subs(dict(zip(xx,xt)),simultaneous=True))
report={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'data':data,'seconds':time.time()-start,'scope':'all variables symbolic; no material specialization and no floating-point arithmetic'}
(OUT/'elementary_independent.json').write_text(json.dumps(report,indent=2),encoding='utf8')
print(json.dumps(report,indent=2))
