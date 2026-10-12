"""Full cylindrical R-R family: exact second-order infinitesimal integrability.

Target uniform scale normalizes R0=r. Nonzero omega chart only; zero twist
is a separate physical endpoint, never discarded from whole classification.
"""
from sage.all import *
from pathlib import Path
import signal,time,json,resource
OUT=Path(__file__).resolve().parent;st=time.time();signal.alarm(170)
vm=int(Path('/proc/self/statm').read_text().split()[0])*4096
resource.setrlimit(resource.RLIMIT_AS,(vm+1800*1024**2,vm+1800*1024**2))
K=FractionField(PolynomialRing(QQ,names=('E','k','omega')));E,k,om=K.gens()
R=LaurentPolynomialRing(K,names=('r','z','lr','A0','A1','A2','B1','W2','Wlog'))
r,z,lr,A0,A1,A2,B1,W2,Wlog=R.gens()
dr=lambda f:R(f).derivative(r)+R(f).derivative(lr)/r
dz=lambda f:R(f).derivative(z)
def integrate(f,var):
 j=R.gens().index(var);ans=R(0)
 for ex,cc in R(f).dict().items():
  ee=list(ex);nn=ee[j]
  if j==0:
   q=ee[2];ee[0]=0;ee[2]=0;rest=R({tuple(ee):cc})
   if nn==-1:ans+=rest*lr**(q+1)/(q+1)
   else:
    ans+=rest*r**(nn+1)*sum(QQ((-1)**h)*factorial(q)/QQ(factorial(q-h))*lr**(q-h)/QQ(nn+1)**(h+1) for h in range(q+1))
  else:
   assert nn!=-1;ee[j]+=1;ans+=R({tuple(ee):cc/(nn+1)})
 return ans
L=1+om**2+k**2;beta=(1-om**2-k**2)/om;gamma=(1+om**2-k**2)/om
C0=matrix(K,[[1+om**2,om*k,0],[om*k,k**2,0],[0,0,E]])
Ci=matrix(K,[[1,-om/k,0],[-om/k,(1+om**2)/k**2,0],[0,0,1/E]])
assert C0*Ci==1
Linv=(E*identity_matrix(K,2)-C0[:2,:2]).inverse().change_ring(R)
C0=C0.change_ring(R);Ci=Ci.change_ring(R)
A=A0+A1*z+A2*z**2
x=A/r**2;y=-beta*A/(2*r**2)+B1*z;w=W2*r**2+Wlog*lr
def linear(x,y,w):
 return matrix(R,[[2*((1+om**2)*x+r*dr(x)+om*r*dr(y)),2*om*k*x+k*r*dr(y),r*(dz(x)+om*dz(y))+E*dr(w)],
                  [2*om*k*x+k*r*dr(y),2*k**2*x,k*r*dz(y)],
                  [r*(dz(x)+om*dz(y))+E*dr(w),k*r*dz(y),2*E*dz(w)]])
H=linear(x,y,w)
Q=matrix(R,3,3,0)
Q[0,0]=(x+r*dr(x))**2+om**2*x**2+4*om*r*x*dr(y)+r**2*dr(y)**2+E*dr(w)**2
Q[0,1]=Q[1,0]=k*(om*x*x+2*r*x*dr(y))
Q[1,1]=k*k*x*x
Q[0,2]=Q[2,0]=r*(x+r*dr(x))*dz(x)+r*r*dr(y)*dz(y)+2*om*r*x*dz(y)+E*dr(w)*dz(w)
Q[1,2]=Q[2,1]=2*k*r*x*dz(y)
Q[2,2]=r*r*(dz(x)**2+dz(y)**2)+E*dz(w)**2
Cross=vector(R,[H[0,2],H[1,2]])
G=-(Q[2,2]+Cross*Linv*Cross)/(2*E)
w2=integrate(G,z)
fJ=x*x+r*dz(x)*dr(w)
x2=integrate(r*(fJ-G),r)/r**2
y2=integrate((-L*x2-r*dr(x2)-(E*G)-Q.trace()/2)/(om*r),r)
T=Q+linear(x2,y2,w2)
checks={'first_trace':H.trace()==0,'first_trace_square':(C0*H).trace()==0,
 'first_det':(Ci*H).trace()==0,
 'w2_z':dz(w2)==G,'second_trace':T.trace()==0,
 'second_trace_square':2*(C0*T).trace()+(H*H).trace()==0,
 'second_det':(Ci*T).trace()-((Ci*H)*(Ci*H)).trace()/2==0}
assert all(checks.values()),checks
Ti=Ci*H*Ci*H*Ci-Ci*T*Ci
def residuals(T):
 vr=dr(T[0,0])+dz(T[0,2])+(T[0,0]-T[1,1])/r
 vt=dr(T[0,1])+dz(T[1,2])+2*T[0,1]/r
 vz=dr(T[0,2])+dz(T[2,2])+T[0,2]/r
 return vector(R,[dr(r*vt),dz(r*vt),dz(vr)-dr(vz)])
checks['full_first_source_curl']=residuals(H)==0
checks['full_first_inverse_curl']=residuals(-Ci*H*Ci)==0
assert all(checks.values()),checks
RC=residuals(T);RI=residuals(Ti)
Bsecond=-RC[0]/(2*k*r)
Athird=-(RC[1]+E*k*k*RI[1])/(k*om)
Wlapderiv=Athird/(E*r)-RI[2]
obs=[RC[0]+E*k*k*RI[0],dr(Bsecond),dr(Athird),
     -k*beta*Athird/2+k*r*r*dz(Bsecond)+RC[1],
     RC[2]+E*RI[2]-om*beta*Athird/(2*r)+om*r*dz(Bsecond),
     dz(Wlapderiv)]
print('source integrability derived',time.time()-st,flush=True)
# Collect by r,z,logr; coefficients are homogeneous quadratics in six amplitudes.
P=PolynomialRing(K,names=('A0','A1','A2','B1','W2','Wlog'));amps=P.gens()
co=[];locations=[]
for ix,f in enumerate(obs):
 dc={}
 for ex,cc in f.dict().items():
  key=tuple(ex[:3]);mon=prod(amps[j]**ex[j+3] for j in range(6))
  dc[key]=dc.get(key,P(0))+cc*mon
 for key,g in sorted(dc.items()):
  if g:co.append(g);locations.append([ix,key])
mons=[prod(amps[j]**ex[j] for j in range(6)) for ex in IntegerVectors(2,6)]
MM=matrix(K,[[g.monomial_coefficient(m) for m in mons] for g in co])
print('coefficient matrix',MM.dimensions(),time.time()-st,flush=True)
save({'ring':R,'C0':C0,'Cinv0':Ci,'x1':x,'y1':y,'w1':w,'H':H,'Q':Q,'x2_particular':x2,
 'y2_particular':y2,'w2_particular':w2,'T':T,'Tinv':Ti,'source_residual':RC,'inverse_residual':RI,
 'obstructions':obs,'coefficient_ring':P,'coefficients':co,'locations':locations,'coefficient_matrix':MM,
 'quadratic_monomials':mons,'checks':checks},str(OUT/'rr_deformation_second_order.sobj'))
rr=MM.echelon_form();basis=[]
for row in rr.rows():
 if row: basis.append(sum(row[j]*mons[j] for j in range(len(mons))))
save({'ring':P,'basis':basis,'source_coefficients':co,'matrix':MM,'checks':checks},str(OUT/'rr_deformation_second_order_relations.sobj'))
rep={'status':'EXACT_NONZERO_TWIST_SECOND_ORDER_COMPATIBILITY_NOT_WHOLE_RR_CLASSIFICATION',
 'checks':checks,'source_terms':[len(f.dict()) for f in obs],
 'coefficient_count':len(co),'quadratic_rank_over_parameter_field':len(basis),
 'basis_over_parameter_field':[str(f.factor()) for f in basis],
 'scope':'Parameter-field echelon divisions are not yet globally protected. All exact original coefficient equations are retained.',
 'seconds':time.time()-st}
(OUT/'rr_deformation_second_order.json').write_text(json.dumps(rep,indent=2),encoding='utf8')
print(json.dumps(rep,indent=2),flush=True)
