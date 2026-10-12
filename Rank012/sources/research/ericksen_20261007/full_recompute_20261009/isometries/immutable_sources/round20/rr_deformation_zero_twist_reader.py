"""Whole zero-twist endpoint: complete first variations and protected order-two obstruction."""
from sage.all import *
from pathlib import Path
import json,signal
OUT=Path(__file__).resolve().parent;signal.alarm(60)
K=FractionField(PolynomialRing(QQ,names=('E','k')));E,k=K.gens()
R=LaurentPolynomialRing(K,names=('r','z','lr','A0','A1','A2','B1','W2','Wlog','Gamma'))
r,z,lr,A0,A1,A2,B,W,H,Gamma=R.gens()
dr=lambda f:R(f).derivative(r)+R(f).derivative(lr)/r
dz=lambda f:R(f).derivative(z)
C0=diagonal_matrix(R,[1,k*k,E]);Ci=diagonal_matrix(R,[1,1/(k*k),1/E])
def linear(x,y,w):
 return matrix(R,[[2*(x+r*dr(x)),k*r*dr(y),r*dz(x)+E*dr(w)],
 [k*r*dr(y),2*k*k*x,k*r*dz(y)],
 [r*dz(x)+E*dr(w),k*r*dz(y),2*E*dz(w)]])
def residual(T):
 vr=dr(T[0,0])+dz(T[0,2])+(T[0,0]-T[1,1])/r
 vt=dr(T[0,1])+dz(T[1,2])+2*T[0,1]/r
 vz=dr(T[0,2])+dz(T[2,2])+T[0,2]/r
 return vector(R,[dr(r*vt),dz(r*vt),dz(vr)-dr(vz)])
A=A0+A1*z+A2*z*z;y=Gamma*lr+B*z+A/r**2;w=W*r*r+H*lr
HH=linear(R(0),y,w)
checks={'first_constant_spectrum':HH.trace()==0 and (C0*HH).trace()==0 and (Ci*HH).trace()==0,
 'complete_first_source_curl':residual(HH)==0,'complete_first_inverse_curl':residual(-Ci*HH*Ci)==0}
x2=-r*r*(dr(y)**2/(k*k-1)+dz(y)**2/(k*k-E))/2
gz=-r*r*dz(y)**2/(2*(E-k*k))-E*dr(w)**2/(2*(E-1))
jac=2*x2+r*dr(x2)+gz
def radial(f,n):return sum(R({ex:cc}) for ex,cc in R(f).dict().items() if ex[0]==n)
checks['pole_squared_protected']=radial(jac,-4)==4*A*A/((k*k-1)*r**4)
j0=R({ex:cc for ex,cc in jac.dict().items() if not any(ex[j] for j in [3,4,5])})
checks['log_shear_protected']=radial(j0,-2)==-E*H*H/(2*(E-1)*r*r)
j1=R({ex:cc for ex,cc in j0.dict().items() if not ex[8]})
checks['fixed_angle_spin_protected']=radial(j1,0)==-Gamma*Gamma/(k*k-1)
P=3*(E-1)*B*B+4*E*(k*k-E)*W*W
checks['remaining_conic']=radial(j1,2)==-r*r*P/(2*(k*k-E)*(E-1))
# Explicit local formal order-two lift of every real point on the remaining conic.
y=B*z;w=W*r*r;HH=linear(R(0),y,w)
Q=matrix(R,[[r*r*dr(y)**2+E*dr(w)**2,0,r*r*dr(y)*dz(y)],
 [0,0,0],[r*r*dr(y)*dz(y),0,r*r*dz(y)**2]])
xx=-B*B*r*r/(2*(k*k-E))
ww=2*B*B*r*r*z/(k*k-E)
yy=E*B*W*(r*r-4*z*z)/(E-1)
TT=Q+linear(xx,yy,ww);TI=Ci*HH*Ci*HH*Ci-Ci*TT*Ci
S=PolynomialRing(K,names=('B','W'));bb,wwa=S.gens();pp=3*(E-1)*bb*bb+4*E*(k*k-E)*wwa*wwa
def modconic(f):
 groups={}
 for ex,cc in R(f).dict().items():
  key=tuple(ex[:3]);groups[key]=groups.get(key,S(0))+cc*bb**ex[6]*wwa**ex[7]
 return all(g.reduce([pp])==0 for g in groups.values())
checks['lift_second_trace']=modconic(TT.trace())
checks['lift_second_trace_square']=modconic(2*(C0*TT).trace()+(HH*HH).trace())
checks['lift_second_det']=modconic((Ci*TT).trace()-((Ci*HH)*(Ci*HH)).trace()/2)
checks['lift_full_second_source_curl']=all(modconic(f) for f in residual(TT))
checks['lift_full_second_inverse_curl']=all(modconic(f) for f in residual(TI))
assert all(checks.values()),{k:v for k,v in checks.items() if not v}
rep={'status':'WHOLE_ZERO_TWIST_FIXED_K_FORMAL_SECOND_ORDER_CLASSIFICATION',
 'checks':checks,'guards':'E>0,k!=0,k²!=1,E!=1,E!=k²; exactly positive simple spectrum at omega=0',
 'remaining_conic':'3(E-1) B1²+4E(k²-E) W2²=0',
 'scope':'A0,A1,A2,Wlog,Gamma vanish. B1,W2 survive precisely on this cone at formal order two; higher integrability and actual deformation are not established.'}
(OUT/'rr_deformation_zero_twist_reader.json').write_text(json.dumps(rep,indent=2),encoding='utf8')
save({'ring':R,'checks':checks,'jacobian_obstruction':jac,'remaining_conic':P,
 'lift_x2':xx,'lift_y2':yy,'lift_w2':ww,'second_metric':TT,'second_inverse_metric':TI},str(OUT/'rr_deformation_zero_twist_reader.sobj'))
print(json.dumps(rep,indent=2),flush=True)
