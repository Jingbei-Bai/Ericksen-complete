"""Linearize about a cylinder with arbitrary functions, before selecting modes."""
from pathlib import Path
import sympy as s,json,time
O=Path(__file__).absolute().parent;start=time.time()
r,z,E,k,w=s.symbols('r z E k omega',nonzero=True)
A=s.Function('A')(z);B=s.Function('B')(z);W=s.Function('W')(r)
gm=(1-w*w-k*k)/w;gp=(1+w*w-k*k)/w
x=A/r**2;y=-gm*A/(2*r*r)+B
C=s.Matrix([[1+w*w,w*k,0],[w*k,k*k,0],[0,0,E]])
CI=s.Matrix([[1,-w/k,0],[-w/k,(1+w*w)/(k*k),0],[0,0,1/E]])
def linear(x,y,W,om):
 return s.Matrix([[2*((1+om*om)*x+r*s.diff(x,r)+om*r*s.diff(y,r)),2*om*k*x+k*r*s.diff(y,r),r*(s.diff(x,z)+om*s.diff(y,z))+E*s.diff(W,r)],[2*om*k*x+k*r*s.diff(y,r),2*k*k*x,k*r*s.diff(y,z)],[r*(s.diff(x,z)+om*s.diff(y,z))+E*s.diff(W,r),k*r*s.diff(y,z),2*E*s.diff(W,z)]])
def residual(H):
 dr=s.diff(H[0,0],r)+s.diff(H[0,2],z)+(H[0,0]-H[1,1])/r
 dt=s.diff(H[0,1],r)+s.diff(H[1,2],z)+2*H[0,1]/r
 dz=s.diff(H[0,2],r)+s.diff(H[2,2],z)+H[0,2]/r
 return [s.factor(s.diff(r*dt,r)),s.factor(s.diff(r*dt,z)),s.factor(s.diff(dr,z)-s.diff(dz,r))],s.factor(r*dt)
H=linear(x,y,W,w);RC,KC=residual(H);RI,KB=residual(-CI*H*CI)
Y=s.Function('Y')(r,z);Hz=linear(0,Y,W,0);Ciz=s.diag(1,1/k**2,1/E)
RZ,KZ=residual(Hz);RIZ,KIZ=residual(-Ciz*Hz*Ciz)
checks={'twist_spectral_trace':s.factor(s.trace(H))==0,'twist_spectral_trace_square':s.factor(s.trace(C*H))==0,'twist_spectral_logdet':s.factor(s.trace(CI*H))==0,
'twist_B_second_derivative':s.factor(RC[0]-2*k*r*s.diff(B,z,2))==0,
'twist_inverse_B_second_derivative':s.factor(RI[0]+2*r*s.diff(B,z,2)/(E*k))==0,
'twist_two_A_third_coefficients':s.factor(gp-gm-2*w)==0}
out={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'twist':{'forward_curls':list(map(str,RC)),'inverse_curls':list(map(str,RI)),'forward_pressure':str(KC),'inverse_pressure':str(KB)},'zero_twist':{'forward_curls':list(map(str,RZ)),'inverse_curls':list(map(str,RIZ)),'forward_pressure':str(KZ),'inverse_pressure':str(KIZ)},'seconds':time.time()-start}
(O/'cylinder_arbitrary_linear.json').write_text(json.dumps(out,indent=2),encoding='utf8');print(json.dumps(out,indent=2))
