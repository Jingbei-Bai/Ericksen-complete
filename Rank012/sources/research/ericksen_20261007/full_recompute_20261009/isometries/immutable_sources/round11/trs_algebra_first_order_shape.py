"""Full four first-order phase equations after harmonic-eta source.

All simple material spectra retained through birational coordinates A,T.
k=1 is a recoverable common source/target Euclidean scale choice.
a1=U!=0 is an explicitly retained phase chart, not global assumption.
"""
from sage.all import *
from pathlib import Path
import time,resource,json
START=time.time();OUT=Path(__file__).resolve().parent
vm=int(Path('/proc/self/statm').read_text().split()[0])*4096
resource.setrlimit(resource.RLIMIT_AS,(vm+550*1024**2,vm+550*1024**2))
R=PolynomialRing(QQ,names=('v','p','j','A','T','B','L'));K=R.fraction_field()
v,p,j,A,T,B,L=map(K,R.gens())
def nk(z):
 z=K(z);n=z.numerator();d=z.denominator();g=n.gcd(d);n=n//g;d=d//g
 lc=d.leading_coefficient();return K(n/lc)/K(d/lc)
Rshape=nk((T*T+j*v*v)/(v*A));s=nk((Rshape+A)/j)
I2=nk(j*s+j*(1+B)/v);P=nk(B*(j-v*A));I1=nk((v**3+v*I2-j-P)/v**2)
AB=nk((j/v-A)/B)
TB=nk((v*(Rshape-A)*AB-j*A)/(2*T))
Pv=nk(3*v*v-2*I1*v+I2)
Av=nk(-j/v**2-(v*Pv-P)/(v*v*B))
Tv=nk((I2*A-A*A-2*j*v+v*(Rshape-A)*Av)/(2*T))
def dv(z):return nk(z.derivative(v)+z.derivative(A)*Av+z.derivative(T)*Tv)
def db(z):return nk(z.derivative(B)+z.derivative(A)*AB+z.derivative(T)*TB)
assert dv(I1)==dv(I2)==db(I1)==db(I2)==0
q2=v-p*p;D=v+L*L
U2=nk(B/(q2*D))
F=nk(A-T*L/v);G=nk(A*L+T)
S=nk(q2*(T*(v-L*L)-v*(Rshape-A)*L)/(v*v*D))
RR=nk((v*Rshape+A*L*L+2*T*L)/(v*D))
Sb,Sl,Sv=db(S),nk(S.derivative(L)),dv(S)
Rb,Rl=db(RR),nk(RR.derivative(L))
Fb,Fl,Fv=db(F),nk(F.derivative(L)),dv(F)
Gb,Gl=db(G),nk(G.derivative(L))
bq0=nk(2*B/q2+2*B/D)
bq=vector(K,[2*B,0,2*B*L/D,0]);bz=vector(K,[0,2*B,0,2*B*L/D])
M=matrix(K,4,4)
M[0,:]=matrix(K,[[q2*L/v,-1,q2/v,0]])
M[1,:]=matrix(K,[[q2,L,0,1]])
M[2,:]=matrix(K,[q2*(vector(K,[F,0,Fl,0])+Fb*bq)+vector(K,[0,G,0,Gl])+Gb*bz])
M[3,:]=matrix(K,[Sb*bq+vector(K,[0,0,Sl,Rl])+Rb*bz])
M=M.apply_map(nk)
# Equation M*(Uq/U,Uz/U,Lq,Lz)+ rhs0+rhsKC*KC+rhsKB*KB+rhsU*U+rhsInvU/U=0.
# Pressure/1-U products are encoded by separate labelled columns.
rhs={
 'constant':vector(K,[2*p*p*L/v**2,2,2*F+q2*(2*Fv+Fb*bq0),2*Sv+Sb*bq0]),
 'KC_over_U':vector(K,[0,0,-1,0]),
 'KB_over_U':vector(K,[0,1,0,0]),
 'one_over_U':vector(K,[2*p/v**2,0,0,0]),
 'KC':vector(K,[0,0,0,p/v]),
 'U':vector(K,[0,0,0,-2*p*q2*F/v**2]),
}
rhs={key:vec.apply_map(nk) for key,vec in rhs.items()}
print('matrix ready',round(time.time()-START,2),flush=True)
# Explicit 4x4 Laplace expansion avoids Sage quotient-field rank heuristics.
det=nk(sum((-1)**i*M[0,i]*M.delete_rows([0]).delete_columns([i]).determinant() for i in range(4)))
wall=nk(L*(2*T*AB+(Rshape-A)*TB)/(v*B)+2*AB**2+2*TB**2/v+j/(v*B)*(2*j/v-A))
ratio=nk(det/wall)
print('determinant ratio',ratio.factor(),flush=True)
save({'field':K,'M':M,'rhs':rhs,'U_squared':U2,'I1':I1,'I2':I2,'P':P,'A_B':AB,'T_B':TB,'A_v':Av,'T_v':Tv,'wall':wall,'det':det,'det_over_wall':ratio,'S':S,'R':RR,'forward_flux_factors':(F,G)},str(OUT/'trs_algebra_first_order_shape.sobj'))
summary={'status':'EXACT_SYSTEM_BUILT_NEEDS_COORDINATE_REPLAY','seconds':time.time()-START,
 'det_factor':str(det.factor()),'wall_factor':str(wall.factor()),'det_over_wall':str(ratio.factor()),
 'chart':'k=1, a1=U!=0; U^2=B/[(v-p^2)(v+L^2)], a2=(v-p^2)L U/v; v>p^2, B>0, A>0, T!=0',
 'unknown_derivatives':['U_q/U','U_zeta/U','L_q','L_zeta'],
 'rows':['phase curl','inverse vertical','forward vertical','harmonic eta reduced by forward vertical'],
 'scope':'All physical constants retained. det=0 wall must remain; no finite phase fiber theorem yet.'}
(OUT/'trs_algebra_first_order_shape.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2),flush=True)
