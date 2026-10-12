"""Formal B-infinity coefficient probe for the finite-angle wall compatibility.
Exploratory until matched to exact rational source/cofactor certificates.
"""
import os,resource,time,json
from pathlib import Path
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from sage.all import *
root=Path(__file__).resolve().parent;start=time.time()
R0=PolynomialRing(QQ,'ii');ii=R0.gen();NF=NumberField(ii**2+1,'I');I=NF.gen()
PR=PolynomialRing(NF,names=('w','p','j','I1','I2','KC','KB'));K=PR.fraction_field()
w,p,j,I1,I2,KC,KB=map(K,PR.gens());v=w*w;q=v-p*p
PREC=10
SR=LaurentSeriesRing(K,'e',default_prec=PREC);e=SR.gen();B=e**(-2)
vm=int(Path('/proc/self/statm').read_text().split()[0])*4096
resource.setrlimit(resource.RLIMIT_AS,(vm+550*1024**2,vm+550*1024**2))
def nk(x):
 x=K(x);n=x.numerator();d=x.denominator();g=n.gcd(d);n=n//g;d=d//g;lc=d.leading_coefficient()
 return K(n/lc)/K(d/lc)
def cl(x):
 x=SR(x);prec=min(PREC,x.precision_absolute())
 if not x:return SR(0).add_bigoh(prec)
 return sum(SR(nk(x[i]))*e**i for i in range(int(x.valuation()),int(prec)) if x[i]).add_bigoh(prec)
def dv(x):
 x=SR(x);prec=min(PREC,x.precision_absolute())
 if not x:return SR(0).add_bigoh(prec)
 return cl(sum(SR(nk(x[i].derivative(w)/(2*w)))*e**i for i in range(int(x.valuation()),int(prec)) if x[i]).add_bigoh(prec))
def db(x):return cl(-e**3*SR(x).derivative()/2)
def sqrtone(x):
 z=cl(x-1);return cl(sum(binomial(QQ(1)/2,k)*z**k for k in range(PREC)))
P=v**3-I1*v*v+I2*v-j
A=cl(j/v-P/(v*B));R=cl(I2-A-j*(1+B)/v)
Theta=cl(v*A*R-j*v*v)
T=cl(I*j/(w*e)*sqrtone(Theta/(-j*j*B/v)))
AB,TB,Av,Tv,Rb,Rv=db(A),db(T),dv(A),dv(T),db(R),dv(R)
alpha=cl((2*T*AB+(R-A)*TB)/(v*B))
beta=cl(2*AB*AB+2*TB*TB/v+j*(2*j/v-A)/(v*B))
L=cl(-beta/alpha);D=cl(v+L*L);U2=cl(B/(q*D))
F=cl(A-T*L/v);G=cl(A*L+T)
Fb,Fl,Fv=cl(AB-TB*L/v),cl(-T/v),cl(Av-Tv*L/v+T*L/v**2)
Gb,Gl=cl(AB*L+TB),A
NN=cl(T*(v-L*L)-v*(R-A)*L)
S=cl(q*NN/(v*v*D))
Sb=cl(q*(TB*(v-L*L)-v*(Rb-AB)*L)/(v*v*D))
Sl=cl(q*((-2*T*L-v*(R-A))*D-2*L*NN)/(v*v*D*D))
Nv=cl(Tv*(v-L*L)+T-(R-A)*L-v*(Rv-Av)*L)
Sv=cl(NN/(v*v*D)+q*Nv/(v*v*D)-q*NN*(2*v*D+v*v)/(v**4*D*D))
Rbph=cl((v*Rb+AB*L*L+2*TB*L)/(v*D))
Rl=cl(((2*A*L+2*T)*D-2*L*(v*R+A*L*L+2*T*L))/(v*D*D))
bq=vector(SR,[2*B,0,cl(2*B*L/D),0]);bz=vector(SR,[0,2*B,0,cl(2*B*L/D)])
bq0=cl(2*B/q+2*B/D)
M=matrix(SR,4,4)
M[0,:]=matrix(SR,[[q*L/v,-1,q/v,0]])
M[1,:]=matrix(SR,[[q,L,0,1]])
M[2,:]=matrix(SR,q*(vector(SR,[F,0,Fl,0])+Fb*bq)+vector(SR,[0,G,0,Gl])+Gb*bz)
M[3,:]=matrix(SR,Sb*bq+vector(SR,[0,0,Sl,Rl])+Rbph*bz)
M=M.apply_map(cl)
rhs={
'constant':vector(SR,[2*p*p*L/v**2,2,2*F+q*(2*Fv+Fb*bq0),2*Sv+Sb*bq0]),
'KC_over_U':vector(SR,[0,0,-1,0]),
'KB_over_U':vector(SR,[0,1,0,0]),
'one_over_U':vector(SR,[2*p/v**2,0,0,0]),
'KC':vector(SR,[0,0,0,p/v]),
'U':vector(SR,[0,0,0,-2*p*q*F/v**2])}
E=M[:2,:2];Einv=matrix(SR,[[E[1,1],-E[0,1]],[-E[1,0],E[0,0]]])/(q*D/v)
N=(M[2:,2:]-M[2:,:2]*Einv*M[:2,2:]).apply_map(cl)
rred={key:(vec[2:]-M[2:,:2]*Einv*vec[:2]).apply_map(cl) for key,vec in rhs.items()}
print('SCHUR ready',time.time()-start,flush=True)
for col in (0,1):
 c={key:cl(N[0,col]*vec[1]-N[1,col]*vec[0]) for key,vec in rred.items()}
 ac=cl(c['constant']+KC*c['KC'])
 dc=cl(KC*c['KC_over_U']+KB*c['KB_over_U']+c['one_over_U']+c['U']*U2)
 norm=cl(dc*dc-U2*ac*ac)
 print('col',col,'precision',norm.precision_absolute(),'valuation',norm.valuation(),flush=True)
 for label,sub in [('all',{}),('KC0',{KC:0}),('KC_KB0',{KC:0,KB:0})]:
  for exponent in range(int(norm.valuation()),int(norm.precision_absolute())):
   coeff=nk(norm[exponent].subs(sub))
   if coeff:
    print('COEFF',col,label,exponent,coeff.factor(),flush=True)
    break
  else:print('ZERO_KNOWN_PART',col,label,flush=True)
save({'N':N,'rhs':rred,'rhs_raw':rhs,'M':M,'L':L,'U2':U2,'A':A,'R':R,'T':T,'Theta':Theta,'B':B,'PREC':PREC},str(root/'trs_angle_wall_asymptotic_probe.sobj'))
