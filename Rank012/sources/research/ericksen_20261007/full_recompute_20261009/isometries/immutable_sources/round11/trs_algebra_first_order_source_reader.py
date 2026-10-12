"""Independent direct tensor-jet replay of all four affine equations."""
from sage.all import *
from pathlib import Path
import json,time
OUT=Path(__file__).resolve().parent;START=time.time()
d=load(str(OUT/'trs_algebra_first_order_shape.sobj'));K=d['field'];v,p,j,A,T,B,L=K.gens()
def nk(z):
 z=K(z);n=z.numerator();de=z.denominator();g=n.gcd(de);n=n//g;de=de//g;lc=de.leading_coefficient()
 return K(n/lc)/K(de/lc)
AB,TB,Av,Tv=(d[x] for x in ('A_B','T_B','A_v','T_v'))
def dv(z):return nk(z.derivative(v)+z.derivative(A)*Av+z.derivative(T)*Tv)
def db(z):return nk(z.derivative(B)+z.derivative(A)*AB+z.derivative(T)*TB)
P=PolynomialRing(K,names=('xi','eta','Lq','Lz','KC','KB'));xi,eta,lq,lz,KC,KB=P.gens();U2=d['U_squared'];q2=v-p*p
def np(f):return P({e:nk(c) for e,c in P(f).dict().items()})
class E:
 def __init__(self,a=0,b=0):self.a=np(a);self.b=np(b)
 @staticmethod
 def c(y):return y if isinstance(y,E) else E(y)
 def __add__(x,y):
  y=E.c(y);return E(x.a+y.a,x.b+y.b)
 __radd__=__add__
 def __neg__(x):return E(-x.a,-x.b)
 def __sub__(x,y):return x+-E.c(y)
 def __mul__(x,y):
  y=E.c(y);return E(x.a*y.a+U2*x.b*y.b,x.a*y.b+x.b*y.a)
 __rmul__=__mul__
 def divU(x):return E(x.b,x.a/U2)
 def zero(x):return not x.a and not x.b
U=E(0,1);a=[U,E(0,q2*L/v)]
aq=[xi*U,E(0,2*p*p*L/v**2+q2*lq/v+q2*L*xi/v)]
az=[eta*U,E(0,q2*lz/v+q2*L*eta/v)]
kdiag=(1/q2,q2/v);kidiag=(q2,v/q2)
ss=nk(d['I2']/j-(1+B)/v)
aa=[[U2,q2*L*U2/v],[q2*L*U2/v,q2*q2*L*L*U2/v**2]]
# r=(-Ki22*a2,Ki11*a1); reconstruct original G, independently of saved S,R.
rr=[E(0,-L),E(0,q2)]
G=[]
for i in range(2):
 row=[]
 for z in range(2):
  cross=a[i]*rr[z]+rr[i]*a[z]
  assert not cross.b
  row.append(nk((A/j*kdiag[i] if i==z else 0)+v/B*(ss-2*A/j)*aa[i][z]-T/(j*B)*K(cross.a)))
 G.append(row)
H=[[j*G[1][1],-j*G[0][1]],[-j*G[1][0],j*G[0][0]]]
checks={'det_H':nk(H[0][0]*H[1][1]-H[0][1]**2-j)==0,
        'saved_S':nk(kdiag[1]*H[0][1]-d['S'])==0,
        'saved_R':nk(kdiag[1]*H[1][1]-d['R'])==0}
bq=np(2*B/q2+2*B/(v+L*L)+2*B*xi+2*B*L/(v+L*L)*lq)
bz=np(2*B*eta+2*B*L/(v+L*L)*lz)
def hq(z):return np(2*dv(z)+db(z)*bq+z.derivative(L)*lq)
def hz(z):return np(db(z)*bz+z.derivative(L)*lz)
Hq=[[hq(z) for z in row] for row in H];Hz=[[hz(z) for z in row] for row in H]
flux=[sum((v*H[i][z]*a[z] for z in range(2)),E()) for i in range(2)]
V=sum((2*H[0][z]*a[z]+v*Hq[0][z]*a[z]+v*H[0][z]*aq[z]+v*Hz[1][z]*a[z]+v*H[1][z]*az[z] for z in range(2)),E())
etacov=[p*a[0],p*a[1]+E(q2/v)]
etaq=[p*aq[0],p*aq[1]+E(2*p*p/v**2)]
etaz=[p*az[0],p*az[1]]
lap=sum((Hq[0][z]*etacov[z]+H[0][z]*etaq[z]+Hz[1][z]*etacov[z]+H[1][z]*etaz[z] for z in range(2)),E())
original=[(aq[1]-az[0]+E(2*p/v**2)).divU(),
          (2*a[0]+q2*aq[0]+v/q2*az[1]+E(KB)).divU(),
          (V-E(KC)).divU(),lap-E(p/v)*(V-E(KC))]
for i in range(4):
 target=E(sum(d['M'][i,z]*w for z,w in enumerate((xi,eta,lq,lz)))+d['rhs']['constant'][i]+KC*d['rhs']['KC'][i],
          d['rhs']['U'][i]+(d['rhs']['one_over_U'][i]+KC*d['rhs']['KC_over_U'][i]+KB*d['rhs']['KB_over_U'][i])/U2)
 res=original[i]-target
 checks['full_affine_original_source_'+str(i)]=res.zero()
assert all(checks.values()),checks
out={'status':'PASS','seconds':time.time()-START,'checks':checks,'scope':'All four full affine source equations independently rebuilt from original G/K/a tensor formulas; all parameters retained.'}
(OUT/'trs_algebra_first_order_source_reader.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
