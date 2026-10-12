"""Regular-angle first-order system and its two exact Frobenius sources."""
from sage.all import *
from pathlib import Path
import resource,time,signal,json
OUT=Path(__file__).resolve().parent;START=time.time()
d=load(str(OUT/'trs_algebra_first_order_shape.sobj'))
vm=int(Path('/proc/self/statm').read_text().split()[0])*4096
resource.setrlimit(resource.RLIMIT_AS,(vm+550*1024**2,vm+550*1024**2));signal.alarm(240)
K=d['field'];v,p,j,A,T,B,L=K.gens()
def nk(z):
 z=K(z);n=z.numerator();de=z.denominator();g=n.gcd(de);n=n//g;de=de//g;lc=de.leading_coefficient()
 return K(n/lc)/K(de/lc)
AB,TB,Av,Tv=(d[x] for x in ('A_B','T_B','A_v','T_v'))
def dv(z):return nk(z.derivative(v)+z.derivative(A)*Av+z.derivative(T)*Tv)
def db(z):return nk(z.derivative(B)+z.derivative(A)*AB+z.derivative(T)*TB)
labels=list(d['rhs']);M=d['M'];aug=[list(M.row(i))+[-d['rhs'][key][i] for key in labels] for i in range(4)]
for c in range(4):
 pivot=next(i for i in range(c,4) if aug[i][c]);aug[c],aug[pivot]=aug[pivot],aug[c]
 den=aug[c][c];aug[c]=[nk(x/den) for x in aug[c]]
 for i in range(4):
  if i==c or not aug[i][c]:continue
  a=aug[i][c];aug[i]=[nk(x-a*y) for x,y in zip(aug[i],aug[c])]
 print('solve pivot',c,round(time.time()-START,2),flush=True)
sol={key:vector(K,[aug[i][4+c] for i in range(4)]) for c,key in enumerate(labels)}
for key,z in sol.items():assert all(nk(w)==0 for w in M*z+d['rhs'][key])
save({'source':d,'solution_coefficients':sol,'labels':labels},str(OUT/'trs_algebra_regular_first_order_solution.sobj'))
print('all six source columns replayed',round(time.time()-START,2),flush=True)
U2=d['U_squared']
def np(f):
 if isinstance(f,dict):return {e:nk(c) for e,c in f.items() if c}
 return {(0,0):nk(f)} if f else {}
def padd(x,y):
 out=np(x)
 for e,c in np(y).items():out[e]=nk(out.get(e,0)+c)
 return {e:c for e,c in out.items() if c}
def pneg(x):return {e:-c for e,c in np(x).items()}
def pmul(x,y):
 out={}
 for e,c in np(x).items():
  for ee,cc in np(y).items():
   key=tuple(a+b for a,b in zip(e,ee));out[key]=nk(out.get(key,0)+nk(c*cc))
 return {e:c for e,c in out.items() if c}
def pd(f,which):
 fn={'v':dv,'B':db,'L':lambda x:nk(x.derivative(L))}[which]
 return {e:fn(c) for e,c in np(f).items()}
class E:
 def __init__(self,a=0,b=0):self.a=np(a);self.b=np(b)
 @staticmethod
 def c(y):return y if isinstance(y,E) else E(y)
 def __add__(x,y):
  y=E.c(y);return E(padd(x.a,y.a),padd(x.b,y.b))
 __radd__=__add__
 def __neg__(x):return E(pneg(x.a),pneg(x.b))
 def __sub__(x,y):return x+-E.c(y)
 def __mul__(x,y):
  y=E.c(y);return E(padd(pmul(x.a,y.a),pmul(U2,pmul(x.b,y.b))),padd(pmul(x.a,y.b),pmul(x.b,y.a)))
 __rmul__=__mul__
 def part(x,key):return E(pd(x.a,key),pd(x.b,key))
 def payload(x):return (x.a,x.b)
U=E(0,1)
xx=[]
for i in range(4):
 xx.append(E({(0,0):sol['constant'][i],(1,0):sol['KC'][i]},
             {(0,0):sol['U'][i]+sol['one_over_U'][i]/U2,(1,0):sol['KC_over_U'][i]/U2,(0,1):sol['KB_over_U'][i]/U2}))
xi,eta,lq,lz=xx
bq=E(2*B/(v-p*p)+2*B/(v+L*L))+2*B*xi+E(2*B*L/(v+L*L))*lq
bz=2*B*eta+E(2*B*L/(v+L*L))*lz
def dx(z):return 2*z.part('v')+bq*z.part('B')+lq*z.part('L')+U*E(z.b)*xi
def dy(z):return bz*z.part('B')+lz*z.part('L')+U*E(z.b)*eta
print('starting commutator 0',round(time.time()-START,2),flush=True)
C0=dx(eta)-dy(xi)
save({'source':d,'pressure_variables':('KC','KB'),'frobenius':[C0.payload()]},str(OUT/'trs_algebra_regular_frobenius_first.sobj'))
print('commutator 0 saved',round(time.time()-START,2),flush=True)
C1=dx(lz)-dy(lq)
save({'source':d,'pressure_variables':('KC','KB'),'U_squared':U2,'solution':[z.payload() for z in xx],'B_q':bq.payload(),'B_zeta':bz.payload(),'frobenius':[C0.payload(),C1.payload()]},str(OUT/'trs_algebra_regular_frobenius.sobj'))
signal.alarm(0)
summary={'status':'EXACT_FROBENIUS_SOURCES','seconds':time.time()-START,'solution_source_checks':24,
 'coefficient_counts':[[len(z.a),len(z.b)] for z in (C0,C1)],
 'base_numerator_terms':[[sum(len(c.numerator().dict()) for c in f.values()) for f in z.payload()] for z in (C0,C1)],
 'scope':'Necessary exact sources only on W!=0, a1!=0. No finite phase fiber or singular-wall exclusion yet.'}
(OUT/'trs_algebra_regular_frobenius.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
