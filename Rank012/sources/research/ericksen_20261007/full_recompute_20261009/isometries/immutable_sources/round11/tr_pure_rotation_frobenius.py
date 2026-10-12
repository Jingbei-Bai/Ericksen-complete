"""Entire p=0 first-order system, outside the explicitly retained det wall.

Pressure constants C=KC,D=KB are not specialized. Source equations are the
phase curl, both genuine vertical equations, and newly proved Delta eta=0.
"""
from sage.all import *
from pathlib import Path
import time,resource,json,sys
START=time.time();OUT=Path(__file__).resolve().parent
dat=load(str(OUT/'trs_algebra_first_order_shape.sobj'))
vm=int(Path('/proc/self/statm').read_text().split()[0])*4096
resource.setrlimit(resource.RLIMIT_AS,(vm+650*1024**2,vm+650*1024**2))
P=PolynomialRing(QQ,names=('v','j','A','T','B','L'));K=P.fraction_field()
v,j,A,T,B,L=map(K,P.gens());old=dat['field'];hom=old.hom([v,0,j,A,T,B,L],K)
def nk(z):
 z=K(z);n=z.numerator();d=z.denominator();g=n.gcd(d);return K(n//g)/K(d//g)
M=matrix(K,4,4,[hom(z) for z in dat['M'].list()])
rhs={key:vector(K,[hom(z) for z in val]) for key,val in dat['rhs'].items()}
AB,TB,Av,Tv=[hom(dat[key]) for key in ['A_B','T_B','A_v','T_v']]
def db(z):return nk(z.derivative(B)+z.derivative(A)*AB+z.derivative(T)*TB)
def dv(z):return nk(z.derivative(v)+z.derivative(A)*Av+z.derivative(T)*Tv)
def dl(z):return nk(z.derivative(L))
E2=v*(v+L*L)/B
det=hom(dat['det']); print('loaded',round(time.time()-START,2),flush=True)
if '--reload' not in sys.argv:
 # Cramer replaces one column; no pivot other than the named complete det.
 X=[]
 for lab in ['constant','KC_over_U','KB_over_U']:
  vals=[]
  for col in range(4):
   W=matrix(M);W.set_column(col,-rhs[lab]);vals.append(nk(W.det()/det))
  X.append(vector(K,vals));print('solved',lab,[(len(x.numerator().dict()),len(x.denominator().dict())) for x in vals],flush=True)
 save({'P':P,'X':X,'det':det,'E_squared':E2,'I1':hom(dat['I1']),'I2':hom(dat['I2']),'A_B':AB,'T_B':TB,'A_v':Av,'T_v':Tv},str(OUT/'tr_pure_rotation_first_order_fields.sobj'))
else:
 X=load(str(OUT/'tr_pure_rotation_first_order_fields.sobj'))['X'];X=[vector(K,z) for z in X]
if '--solve-only' in sys.argv:sys.exit(0)
Bq=[nk(2*B/v+2*B/(v+L*L)+2*B*X[0][0]+2*B*L/(v+L*L)*X[0][2])]+[nk(2*B*z[0]+2*B*L/(v+L*L)*z[2]) for z in X[1:]]
Bz=[nk(2*B*z[1]+2*B*L/(v+L*L)*z[3]) for z in X]
save({'Bq':Bq,'Bz':Bz},str(OUT/'tr_pure_rotation_B_fields.sobj'))
def act(f,direction,pressureindex):
 bfield,lfield=(Bq,[z[2] for z in X]) if direction==0 else (Bz,[z[3] for z in X])
 return nk((2*dv(f) if direction==0 and pressureindex==0 else 0)+bfield[pressureindex]*db(f)+lfield[pressureindex]*dl(f))
answers={}
for label,iq,iz in [('logU',1,0),('L',3,2)]:
 out={}
 out['1']=nk(act(X[0][iq],0,0)-act(X[0][iz],1,0))
 print(label,'constant',len(out['1'].numerator().dict()),round(time.time()-START,2),flush=True)
 for cidx,clab in [(1,'C'),(2,'D')]:
  val=act(X[0][iq],0,cidx)+act(X[cidx][iq],0,0)-X[0][0]*X[cidx][iq]
  val-=act(X[0][iz],1,cidx)+act(X[cidx][iz],1,0)-X[0][1]*X[cidx][iz]
  out['E'+clab]=nk(val);print(label,'E'+clab,len(out['E'+clab].numerator().dict()),round(time.time()-START,2),flush=True)
 for aa,bb,lab in [(1,1,'CC'),(1,2,'CD'),(2,2,'DD')]:
  pairs=[(aa,bb)] if aa==bb else [(aa,bb),(bb,aa)]
  val=sum(act(X[b][iq],0,a)-X[a][0]*X[b][iq]-act(X[b][iz],1,a)+X[a][1]*X[b][iz] for a,b in pairs)
  out[lab]=nk(E2*val);print(label,lab,len(out[lab].numerator().dict()),round(time.time()-START,2),flush=True)
 answers[label]=out
 save({'P':P,'source':label,'coefficients':out,'det':det},str(OUT/('tr_pure_rotation_'+label+'_frobenius.sobj')))
 print('saved',label,flush=True)
save({'P':P,'answers':answers,'det':det},str(OUT/'tr_pure_rotation_frobenius.sobj'))
print('COMPLETE',round(time.time()-START,2),flush=True)
