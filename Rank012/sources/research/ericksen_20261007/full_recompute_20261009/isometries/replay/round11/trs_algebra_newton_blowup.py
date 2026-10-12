"""Exact local Newton blowup at the two common complex B-infinity directions.
All parameters retained. A formal certificate computation, not an actual real limit.
"""
from sage.all import *
from pathlib import Path
import time,resource,signal,json,sys
OUT=Path(__file__).resolve().parent;START=time.time()
sd=load(str(OUT/'trs_algebra_regular_first_order_solution.sobj'))
vm=int(Path('/proc/self/statm').read_text().split()[0])*4096
resource.setrlimit(resource.RLIMIT_AS,(vm+550*1024**2,vm+550*1024**2));signal.alarm(240)
NF=QuadraticField(-1,'ii');ii=NF.gen()
PR=PolynomialRing(NF,names=('v','p','j','I1','I2','KC','KB','h'))
K=PR.fraction_field();v,p,j,I1,I2,KC,KB,h=map(K,PR.gens())
PREC=int(sys.argv[1]) if len(sys.argv)>1 else 10
POWER=int(sys.argv[2]) if len(sys.argv)>2 else 2
SR=LaurentSeriesRing(K,'e',default_prec=PREC);e=SR.gen()
def nk(x):
 x=K(x);n=x.numerator();d=x.denominator();g=n.gcd(d);n=n//g;d=d//g;lc=d.leading_coefficient()
 return K(n/lc)/K(d/lc)
def cl(x):
 x=SR(x);prec=min(PREC,x.precision_absolute())
 if not x:return SR(0).add_bigoh(prec)
 return sum(SR(nk(x[i]))*e**i for i in range(int(x.valuation()),int(prec)) if x[i]).add_bigoh(prec)
def cp(x,var):
 x=SR(x);prec=min(PREC,x.precision_absolute())
 return cl(sum(SR(nk(x[i].derivative(var)))*e**i for i in range(int(x.valuation()),int(prec)) if x[i]))
def sqrtone(x):
 z=cl(x-1);return cl(sum(binomial(QQ(1)/2,k)*z**k for k in range(PREC)))
q=v-p*p;c=-j*j/v;P=v**3-I1*v*v+I2*v-j
l0=ii*v/p;D0=nk(v+l0*l0);l=cl(l0+h*e)
areg=cl(j/v+P*e**(2*POWER)/j**2)
breg=c
f1=-v*I2/j+2-P/j+v**3/j
f2=P*(v*I2-3*j)/j**2
f3=P**2/j**2
treg=cl(c*sqrtone(1+f1*e**(2*POWER)/c+f2*e**(4*POWER)/c**2+f3*e**(6*POWER)/c**3))
# U=y/e^2; the common leading radical branch has y0=-jp/[v(v-p^2)].
y0=-j*p/(v*q)
y=cl(y0*sqrtone(D0/(v+l*l)))
U=cl(y/e**POWER);Ui=cl(e**POWER/y)
assert nk(y0*y0-c/(q*D0))==0
pcache={}
def ppow(which,n):
 key=(which,n)
 if key not in pcache:
  xx={'a':areg,'t':treg,'l':l}[which]
  pcache[key]=cl(xx**n)
 return pcache[key]
evalcache={}
def pseries(poly):
 # Return polynomial as e^m times a power series, preserving relative precision.
 key=poly
 if key in evalcache:return evalcache[key]
 groups={};m=min(-POWER*int(ex[4])-2*POWER*int(ex[5]) for ex in poly.dict()) if poly else 0
 for ex,coef in poly.dict().items():
  ev,ep,ej,ea,et,eb,el=map(int,ex)
  key=(ea,et,eb,el)
  groups[key]=groups.get(key,K(0))+K(coef)*v**ev*p**ep*j**ej
 out=SR(0)
 for (ea,et,eb,el),coef in groups.items():
  shift=-POWER*et-2*POWER*eb-m
  if shift>=PREC:continue
  out+=coef*c**eb*e**shift*ppow('a',ea)*ppow('t',et)*ppow('l',el)
 out=cl(out)
 evalcache[poly]=(m,out)
 return m,out
def evrat(x):
 nn,dd=x.numerator(),x.denominator()
 if not nn:return SR(0).add_bigoh(PREC)
 mn,sn=pseries(nn);md,ss=pseries(dd)
 return cl(e**(mn-md)*sn/ss)
cols={}
for key,vec in sd['solution_coefficients'].items():
 cols[key]=[evrat(x) for x in vec]
 print('column',key,'seconds',round(time.time()-START,2),flush=True)
save({'columns':cols,'base':K,'series':SR,'U':U,'Ui':Ui,'L':l,'L0':l0,'B':c/e**(2*POWER),'precision':PREC},str(OUT/('trs_algebra_newton_blowup_columns_'+str(POWER)+'.sobj')))
xx=[]
for row in range(4):
 x=cl(cols['constant'][row]+KC*cols['KC'][row]+U*cols['U'][row]+Ui*(KC*cols['KC_over_U'][row]+KB*cols['KB_over_U'][row]+cols['one_over_U'][row]))
 xx.append(x);print('row',row,'valuation',x.valuation(),'precision',x.precision_absolute(),flush=True)
xi,eta,lq,lz=xx
D=cl(v+l*l)
bqr=cl(2/q+2/D+2*xi+2*l*lq/D)
bzr=cl(2*eta+2*l*lz/D)
eq=cl(-e/(POWER*v)-e*bqr/(2*POWER));ez=cl(-e*bzr/(2*POWER))
hq=cl((lq-2*ii/p-h*eq)/e);hz=cl((lz-h*ez)/e)
def dx(x):return cl(2*cp(x,v)+eq*x.derivative()+hq*cp(x,h))
def dz(x):return cl(ez*x.derivative()+hz*cp(x,h))
C0=cl(dx(eta)-dz(xi));print('C0 done',round(time.time()-START,2),flush=True)
C1=cl(dx(lz)-dz(lq));print('C1 done',round(time.time()-START,2),flush=True)
CB=cl(C0+l*C1/D)
res={'C0':C0,'C1':C1,'CB':CB}
for name,s in res.items():
 print(name,'valuation',s.valuation(),'precision',s.precision_absolute(),flush=True)
 for n in range(int(s.valuation()),min(int(s.precision_absolute()),int(s.valuation())+4)):
  print(name,n,nk(s[n]).factor(),flush=True)
save({'base':K,'series':SR,'precision':PREC,'solution':xx,'vector_fields':(eq,ez,hq,hz),'residuals':res},str(OUT/('trs_algebra_newton_blowup_'+str(POWER)+'.sobj')))
(OUT/('trs_algebra_newton_blowup_'+str(POWER)+'.json')).write_text(json.dumps({'status':'EXACT_LOCAL_SERIES_ONLY','precision':PREC,'seconds':time.time()-START,'residuals':{name:{'valuation':int(s.valuation()),'precision':int(s.precision_absolute()),'coefficients':{str(n):str(nk(s[n]).factor()) for n in range(int(s.valuation()),min(int(s.precision_absolute()),int(s.valuation())+4))}} for name,s in res.items()}},indent=2))
signal.alarm(0)

