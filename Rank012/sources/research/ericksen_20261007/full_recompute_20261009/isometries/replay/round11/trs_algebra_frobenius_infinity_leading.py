"""Exact highest B powers of regular Frobenius sources, all constants free.

Finite polynomial coefficient extraction, not numerical asymptotics.
The retained leading symbols are not asserted to give a nonzero resultant.
"""
from sage.all import *
from pathlib import Path
import resource,time,signal,json
OUT=Path(__file__).resolve().parent;START=time.time()
d=load(str(OUT/'trs_algebra_regular_frobenius.sobj'));K=d['source']['field']
vm=int(Path('/proc/self/statm').read_text().split()[0])*4096
resource.setrlimit(resource.RLIMIT_AS,(vm+550*1024**2,vm+550*1024**2));signal.alarm(240)
R=PolynomialRing(QQ,names=('v','p','j','I2','P','L'));F=R.fraction_field();v,p,j,I2,P,L=map(F,R.gens())
DEPTH=10
def nf(x):
 x=F(x);n=x.numerator();de=x.denominator();g=n.gcd(de);n=n//g;de=de//g;lc=de.leading_coefficient()
 return F(n/lc)/F(de/lc)
def conv(x,y):
 out=[F(0)]*min(DEPTH+1,len(x)+len(y)-1)
 for a,xx in enumerate(x):
  if not xx:continue
  for b,yy in enumerate(y):
   if a+b<len(out):out[a+b]+=xx*yy
 return [nf(z) for z in out]
# T^2 = (-j^2/v)*B*(1+f1/B+f2/B^2+f3/B^3).
fn=[F(1),nf(-v*I2/j+2-P/j+v**3/j),nf(P*(v*I2-3*j)/j**2),nf(P**2/j**2)]
fpowers=[[F(1)]]
def fpow(n):
 while len(fpowers)<=n:fpowers.append(conv(fpowers[-1],fn))
 return fpowers[n]
acache={};prod_cache={}
def coeffs(a,n):
 key=(a,n)
 if key not in prod_cache:
  if a not in acache:acache[a]=[nf(binomial(a,h)*j**(a-h)*(-P)**h/v**a) for h in range(min(a,DEPTH)+1)]
  prod_cache[key]=[nf(z*(-j*j/v)**n) for z in conv(acache[a],fpow(n))]
 return prod_cache[key]
def poly_lead(poly):
 groups={};tops=[-10**9,-10**9]
 for e,c in poly.dict().items():
  ev,ep,ej,ea,et,eb,el=map(int,e);key=(ea,et,eb)
  mon=R.monomial(ev,ep,ej,0,0,el)*QQ(c)
  groups[key]=groups.get(key,R(0))+mon
  tops[et%2]=max(tops[et%2],eb+et//2)
 totals=[{},{}]
 for (ea,et,eb),coef in groups.items():
  parity=et%2;n=et//2;top=eb+n
  upto=min(DEPTH,top-tops[parity]+DEPTH)
  if upto<0:continue
  vals=coeffs(ea,n)
  for drop,z in enumerate(vals[:upto+1]):
   deg=top-drop;totals[parity][deg]=totals[parity].get(deg,F(0))+F(coef)*z
 best=[]
 for parity in range(2):
  for deg in sorted(totals[parity],reverse=True):
   c=nf(totals[parity][deg])
   if c:
    best.append((2*deg+parity,c,parity));break
 if not best:raise RuntimeError('No leading coefficient within exact depth '+str(DEPTH))
 return max(best,key=lambda x:x[0])
def rational_lead(x):
 nn=poly_lead(x.numerator());dd=poly_lead(x.denominator());power=nn[0]-dd[0];parity=nn[2]-dd[2]
 c=nf(nn[1]/dd[1])
 if parity<0:c=nf(c/(-j*j/v));parity+=2
 return power,c,parity
items=[]
for row,C in enumerate(d['frobenius']):
 for odd,pol in enumerate(C):
  for mon,x in sorted(pol.items()):
   if not x:continue
   power,c,parity=rational_lead(x);power+=odd
   item={'row':row,'U_parity':odd,'pressure_monomial':tuple(mon),'twice_B_power':power,'tau_parity':parity,'coefficient':c}
   items.append(item)
   print('leading',row,odd,mon,'2degree',power,'tau',parity,'terms',len(c.numerator().dict()),round(time.time()-START,2),flush=True)
   save({'base':F,'tau_squared':-j*j/v,'u_squared':1/((v-p*p)*(v+L*L)),'items':items},str(OUT/'trs_algebra_frobenius_infinity_partial.sobj'))
save({'base':F,'tau_squared':-j*j/v,'u_squared':1/((v-p*p)*(v+L*L)),'items':items},str(OUT/'trs_algebra_frobenius_infinity_leading.sobj'))
summary={'status':'EXACT_COEFFICIENT_LEADING_ONLY','seconds':time.time()-START,'items':[{key:(str(value.factor()) if key=='coefficient' else value) for key,value in item.items()} for item in items],
 'scope':'The algebraic roots tau^2=-j^2/v, u^2=1/[(v-p^2)(v+L^2)] remain. Terms with matching order must be combined before any resultant claim.'}
(OUT/'trs_algebra_frobenius_infinity_leading.json').write_text(json.dumps(summary,indent=2));signal.alarm(0);print(json.dumps({'seconds':time.time()-START,'items':len(items)}),flush=True)
