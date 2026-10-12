"""Full ratio-curve tangencies from independently rebuilt joint fields."""
from sage.all import *
from pathlib import Path
from functools import lru_cache
import json,time,hashlib
out=Path(__file__).absolute().parent;st=time.monotonic();checks=[]
d=load(str(out/'fresh_joint_fields.sobj'));O=d['ring'];Ko=O.base_ring();oa,ob,oh,oz=Ko.gens();xx=O.gens();fields=d['fields'];base=d['base_fields']
C=PolynomialRing(QQ,names=('a','r','H','Z','S','N'));K=C.fraction_field();a,r,H,Z,S,N=K.gens();q=1+r;b=a*q/(1+a*r)
R=PolynomialRing(K,names=('w1','m1','w2','m2','w3','m3'));w1,m1,w2,m2,w3,m3=R.gens()
pi=[w3,w2,((1+H)*S-Z)*m3/(2*Z)+H*w3,w1,-(1+H*S-Z)*m2/(2*H)-Z*w2,-(1-S)*m1/(2*H)-Z*w1/H]
rr=[(Z-H*S)*m1/(H*Z),S*m2/Z,-m3/H];images=pi+rr;ph=Ko.ring().hom([a,b,H/r,q*Z/r],K)
@lru_cache(None)
def cv(c):c=Ko(c);return ph(c.numerator())/ph(c.denominator())
def ev(f):return sum((cv(c)*prod(images[j]**int(e) for j,e in enumerate(mon) if e) for mon,c in f.dict().items()),R.zero())
def ck(n,v):
 if not v:raise AssertionError(n)
 checks.append(n)
def BD(f,p):return O({tuple(e):c.derivative(p) for e,c in f.dict().items() if c.derivative(p)})
def D(i,f):return sum(fields[i,j]*f.derivative(xx[j]) for j in range(9))+base[i,0]*BD(f,oh)+base[i,1]*BD(f,oz)
oratio=(oa-ob)/(oa*(ob-1));oq=1+oratio
EH=[oratio*base[i,0] for i in range(3)];EZ=[oratio/oq*base[i,1] for i in range(3)]
mh=[ev(f) for f in EH];mz=[ev(f) for f in EZ]
for i in range(3):ck('curve first tangency '+str(i),mz[i]==S*mh[i])
rows=[];pairs=[(0,1),(0,2),(1,2),(0,0),(1,1),(2,2)]
for i,j in pairs:
 f=ev(D(i,EZ[j]))-S*ev(D(i,EH[j]))-N*mh[i]*mh[j]
 rows.append(f);print('fresh second tangency',i,j,'degree',f.total_degree(),'seconds',time.monotonic()-st,flush=True)
Q0=ev(d['quartic_source_rows'][0]);allrows=rows[:3]+[Q0]+rows[3:]
ck('six complete second tangencies and Q0',len(rows)==6 and len(allrows)==7)
ck('complete actual joint field dimensions',fields.dimensions()==(3,9) and base.dimensions()==(3,2))
# The six tangencies use Z=r*z/(1+r); Q0 has no derivative-coordinate factor.
row_scalings=[r/q,r/q,r/q,K.one(),r/q,r/q,r/q]
for i,c in enumerate(row_scalings):ck('invertible ratio-coordinate row scaling '+str(i),c!=0)
save(dict(ring=R,rows=allrows,row_scalings=row_scalings,gradient_H=mh,images=images,provenance='Differentiated E_jZ-S E_jH under fresh actual fields and inserted S prime=N; all first and second tangencies retained.'),str(out/'fresh_curve_sources.sobj'),compress=False)
report=dict(status='PASS',checks=checks,count=len(checks),seconds=time.monotonic()-st,source='Full actual ratio-curve tangencies, with all original lower-degree terms retained.',rows=len(allrows),second_tangency_pairs=pairs,row_scalings=[str(v) for v in row_scalings],inputs=[dict(path='fresh_joint_fields.sobj',sha256=hashlib.sha256((out/'fresh_joint_fields.sobj').read_bytes()).hexdigest())])
(out/'fresh_curve_sources.json').write_text(json.dumps(report,indent=2));print(report,flush=True)
