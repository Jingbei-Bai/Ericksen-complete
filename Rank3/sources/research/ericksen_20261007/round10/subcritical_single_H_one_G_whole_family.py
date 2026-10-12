"""Direct full derivative source on the H=1 conic kappa family; moving Z."""
from sage.all import *
from pathlib import Path
import os,resource,time
p=Path(__file__).resolve().parent;st=time.monotonic()
g=load(str(p/'independent_partial_D_Newton_gauges_reader.sobj'));ev=load(str(p/'independent_delta_subcritical_even_full_source.sobj'))
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+500*1024**2,)*2)
K0=FractionField(PolynomialRing(QQ,names=('k','r','Z')));k,r,Z=K0.gens();q=1+r;m=1+k
K=FractionField(PolynomialRing(K0,'w'));w=K.gen();T=PolynomialRing(K,'x');E=T.quotient(T.gen()**2-2*m/k*(k+w*w),'x');x=E.gen()
old=g['ring'].base_ring();OK=ev['ring'].base_ring();oa,orr,H0,Z0=OK.gens();ph=old.hom([k/r,k*q/(r*m),1/r,q*Z/r],K0);po=OK.hom([k/r,r,1,Z],K0)
R=PolynomialRing(K0,names=ev['ring'].variable_names());vs=R.gens();cache={}
def cp(f):
 field=f.parent().base_ring();hh=ph if field==old else po
 out={}
 for ee,c in f.dict().items():
  ky=(str(field),c)
  if ky not in cache:cache[ky]=hh(c)
  out[tuple(map(int,ee))]=cache[ky]
 return R(out)
items=[]
for gi,qi in [(1,6)]:
 it=g['items'][gi];fv=list(map(cp,it['field_numerators']));bv=list(map(cp,it['base_numerators']))
 OQ=next(t['full_source'] for t in ev['items'] if int(t['Q_index'])==qi);Q=cp(OQ)
 G=sum((fv[j]*Q.derivative(vs[j]) for j in range(9)),R.zero())
 for j,param in enumerate((H0,Z0)):
  scale=orr if j==0 else orr/(1+orr)
  df=OQ.parent()({tuple(map(int,e)):scale*c.derivative(param) for e,c in OQ.dict().items() if c.derivative(param)})
  G+=bv[j]*cp(df)
 out={};images=(x,E.one(),-x,E(w),E.zero(),E.zero(),E(-2*w),E(-2),E.one())
 for ee,c in G.dict().items():
  e=tuple(map(int,ee))
  if sum(e)!=6 or e[4] or e[5]:continue
  power=13-2*sum(e)-sum(e[6:])+e[6]+e[7]
  val=E(c)*prod(images[j]**e[j] for j in range(9) if e[j])
  out[power]=out.get(power,E.zero())+val
 out={j:v for j,v in out.items() if v}
 print('G WHOLE FAMILY',gi,qi,'REMAINDER',out,'SECONDS',time.monotonic()-st,flush=True)
 items.append({'gauge':gi,'Q_index':qi,'source':G,'whole_family_remainder':out})
save({'ring':E,'source_ring':R,'items':items},str(p/'subcritical_single_H_one_G_whole_family.sobj'),compress=False)
print('DONE',time.monotonic()-st,flush=True)
