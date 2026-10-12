"""All nine complete own-active derivative sources with all positive powers."""
from sage.all import *
from pathlib import Path
import time
p=Path(__file__).resolve().parent;st=time.monotonic();ch=load(str(p/'independent_subcritical_single_RR_full_regular_chart_H2.sobj'));P=ch['ring'];K=P.base_ring();k,r=K.gens();de=P.gen(6);kap=P.gen(8);g=load(str(p/'independent_partial_D_Newton_gauges_reader.sobj'));ev=load(str(p/'independent_delta_subcritical_even_full_source.sobj'));old=g['ring'].base_ring();OK=ev['ring'].base_ring();a0,r0,H0,Z0=OK.gens();m=1+k;q=1+r;ph=old.ring().hom([k/r,k*q/(r*m),2/r,-q/r],K);po=OK.ring().hom([k/r,r,2,-1],K);cache={};R=PolynomialRing(K,names=ev['ring'].variable_names());vs=R.gens();hm=R.hom(list(ch['chart_p'])+list(ch['chart_Rbar']),P)
def cp(f):
 field=f.parent().base_ring();hh=ph if field==old else po;out={}
 for ee,c in f.dict().items():
  key=(str(field),c)
  if key not in cache:cache[key]=hh(c.numerator())/hh(c.denominator())
  out[tuple(map(int,ee))]=cache[key]
 return R(out)
items=[]
for gi in (1,2,4):
 it=g['items'][gi];fv=list(map(cp,it['field_numerators']));bv=list(map(cp,it['base_numerators']))
 for qi in (1,5,6):
  OQ=next(it['full_source'] for it in ev['items'] if int(it['Q_index'])==qi);Q=cp(OQ);parts=[fv[j]*Q.derivative(vs[j]) for j in range(9)]
  for j,param in enumerate((H0,Z0)):
   scale=r0 if j==0 else r0/(1+r0);dd=OQ.parent()({tuple(map(int,ee)):scale*c.derivative(param) for ee,c in OQ.dict().items() if c.derivative(param)});parts.append(bv[j]*cp(dd))
  G=sum(parts,R.zero());layers={}
  for ee,c in G.dict().items():
   deg=sum(ee);rd=sum(ee[6:]);layers.setdefault((deg,rd),{})[tuple(map(int,ee))]=c
  ans=P.zero();regular=[]
  for (deg,rd),co in layers.items():
   power=13-2*deg-rd;dp=6-deg;assert dp>=0;f=hm(R(co))
   if power<0:
    assert all(int(ee[8])>=-power for ee in f.dict()),(gi,qi,deg,rd)
    f=P({tuple(int(n)+(power if j==8 else 0) for j,n in enumerate(ee)):c for ee,c in f.dict().items()});regular.append((deg,rd,power));power=0
   ans+=kap**power*de**dp*f
  print('FULL G',gi,qi,'terms',len(ans.dict()),'sec',time.monotonic()-st,flush=True);items.append({'gauge':gi,'direction':it['direction'],'Q_index':qi,'source':ans,'regularity':regular})
  save({'ring':P,'items':items,'scope':__doc__},str(p/'independent_subcritical_single_RR_full_complete_E_chart_H2.sobj'),compress=False)
print('PASS full nine G sources',len(items),time.monotonic()-st,flush=True)
