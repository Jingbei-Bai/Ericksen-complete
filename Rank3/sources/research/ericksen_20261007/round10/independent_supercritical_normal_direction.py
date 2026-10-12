"""Complete direction consequences of saved true normal minors off lambda=0."""
from sage.all import *
from pathlib import Path
import json,time
p=Path(__file__).resolve().parent;st=time.monotonic()
d=load(str(p/'independent_supercritical_normal_minors.sobj'));K=d['field'];a,r,H,Z,s,z=K.gens();BR=K.ring();lam=s*z+s+z
P=PolynomialRing(QQ,names=('ds','dz'));ss,zz=P.gens();ph=BR.hom([1,1,1,1,ss,zz],P)
gens=[];items=[];checks=[]
def ck(n,v):assert v,n;checks.append(n)
for it in d['items']:
 f=it['determinant']
 if not f:continue
 f=f/lam;fac=f.numerator().factor();core=BR.one();unit=K(fac.unit())
 for g,m in fac:
  hasdir=(g.degree(BR.gen(4))>0 or g.degree(BR.gen(5))>0)
  if hasdir and g not in (BR.gen(4),BR.gen(5)):core*=g**m
  else:unit*=K(g)**m
 ck('factor replay '+str(it['source_rows']),K(core)*unit/f.denominator()==f)
 # Every removed factor is exactly a protected support/base/material unit.
 allowed=[s,z,r,r+1,H,Z,H+1,Z-H,Z-1,H*Z-H+Z]
 def const(v):return v.numerator().total_degree()==0 and v.denominator().total_degree()==0
 for g,m in unit.factor():
  ck('removed unit '+str(g),any(const(K(g)/v) for v in allowed) or const(K(g)))
 groups={}
 for ee,co in core.dict().items():
  e=tuple(map(int,ee));ck('material-free residual coefficient',e[0]==e[1]==0)
  key=e[2:4];groups[key]=groups.get(key,P.zero())+QQ(co)*ss**e[4]*zz**e[5]
 rows=tuple(groups.values());gens.extend(rows)
 items.append({'source_rows':it['source_rows'],'residual_core':core,'removed_unit':unit/f.denominator(),'base_coefficients':groups})
I=P.ideal(gens);gb=I.groebner_basis();print('GB',gb,'seconds',time.monotonic()-st,flush=True)
targets=[(ss-1)**2*(2*ss*ss-13*ss+2),(zz-1)**2*(2*zz*zz-13*zz+2),(ss-1)*(zz-1)*(2*zz*zz-13*zz+2),(ss-1)*(zz-1)*(ss-zz)]
lifts=[]
for f in targets:
 if f in I:
  co=f.lift(I);ck('original coefficient lift',sum((g*c for g,c in zip(gens,co)),P.zero())==f);lifts.append((f,tuple(co)))
 else:print('TARGET NOT MEMBER',f,flush=True)
save({'ring':P,'sources':gens,'items':items,'groebner_basis':gb,'target_lifts':lifts,'checks':checks,'scope':__doc__},str(p/'independent_supercritical_normal_direction.sobj'),compress=False)
print('PASS',len(checks),'LIFTS',len(lifts),'seconds',time.monotonic()-st,flush=True)
