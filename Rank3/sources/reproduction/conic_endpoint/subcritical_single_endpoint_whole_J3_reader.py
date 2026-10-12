"""Whole-base, ALL t/kappa J3 identities for actual range and B sources.

J3=(w,e0,e1,delta). No claim is made for the other H sources.
These six exact identities suffice for B's w factor after eliminating
the three genuine range rows R2,R3,Q0, before solving any other row.
"""
from sage.all import *
from pathlib import Path
import os,resource,time,json,sys
out=Path(__file__).absolute().parent;p=out/'source_inputs';st=time.monotonic()
qf=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'))
full=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'))
old=full['coefficient_field'];OK=qf['ring'].base_ring()
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+500*1024**2,)*2)
K=PolynomialRing(QQ,names=('k','r','H','Z')).fraction_field();k,r,H,Z=K.gens();q=1+r;m=1+k
P=PolynomialRing(K,names=('x','z','e','t','kap'));x,z,e,t,kap=P.gens();beta=H*(Z-1)/(Z*(H+1))
images=(x,P.one(),z,P.zero(),(t*z-kap*e)/beta,P.zero(),P.zero(),t,P.one())
ph=old.ring().hom([k/r,k*q/(r*m),H/r,q*Z/r],K);po=OK.ring().hom([k/r,r,H,Z],K);cache={};checks=[]
def cv(c,F):
 ky=(str(F),c)
 if ky not in cache:
  hh=ph if F==old else po;cache[ky]=hh(c.numerator())/hh(c.denominator())
 return cache[ky]
def fmap(dd,F):
 ans={}
 for ee,c in dd.items():
  if any(n and not image for n,image in zip(ee,images)):continue
  ky=(sum(ee),sum(ee[6:]));v=cv(c,F)*prod(image**int(n) for image,n in zip(images,ee) if n)
  ans[ky]=ans.get(ky,P.zero())+v
 return ans
def normalized(mapped,off,dg):
 return sum((kap**(off-2*dd-rr+10)*v for (dd,rr),v in mapped.items() if dd==dg),P.zero())
purey='--purey' in sys.argv
for i in ((4,5) if purey else (2,3)):
 f=qf['full_R_sources'][i];ans=normalized(fmap(f.dict(),f.parent().base_ring()),8,4)
 assert not ans,('R',i);checks.append('complete R'+str(i)+' whole J3')
f=qf['full_Q0'];ans=normalized(fmap(f.dict(),f.parent().base_ring()),7,3)
assert not ans,'Q0';checks.append('complete Q0 whole J3')
pool=full['polynomial_pool'];pc={}
def pm(j):
 if j not in pc:pc[j]=fmap(pool[j],old)
 return pc[j]
for j in ((16,) if purey else (6,14,16)):
 ans=P.zero()
 for sign,li,ri in full['records'][j]['full_signed_product_DAG']:
  for (dg,rd),a in pm(li).items():
   for (dg2,rd2),b in pm(ri).items():
    if dg+dg2!=5:continue
    pw=21-2*(dg+dg2)-rd-rd2;assert pw>=0
    ans+=sign*kap**pw*a*b
 assert not ans,('H',j)
 checks.append('complete H'+str(j)+' whole J3');print('H',j,'PASS seconds',time.monotonic()-st,flush=True)
stem='subcritical_single_pure_y_whole_J3_reader' if purey else 'subcritical_single_endpoint_whole_J3_reader'
save({'ring':P,'checks':checks,'scope':__doc__},str(out/(stem+'.sobj')),compress=False)
(out/(stem+'.json')).write_text(json.dumps({'status':'PASS','checks':len(checks),'identities':checks,'seconds':time.monotonic()-st},indent=2)+'\n')
print('PASS',len(checks),'seconds',time.monotonic()-st,flush=True)
