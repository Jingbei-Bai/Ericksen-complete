"""Whole-base FULL pair sources on retained-incidence p12 chart.

For every apparent kappa divisor, subtract a polynomial multiple of
the ACTUAL incidence and divide only the exact kappa factor. All
resulting polynomials are necessary on incidence=0,kappa!=0.
"""
from sage.all import *
from pathlib import Path
import time,os,resource
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_ratio_p_boundary_global_pair_syzygies.sobj'))
K=PolynomialRing(QQ,names=('k','r','H','Z')).fraction_field();k,r,H,Z=K.gens();m=1+k;q=1+r
P=PolynomialRing(K,names=('y','z','w','e0','e1','e2','delta','t','v','kap'));y,z,w,e0,e1,e2,de,t,v,kap=P.gens()
alpha=(Z-H)/(Z*(H+1));beta=H*(Z-1)/(Z*(H+1));imgs=(P.one(),y,z,w,(t*z-kap*e2)/beta,(v-kap*e0)/alpha,v,t,P.one());inc0=v*y-alpha/beta*t*w;inc=inc0-kap*e1
old=d['ring'].base_ring();ph=old.ring().hom([k/r,k*q/(r*m),H/r,q*Z/r],K);cache={};items=[]
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+500*1024**2,)*2)
def cv(c):
 if c not in cache:cache[c]=ph(c.numerator())/ph(c.denominator())
 return cache[c]
def om(dd):return sum((cv(c)*prod(im**int(n) for im,n in zip(imgs,ee) if n) for ee,c in dd.items()),P.zero())
def exactdiv(f,g):
 ee,cc=next(iter(g.lt().dict().items()));quo=P.zero();rem=f
 while rem:
  le,lc=next(iter(rem.lt().dict().items()));diff=tuple(int(a)-int(b) for a,b in zip(le,ee));assert min(diff)>=0,'not an incidence multiple'
  pp=P({diff:lc/cc});quo+=pp;rem-=pp*g
 return quo
for it in d['items']:
 layers={};cert=[];ans=P.zero()
 for ee,c in it['full'].dict().items():layers.setdefault((sum(ee),sum(ee[6:])),{})[tuple(map(int,ee))]=c
 for (deg,rd),co in layers.items():
  pw=9-2*deg-rd;dp=4-deg;assert dp>=0;ff=om(co)
  if pw<0:
   for j in range(-pw):
    f0=P.hom(list(P.gens()[:9])+[P.zero()],P)(ff);quo=P.zero()
    if f0:
     quo=exactdiv(f0,inc0)
     ff-=quo*inc
    assert all(ee[9]>=1 for ee in ff.dict())
    ff=P({tuple(int(nn)-(1 if a==9 else 0) for a,nn in enumerate(ee)):c for ee,c in ff.dict().items()});cert.append((deg,rd,j,quo))
  else:ff*=kap**pw
  ans+=de**dp*ff
 items.append({'pair':it['pair'],'source':ans,'division_certificates':cert});print('PAIR',it['pair'],'terms',len(ans.dict()),'sec',time.monotonic()-st,flush=True)
save({'ring':P,'items':items,'incidence':inc,'chart_p':imgs[:6],'chart_Rbar':imgs[6:],'scope':__doc__},str(p/'subcritical_single_p12_pair_generic.sobj'),compress=False)
bp=P.hom([P.zero(),z,P.zero(),e0,P.zero(),e2,P.zero(),P.zero(),P.zero(),P.zero()],P)
for it in items:print('XZ BOUNDARY PAIR',it['pair'],bp(it['source']),flush=True)
print('DONE',time.monotonic()-st,flush=True)
