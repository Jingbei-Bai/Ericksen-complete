"""COMPLETE single RR chart with p12=1 and retained incidence, no p13 inverse.

The normal forms equal the original normalized sources on incidence=0
and kappa nonzero; they are regular polynomials on the ambient chart.
"""
from sage.all import *
from pathlib import Path
import time,sys,resource,os
p=Path(__file__).resolve().parent;st=time.monotonic();hfix=int(sys.argv[1]) if len(sys.argv)>1 else 2
qf=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'));R=qf['ring'];OK=R.base_ring()
hs=load(str(p/'independent_delta_subcritical_H_initials.sobj'));full=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'));old=load(str(p/'independent_delta_H_global_L_module.sobj'))['original_field']
ev=load(str(p/'independent_delta_subcritical_even_full_source.sobj'));pairdata=load(str(p/'independent_ratio_p_boundary_global_pair_syzygies.sobj'))
K=PolynomialRing(QQ,names=('k','r')).fraction_field();k,r=K.gens();H=K(hfix);Z=K(-1);m=1+k;q=1+r
ph=OK.ring().hom([k/r,r,H,Z],K);cache={}
P=PolynomialRing(K,names=('y','z','w','e0','e1','e2','delta','t','v','kap'));y,z,w,e0,e1,e2,de,t,v,kap=P.gens()
alpha=(Z-H)/(Z*(H+1));beta=H*(Z-1)/(Z*(H+1));imgs=(P.one(),y,z,w,(t*z-kap*e2)/beta,(v-kap*e0)/alpha,v,t,P.one());es=(e0,e1,e2)
inc0=v*y-alpha/beta*t*w;inc=inc0-kap*e1
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+550*1024**2,)*2)
def cv(c):
 c=OK(c)
 if c not in cache:cache[c]=ph(c.numerator())/ph(c.denominator())
 return cache[c]
def cp(f):return R({tuple(map(int,ee))+(0,)*(9-len(ee)):OK(c) for ee,c in f.dict().items()})
def mp(f):return sum((cv(c)*prod(im**int(n) for im,n in zip(imgs,ee) if n) for ee,c in cp(f).dict().items()),P.zero())
def norm(f,offset,discard):
 ans=P.zero()
 for ee,c in cp(f).dict().items():
  pd=sum(ee[:6]);rd=sum(ee[6:])
  if (pd,rd) in discard:continue
  ke,dd=offset[0]-2*pd-3*rd,offset[1]-pd-rd;assert ke>=0 and dd>=0,(pd,rd,ke,dd)
  ans+=cv(c)*kap**ke*de**dd*prod(im**int(n) for im,n in zip(imgs,ee) if n)
 return ans
DC=matrix(R,qf['D'])*matrix(R,qf['C']);sources=[];labels=[]
for j,rf in enumerate(qf['full_R_sources']):
 sources.append(norm(rf,(8,4),{(3,1)})+sum((mp(DC[j,i])*es[i] for i in range(3)),P.zero()));labels.append(('R',j))
sources.append(norm(qf['full_Q0'],(7,3),{(1,2)})+sum((mp(qf['Q0_L_cofactors'][i])*es[i] for i in range(3)),P.zero()));labels.append(('Q0',0))
print('R/Q built',time.monotonic()-st,flush=True)
oph=old.ring().hom([k/r,k*q/(r*m),H/r,q*Z/r],K);oldcache={};pool=full['polynomial_pool']
def oldcv(c):
 c=old(c)
 if c not in oldcache:oldcache[c]=oph(c.numerator())/oph(c.denominator())
 return oldcache[c]
def om(dd):return sum((oldcv(c)*prod(im**int(n) for im,n in zip(imgs,ee) if n) for ee,c in dd.items()),P.zero())
for j,it in enumerate(hs['items']):
 coeffs={}
 for sign,lr,rr in full['records'][j]['full_signed_product_DAG']:
  for le,lc in pool[lr].items():
   for re,rc in pool[rr].items():
    ee=tuple(int(a)+int(b) for a,b in zip(le,re))
    if (sum(ee[:6]),sum(ee[6:]))==(3,2):continue
    coeffs[ee]=coeffs.get(ee,K.zero())+sign*oldcv(lc)*oldcv(rc)
 ans=sum((om(it['L_cofactors_original'][i])*es[i] for i in range(3)),P.zero())
 for ee,c in coeffs.items():
  if not c:continue
  pd=sum(ee[:6]);rd=sum(ee[6:]);ke,dd=11-2*pd-3*rd,5-pd-rd;assert ke>=0 and dd>=0,(j,pd,rd,ke,dd)
  ans+=c*kap**ke*de**dd*prod(im**int(n) for im,n in zip(imgs,ee) if n)
 sources.append(ans);labels.append(('H',j));print('H',j,len(ans.dict()),'sec',time.monotonic()-st,flush=True)
for it in ev['items']:
 src=cp(it['full_source']);ell=matrix(R,qf['L'])*vector(R,R.gens()[6:]);rem=src-OK(it['square_unit'])*ell[int(it['square_component'])]**2-sum((cp(it['A_p_quadratic'][i])*ell[i] for i in range(3)),R.zero())
 ans=norm(rem,(8,4),set())+cv(it['square_unit'])*es[int(it['square_component'])]**2+sum((mp(it['A_p_quadratic'][i])*es[i] for i in range(3)),P.zero())
 sources.append(ans);labels.append(('even',it['Q_index']))
save({'ring':P,'sources':sources,'labels':labels,'chart_p':imgs[:6],'chart_Rbar':imgs[6:],'incidence':inc,'H':hfix,'scope':__doc__},str(p/('subcritical_single_p12_full_chart_H'+str(hfix)+'.partial.sobj')),compress=False)
pairregs=[]
def exactdiv(f,g):
 ee,cc=next(iter(g.lt().dict().items()));quo=P.zero();rem=f
 while rem:
  le,lc=next(iter(rem.lt().dict().items()));diff=tuple(int(a)-int(b) for a,b in zip(le,ee));assert min(diff)>=0,'not an incidence multiple'
  pp=P({diff:lc/cc});quo+=pp;rem-=pp*g
 return quo
def kappadiv(f,n):
 for unused in range(n):
  f0=P.hom(list(P.gens()[:9])+[P.zero()],P)(f)
  if f0:
   co=exactdiv(f0,inc0)
   f-=co*inc
  assert all(ee[9]>=1 for ee in f.dict())
  f=P({tuple(int(nn)-(1 if j==9 else 0) for j,nn in enumerate(ee)):c for ee,c in f.dict().items()})
 return f
for it in pairdata['items']:
 layers={}
 for ee,c in it['full'].dict().items():
  ky=(sum(ee),sum(ee[6:]));layers.setdefault(ky,{})[tuple(map(int,ee))]=c
 ans=P.zero()
 for (deg,rd),co in layers.items():
  pw=9-2*deg-rd;dp=4-deg;assert dp>=0;ff=om(co)
  ff=kappadiv(ff,-pw) if pw<0 else ff*kap**pw
  ans+=de**dp*ff
 sources.append(ans);labels.append(('pair',it['pair']));print('PAIR',it['pair'],len(ans.dict()),flush=True)
sources.append(inc);labels.append(('incidence',0))
stem='subcritical_single_p12_full_chart_H'+str(hfix)
save({'ring':P,'sources':sources,'labels':labels,'chart_p':imgs[:6],'chart_Rbar':imgs[6:],'incidence':inc,'H':hfix,'scope':__doc__},str(p/(stem+'.sobj')),compress=False)
print('DONE',len(sources),'sec',time.monotonic()-st,flush=True)
