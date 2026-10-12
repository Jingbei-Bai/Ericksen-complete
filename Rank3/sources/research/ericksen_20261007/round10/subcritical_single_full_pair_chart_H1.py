"""The three FULL global pair equations, regularized on single RR chart.

The normalization is kappa*S_h/tau^2. No layer is discarded.
"""
from sage.all import *
from pathlib import Path
import time,os,resource,json
p=Path(__file__).resolve().parent;st=time.monotonic()
ch=load(str(p/'independent_subcritical_single_RR_full_regular_chart_H1.sobj'))
d=load(str(p/'independent_ratio_p_boundary_global_pair_syzygies.sobj'))
P=ch['ring'];K=P.base_ring();k,r=K.gens();m=1+k;q=1+r;de=P.gen(6);kap=P.gen(8)
old=d['ring'].base_ring();ph=old.ring().hom([k/r,k*q/(r*m),1/r,-q/r],K)
R=PolynomialRing(K,names=d['ring'].variable_names());hm=R.hom(list(ch['chart_p'])+list(ch['chart_Rbar']),P)
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+450*1024**2,)*2)
cache={};items=[];checks=[]
for it in d['items']:
 layers={}
 for ee,c in it['full'].dict().items():
  if c not in cache:cache[c]=ph(c.numerator())/ph(c.denominator())
  deg=sum(ee);rd=sum(ee[6:]);layers.setdefault((deg,rd),{})[tuple(map(int,ee))]=cache[c]
 ans=P.zero();reg=[]
 for (deg,rd),co in layers.items():
  pw=9-2*deg-rd;dp=4-deg;assert dp>=0
  f=hm(R(co))
  if pw<0:
   assert all(int(ee[8])>=-pw for ee in f.dict()),(it['pair'],deg,rd,pw)
   checks.append((it['pair'],deg,rd,pw));reg.append((deg,rd,pw))
   f=P({tuple(int(n)+(pw if j==8 else 0) for j,n in enumerate(ee)):c for ee,c in f.dict().items()});pw=0
  ans+=kap**pw*de**dp*f
 items.append({'pair':it['pair'],'source':ans,'regularity':reg});print('PAIR',it['pair'],'terms',len(ans.dict()),'sec',time.monotonic()-st,flush=True)
save({'ring':P,'items':items,'checks':checks,'scope':__doc__},str(p/'independent_subcritical_single_RR_full_pair_chart_H1.sobj'),compress=False)
print('PASS',len(checks),'sec',time.monotonic()-st,flush=True)
