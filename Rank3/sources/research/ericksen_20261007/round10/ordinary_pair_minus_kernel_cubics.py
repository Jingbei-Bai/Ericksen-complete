"""Literal pure-k cubic coefficients from original full Q/R/H sources.
At v=-2 the incoming family is exactly RR=(-2*l,u,-2). The pure kernel
ray is p=(1,-k,0,0,k,l), eta=tau=0. Full signed H DAGs are retained.
These are raw cubic coefficients, not by themselves cokernel equations.
"""
from sage.all import *
from pathlib import Path
import time,json,resource,os
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
g=load(str(p/'independent_ratio_p_boundary_second_gauge.sobj'));hd=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'));oq=load(str(p.parent/'round8/rank3_joint_fiber_leading.sobj'))
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+650*1024**2,)*2)
OK=g['coefficient_field'];O=PolynomialRing(OK,names=('l','u'));l,u=O.gens()
def factor(f):
 dd={}
 for ee,c in (f if isinstance(f,dict) else f.dict()).items():
  e=tuple(map(int,ee))
  if e[2] or e[3]:continue
  kd=e[1]+e[4]
  if kd>3:continue
  key=(sum(e),kd)
  dd[key]=dd.get(key,O.zero())+c*(-1)**e[1]*(-2)**(e[6]+e[8])*l**(e[5]+e[6])*u**e[7]
 return dd
rows=[];records=[]
for qi in (0,1,5,6):
 f=oq['quartic_source_rows'][qi];D=int(f.total_degree());dd=factor(f);rows.append(dd.get((D,3),O.zero()));records.append(('Q',qi))
for j,f in enumerate(g['six_full_sources']):
 rows.append(factor(f).get((4,3),O.zero()));records.append(('R',j))
cache={}
def pool(i):
 if i not in cache:cache[i]=factor(hd['polynomial_pool'][i])
 return cache[i]
for j,it in enumerate(hd['records']):
 ans=O.zero()
 for sign,li,ri in it['full_signed_product_DAG']:
  A=pool(li);B=pool(ri)
  for (dl,kl),fl in A.items():
   if (5-dl,3-kl) in B:ans+=sign*fl*B[(5-dl,3-kl)]
 rows.append(ans);records.append(('H',j));print('H',j,'rawcube terms',len(ans.dict()),flush=True)
K=PolynomialRing(QQ,names=('a','r','H','Z')).fraction_field();a,r,H,Z=K.gens();q=1+r;m=1+a*r;ph=OK.hom([a,a*q/m,H/r,q*Z/r],K);U=PolynomialRing(K,names=('l','u'));ll,uu=U.gens();co={}
def cv(c):
 c=OK(c)
 if c not in co:co[c]=ph(c)
 return co[c]
mapped=[U({tuple(map(int,e)):cv(c) for e,c in f.dict().items()}) for f in rows]
pair=load(str(p/'ordinary_pair_incoming_full_twojet.sobj'))
for it in pair['items']:
 f=it['full_homogeneous_pullback'];ans=U.zero()
 for ee,c in f.dict().items():
  e=tuple(map(int,ee))
  if e[4] or e[5] or e[7] or e[8]:continue
  if e[3]+e[6]!=3:continue
  ans+=K(c)*(-1)**e[3]*(-2)**e[2]*ll**e[0]*uu**e[1]
 mapped.append(ans);records.append(('pair',it['pair']))
for i,z in enumerate(mapped):print('ROW',i,records[i],'value',z.factor() if z else 0,flush=True)
save({'ring':U,'raw_ring':O,'raw_cubic_rows':rows,'mapped_cubic_rows':mapped,'records':records,'scope':__doc__},str(p/'ordinary_pair_minus_kernel_cubics.sobj'),compress=False)
print('DONE',time.monotonic()-st,flush=True)
