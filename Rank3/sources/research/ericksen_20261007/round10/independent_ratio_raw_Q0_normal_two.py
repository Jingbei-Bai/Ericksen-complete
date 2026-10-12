"""The three UNMULTIPLIED full hatFi(Q0) sources: complete family two-jets.

Their polynomial fields are recovered by exact monomial division of the
stored p_s*hatFi factors. This is a coefficient identity, not localization
at p_s. All three sources remain necessary on inactive-D supports.
"""
from sage.all import *
from pathlib import Path
from itertools import combinations_with_replacement
import time,os,resource,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_ratio_global_Q0_active_degreefall.sobj'));U=d['family_ring'];K=U.base_ring();l,u,v=U.gens();A=d['lambda']*l
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+650*1024**2,)*2)
pool=dict(d['polynomial_pool']);exps=[(0,)*6]
for j in range(6):
 e=[0]*6;e[j]=1;exps.append(tuple(e))
for i,j in combinations_with_replacement(range(6),2):
 e=[0]*6;e[i]+=1;e[j]+=1;exps.append(tuple(e))
slots={e:i for i,e in enumerate(exps)};cache={}
def factor(ref):
 if ref in cache:return cache[ref]
 out={}
 for ee,c in pool[ref].items():
  e=tuple(map(int,ee));ne=e[1:5];sz=sum(ne)
  if sz>2:continue
  for j in range(min(e[6],2-sz)+1):
   key=(sum(e),ne+(j,));out[key]=out.get(key,U.zero())+K(c)*binomial(e[6],j)*A**(e[6]-j)*l**e[5]*u**e[7]*v**e[8]
 cache[ref]={e:f for e,f in out.items() if f};return cache[ref]
def prod2(ff,gg):
 out=[U.zero() for _ in exps]
 for (dl,nl),fl in ff.items():
  for (dr,nr),fr in gg.items():
   ne=tuple(nl[j]+nr[j] for j in range(5));tp=5-dl-dr
   if tp<0 or sum(ne)+tp>2:continue
   out[slots[ne+(tp,)]]+=fl*fr
 return vector(U,out)
def remove_monomial(ref,pcol,newref):
 dd={}
 for ee,c in pool[ref].items():
  e=list(map(int,ee));assert e[pcol]>=1;e[pcol]-=1;dd[tuple(e)]=c
 pool[newref]=dd
items=[]
for i,pcol in enumerate((0,2,4)):
 dag=[]
 for j in range(9):
  name='rawF%d_%d'%(i,j);remove_monomial('p%dF%d_%d'%(pcol,i,j),pcol,name);dag.append((name,'Qd%d'%j))
 for j in range(2):
  name='rawB%d_%d'%(i,j);remove_monomial('p%dB%d_%d'%(pcol,i,j),pcol,name);dag.append((name,'Qb%d'%j))
 value=sum((prod2(factor(aa),factor(bb)) for aa,bb in dag),vector(U,[0]*len(exps)))
 assert value[0]==0
 items.append({'direction':i,'full_product_DAG':dag,'coefficients':value[1:]})
 print(i,'first-normal counts',[len(q.dict()) for q in value[1:7]],'quadratic nonzero',sum(bool(q) for q in value[7:]),'seconds',round(time.monotonic()-st,1),flush=True)
 save({'ring':U,'coefficient_field':K,'normal_exponents':exps[1:],'polynomial_pool':pool,'items':items,'source':'independent_ratio_global_Q0_active_degreefall.sobj','scope':__doc__},str(p/'independent_ratio_raw_Q0_normal_two.sobj'),compress=False)
(p/'independent_ratio_raw_Q0_normal_two.json').write_text(json.dumps({'sources':len(items),'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n')
print('DONE',time.monotonic()-st,flush=True)
