"""Primitive polynomial row elimination of the COMPLETE second cokernel.

No fraction-field matrix solve. Every row content and elimination multiplier
is recorded. Results are coefficient-generic until their pivots are audited.
"""
from sage.all import *
from pathlib import Path
import time,os,resource,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_ratio_incoming12_generic_nullwall_all_cokernel.sobj'));U=d['ring'];K=U.base_ring();R=K.ring();l=U.gen();mod=d['relation']
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+650*1024**2,)*2)
powers=[2,1,1,0,0,0,1,0,0,2,1,1,2,1,1];rows={};initial_scalings={};pool={};operations=[];nextid=0
def add(v,label):
 global nextid
 den=lcm([cc.denominator() for cc in v]);pv=[R(cc*den) for cc in v];g=gcd(pv);pv=[q//g for q in pv] if g else pv
 if any(pv):
  lc=next(q.leading_coefficient() for q in pv if q);pv=[R(q/lc) for q in pv];g*=lc
 ident=nextid;nextid+=1;pool[ident]=pv
 return ident,den,g
for it in d['items']:
 row=int(it['source_row']);vv=[(f*l**j)%mod for f,j in zip(it['rhs'],powers)];nz=[f for f in vv if f]
 if not nz:continue
 parity=0 if nz[0][0] else 1;assert all(f[1-parity]==0 for f in vv)
 cv=[f[parity] for f in vv];den=lcm([cc.denominator() for cc in cv]);pv=[R(cc*den) for cc in cv];g=gcd(pv);rows[row]=[f//g for f in pv];initial_scalings[row]={'denominator':den,'content':g,'l_parity':parity}
groups={'even':[0,2,3,5,11],'odd':[1,4,8,10,12,14]};modules={}
def complexity(v):return sum(len(f.dict()) for f in v),max([int(f.degree()) for f in v if f]+[0])
for name,cols in groups.items():
 work=[]
 for j,row in rows.items():
  v=[row[c] for c in cols]
  if any(v):
   ident,den,g=add(v,('source',j));operations.append({'kind':'source','id':ident,'source_row':j,'columns':cols,'denominator':den,'content':g});work.append(ident)
 basis=[]
 for col in range(len(cols)):
  choices=[i for i in work if pool[i][col]]
  if not choices:continue
  pivot=min(choices,key=lambda i:complexity(pool[i]));work.remove(pivot);pv=pool[pivot];basis.append((col,pivot));fresh=[]
  for old in work:
   row=pool[old]
   if not row[col]:fresh.append(old);continue
   mul1,mul2=pv[col],row[col];v=[mul1*a-mul2*b for a,b in zip(row,pv)]
   ident,den,g=add(v,('elimination',old,pivot))
   operations.append({'kind':'eliminate','id':ident,'left':old,'pivot':pivot,'left_multiplier':mul1,'pivot_multiplier':mul2,'denominator':den,'content':g})
   if any(pool[ident]):fresh.append(ident)
  work=fresh
  print(name,'pivot',col,'size',complexity(pv),'remaining',len(work),'sec',time.monotonic()-st,flush=True)
  save({'ring':R,'modules':modules,'current':name,'basis':basis,'pool':pool,'operations':operations,'initial_scalings':initial_scalings},str(p/'independent_ratio_incoming12_second_module_fractionfree.partial.sobj'),compress=False)
 modules[name]={'columns':cols,'basis':basis,'rank':len(basis),'echelon_rows':[pool[i] for c,i in basis]}
save({'ring':R,'modules':modules,'pool':pool,'operations':operations,'initial_scalings':initial_scalings,'source_rows':rows,'scope':__doc__},str(p/'independent_ratio_incoming12_second_module_fractionfree.sobj'),compress=False)
meta={'seconds':time.monotonic()-st,'modules':{n:{'rank':v['rank'],'supports':[[v['columns'][i] for i,f in enumerate(row) if f] for row in v['echelon_rows']]} for n,v in modules.items()},'operations':len(operations),'scope':__doc__}
(p/'independent_ratio_incoming12_second_module_fractionfree.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta,indent=2),flush=True)
