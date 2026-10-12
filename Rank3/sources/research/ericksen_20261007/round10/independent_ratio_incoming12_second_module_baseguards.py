"""Inspect coefficient-generic row pivots without specializing h,z or spectra."""
from sage.all import *
from pathlib import Path
import time,json,resource,os
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_ratio_incoming12_second_module_fractionfree.sobj'));R=d['ring'];A=PolynomialRing(QQ,names=('a','b'));a,b=A.gens()
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+500*1024**2,)*2)
items=[]
def check(label,f):
 f=R(f)
 if not f:return
 groups={}
 for ee,c in f.dict().items():
  e=tuple(map(int,ee));groups.setdefault(e[2:],{})[e[:2]]=c
 key=max(groups);lead=A(groups[key]);fac=lead.factor()
 old=[a,b,a-1,b-1,a-b];res=lead
 for v in old:
  while res and res%v==0:res=res//v
 items.append({'label':label,'base_exponent':key,'coefficient':lead,'coefficient_factors':fac,'old_guard_only':res.degree()==0})
for j,it in enumerate(d['operations']):
 check('operation%d content'%j,it['content'])
 check('operation%d denominator'%j,it['denominator'])
 if it['kind']=='eliminate':check('operation%d pivot'%j,it['left_multiplier'])
for j,it in d['initial_scalings'].items():
 check('source%d content'%j,it['content']);check('source%d denominator'%j,it['denominator'])
c=load(str(p/'independent_ratio_incoming12_generic_nullwall_all_cokernel.sobj'))
for key,v in c['pivot_norms'].items():
 check(key+' numerator',v.numerator());check(key+' denominator',v.denominator())
j=load(str(p/'independent_ratio_incoming12_nullwall_second_compatibility.sobj'))
for key in ('range_denominator_norm','reference_leading_norm'):
 v=j[key];check(key+' numerator',v.numerator());check(key+' denominator',v.denominator())
bad=[it for it in items if not it['old_guard_only']]
save({'items':items,'bad':bad,'source':'independent_ratio_incoming12_second_module_fractionfree.sobj'},str(p/'independent_ratio_incoming12_second_module_baseguards.sobj'),compress=False)
meta={'count':len(items),'bad':[{'label':it['label'],'exponent':list(it['base_exponent']),'factors':str(it['coefficient_factors'])} for it in bad],'seconds':time.monotonic()-st,'scope':__doc__}
(p/'independent_ratio_incoming12_second_module_baseguards.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta,indent=2),flush=True)
