"""Exact WHOLE-family linear-normal factorization for the relative Rees bridge.

This is stronger than a first jet at the extra-tau wall. No state division
is used. The range IFT is a separate local argument using saved pivots.
"""
from sage.all import *
from pathlib import Path
import os,resource,time,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_ratio_incoming12_generic_normal.sobj'));U=d['ring'];K=d['coefficient_field'];l,u,v=U.gens();a,b,h,z=K.gens();M=d['matrix'];r=(a-b)/(a*(b-1));q=1+r;Z=r*z/q;kap=(Z-1)/Z;dd=-(a-b)*z*z/(b*(a-1));g=l*l-dd;w=v+2;items=[];checks=0
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+600*1024**2,)*2)
def ck(test):
 global checks
 assert test;checks+=1
def div(f,s):
 # Avoid libSingular // over a nested rational function coefficient field.
 fd={tuple(map(int,e)):K(c) for e,c in f.dict().items()};sd={tuple(map(int,e)):K(c) for e,c in s.dict().items()};lead=max(sd);out={}
 while fd:
  e=max(fd);ck(all(e[i]>=lead[i] for i in range(3)));ee=tuple(e[i]-lead[i] for i in range(3));cc=fd[e]/sd[lead];out[ee]=out.get(ee,K(0))+cc
  for se,sc in sd.items():
   key=tuple(ee[i]+se[i] for i in range(3));fd[key]=fd.get(key,K(0))-cc*sc
   if not fd[key]:del fd[key]
 z=U(out);ck(z*s==f);return z
for j in range(M.nrows()):
 eta=U(M[j,4]);kk=U(M[j,3]-M[j,0]);tt=U(M[j,5]+kap*l*M[j,0])
 ce=div(eta,u);ckk=div(kk,w)
 t0=U(tt.subs({u:0}));t00=U(t0.subs({v:-2}))
 tu=div(tt-t0,u);tw=div(t0-t00,w);tg=div(t00,g)
 ck(eta==u*ce);ck(kk==w*ckk);ck(tt==u*tu+w*tw+g*tg)
 items.append({'source_row':j,'five_coefficients':[ce,ckk,tu,tw,tg],'transformed_normal_columns':[kk,eta,tt]})
save({'ring':U,'coefficient_field':K,'wall_polynomial':g,'normal_transform':'p13=range0-k+kap*l*tau, p31=k; eta unchanged','kappa':kap,'module_generators':['u*eta','(v+2)*k','u*tau','(v+2)*tau','g*tau'],'items':items,'range_columns':[0,1,2],'range_source_rows':[1,17,5],'source':'independent_ratio_incoming12_generic_normal.sobj','scope':__doc__},str(p/'independent_ratio_incoming12_relative_linear_module.sobj'),compress=False)
# P must be simple over the generic base for every fixed allowed spectrum.
j=load(str(p/'independent_ratio_incoming12_nullwall_second_compatibility.sobj'));pp=list(j['primitive_quadratic']);disc=pp[1]**2-4*pp[0]*pp[2];ck(bool(disc))
A=PolynomialRing(QQ,names=('a','b'));aa,bb=A.gens();metadata=[]
for label,f in [('P_leading_numerator',pp[2].numerator()),('P_discriminant_numerator',disc.numerator()),('P_discriminant_denominator',disc.denominator())]:
 groups={}
 for ee,cc in f.dict().items():
  e=tuple(map(int,ee));groups.setdefault(e[2:],{})[e[:2]]=cc
 key=max(groups);lead=A(groups[key]);res=lead
 for guard in (aa,bb,aa-1,bb-1,aa-bb):
  while res and res%guard==0:res=res//guard
 metadata.append({'label':label,'base_exponent':key,'leading_coefficient':lead,'old_guard_only':res.degree()==0})
print('discriminant',[{'label':it['label'],'base_exponent':it['base_exponent'],'factor':str(it['leading_coefficient'].factor()),'old_guard_only':it['old_guard_only']} for it in metadata],flush=True)
save({'coefficient_field':K,'P':pp,'discriminant':disc,'base_leading_checks':metadata},str(p/'independent_ratio_incoming12_relative_quadratic_discriminant.sobj'),compress=False)
meta={'checks':checks,'sources':len(items),'seconds':time.monotonic()-st,'P_checks':[{'label':it['label'],'base_exponent':list(it['base_exponent']),'factors':str(it['leading_coefficient'].factor()),'old_guard_only':it['old_guard_only']} for it in metadata],'scope':__doc__}
(p/'independent_ratio_incoming12_relative_linear_module.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta,indent=2),flush=True)
