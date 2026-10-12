"""Division-free true second cokernel of all 37 full sources.

All range corrections are removed by explicit polynomial row combinations.
No tangent correction dl,du,dv and no tau kernel coordinate is suppressed.
"""
from sage.all import *
from pathlib import Path
import os,resource,time,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_ratio_incoming12_generic_nullwall_all.sobj'));assert d['complete'];U=d['ring'];K=U.base_ring();l=U.gen();mod=d['relation'];items={int(it['source_row']):it for it in d['items']}
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+700*1024**2,)*2)
def mul(f,g):return (f*g)%mod
A=[[items[j]['linear_row'][c] for c in (1,2)] for j in (1,17)]
det=(A[0][0]*A[1][1]-A[0][1]*A[1][0])%mod;adj=[[A[1][1],-A[0][1]],[-A[1][0],A[0][0]]];c5=items[5]['linear_row'][0];out=[];checks=0
for row in range(37):
 it=items[row];cf=[sum(mul(it['linear_row'][j+1],adj[j][i]) for j in (0,1))%mod for i in (0,1)]
 co=[mul(det,c5),-mul(c5,cf[0]),-mul(c5,cf[1]),-mul(det,it['linear_row'][0])]
 rows=[row,1,17,5]
 for j in range(6):
  assert sum(mul(c,items[rr]['linear_row'][j]) for c,rr in zip(co,rows))%mod==0;checks+=1
 rhs=[sum(mul(c,items[rr]['rhs'][j]) for c,rr in zip(co,rows))%mod for j in range(15)]
 out.append({'source_row':row,'rhs':rhs,'cofactors':co,'range_rows':[1,17,5]})
 print('row',row,'support',[j for j,f in enumerate(rhs) if f],'sec',round(time.monotonic()-st,2),flush=True)
dd=-mod[0]
def norm(c):
 c=U(c)%mod;return K(c[0]**2-dd*c[1]**2)
norms={key:norm(c) for key,c in [('two_range_determinant',det),('odd_range_pivot',c5)]}
facts={key:f.factor() for key,f in norms.items()}
save({'ring':U,'relation':mod,'quadratic_exponents':d['quadratic_exponents'],'tangent_order':d['tangent_order'],'items':out,'range_rows':[1,17,5],'two_range_determinant':det,'odd_range_pivot':c5,'pivot_norms':norms,'pivot_norm_factors':facts,'scope':__doc__},str(p/'independent_ratio_incoming12_generic_nullwall_all_cokernel.sobj'),compress=False)
(p/'independent_ratio_incoming12_generic_nullwall_all_cokernel.json').write_text(json.dumps({'checks':checks,'supports':{str(it['source_row']):[j for j,f in enumerate(it['rhs']) if f] for it in out},'pivot_norm_factors':{key:str(v) for key,v in facts.items()},'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n')
print('DONE',time.monotonic()-st,flush=True)
