"""Every cyclic denominator used by the raw RR2 rankone reader is protected."""
from sage.all import *
from pathlib import Path
import os,resource,time,json
p=Path(__file__).resolve().parent;st=time.time()
d=load(str(p/'independent_subcritical_two_RR_inactive_rankone_raw_reader.sobj'))
vm=int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE');resource.setrlimit(resource.RLIMIT_AS,(vm+350*1024**2,vm+350*1024**2))
P=d['minors'][0]['minor'].parent().ring();K=P.base_ring();a,r,H,Z=K.gens();KR=K.ring();ph=KR.hom([a,r,1,-1],KR);guards=[KR(a),KR(r),KR(1+r),KR(1+a*r),KR(a-1)];checked={};count=0
def guard(c):
 global count
 c=K(c);dd=ph(c.denominator())
 if dd in checked:return
 if not dd:raise AssertionError('cyclic denominator zero')
 rem=dd
 for f in guards:
  while True:
   quo,res=rem.quo_rem(f)
   if res:break
   rem=quo
 if rem.total_degree()!=0:raise AssertionError(('unprotected material denominator',rem))
 checked[dd]=True;count+=1
for row in d['original_rows'].values():
 for f in row:
  for pp in [f.numerator(),f.denominator()]:
   for c in pp.coefficients():guard(c)
for f in d['raw_H15_layers'].values():
 for c in f.coefficients():guard(c)
for f in d['complete_even_equations']:
 for c in f.coefficients():guard(c)
for item in d['minors']:
 for key in ['wco','tco']:guard(item[key])
 for c in item['den'].coefficients():guard(c)
for f in [d['range_base_determinant'],d['even_resultant']]:guard(f)
save({'cyclic_denominators':list(checked),'allowed_material_guards':guards,'scope':'a,r,1+r,1+ar,a-1 are nonzero for each normalized positive distinct spectrum'},str(p/'independent_subcritical_two_RR_inactive_rankone_material_units_reader.sobj'),compress=False)
rep={'PASS':True,'distinct_cyclic_denominators':count,'seconds':time.time()-st,'guards':['a','r','q=1+r','m=1+ar','a-1']};(p/'independent_subcritical_two_RR_inactive_rankone_material_units_reader.json').write_text(json.dumps(rep,indent=2));print('PASS',rep,flush=True)
