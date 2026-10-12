"""Protected-denominator audit of the global even p-cubic/RR-linear module."""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_even_cubic_RR_linear_module.sobj'));src=load(str(p/'independent_ratio_weighted_degreefall.sobj'));guard=src['oldguard_numerator'];K=d['p_ring'].base_ring();checks=0;items=[]
for it in d['items']:
 assert it['member'];factors=[]
 for j,f in enumerate(it['polynomial_cofactors']):
  for c in f.coefficients():
   c=K(c)
   for g,n in c.denominator().factor():
    assert not guard.quo_rem(g)[1];checks+=1;factors.append((j,str(g),int(n)))
 print('SOURCE',it['source'],'protectedcofactor_denominators',sorted(set(factors)),flush=True);items.append({'source':it['source'],'factors':factors})
save({'source_file':'independent_even_cubic_RR_linear_module.sobj','items':items,'checks':checks,'old_guard':guard,'scope':__doc__},str(p/'independent_even_cubic_module_unit_reader.sobj'),compress=False)
report={'status':'PASS','checks':checks,'seconds':time.monotonic()-st,'scope':__doc__};(p/'independent_even_cubic_module_unit_reader.json').write_text(json.dumps(report,indent=2)+'\n');print('DONE',report,flush=True)
