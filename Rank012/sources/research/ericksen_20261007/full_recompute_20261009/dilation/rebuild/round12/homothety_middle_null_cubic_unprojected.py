"""Check whether the generic cubic source vanishes before forcing projection."""
from sage.all import *
from pathlib import Path
import json
OUT=Path(__file__).absolute().parent;d=load(str(OUT/'homothety_middle_null_cubic_solution.sobj'));src=d['source'];S=src['ring'];E,t2,h2,z2=S.gens();sol=d['solution'];sub=dict(zip(S.gens(),sol));K=S.base_ring()
rem=[K(f.subs(sub)) for f in src['raw_cubic']]
record={'status':'EXACT_UNPROJECTED_CUBIC_CHECK','nonzero':[(i,str(z.factor())) for i,z in enumerate(rem) if z],
 'zero_count':sum(z==0 for z in rem)}
save({'field':K,'residuals':rem},str(OUT/'homothety_middle_null_cubic_unprojected.sobj'))
(OUT/'homothety_middle_null_cubic_unprojected.json').write_text(json.dumps(record,indent=2));print(json.dumps(record))
