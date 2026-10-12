from sage.all import *
from pathlib import Path
from itertools import combinations
import time,json
OUT=Path(__file__).absolute().parent;d=load(str(OUT/'homothety_first_jet_reduction.sobj'));M=d['matrix'];K=M.base_ring();a,b,s,t=K.gens();START=time.time()
records=[];factors={}
for rows in combinations(range(9),6):
 det=M.matrix_from_rows(rows).det()
 if not det:continue
 ff=det.factor();records.append((rows,det))
 extra=[]
 for fac,power in ff:
  if fac not in [s,t,a,b,a-1,b-1,a-b]:
   if any(m[2] or m[3] for m in fac.numerator().dict()):extra.append((fac,power))
 key=str(extra)
 if key not in factors:
  factors[key]=rows;print('rows',rows,'factor',ff,flush=True)
save({'minors':records,'representatives':factors},str(OUT/'homothety_first_jet_minor_atlas.sobj'))
(OUT/'homothety_first_jet_minor_atlas.json').write_text(json.dumps({'nonzero_minors':len(records),'representatives':[(str(rows),key) for key,rows in factors.items()],'seconds':time.time()-START},indent=2))
print('COMPLETE',len(records),'unique',len(factors),time.time()-START,flush=True)
