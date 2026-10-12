from sage.all import *
from pathlib import Path
import json
from itertools import combinations
OUT=Path(__file__).absolute().parent;d=load(str(OUT/'homothety_invariant_linear_chart.sobj'));M=d['matrix'];d1=d['source'];Mt=d1['matrix'];K=M.base_ring();a,b,s,t=K.gens();out={}
for name,sub in [('s0',{s:0}),('t0',{t:0}),('axis',{s:0,t:0})]:
 MM=M.apply_map(lambda z:z.subs(sub));TT=Mt.apply_map(lambda z:z.subs(sub));out[name]={'observation_rank':MM.rank(),'t_symbol_rank':TT.rank()}
 print(name,out[name],flush=True)
 if name=='axis':
  out[name]['observation_kernel']=[[str(z) for z in row] for row in MM.right_kernel().basis()]
  print('axis_obs_kernel',out[name]['observation_kernel'],flush=True)
 else:
  detobs=MM.matrix_from_columns([i for i in range(9) if i!=0]).det()
  out[name]['observation_protected_minor']=str(detobs.factor())
  for rows in combinations(range(9),6):
   dd=TT.matrix_from_rows(rows).det()
   if not dd:continue
   # Require a single monomial in the surviving real direction variable.
   surviving=3 if name=='s0' else 2
   num=dd.numerator();den=dd.denominator()
   if len({m[surviving] for m in num.dict()})==1 and len({m[surviving] for m in den.dict()})==1:
    out[name]['t_protected_rows']=rows;out[name]['t_protected_minor']=str(dd.factor());break
  print('protected',out[name],flush=True)
save({'source':d,'boundary_ranks':out},str(OUT/'homothety_observation_boundary_rank.sobj'))
(OUT/'homothety_observation_boundary_rank.json').write_text(json.dumps(out,indent=2))
