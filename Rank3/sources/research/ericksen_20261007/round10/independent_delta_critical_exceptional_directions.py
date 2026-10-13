"""Exact seven-direction rank-drop reduction, with both algebraic conjugates."""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_delta_critical_ranktwo_structure.sobj'));T=d['ring'];s,z=T.gens();U=PolynomialRing(QQ,'ww');K=NumberField(U.gen()**2-QQ(13)/2*U.gen()+1,'omega');w=K.gen();cases=[(K.one(),K.one()),(K.one(),w),(w,K.one()),(w,w)];items=[]
for ix,(sv,zv) in enumerate(cases):
 ph=T.hom([sv,zv],K);M=matrix(K,[[ph(v) for v in row] for row in d['distinct_rows']]);rr=M.row_space().basis_matrix();items.append({'case':ix,'s':sv,'z':zv,'matrix':M,'row_basis':rr,'rank':M.rank()});print('CASE',ix,'s,z',sv,zv,'rank',M.rank(),'rowbasis',list(rr.rows()),flush=True)
save({'number_field':K,'minimal_polynomial':K.polynomial(),'cases':items,'all_conjugates_included':True,'scope':__doc__},str(p/'independent_delta_critical_exceptional_directions.sobj'),compress=False)
(p/'independent_delta_critical_exceptional_directions.json').write_text(json.dumps({'ranks':[int(v['rank']) for v in items],'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n')
