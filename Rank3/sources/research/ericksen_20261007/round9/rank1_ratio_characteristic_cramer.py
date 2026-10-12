"""Cramer's rule with its protected determinant kept factored."""
from pathlib import Path
import time,json
from copy import copy
from sage.all import load,save
out=Path(__file__).resolve().parent;start=time.monotonic()
d=load(str(out/'rank1_ratio_characteristic_source.sobj'));R=d['ring'];M=d['t_coefficient_matrix'];b=d['constant_terms'];M0=M.matrix_from_rows(range(6));det=d['selected_minor']
N=[];records=[]
for j in range(6):
 A=copy(M0)
 for i in range(6):A[i,j]=-b[i]
 n=A.det();N.append(n)
 records.append({'t_index':j,'degree':int(n.total_degree()),'terms':len(n.dict())})
 save({'ring':R,'parameters':d['parameters'],'numerators':N,'denominator':det,'records':records},str(out/'rank1_ratio_characteristic_cramer_partial.sobj'))
 print('Cramer',j,'degree',n.total_degree(),'terms',len(n.dict()),'seconds',round(time.monotonic()-start,2),flush=True)
payload=dict(d);payload.update(t_numerators=N,t_denominator=det,records=records,
 proof_rule='Cramer replacement-column determinants for the six selected source equations; determinant was independently checked to be protected_scalar*p23.')
save(payload,str(out/'rank1_ratio_characteristic_cramer.sobj'))
print('saved_t_graph',round(time.monotonic()-start,2),flush=True)
res=[]
for i in (6,7):
 rr=sum((M[i,j]*N[j] for j in range(6)),R.zero())+b[i]*det
 res.append(rr)
 print('residual',i,'degree',rr.total_degree(),'terms',len(rr.dict()),'seconds',round(time.monotonic()-start,2),flush=True)
payload.update(residual_numerators=res)
save(payload,str(out/'rank1_ratio_characteristic_cramer.sobj'))
report={'status':'EXACT_CRAMER_GRAPH_AND_RESIDUALS','records':records,
        'residuals':[{'source_row':i+6,'degree':int(p.total_degree()),'terms':len(p.dict())} for i,p in enumerate(res)],
        'elapsed_seconds':time.monotonic()-start,
        'scope':'Source-level Cramer graph with protected determinant retained. No fiber dimension or exclusion claim.'}
(out/'rank1_ratio_characteristic_cramer.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2),flush=True)
