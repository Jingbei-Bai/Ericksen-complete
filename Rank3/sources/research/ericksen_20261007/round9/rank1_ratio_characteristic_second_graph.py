"""Protected u reconstruction from the two diagonal L tangencies."""
from pathlib import Path
import time,json
from sage.all import load,save,vector
out=Path(__file__).resolve().parent;start=time.monotonic()
d=load(str(out/'rank1_ratio_characteristic_second_tangency.sobj'));R=d['ring'];K=R.base_ring();a,b,c,h=d['parameters']
r=(a-b)/(a*(b-1));q=1+r;z=q*h+c;s=R.gen(3)
U=d['u_coefficient_matrix'];low=d['scaled_lower_rows'];M=U.matrix_from_rows([4,8]);det=M.det()
checks={'protected_u_determinant':det==-2*q*q*(1-r*c/q)**2}
un=-M.inverse()*vector(R,[low[4],low[8]])
checks['diagonal_22']=low[4]+sum(U[4,j]*un[j] for j in range(2))==0
checks['diagonal_33']=low[8]+sum(U[8,j]*un[j] for j in range(2))==0
assert all(checks.values()),checks
res=[low[i]+sum(U[i,j]*un[j] for j in range(2)) for i in range(9)]
payload=dict(d);payload.update(u_numerators_over_p23=un,remaining_second_tangencies=res,checks=checks,
 protected_u_determinant=det)
save(payload,str(out/'rank1_ratio_characteristic_second_graph.sobj'))
print('u_graph_saved',round(time.monotonic()-start,2),flush=True)
records=[]
for i in [1,2,5]:
 f=res[i];rec={'row':i,'degree':int(f.total_degree()),'terms':len(f.dict()),
 'support':[str(R.monomial(*m)) for m in f.dict()]}
 records.append(rec);print(rec,flush=True)
ql=load(str(out/'rank1_ratio_characteristic_quartics_partial.sobj'))
Q=[];qrec=[]
for i in [0,10,11,13,15,18,19]:
 f=ql['Q_scaled_lower'][i]+s*sum(ql['Q_u_coefficients'][i][j]*un[j] for j in range(2))
 Q.append(f);rec={'row':i,'degree':int(f.total_degree()) if f else -1,'terms':len(f.dict())}
 qrec.append(rec);print('quartic',rec,round(time.monotonic()-start,2),flush=True)
 payload.update(selected_quartic_indices=[0,10,11,13,15,18,19][:len(Q)],selected_quartic_rows=Q)
 temp=out/'rank1_ratio_characteristic_second_graph_next.sobj';save(payload,str(temp));temp.replace(out/'rank1_ratio_characteristic_second_graph.sobj')
report={'status':'PASS_PROTECTED_ACTUAL_U_GRAPH','checks':checks,'mixed_records':records,
 'selected_quartics':qrec,'u_scale':'u_j=un_j/p23','quartic_scale':'p23^2',
 'scope':'Actual E1 pencil necessary source system. No fiber dimension conclusion.',
 'elapsed_seconds':time.monotonic()-start}
(out/'rank1_ratio_characteristic_second_graph.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2),flush=True)
