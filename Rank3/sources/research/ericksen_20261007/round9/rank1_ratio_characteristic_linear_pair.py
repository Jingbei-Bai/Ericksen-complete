"""Two actual mixed equations as a 2x2 affine system in p21,p31.

The determinant's zero stratum is retained, not implicitly localized away.
"""
from pathlib import Path
import time,json
from sage.all import load,save,matrix,vector
out=Path(__file__).resolve().parent;start=time.monotonic()
d=load(str(out/'rank1_ratio_characteristic_second_graph.sobj'));R=d['ring'];K=R.base_ring();a,b,c,h=d['parameters'];u,v,x,s,y=R.gens()
f=d['remaining_second_tangencies'][1];g=d['remaining_second_tangencies'][2]
scales=[f.leading_coefficient(),g.leading_coefficient()];f=f/scales[0];g=g/scales[1]
M=matrix(R,[[f.derivative(x),f.derivative(y)],[g.derivative(x),g.derivative(y)]])
bb=vector(R,[f.subs({x:0,y:0}),g.subs({x:0,y:0})])
assert vector(R,[f,g])==M*vector(R,[x,y])+bb
assert all(not e.derivative(x) and not e.derivative(y) for e in M.list()+list(bb))
det=M.det();xn=M[0,1]*bb[1]-M[1,1]*bb[0];yn=M[1,0]*bb[0]-M[0,0]*bb[1]
checks={'cramer_row0':M[0,0]*xn+M[0,1]*yn+bb[0]*det==0,
 'cramer_row1':M[1,0]*xn+M[1,1]*yn+bb[1]*det==0}
assert all(checks.values())
support=[str(R.monomial(*m)) for m in det.dict()]
s_power=min(m[3] for m in det.dict())
det_reduced=R({tuple(e-(s_power if j==3 else 0) for j,e in enumerate(m)):co for m,co in det.dict().items()})
assert s**s_power*det_reduced==det
payload={'ring':R,'parameters':d['parameters'],'linear_source_rows':[f,g],
 'row_normalization_factors':scales,'matrix':M,'constant_vector':bb,'determinant':det,
 'determinant_s_power':int(s_power),'reduced_determinant':det_reduced,
 'cramer_numerators':[xn,yn],'solved_variables':[x,y],
 'remaining_mixed_source':d['remaining_second_tangencies'][5],
 'quartic_source_Q0':d['selected_quartic_rows'][0],
 'checks':checks,'coverage':'det!=0 Cramer graph and det=0 retained as separate necessary affine consistency stratum.'}
save(payload,str(out/'rank1_ratio_characteristic_linear_pair.sobj'))
report={'status':'PASS_EXACT_AFFINE_PAIR','checks':checks,'s_power':int(s_power),
 'determinant_terms':len(det.dict()),'determinant_support':support,
 'reduced_determinant_degree':int(det_reduced.total_degree()),
 'cramer_terms':[len(xn.dict()),len(yn.dict())],'elapsed_seconds':time.monotonic()-start,
 'scope':payload['coverage']}
(out/'rank1_ratio_characteristic_linear_pair.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2),flush=True)
