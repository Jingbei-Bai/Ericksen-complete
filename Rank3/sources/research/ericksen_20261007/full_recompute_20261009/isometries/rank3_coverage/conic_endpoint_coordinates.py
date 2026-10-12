from pathlib import Path
from sage.all import *
import json
h=Path(__file__).absolute().parent;o=h.parents[2]/'round10'
d=load(str(o/'subcritical_single_cyclic_three_p_conic_quadratic.sobj'));M=d['matrix'];R=d['range_columns'];N=d['tangent_columns'];E=d['ring'];K=E.base_ring();k,r=K.gens();x=E.gen();checks=[]
def ck(n,v):assert v,n;checks.append(n)
ck('full37 actual normals annihilate allthree proposed tangent columns',not bool(M*N))
J=R.augment(N);dt=J.det();ck('proposed range and tangent coordinates form true nine-dimensional basis',bool(dt))
# Norm of this coordinate determinant must have only protected spectral factors.
nm=matrix(K,2,2,lambda i,j:(dt*x**j)[i]).det();guard=K.ring()(k*r*(1+k)*(1+r))
for part in [nm.numerator(),nm.denominator()]:
 for fac,pow in part.factor():
  ck('coordinate-unit factor '+str(fac),fac in QQ or guard.quo_rem(fac)[1]==0)
report=dict(status='PASS',checks=checks,coordinate_determinant=str(dt),coordinate_norm=str(nm),source_evidence='rank3/source_replays/subcritical_single_cyclic_three_p_conic_reader/replay_record.json verifies every entry of M against newly rebuilt whole H1 chart/G9 sources')
(h/'conic_endpoint_coordinates.json').write_text(json.dumps(report,indent=2));print(report,flush=True)
