"""Original six N-free direction equations centered at the known graph."""
from pathlib import Path
import json,time
from sage.all import load,save,PolynomialRing,QQ,matrix
from itertools import combinations
out=Path(__file__).resolve().parent;start=time.monotonic()
d=load(str(out/'general_curve_unit_gradient_directions.sobj'));oldR=d['ring'];oldK=oldR.base_ring()
C=PolynomialRing(QQ,names=('H','Z','S'));K=C.fraction_field();H,Z,S=K.gens();L=H*S-Z
R=PolynomialRing(K,names=('X','Y','T'));X,Y,T=R.gens();images=[X+S/(2*Z),Y-L/(2*H*Z),T-L/(2*H*Z)]
def transfer(f):
 p=R.zero()
 for m,c in f.dict().items():
  assert not c.derivative(oldK.gen(0)) and not c.derivative(oldK.gen(1)) and not c.derivative(oldK.gen(5))
  co=K(c.numerator()(0,0,H,Z,S,0))/K(c.denominator()(0,0,H,Z,S,0))
  p+=co*images[0]**int(m[0])*images[1]**int(m[1])*images[2]**int(m[2])
 return p
sources=[transfer(f) for f in d['mixed_differences']+d['old_cubic_rows']];assert all(not f.constant_coefficient() for f in sources)
linear=matrix(K,[[f.monomial_coefficient(v) for v in R.gens()] for f in sources]);records=[]
for ii in combinations(range(6),3):
 det=linear.matrix_from_rows(ii).det();facts=det.factor() if det else None
 rec={'rows':ii,'zero':not bool(det),'factorization':str(facts) if facts is not None else None};records.append(rec);print(rec,flush=True)
tmp=out/'general_curve_centered_six_sources_next.sobj';save({'ring':R,'sources':sources,'linear_matrix':linear,'simple_center':images,'records':records,'scope':'Exact centered source equations and infinitesimal Jacobian minors. Full rank only gives local simplicity, not absence of other components.'},str(tmp));tmp.replace(out/'general_curve_centered_six_sources.sobj')
(out/'general_curve_centered_six_sources.json').write_text(json.dumps({'records':records,'elapsed':time.monotonic()-start},indent=2)+'\n')
