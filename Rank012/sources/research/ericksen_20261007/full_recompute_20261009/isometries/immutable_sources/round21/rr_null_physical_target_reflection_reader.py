"""Full-parameter target reflection identifies the two null-axis sheets."""
from sage.all import *
from pathlib import Path
import time,signal,json
O=Path(__file__).resolve().parent;st=time.time();signal.alarm(90)
p=load(str(O/'rr_null_physical_quadratics_plus.sobj'));m=load(str(O/'rr_null_physical_quadratics_minus.sobj'))
K=p['field'];T=p['ring'];s,t,pp,q,c=K.gens();cf=K.hom([s,-t,pp,q,c],K)
def reflect(f):return T({ex:cf(co)*((-1)**ex[1]) for ex,co in f.dict().items()})
checks={f'conic_{j}_full_reflection':reflect(p['raw_quadratics'][j])==m['raw_quadratics'][j] for j in range(3)}
assert all(checks.values()),checks
rep={'status':'PASS_FULL_PARAMETER_TWO_SHEET_REFLECTION','checks':checks,'seconds':time.time()-st,
 'map':'Q_minus(s,t,p,q,c;alpha,beta,chi)=Q_plus(s,-t,p,q,c;alpha,-beta,chi)',
 'geometry':'target O=diag(1,-1,1); angular pseudovector transforms as det(O)*O*Omega'}
(O/'rr_null_physical_target_reflection_reader.json').write_text(json.dumps(rep,indent=2),encoding='utf8');print(json.dumps(rep,indent=2),flush=True)
