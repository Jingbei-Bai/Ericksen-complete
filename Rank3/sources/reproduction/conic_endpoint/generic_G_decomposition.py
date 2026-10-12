"""Exact complete-source monomial decomposition, valid over Q(k,r,eta,zeta)."""
from sage.all import *
import json,hashlib
def verify_and_export(P,G,out):
    parts=[{},{}]
    for ex,c in G.dict().items():
        e=list(map(int,ex))
        if e[7]:
            e[7]-=1;parts[0][tuple(e)]=c
        else:
            assert e[8]>0,'G has a monomial outside (t,kappa)'
            e[8]-=1;parts[1][tuple(e)]=c
    Gt,Gk=(P(v) for v in parts)
    assert G==P.gen(7)*Gt+P.gen(8)*Gk
    record={'status':'PASS_COMPLETE_GENERIC_G_IN_T_KAPPA',
            'coefficient_field':'QQ(k,r,eta,zeta)','variables':list(P.variable_names()),
            'identity':'G = t*Gt + kappa*Gkappa','priority':['t','kappa'],
            'G_terms':len(G.dict()),'Gt_terms':len(Gt.dict()),'Gkappa_terms':len(Gk.dict()),
            'delta_square_cofactor_needed':False,'exact_polynomial_identity_verified':True}
    sparse=out/'generic_G_ideal_decomposition_sparse.json'
    if not sparse.exists():
        sparse.write_text(json.dumps({**record,
            'Gt':[{'exponents':list(map(int,e)),'coefficient':str(c)} for e,c in Gt.dict().items()],
            'Gkappa':[{'exponents':list(map(int,e)),'coefficient':str(c)} for e,c in Gk.dict().items()]},indent=2)+'\n')
    record['sparse_cofactors_sha256']=hashlib.sha256(sparse.read_bytes()).hexdigest()
    (out/'generic_G_ideal_decomposition.json').write_text(json.dumps(record,indent=2)+'\n')
    return record
