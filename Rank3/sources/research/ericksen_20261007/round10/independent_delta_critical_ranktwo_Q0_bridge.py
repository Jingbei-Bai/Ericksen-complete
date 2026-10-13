"""Exact source link from the full critical canonical-amplitude graph to F8.

No inverse of C is used: the determinant-cleared Q0 compatibility is an
explicit combination of its row and the three R rows.
"""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic()
d=load(str(p/'independent_delta_critical_canonical_amplitude.sobj'))
q0=load(str(p/'independent_delta_subcritical_Q0_direction.sobj'))
K=d['field'];a,r,H,Z,s,z=K.gens();q=1+r;M=d['material_matrix']
C=M[:3,:3];off=vector(K,[M[i,4] for i in range(3)]);qe=vector(K,[M[3,j] for j in range(3)])
det=C.det();adj=C.adjugate();compat=det*M[3,4]-qe*adj*off
assert all(M[i,3]==0 for i in range(M.nrows()))
ph=q0['field'].ring().hom([a,r,H,Z,1/s,z],K)
f=q0['condition'];expected=ph(f.numerator())/ph(f.denominator())
ratio=compat/expected
assert compat==ratio*expected
print('RATIO',ratio.factor(),'seconds',time.monotonic()-st,flush=True)
# This is a polynomial row identity after clearing only the recorded source
# denominators. It does not use det(C) != 0.
row=det*M.row(3)-sum(((qe*adj)[i]*M.row(i) for i in range(3)),vector(K,5))
assert tuple(row[:4])==(0,0,0,0) and row[4]==compat
save({'source_files':['independent_delta_critical_canonical_amplitude.sobj','independent_delta_subcritical_Q0_direction.sobj'],
      'field':K,'C':C,'C_determinant':det,'Q0_compatibility':compat,
      'direction_condition_image':expected,'protected_ratio':ratio,
      'row_multipliers':(det,tuple(-v for v in qe*adj)),
      'row_identity':row,'scope':__doc__},str(p/'independent_delta_critical_ranktwo_Q0_bridge.sobj'),compress=False)
(p/'independent_delta_critical_ranktwo_Q0_bridge.json').write_text(json.dumps({'ratio':str(ratio.factor()),'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n')
