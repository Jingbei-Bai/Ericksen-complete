"""Exact generic solve of the full middle-null cubic face."""
from sage.all import *
from pathlib import Path
import json,time
OUT=Path(__file__).absolute().parent;START=time.time();d=load(str(OUT/'homothety_middle_null_cubic_face.sobj'))
S=d['ring'];K=S.base_ring();a,p,q,r,A2,B2,k1=K.gens();vv=S.gens();src=d['forcing_free_cubic']
M=matrix(K,[[f.monomial_coefficient(v) for v in vv] for f in src]);rhs=vector(K,[-f.constant_coefficient() for f in src])
def sample(c):return QQ(c.numerator()(QQ(1)/2,1,3,2,3,4,5))/QQ(c.denominator()(QQ(1)/2,1,3,2,3,4,5))
ms=M.apply_map(sample);rows=[];rank=0
order=sorted(range(M.nrows()),key=lambda i:sum(len(c.numerator().dict())+len(c.denominator().dict()) for c in M.row(i) if c))
for i in order:
 rk=ms.matrix_from_rows(rows+[i]).rank()
 if rk>rank:rank=rk;rows.append(i)
cols=list(ms.matrix_from_rows(rows).pivots());assert len(cols)==4
A=M.matrix_from_rows_and_columns(rows,cols);sol=A.inverse()*vector(K,[rhs[i] for i in rows]);res=M*sol-rhs
summary={'status':'EXACT_FACE_ONLY','rows':rows,'minor':str(A.det().factor()),
 'solution':{str(v):str(z.factor()) if z else '0' for v,z in zip(vv,sol)},
 'nonzero_compatibility':[(i,str(z.factor())) for i,z in enumerate(res) if z],
 'seconds':float(time.time()-START)}
save({'source':d,'matrix':M,'rhs':rhs,'solution':sol,'residuals':res,'summary':summary},str(OUT/'homothety_middle_null_cubic_solution.sobj'))
(OUT/'homothety_middle_null_cubic_solution.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
