"""Exact n=lambda exceptional-wall audit, retaining every allowed material."""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_delta_critical_line_native_resultants.sobj'));B=PolynomialRing(QQ,'wn');n=B.gen();T=PolynomialRing(B,'wt');t=T.gen();ph=d['ring'].hom([t,T(n),T(n)],T);polys=[]
for f in d['flat_polynomials']:
 g=ph(f);qq,rr=g.quo_rem(t);assert rr==0;polys.append(qq)
print('WALL_POLYS',[f.factor() for f in polys],flush=True)
f,g=polys[:2];m=int(f.degree());k=int(g.degree());fc=list(reversed(f.list()));gc=list(reversed(g.list()));rows=[]
for i in range(k):rows.append([B.zero()]*i+fc+[B.zero()]*(k-1-i))
for i in range(m):rows.append([B.zero()]*i+gc+[B.zero()]*(m-1-i))
S=matrix(B,rows);result=S.det();factor=result.factor();print('WALL_RESULTANT',factor,'seconds',time.monotonic()-st,flush=True)
bound=sum(max(f.degree() for f in S.row(i)) for i in range(S.nrows()));assert result.degree()<=bound
for j in range(int(bound)+1):assert matrix(ZZ,S.nrows(),S.ncols(),[f(j) for f in S.list()]).det()==result(j)
save({'ring':T,'polynomials_after_t_unit':polys,'sylvester':S,'resultant':result,'factorization':factor,'strict_degree_bound':bound,'exact_integer_nodes':int(bound)+1,'scope':__doc__},str(p/'independent_delta_critical_line_equal_material_wall.sobj'),compress=False)
(p/'independent_delta_critical_line_equal_material_wall.json').write_text(json.dumps({'resultant_factors':str(factor),'degree_bound':int(bound),'integer_nodes':int(bound)+1,'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n')
