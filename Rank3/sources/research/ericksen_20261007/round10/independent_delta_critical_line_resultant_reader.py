"""Exact polynomial factor replay and integer-grid audit of three resultants."""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_delta_critical_line_native_resultants.sobj'));R=d['ring'];t,n,l=R.gens();B=PolynomialRing(QQ,names=('nn','ll'));nn,ll=B.gens();ph=R.hom([B.zero(),nn,ll],B);fs=d['flat_polynomials'];items=[];checks=[]
def syl(f,g):
 m=int(f.degree(t));k=int(g.degree(t));fc=[ph(f.coefficient({t:i})) for i in range(m,-1,-1)];gc=[ph(g.coefficient({t:i})) for i in range(k,-1,-1)];rows=[]
 for i in range(k):rows.append([B.zero()]*i+fc+[B.zero()]*(k-1-i))
 for i in range(m):rows.append([B.zero()]*i+gc+[B.zero()]*(m-1-i))
 return matrix(B,rows)
for i,j in ((0,1),(0,2),(1,2)):
 text=(p/f'independent_delta_critical_line_resultant_{i}{j}.txt').read_text().strip();res=B(text.replace('\n',''));fac=res.factor();assert fac.value()==res;S=syl(fs[i],fs[j]);bounds=[sum(max(int(f.degree(v)) for f in S.row(k)) for k in range(S.nrows())) for v in (nn,ll)];assert res.degree(nn)<=bounds[0] and res.degree(ll)<=bounds[1]
 for v1 in range(1,bounds[0]+2):
  for v2 in range(1,bounds[1]+2):
   mm=matrix(ZZ,S.nrows(),S.ncols(),[f(v1,v2) for f in S.list()]);assert mm.det()==res(v1,v2),(i,j,v1,v2)
 checks.append({'pair':(i,j),'bounds':bounds,'integer_nodes':(bounds[0]+1)*(bounds[1]+1)})
 items.append({'pair':(i,j),'resultant':res,'factorization':fac,'sylvester':S,'bounds':bounds});print('RESULTANT',i,j,'factor',fac,'grid',checks[-1],'seconds',time.monotonic()-st,flush=True)
save({'ring':B,'items':items,'checks':checks,'scope':__doc__},str(p/'independent_delta_critical_line_resultant_reader.sobj'),compress=False)
(p/'independent_delta_critical_line_resultant_reader.json').write_text(json.dumps({'checks':checks,'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n')
