"""Literal complete beta normal matrices on every four-p full-RR support."""
from sage.all import *
from pathlib import Path
import time
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_RR_torus_missing_p_sources.sobj'));old=d['field'];a0,r0,H0,Z0,x0,b0,c0=old.gens()
K=PolynomialRing(QQ,names=('na','nr','nH','nZ','ny')).fraction_field();a,r,H,Z,y=K.gens();q=1+r;P=PolynomialRing(K,names=('nb','nc'));b,c=P.gens()
def red(f):
 dd={}
 for ex,co in P(f).dict().items():
  e=tuple(map(int,ex));key=(e[0]%2,e[1]%2);dd[key]=dd.get(key,K.zero())+co*(-q)**(e[0]//2)*(-q/r)**(e[1]//2)
 return P(dd)
items=[]
for chart,it in enumerate(d['items']):
 xf=[-y*b/q,P(y),y*b][chart];ph=old.ring().hom([P(a),P(r),P(H),P(Z),xf,b,c],P)
 def cv(v):
  v=old(v);nn=red(ph(v.numerator()));dd=red(ph(v.denominator()));assert len(dd.dict())==1,(it['label'],dd)
  ex,co=next(iter(dd.dict().items()));i,j=map(int,ex);return red(nn*b**i*c**j/(co*(-q)**i*(-q/r)**j))
 M=it['matrix'];sgn=[-1,0,1][chart];rows=[]
 for i in range(M.nrows()):
  f=M[i,3];assert cv(f)==0,(it['label'],i)
  jb=cv(b0*f.derivative(b0)+sgn*x0*f.derivative(x0));jc=cv(c0*f.derivative(c0))
  rows.append([jb,jc]+[cv(M[i,j]) for j in (0,1,2,4)])
 MM=matrix(P,rows);eqs=[];rhs=[[],[]]
 for i in range(MM.nrows()):
  for j in range(6):
   f=MM[i,j]
   if not f:continue
   assert len(f.dict())==1,(it['label'],i,j,f)
   ex=tuple(map(int,next(iter(f.dict()))));row=[0]*(MM.nrows()+5);row[i]=1
   if j<5:row[MM.nrows()+j]=1
   eqs.append(row)
   for k in range(2):rhs[k].append(ex[k])
 AA=matrix(GF(2),eqs);sol=[AA.solve_right(vector(GF(2),v)) for v in rhs];rp=[tuple(int(sol[k][i]) for k in range(2)) for i in range(MM.nrows())];cp=[tuple(int(sol[k][MM.nrows()+j]) for k in range(2)) for j in range(5)]+[(0,0)]
 N=matrix(K,MM.nrows(),6)
 for i in range(MM.nrows()):
  ri,rj=rp[i];iv=b**ri*c**rj/((-q)**ri*(-q/r)**rj)
  for j in range(6):
   ci,cj=cp[j];f=red(MM[i,j]*b**ci*c**cj*iv);assert f==P(f.constant_coefficient());N[i,j]=K(f.constant_coefficient())
 basis=[];pivs=[];rs=[];pivots=[]
 for i in range(N.nrows()):
  v=N.row(i)
  for j,w in zip(pivs,basis):v=v-v[j]*w
  if v:
   j=next(j for j in range(6) if v[j]);pv=v[j];basis.append(v/pv);pivs.append(j);rs.append(i);pivots.append(pv)
   print('NORMAL PIVOT',it['label'],i,j,pv.factor(),'seconds',time.monotonic()-st,flush=True)
  if len(rs)==6:break
 items.append({'label':it['label'],'matrix':N,'extension_matrix':MM,'row_phases':rp,'column_phases':cp,'pivot_rows':rs,'pivot_columns':pivs,'pivots':pivots})
 save({'field':K,'root_ring':P,'items':items,'scope':__doc__},str(p/'independent_RR_torus_missing_p_super_normal.sobj'),compress=False)
print('DONE','seconds',time.monotonic()-st,flush=True)
