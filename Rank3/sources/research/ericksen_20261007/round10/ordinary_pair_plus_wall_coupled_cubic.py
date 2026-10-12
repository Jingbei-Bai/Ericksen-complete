"""True eight-generator coupled cubic module on m*l^2=k, u!=0, v=-2.
Generators: range5,(v+2)*kernel,(m*l^2-k)*kernel^2,kernel^3.
No pure cubic is used without the full range columns.
"""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
d=load(str(p/'ordinary_pair_full_normal_blocks.sobj'));c=load(str(p/'ordinary_pair_minus_kernel_cubics.sobj'));P0=d['ring'];K=P0.base_ring();a,r,H,Z=K.gens();q=1+r;k=a*r;m=1+k;ll,uu,vv=P0.gens();P=PolynomialRing(K,names=('l','u'));l,u=P.gens();labels=d['labels'];M=d['matrix'];cm=c['ring'].base_ring().hom(list(K.gens()),K)
def ck(n,x):assert bool(x),n;checks.append(n)
def pull(f):
 dd={}
 for e,z in f.dict().items():
  key=(int(e[0]),int(e[1]));dd[key]=dd.get(key,K.zero())+K(z)*(-2)**int(e[2])
 return P(dd)
def cu(lab):
 typ,i=lab
 j=i if typ=='old' else 28+i
 return P({tuple(map(int,e)):cm(z) for e,z in c['mapped_cubic_rows'][j].dict().items()})
def row(lab):
 z=list(M.row(labels.index(lab)));wk=pull((z[3]-z[0]).derivative(vv));ck('whole k kernel '+str(lab),pull(z[3]-z[0])==0)
 return [pull(z[j]) for j in (0,1,2,4,5)]+[wk,P(8) if lab==('old',2) else P.zero(),cu(lab)]
allrs=[('old',i) for i in (0,1,2,3,5,7,11,13,17,19)]+[('pair',0),('pair',2)];A=matrix(P,[row(z) for z in allrs])
def det(B):
 if not B.nrows():return P.one()
 i=min(range(B.nrows()),key=lambda i:sum(bool(z) for z in B.row(i)))
 return sum(((-1)**(i+j)*B[i,j]*det(B.delete_rows([i]).delete_columns([j])) for j in range(B.ncols()) if B[i,j]),P.zero())
sets=[(('old',1),('old',5),('old',7),('old',11),('old',13),('pair',2)),(('old',1),('old',5),('old',7),('old',11),('pair',0),('pair',2)),(('old',17),('old',5),('old',7),('old',11),('old',13),('pair',2))]
items=[]
def rem(f):
 ans=P.zero()
 for e,z in f.dict().items():ans+=z*(k/m)**(int(e[0])//2)*l**(int(e[0])%2)*u**int(e[1])
 return ans
for rs in sets:
 z=det(A.matrix_from_rows_and_columns([allrs.index(x) for x in rs],(0,1,2,4,5,7)));rr=rem(z)
 print('ROWS',rs,'terms',len(z.dict()),'wallterms',len(rr.dict()),'wallfactors',rr.factor() if rr else 0,'sec',time.monotonic()-st,flush=True)
 items.append({'rows':rs,'six_minor':z,'wall_remainder':rr})
 save({'ring':P,'matrix':A,'labels':allrs,'items':items,'checks':checks,'scope':__doc__},str(p/'ordinary_pair_plus_wall_coupled_cubic.sobj'),compress=False)
print('DONE',time.monotonic()-st,flush=True)
