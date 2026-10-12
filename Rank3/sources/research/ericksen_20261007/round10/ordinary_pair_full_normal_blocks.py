"""Genuine full incoming normal matrix and exact eta-row elimination on u!=0."""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
old=load(str(p/'independent_ratio_incoming12_generic_normal.sobj'));pairs=load(str(p/'ordinary_pair_incoming_full_twojet.sobj'))
K=pairs['coefficient_field'];a,r,H,Z=K.gens();q=1+r;m=1+a*r
O=old['ring'];OK=O.base_ring();ph=OK.hom([a,a*q/m,H/r,q*Z/r],K)
P=PolynomialRing(K,names=('l','u','v'));l,u,v=P.gens();cache={}
def ck(n,z):assert bool(z),n;checks.append(n)
def cv(c):
 c=OK(c)
 if c not in cache:cache[c]=ph(c)
 return cache[c]
def mp(f):return P({tuple(map(int,e)):cv(c) for e,c in O(f).dict().items()})
ids=[0,1,2,3,4,5,6,7,8,9,11,13,17,19,29,32,35]
rows=[];labels=[]
for i in ids:
 rows.append([mp(old['matrix'][i,j]) for j in range(6)]);labels.append(('old',i))
for it in pairs['items']:
 rows.append([P({tuple(map(int,e)):K(c) for e,c in z.dict().items()}) for z in it['jet_coefficients'][1:7]]);labels.append(('pair',it['pair']))
M=matrix(P,rows);a0=P(M[0,4]/u);ck('eta pivot',a0==Z**2*(H*Z-H+Z)/(7*(Z-1)**2*(Z-H)))
N=[]
for i,row in enumerate(rows):
 t=P(row[4]/u)
 nr=[P(row[j]-t/a0*rows[0][j]) for j in (0,1,2,3,5)]
 N.append(nr)
 print('ROW',labels[i],'eta/u',str(t) if len(t.dict())<=4 else str(len(t.dict()))+'terms','supports',[j for j,c in enumerate(nr) if c],'udeg',[int(c.degree(u)) for c in nr],flush=True)
N=matrix(P,N)
save({'ring':P,'matrix':M,'reduced_matrix':N,'labels':labels,'eta_pivot_unit':a0,'noneta_columns':(0,1,2,3,5),'checks':checks,'scope':__doc__},str(p/'ordinary_pair_full_normal_blocks.sobj'),compress=False)
print('DONE',len(checks),time.monotonic()-st,flush=True)
