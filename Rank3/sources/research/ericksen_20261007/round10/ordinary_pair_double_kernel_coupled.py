"""Full relative two-kernel module at u=0,v=-2, including tangent couplings.
J=(range4,u*eta,(v+2)*k,eta^2,eta*k,k^2), k=n3,n0=range0-k.
All raw quadratic columns are kept together with complete range columns.
"""
from sage.all import *
from pathlib import Path
import time,json,resource,os
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
d=load(str(p/'independent_ratio_incoming12_generic_normal.sobj'));qq=load(str(p/'independent_ratio_incoming12_generic_normal_two.sobj'));pairs=load(str(p/'ordinary_pair_incoming_full_twojet.sobj'))
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+750*1024**2,)*2)
O=d['ring'];OK=O.base_ring();ol,ou,ov=O.gens();K=pairs['coefficient_field'];a,r,H,Z=K.gens();q=1+r;m=1+a*r;ph=OK.hom([a,a*q/m,H/r,q*Z/r],K)
P=PolynomialRing(K,'l');l=P.gen();cache={};exps=[tuple(map(int,e)) for e in qq['normal_exponents']]
def ck(n,z):assert bool(z),n;checks.append(n)
def cv(c):
 c=OK(c)
 if c not in cache:cache[c]=ph(c)
 return cache[c]
def mp(f):
 dd={}
 for e,c in O(f).dict().items():
  if e[1]:continue
  j=int(e[0]);dd[j]=dd.get(j,K.zero())+cv(c)*(-2)**int(e[2])
 return P(dd)
def ex(i,j):
 z=[0]*6;z[i]+=1;z[j]+=1;return tuple(z)
def rowit(i):
 if i<qq['coefficient_matrix'].nrows():return list(qq['coefficient_matrix'].row(i))
 return list(load(str(p/('independent_ratio_incoming12_generic_normal_two.row%d.sobj'%i)))['coefficients'])
rows=[];labels=[]
for i in range(d['matrix'].nrows()):
 M=list(d['matrix'].row(i));two=rowit(i)
 ck('eta u divisibility at boundary '+str(i),not mp(M[4]))
 ck('kernel k at minus two '+str(i),not mp(M[3]-M[0]))
 e2=two[exps.index(ex(4,4))];ek=two[exps.index(ex(3,4))]-two[exps.index(ex(0,4))]
 k2=two[exps.index(ex(0,0))]-two[exps.index(ex(0,3))]+two[exps.index(ex(3,3))]
 row=[mp(M[j]) for j in (0,1,2,5)]+[mp(M[4].derivative(ou)),mp((M[3]-M[0]).derivative(ov)),mp(e2),mp(ek),mp(k2)]
 rows.append(row);labels.append(('old',i));print('ROW',i,'SUP',[j for j,c in enumerate(row) if c],flush=True)
PU=pairs['jet_ring'];pl,pu,pv=PU.gens();pes=[tuple(e) for e in pairs['normal_exponents']]
def pm(f):
 dd={}
 for e,c in f.dict().items():
  if e[1]:continue
  j=int(e[0]);dd[j]=dd.get(j,K.zero())+K(c)*(-2)**int(e[2])
 return P(dd)
for it in pairs['items']:
 val=it['jet_coefficients'];M=val[1:7]
 e2=val[pes.index(ex(4,4))];ek=val[pes.index(ex(3,4))]-val[pes.index(ex(0,4))];k2=val[pes.index(ex(0,0))]-val[pes.index(ex(0,3))]+val[pes.index(ex(3,3))]
 row=[pm(M[j]) for j in (0,1,2,5)]+[pm(M[4].derivative(pu)),pm((M[3]-M[0]).derivative(pv)),pm(e2),pm(ek),pm(k2)]
 rows.append(row);labels.append(('pair',it['pair']));print('PAIR',it['pair'],'SUP',[j for j,c in enumerate(row) if c],flush=True)
MM=matrix(P,rows)
save({'ring':P,'matrix':MM,'labels':labels,'columns':('range0','n1','n2','tau','u_eta','w_k','eta_square','eta_k','k_square'),'checks':checks,'scope':__doc__},str(p/'ordinary_pair_double_kernel_coupled.sobj'),compress=False)
# A finite-field rank only selects a compact exact minor; it is not a proof.
F=GF(32003);bp={a:F(2),r:F(-1)/3,H:F(2),Z:F(-1)}
br=K.ring()
def hom(z):return sum((F(c)*prod(bp[K.gens()[j]]**int(e[j]) for j in range(4)) for e,c in z.dict().items()),F.zero())
def evc(c):
 c=K(c);return hom(c.numerator())/hom(c.denominator())
FM=matrix(F,[[sum(evc(c)*F(3)**int(e) for e,c in z.dict().items()) for z in row] for row in MM.rows()])
pivs=FM.transpose().pivots();print('ALGORITHM_ONLY_RANK',FM.rank(),'PIVOT_ROWS',pivs,flush=True)
save({'ring':P,'matrix':MM,'labels':labels,'columns':('range0','n1','n2','tau','u_eta','w_k','eta_square','eta_k','k_square'),'algorithm_only_rows':pivs,'checks':checks,'scope':__doc__},str(p/'ordinary_pair_double_kernel_coupled.sobj'),compress=False)
print('DONE',len(checks),time.monotonic()-st,flush=True)
