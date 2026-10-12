"""Genuine seven-generator module on u=0, retaining the v=-2 intersection.

Generators: five range normals, u*eta, eta^2. Source Taylor divisibility
follows from exact whole-family vanishing and the eta kernel at u=0.
No raw quadratic is individually treated as a cokernel equation.
"""
from sage.all import *
from pathlib import Path
import time,resource,os,json
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
d=load(str(p/'independent_ratio_incoming12_generic_normal.sobj'));q=load(str(p/'independent_incoming_simple_wall_quadratics.sobj'));pairs=load(str(p/'ordinary_pair_incoming_full_twojet.sobj'))
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+650*1024**2,)*2)
U=d['ring'];OK=U.base_ring();ol,ou,ov=U.gens();K=pairs['coefficient_field'];a,r,H,Z=K.gens();qq=1+r;mm=1+a*r
ph=OK.hom([a,a*qq/mm,H/r,qq*Z/r],K);P=PolynomialRing(K,names=('oc_l','oc_v'));l,v=P.gens();cache={}
def ck(n,b):assert bool(b),n;checks.append(n)
def cv(c):
 c=OK(c)
 if c not in cache:cache[c]=ph(c)
 return cache[c]
def mp(f):
 f=U(f);dd={}
 for e,c in f.dict().items():
  e=tuple(map(int,e))
  if e[1]:continue
  dd[(e[0],e[2])]=cv(c)
 return P(dd)
et=next(x for x in q['items'] if x['name']=='eta');eta2={int(x['row']):x['second'] for x in et['values']}
rows=[];labels=[];sel=list(range(10))+[11,13,17,19,29,32,35]
for i in sel:
 ck('eta kernel at u0 '+str(i),not d['matrix'][i,4].subs({ou:0}))
 row=[mp(d['matrix'][i,j]) for j in (0,1,2,3,5)]+[mp(d['matrix'][i,4].derivative(ou)),mp(eta2.get(i,U.zero()))]
 rows.append(row);labels.append(('old',i));print('ROW',i,'SUPPORT',[j for j,c in enumerate(row) if c],'sec',time.monotonic()-st,flush=True)
PU=pairs['jet_ring'];pl,pu,pv=PU.gens();es=pairs['normal_exponents'];etaexp=(0,0,0,0,2,0)
def pmap(f):return P({(int(e[0]),int(e[2])):c for e,c in f.dict().items() if not e[1]})
for it in pairs['items']:
 vals=it['jet_coefficients'];first=vals[1:7]
 ck('pair eta kernel '+str(it['pair']),not first[4])
 row=[pmap(first[j]) for j in (0,1,2,3,5)]+[pmap(first[4].derivative(pu)),pmap(vals[es.index(etaexp)])]
 rows.append(row);labels.append(('pair',it['pair']))
M=matrix(P,rows)
save({'ring':P,'matrix':M,'labels':labels,'columns':('n0','n1','n2','n3','tau','u_eta','eta_square'),'checks':checks,'scope':__doc__},str(p/'ordinary_pair_eta_coupled_matrix.sobj'),compress=False)
print('DONE',len(checks),time.monotonic()-st,flush=True)
