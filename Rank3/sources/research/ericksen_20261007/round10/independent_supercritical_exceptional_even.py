"""Original even boundary sources on the canonical supercritical finite cover."""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic()
d=load(str(p/'independent_delta_critical_linear_sources.sobj'));can=load(str(p/'independent_delta_critical_canonical_amplitude.sobj'))
K=can['field'];a,r,H,Z,s,z=K.gens();q=1+r;P=can['extension_ring'];b,c=P.gens();old=d['field']
ph=old.ring().hom([P(a),P(r),P(H),P(Z),P(s),c*z,s*b,c],P)
def red(f):
 out={}
 for ee,v in P(f).dict().items():
  e=tuple(map(int,ee));key=(e[0]%2,e[1]%2)
  out[key]=out.get(key,K.zero())+v*(-q)**(e[0]//2)*(-q/r)**(e[1]//2)
 return P(out)
def cv(v):
 v=old(v);nn=red(ph(v.numerator()));dd=red(ph(v.denominator()));assert len(dd.dict())==1
 ee,co=next(iter(dd.dict().items()));i,j=map(int,ee)
 return red(nn*b**i*c**j/(co*(-q)**i*(-q/r)**j))
E=PolynomialRing(K,names=('se0','se1','se2'));es=E.gens();even=[]
for f in d['even_sources']:
 dd={}
 for ee,co in f.dict().items():
  e=tuple(map(int,ee));powers=[sum(e[j]*can['column_phases'][j][k] for j in range(3)) for k in range(2)]
  v=red(cv(co)*b**powers[0]*c**powers[1]);assert v==P(v.constant_coefficient())
  dd[e]=K(v.constant_coefficient())
 even.append(E(dd))
save({'field':K,'e_ring':E,'even':even,'root_e_phases':can['column_phases'][:3],'scope':__doc__},str(p/'independent_supercritical_exceptional_even.sobj'),compress=False)
V=load(str(p/'independent_supercritical_normal_minors.sobj'))['reduced_matrix']
# Line lambda=0, retaining z != 0,-1. No additional direction is divided silently.
lk=PolynomialRing(QQ,names=('la','lr','lH','lZ','lt')).fraction_field();aa,rr,hh,zz,t=lk.gens();lh=K.ring().hom([aa,rr,hh,zz,-t/(t+1),t],lk)
def lc(v):v=K(v);return lh(v.numerator())/lh(v.denominator())
LV=matrix(lk,V.nrows(),4,[lc(v) for v in V.list()]);A=LV[:3,:3];det=A.det();print('LINE DET',det.factor(),flush=True)
graph=-A.adjugate()*LV[:3,3]/det
print('LINE GRAPH',[v.factor() for v in graph.list()],flush=True)
assert LV*vector(lk,list(graph.list())+[1])==0
ev=[]
for f in even:
 v=sum((lc(co)*prod(graph[j,0]**int(ex[j]) for j in range(3)) for ex,co in f.dict().items()),lk.zero());ev.append(v)
 print('LINE EVEN',len(ev)-1,v.factor(),'seconds',time.monotonic()-st,flush=True)
save({'field':lk,'normal_matrix':LV,'determinant':det,'e_graph':graph,'even_pullbacks':ev,'scope':'Full supercritical lambda=0 first-normal candidate, three actual old even sources.'},str(p/'independent_supercritical_line_even.sobj'),compress=False)
print('DONE','seconds',time.monotonic()-st,flush=True)
