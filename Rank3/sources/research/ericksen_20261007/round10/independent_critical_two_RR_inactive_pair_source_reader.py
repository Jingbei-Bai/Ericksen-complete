from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(s,b):
 if not bool(b):raise AssertionError(s)
 checks.append(s)
src=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'));O=src['ring'];K=O.base_ring();a,r,H,Z=K.gens();q=1+r;al=(Z-H)/(Z*(H+1));be=H*(Z-1)/(Z*(H+1));P=PolynomialRing(K,'T');T=P.gen();F=P.fraction_field();vals=[0,0,1,0,T,0,0,be*T,1];flows=[{5:-1/al},{3:-1/(al*T)},{7:1}]
full=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'));old=full['coefficient_field'];pool=full['polynomial_pool'];ph=old.ring().hom([a,a*q/(1+a*r),H/r,q*Z/r],K);cache={}
def cv(c):
 if c not in cache:cache[c]=ph(c.numerator())/ph(c.denominator())
 return cache[c]
def face(f,i,j):return {tuple(map(int,ex)):co for ex,co in f.dict().items() if (sum(ex[:6]),sum(ex[6:]))==(i,j)}
def pull(dd,flow=None,orig=False):
 z=F.zero()
 for ex,co in dd.items():
  ee=list(map(int,ex));cc=cv(co) if orig else K(co)
  if flow is None:z+=cc*prod(F(vals[k])**ee[k] for k in range(9))
  else:
   for k,v in flow.items():
    if not ee[k]:continue
    xx=ee.copy();xx[k]-=1;z+=ee[k]*F(v)*cc*prod(F(vals[h])**xx[h] for h in range(9))
 return P(z)
rows=[]
for j,f in enumerate(src['full_R_sources']):
 top=face(f,3,1);coef=[pull(top,x) for x in flows];ck('R constant zero '+str(j),pull(face(f,4,0))==0);ck('R free e2 coefficient zero '+str(j),coef[2]==0);rows.append(coef[:2]+[pull(face(f,1,2))])
f=src['full_Q0'];coef=[pull(face(f,1,2),x) for x in flows];ck('Q0 constant zero',pull(face(f,2,1))==0);ck('Q0 free e2 coefficient zero',coef[2]==0);rows.append(coef[:2]+[P.zero()])
lay={(3,2):{},(4,1):{},(1,3):{}}
for sign,ii,jj in full['records'][6]['full_signed_product_DAG']:
 for ex,c in pool[ii].items():
  for ey,b in pool[jj].items():
   ee=tuple(int(x)+int(y) for x,y in zip(ex,ey));key=(sum(ee[:6]),sum(ee[6:]))
   if key in lay:lay[key][ee]=lay[key].get(ee,old.zero())+sign*c*b
coef=[pull(lay[(3,2)],x,True) for x in flows];ck('complete H6 constant zero',pull(lay[(4,1)],orig=True)==0);ck('complete H6 free e2 coefficient zero',coef[2]==0);rows.append(coef[:2]+[pull(lay[(1,3)],orig=True)])
M=matrix(P,rows);saved=load(str(p/'independent_critical_two_RR_inactive_pair_certificate.sobj'));sp=saved['ring'];sg=sp.gen();mapc=sp.hom([T],P);ck('original R/Q0/H6 multiplyback selected full boundary rows',M==matrix(P,[list(saved['matrix'].row(j)) for j in list(range(7))+[13]]))
first=M.matrix_from_rows([0,2,6]).det();second=M.matrix_from_rows([0,2,7]).det();AA=(-5*H+4*Z)*T*T+q*(-4*H+5*Z);BB=5*T*T+4*q
uA=K(F(first)/(T*T*AA));uB=K(F(second)/(T**4*BB));ck('true first three-minor strip',first==uA*T*T*AA);ck('true second three-minor strip',second==uB*T**4*BB);ck('exact protected cofactor identity',5*AA-(-5*H+4*Z)*BB==9*q*Z)
for nm,uu in [('A',uA),('B',uB)]:
 ck('nonzero original minor unit '+nm,uu!=0)
 for ff,nn in uu.factor():
  if ff in QQ:continue
  ck('unit only pure base factors '+nm,ff.degree(K.ring().gen(0))==0 and ff.degree(K.ring().gen(1))==0)
 print(nm,'UNIT',uu.factor(),flush=True)
rec={'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st,'scope':'critical two-RR inactive rankD2 paired support, all materials, arbitrary RR ratio and all finite e, all ramification; original full-source proof'}
save({'checks':checks,'matrix':M,'A':AA,'B':BB,'unit_A':uA,'unit_B':uB,'record':rec},str(p/'independent_critical_two_RR_inactive_pair_source_reader.sobj'),compress=False);(p/'independent_critical_two_RR_inactive_pair_source_reader.json').write_text(json.dumps(rec,indent=2)+'\n');print('PASS',rec,flush=True)
