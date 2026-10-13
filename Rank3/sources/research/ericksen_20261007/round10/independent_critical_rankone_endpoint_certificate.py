from sage.all import *
from pathlib import Path
import json,time
p=Path(__file__).resolve().parent;st=time.time();checks=[]
def ck(n,b):
 if not b:raise AssertionError(n)
 checks.append(n)
s=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'));O=s['ring'];K=O.base_ring();a,r,H,Z=K.gens();q=1+r;m=1+a*r;KR=K.ring();P=PolynomialRing(K,'t');t=P.gen();ev=O.hom([0,1,0,0,0,0,0,t,1],P);ph=KR.hom([a,r,1,-1],K)
def cp(f):return P({int(e):ph(c.numerator())/ph(c.denominator()) for e,c in f.dict().items()})
def face(f,i,j):return O({ex:c for ex,c in f.dict().items() if (sum(ex[:6]),sum(ex[6:]))==(i,j)})
f=ev(face(s['full_R_sources'][4],1,2));rd=load(str(p/'independent_critical_rankone_original_reader.sobj'));g=ev(rd['H_raw_layers'][16][(1,3)])
for j,rr in enumerate(s['full_R_sources']):ck('endpoint R40 '+str(j),ev(face(rr,4,0))==0)
for j,it in rd['H_raw_layers'].items():ck('endpoint H41 '+str(j),ev(it[(4,1)])==0)
ck('cyclic F12 source exact',cp(f)==(t**2+q)/(16*r));ck('cyclic H13 source exact',cp(g)==-(17*t**2+13*q)/(96*r))
Syl=f.sylvester_matrix(g);Sc=Syl.apply_map(lambda x:ph(K(x).numerator())/ph(K(x).denominator()));Rc=Sc.det()
ck('generic resultant specialization',Rc==q**2/(147456*r**4))
# Every defining coefficient denominator and each even unit is valid at the test
# base for every material. Only these protected material factors may occur.
def monic(x):
 x=KR(x);return x/x.leading_coefficient()
allowed={monic(z) for z in [a,r,q,m,a-1]};seen=set()
for lab,ff in [('F12',f),('H13',g)]:
 for j,c in enumerate(ff.coefficients()):
  de=ph(c.denominator());ck('coefficient cyclic denominator '+lab+str(j),de!=0)
  for z,n in de.factor():
   if KR(z).is_constant():continue
   z=monic(z);ck('protected coefficient denominator '+lab+str(j)+'/'+str(z),z in allowed);seen.add(str(z))
even=load(str(p/'independent_delta_subcritical_even_full_source.sobj'));units=[]
for j,it in enumerate(even['items']):
 ck('endpoint full even p4 '+str(j),ev(O(it['pure_p4']))==0)
 for h,co in enumerate(it['A_p_quadratic']):ck('endpoint full even linear '+str(j)+'/'+str(h),ev(O(co))==0)
 u=K(it['square_unit']);uc=ph(u.numerator())/ph(u.denominator());ck('endpoint square unit cyclic '+str(j),uc!=0);units.append(u)
save({'F12':f,'H13':g,'source_indices':{'R':4,'H':16},'Sylvester':Syl,'cyclic_resultant':Rc,'full_even_square_units':units,'checks':checks},str(p/'independent_critical_rankone_endpoint_certificate.sobj'),compress=False)
rep={'PASS':True,'checks':len(checks),'seconds':time.time()-st,'Sylvester_shape':[Syl.nrows(),Syl.ncols()],'cyclic_resultant':str(Rc.factor()),'protected_material_denominators':sorted(seen),'no_F12_separation_assumed':'complete even first forces e=0, original R40/H41 vanish, then delta unit gives literal F12/H13'};(p/'independent_critical_rankone_endpoint_certificate.json').write_text(json.dumps(rep,indent=2));print(rep,flush=True)
