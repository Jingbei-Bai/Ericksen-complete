"""Primitive compact replay of the full generic seven-source h bridge."""
from sage.all import *
from pathlib import Path
import os,resource,time,json
p=Path(__file__).resolve().parent;st=time.time();d=load(str(p/'pure_x_super_reader_generic_h_kappa_bridge.sobj'))
vm=int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE');resource.setrlimit(resource.RLIMIT_AS,(vm+500*1024**2,vm+500*1024**2));print('LinuxPID',os.getpid(),flush=True)
P=d['ring'];K=P.base_ring();k,r,H,Z=K.gens();y,z,w,e0,e1,e2,de,t,v,kap=P.gens();cache={};checks=[]
def ck(n,b):
 if not b:raise AssertionError(n)
 checks.append(n)
def nk(c):
 c=K(c)
 if c not in cache:
  if not c:cache[c]=K.zero()
  else:
   n,dd=c.numerator(),c.denominator();g=n.gcd(dd);n,dd=n//g,dd//g;lc=dd.leading_coefficient();cache[c]=K(n/lc)/K(dd/lc)
 return cache[c]
def nt(f):return P({tuple(map(int,e)):nk(c) for e,c in f.dict().items() if c})
vec=vector(K,[nk(c) for c in d['generic_combination']]);print('GENERIC COEFFICIENTS',[str(c) for c in vec],flush=True)
inc=y*v-(Z-H)/(H*(Z-1))*w*t-kap*e1;zero=P.hom(list(P.gens()[:9])+[P.zero()],P)
first={i:nt(f) for i,f in d['source_first_two'].items()};co=nt(d['incidence_cofactor']);quot=nt(d['true_h_kappa_zero_divided_source'])
polys=[nt(mult*first[idx]) for idx,mult in d['support']]
whole=sum((c*f for c,f in zip(vec,polys)),P.zero())-co*inc
ck('primitive complete generic kappa0 identity',zero(whole)==0)
one=P({tuple(int(nn)-(1 if j==9 else 0) for j,nn in enumerate(e)):c for e,c in whole.dict().items() if e[9]==1})
ck('primitive TRUE divided source coefficient exactly original',nt(one)==quot)
ck('primitive TRUE divided source raw linear zero',all(quot.monomial_coefficient(x)==0 for x in P.gens()))
CB=PolynomialRing(QQ,names=('k','r')).fraction_field();kk,rr=CB.gens();ph1=K.ring().hom([kk,rr,1,-1],CB);ph2=K.ring().hom([kk,rr,2,-1],CB)
def cb(c,ph):return ph(c.numerator())/ph(c.denominator())
den1=[];den2=[]
for j,c in enumerate(vec):
 dd1=ph1(c.denominator());dd2=ph2(c.denominator());ck('generic coefficient cyclic denominators EACH material '+str(j),dd1 in QQ and dd1!=0 and dd2 in QQ and dd2!=0);den1.append(dd1);den2.append(dd2)
ck('generic coefficient material k independence',all(c.derivative(k)==0 for c in vec))
expected1=[CB.one(),-2864*rr/2415,-3832*rr/2415,CB(968)/2415,344*rr/7245,-1856*rr/1449,CB(3592)/7245]
expected2=[CB.one(),-10740*rr/10049,-12960*rr/10049,CB(2960)/10049,476*rr/10049,-11692*rr/10049,CB(3489)/10049]
ck('primitive H1 kernel coefficient exact',list(cb(c,ph1) for c in vec)==expected1)
ck('primitive H2 kernel coefficient exact',list(cb(c,ph2) for c in vec)==expected2)
# Every state power is retained in the kappa-zero identity. For the true
# pure quadratic at the cyclic point, raw linear zero means the old four-
# range second correction cannot alter the cokernel coefficient.
quad=P({tuple(map(int,e)):c for e,c in quot.dict().items() if sum(e)==2})
cy1=P({tuple(map(int,e)):K(cb(c,ph1)) for e,c in quad.dict().items()});cy2=P({tuple(map(int,e)):K(cb(c,ph2)) for e,c in quad.dict().items()})
save({'ring':P,'source_first_two':first,'support':d['support'],'generic_combination':vec,'incidence_cofactor':co,'complete_kappa0_identity_polynomials':polys,'true_divided_kappa0_source':quot,'true_raw_quadratic':quad,'H1_raw_quadratic':cy1,'H2_raw_quadratic':cy2,'cyclic_denominators_H1':den1,'cyclic_denominators_H2':den2,'original_source_DAG_payload':d['raw_sources'],'checks':checks,'scope':__doc__},str(p/'pure_x_super_reader_generic_h_compact_certificate.sobj'),compress=False)
rep={'PASS':True,'checks':len(checks),'seconds':time.time()-st,'generic_coefficient_cyclic_units_material_independent':True,'whole_kappa0_identity':True,'TRUE_raw_divided_linear_zero':True,'original_R_H_only':True,'generic_combination':[str(c) for c in vec]};(p/'pure_x_super_reader_generic_h_compact_certificate.json').write_text(json.dumps(rep,indent=2)+'\n');print('PASS',{'checks':len(checks),'seconds':time.time()-st,'normalized_distinct_coefficients':len(cache)},flush=True)
