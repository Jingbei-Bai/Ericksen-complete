"""Exact displayed Laurent coefficients, material signs, and critical-line input bridge.

The two archived ray coefficients are comparison targets only. The parent
ordinary audit independently reconstructs the same values from fresh raw
third rows. All identities below retain both material parameters.
"""
from sage.all import *
from pathlib import Path
import json,hashlib
out=Path(__file__).absolute().parent; research=out.parent.parent
checks=[];imports=[]
def ck(n,v):
 assert bool(v),n
 checks.append(n)
def old(n):
 p=research/'round10'/(n+'.sobj');imports.append(dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
 return load(str(p))
R=PolynomialRing(QQ,names=('A','B'));A,B=R.gens();K=R.fraction_field()
cm=2352*A**3*B**2+961*A**2*B**3-672*A**3*B-3970*A**2*B**2-1984*A*B**3+48*A**3+577*A**2*B+1664*A*B**2+1024*B**3
cp=1200*A**3*B**2+1089*A**2*B**3+480*A**3*B-2946*A**2*B**2-2112*A*B**3+48*A**3-575*A**2*B+1792*A*B**2+1024*B**3
minus=K(QQ(9216)/49*A**2*(B-1)**4*(B-A)**14*(7*A*B-A-6*B)**2*(11*A*B-3*A-8*B)**4*cm)/((A-1)**18*B**17*(5*A*B-A-4*B)**6*(2*A*B-A-B)**2)
plus=K(QQ(9216)/49*A**2*(B-1)**4*(B-A)**12*(5*A*B+A-6*B)**2*(5*A*B+3*A-8*B)**4*cp)/((A-1)**18*B**17*(3*A*B+A-4*B)**6)
for name,f in [('independent_incoming_third_base_infinity_probe',minus),('independent_incoming_third_positive_ray_probe',plus)]:
 d=old(name);v=d['leading_coefficient'];h=v.parent().ring().hom([A,B],K)
 ck('printed coefficient equals exact ray comparison '+name,f==h(v.numerator())/h(v.denominator()))
S=PolynomialRing(QQ,names=('s','t'));s,t=S.gens();sf=S.fraction_field();a=1+s;b=1+s+t
positive0=3313*s**2*(s+1)**3+s*t*(s+1)**2*(7587*s+610)+t**2*(s+1)*(5235*s**2+548*s+49)+(31*s-1)**2*t**3
positive2=1024*s**5+2048*s**4*t+1024*s**3*t**2+4800*s**4+9536*s**3*t+4736*s**2*t**2+6528*s**3+13504*s**2*t+6977*s*t**2+2752*s**2+6016*s*t+3313*t**2
positive1=48*s**5+144*s**4*t+144*s**3*t**2+48*s**2*t**3+145*s**4+866*s**3*t+1297*s**2*t**2+576*s*t**3+146*s**3+1234*s**2*t+2816*s*t**2+1728*t**3+49*s**2+512*s*t+2752*t**2
ck('negative identity positive decomposition',cm(a,b)==positive0)
ck('negative cycle2 positive decomposition',b**5*cm(1/sf(b),sf(a)/b)==positive2)
ck('positive cycle1 positive decomposition',a**5*cp(sf(b)/a,1/sf(a))==positive1)
ck('cycle2 all positive coefficients',all(c>0 for c in positive2.coefficients()))
ck('cycle1 all positive coefficients',all(c>0 for c in positive1.coefficients()))
for n,c in [(2,1),(5,1),(7,1),(11,3),(53,13),(65,17)]:
 f=n*A*B-c*A-(n-c)*B
 ck('negative identity factor '+str((n,c)),f(a,b)==n*s*(s+1+t)+c*t)
 ck('negative cycle2 factor '+str((n,c)),b*b*f(1/sf(b),sf(a)/b)==-(n-c)*a*(b-1)-c*(b-a))
for n,c in [(3,1),(5,1),(5,3),(27,13),(31,17)]:
 f=n*A*B+c*A-(n+c)*B
 ck('positive cycle1 factor '+str((n,c)),a*a*f(sf(b)/a,1/sf(a))==c*a*(b-1)+n*(b-a))
# Material interpretation of the omega graph: theta may be complex, but x=theta²
# is a real material expression, and no sign of theta is used.
V=PolynomialRing(QQ,names=('a','r'));aa,rr=V.gens();VF=V.fraction_field();q=1+rr;m=1+aa*rr
xx=rr;ck('one omega u=ar squared',aa*xx**2==aa*rr**2)
xx=-q;ck('omega one u=mq',xx*(aa*(xx+1)-1)==m*q)
xx=-VF(q)/rr;ck('omega omega negative product',aa*xx*(aa-xx-1)==-aa*m*q/rr**2)
ck('omega omega exceptional factor',aa-xx==VF(aa*rr+q)/rr)
bb=aa*q/m;ck('omega omega wall is b=a²',bb-aa**2==-aa*(aa-1)*(aa*rr+q)/m)
ck('omega c exact lower bound',232135214-27335384*QQ(13)/2==54455218)
# The critical-line original forms are source-validated by the current even-line
# reader. Rebuild both substitution and flattening, including removed units.
d=old('independent_subcritical_line_material_free');saved=old('independent_delta_critical_line_conics');flat=old('independent_delta_critical_line_native_resultants')
T=saved['ring'];tt=T.gen();F=T.fraction_field();n,la=T.base_ring().gens()
ph=d['material_free_ring'].hom([F(tt),-F(n)/(la*(tt+1)**2),F(n*tt**2)/(tt+1)**2,F.one()],F)
polys=[]
for j,f in enumerate(d['material_free_forms']):
 p=T(ph(f)*la**2*(tt+1)**4);ck('source critical line polynomial '+str(j),p==saved['polynomials'][j]);polys.append(p)
RF=flat['ring'];x,N,L=RF.gens()
for j,(p,unit) in enumerate(zip(polys,[tt+1,tt,T.one()])):
 f,rem=p.quo_rem(unit);ck('critical line removed support unit '+str(j),rem==0)
 fresh=RF.zero()
 for e,c in f.dict().items():
  ck('critical line coefficient integral '+str((j,e)),c.denominator()==1)
  for ex,co in c.numerator().dict().items():fresh+=co*x**int(e)*N**int(ex[0])*L**int(ex[1])
 ck('critical line exact flat resultant input '+str(j),fresh==flat['flat_polynomials'][j])
report=dict(status='PASS',count=len(checks),checks=checks,imports=imports,
 scope='Displayed two Laurent coefficients and every printed positivity identity; omega material interpretation; complete source-to-critical-line/resultant-polynomial bridge. Ray comparison inputs are not new source evidence; root independently rebuilds raw third sources.')
(out/'extra_tau_and_critical_line.json').write_text(json.dumps(report,indent=2))
print('PASS',len(checks))
