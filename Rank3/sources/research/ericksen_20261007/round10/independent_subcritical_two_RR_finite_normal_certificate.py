"""Explicit exact source certificates for the G6 and J19 normal-factor layers.

Scope: representative missing-RR3, t*w*s != 0.  D=tw+ws+s is retained
as a separate whole curve; no claim about other support transports here.
All material parameters are retained via N=a*r^2 in (0,1).
"""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent; start=time.monotonic();checks=[]
def ck(label,v):
 assert v,label
 checks.append(label)
d=load(str(p/'independent_subcritical_two_RR_normal_factor_constants.sobj'))
R=d['ring'];x,W,y=R.gens();items=[]
for idx,it in enumerate(d['items']):
 fs=list(it['HZ_coefficients'].values());I=R.ideal(fs);cert=[]
 for j,g in enumerate(it['basis']):
  vv=list(g.lift(I));ck('basis_source_%s_%s'%(idx,j),sum((a*b for a,b in zip(vv,fs)),R.zero())==g)
  cert.append({'target':g,'cofactors':vv})
 items.append({'sources':fs,'basis_certificates':cert})
# The two safe pivots can be read from the actual three source coefficients.
gZ=x*W+2*x*y-2*W-y;g0=-4*x*W-5*x+2*W+4;gH=-2*x*y-4*x+4*y+5
ck('G_coefficients',set(items[0]['sources'])==set([gZ,g0,gH]))
ck('G_pivot_x_two',gH.subs({x:2})==-3)
ck('G_pivot_x_half',g0.subs({x:QQ(1)/2})==QQ(3)/2)
K=R.fraction_field();wx=(5*x-4)/(2-4*x);yx=(4*x-5)/(4-2*x)
sub=R.hom([K(x),K(wx),K(yx)],K)
ck('G_solve_g0',sub(g0)==0);ck('G_solve_gH',sub(gH)==0)
p3=(x-1)*(7*x*x-13*x+7);ratio=K(sub(gZ)/p3)
ck('G_cubic_source',sub(gZ)==ratio*p3)
print('G_RATIO',ratio,flush=True)
# J19's simple reduced set follows without dividing y-1 globally.
J=R.ideal(items[1]['sources'])
targets=[(y-1)*(y*y-QQ(13)/2*y+1),(W-1)*(W-y),(x-1)*(W-1),x*(y-1)-(W-1)]
jcert=[]
for j,g in enumerate(targets):
 vv=list(g.lift(J));ck('J_classification_source_%s'%j,sum((a*b for a,b in zip(vv,items[1]['sources'])),R.zero())==g)
 jcert.append({'target':g,'cofactors':vv})
# Both omega conjugates are positive.  Retain omega algebraically and check
# the original compatibility source against a manifestly positive formula.
old=load(str(p/'independent_subcritical_two_RR_delta_ratio_coordinates.sobj'))
B=PolynomialRing(QQ,names=('positive_N','omega'));N,om=B.gens();BK=B.fraction_field()
f=old['ring'].hom([BK(N),1/N,BK(om),BK(om)],BK)(old['sources'][-1])
A=om*om+om+1;positive=A*(N-1)**2+2*N*om*(2*om+1)
expected=A*(N*N+1)+2*(om*om-1)*N
ck('J_positive_identity',positive==expected)
fr=BK(f/expected)
print('J_COMPATIBILITY_RATIO',fr,flush=True)
# The ratio is an explicit monomial unit on this branch, no unrecorded
# denominator or omega-dependent zero may be removed.
ck('J_omega_source',f==fr*positive)
ck('J_ratio_monomials',len(fr.numerator().dict())==1 and len(fr.denominator().dict())==1)
ck('J_omega_minpoly_positive_roots',QQ(13)/2>0 and 1>0 and (QQ(13)/2)**2-4>0)
# Verify the earlier full integer Sylvester and point-component certificate.
finite=load(str(p/'independent_subcritical_two_RR_finite_normal_layers.sobj'))
ck('G_resultant_factor_product',finite['resultant_factorization'].value()==finite['resultant'])
ck('J_point_unit',sum((c*f for c,f in zip(finite['J_point_unit_cofactors'],finite['J_point_sources'])),finite['J_point_ring'].zero())==1)
C=PolynomialRing(QQ,'norm_N');nn=C.gen()
ck('G_positive_quadratic',133*nn*nn-218*nn+133==133*(nn-1)**2+48*nn)
payload={'source_file':'independent_subcritical_two_RR_normal_factor_constants.sobj','ring':R,'items':items,'G_cubic':p3,'G_cubic_ratio':ratio,'G_graph':(wx,yx),'J_classification_certificates':jcert,'J_omega_ring':B,'J_omega_source':f,'J_omega_unit':fr,'J_positive_polynomial':positive,'J_omega_minpoly':om*om-QQ(13)/2*om+1,'finite_source_file':'independent_subcritical_two_RR_finite_normal_layers.sobj','checks':checks,'scope':__doc__}
save(payload,str(p/'independent_subcritical_two_RR_finite_normal_certificate.sobj'),compress=False)
(p/'independent_subcritical_two_RR_finite_normal_certificate.json').write_text(json.dumps({'status':'PASS','checks':len(checks),'seconds':time.monotonic()-start,'scope':__doc__,'G_cubic_ratio':str(ratio),'J_omega_unit':str(fr)},indent=2)+'\n')
print('PASS',len(checks),'SECONDS',time.monotonic()-start,flush=True)
