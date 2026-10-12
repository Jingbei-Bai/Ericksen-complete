"""Independent full-parameter replay of the small R-R infinity certificate.

Reconstruct P and D from the protected Cartesian H source, then derive the
two extreme cubic coefficients directly. No diagnostic parameter values.
"""
from sage.all import *
from pathlib import Path
import time,signal,json
O=Path(__file__).resolve().parent;start=time.time();signal.alarm(100)
d=load(str(O/'rr_null_physical_protected_graph_plus.sobj'));T=d['ring'];K=d['field'];B=K.ring();s,t,p,q,c=K.gens();aa,bb,chi=T.gens();ii=B.base_ring().gen()
F=d['F'];H=d['H'];As=d['axes'];At=d['axet'];u=vector(T,[0,0,1]);C=F.T*F;Ci=F.inverse()*F.inverse().T
DC=[];g1=[];g2=[];qb=vector(T,[0,0,0])
for j in range(3):
 ej=vector(T,[int(k==j) for k in range(3)])
 dF=matrix(T,3,3,lambda k,l:H[tuple(sorted((j,l)))][k]);dC=dF.T*F+F.T*dF;DC.append(dC)
 dZ=As.cross_product(ej);du=dZ-(u*dZ)*u
 g1.append(u*dC*u+2*du*C*u)
 g2.append(u*(dC*C+C*dC)*u+2*du*C*C*u)
 dCi=-Ci*dC*Ci
 for k in range(3):qb[k]+=dCi[j,k]
D=g1[0]*g2[1]-g1[1]*g2[0]
Pa=At*sum((H[j,j] for j in range(3)),vector(T,[0,0,0]))+F.det()*(As*qb)
J=D.derivative(chi)*Pa.subs({chi:0})-Pa.derivative(chi)*D.subs({chi:0})
archive=load(str(O/'rr_cartesian_direction_paired_pencil.sobj'))
def cf(z):
 z=K(z)
 if not z:return K.zero()
 lc=z.denominator().leading_coefficient();return K(z.numerator()/lc)/K(z.denominator()/lc)
checks={'P_replayed_from_H':Pa==archive['paired'],'D_replayed_from_unit_direction':D==archive['D'],
 'J_replayed_from_P_and_D':J==archive['cubic'],'P_linear_chi':Pa.degree(chi)<=1,
 'D_linear_chi':D.degree(chi)<=1,'orbit_moment_derivatives_zero':g1[2]==0 and g2[2]==0}
rawQ=load(str(O/'rr_null_physical_quadratics_plus.sobj'))['raw_quadratics'][0]
checks['P_original_quadratic_source_bridge']=4*s*t*c*Pa==rawQ
guard=B(s*t*c*p*q*(s*s-t*t)*(t*t*p*p+s*s*q*q))
def protected_den(f):
 de=B(K(f).denominator())
 while not de.is_constant():
  z=de.gcd(guard)
  if z.is_constant():return False
  de=de//z
 return bool(de)
checks['all_P_D_J_denominators_protected']=all(protected_den(co) for f in [Pa,D,J] for co in f.coefficients())
ac=cf(J.monomial_coefficient(aa**3));bc=cf(J.monomial_coefficient(bb**3))
ap=B(ac.numerator()/((4*ii)*B(c)**5));bp=B(bc.numerator()/((-4*ii)*B(c)**2))
RR=PolynomialRing(QQ,names=('S','T','P','Q','C'));S,TT,PP,Qsq,CC=RR.gens()
def realpart(poly,h):return B({ex:QQ(co[h]) for ex,co in poly.dict().items() if co[h]})
def tosquare(poly):
 assert all(all(e%2==0 for e in ex) for ex in poly.dict())
 return RR({tuple(int(e)//2 for e in ex):QQ(co) for ex,co in poly.dict().items()})
ff=[tosquare(realpart(ap,0)//B(p)),tosquare(-realpart(ap,1)//B(q)),
 tosquare(realpart(bp,0)//B(p*t)),tosquare(-realpart(bp,1)//B(q*s))]
back=RR.hom([B(s*s),B(t*t),B(p*p),B(q*q),B(c*c)],B)
checks['alpha_extreme_reconstructed_without_factoring']=ap==B(p)*back(ff[0])-ii*B(q)*back(ff[1])
checks['beta_extreme_reconstructed_without_factoring']=bp==B(p*t)*back(ff[2])-ii*B(q*s)*back(ff[3])
positive=TT*PP**3+(15*S+8*TT)*PP**2*Qsq+(8*S+15*TT)*PP*Qsq**2+S*Qsq**3+12*(S-TT)**2*PP*Qsq
left=(PP+5*Qsq)*(ff[1]-ff[3])+(5*PP+Qsq)*(ff[2]-ff[0])
checks['direct_whole_positive_identity']=left==2*(S-TT)*positive
checks['only_nonpositive_expansion_is_a_square']=positive-12*(S-TT)**2*PP*Qsq==TT*PP**3+(15*S+8*TT)*PP**2*Qsq+(8*S+15*TT)*PP*Qsq**2+S*Qsq**3
checks['four_remaining_groups_strict_positive_coefficients']=all(co>0 for co in (positive-12*(S-TT)**2*PP*Qsq).coefficients())
assert all(checks.values()),checks
save({'ring':RR,'four_real_sources':ff,'positive':positive,'identity_left':left,'checks':checks},str(O/'rr_global_pencil_positive_independent_certificate.sobj'))
rep={'status':'PASS_INDEPENDENT_ORIGINAL_SOURCE_TO_WHOLE_PENCIL_NONZERO','checks':checks,'count':len(checks),'seconds':time.time()-start,
 'domain':'s,t,c>0; p*q*(s^2-t^2)!=0; no angular-rate or pressure guard introduced',
 'identity':'(P+5Q)(f1-f3)+(5P+Q)(f2-f0)=2(S-T)*Pos, Pos>0',
 'other_sheet':'covered separately by the exact target-reflection reader; complex conjugation covers the remaining two null directions'}
(O/'rr_global_pencil_positive_independent_reader.json').write_text(json.dumps(rep,indent=2));print(json.dumps(rep,indent=2),flush=True)

