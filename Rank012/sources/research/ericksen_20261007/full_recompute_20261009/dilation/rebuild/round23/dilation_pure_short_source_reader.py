"""Whole positive-simple pure D/D: short original-source certificate."""
from sage.all import *
from pathlib import Path
import signal,time,json
signal.alarm(100);O=Path(__file__).absolute().parent;st=time.time()
d=load(str(O.parent/'round12/homothety_first_jet_reduction.sobj'))
R0=d['ring']; K=R0.base_ring();a,b,s,t=K.gens();kap=a*(b-1)/(a-b);la=kap+1
pp=R0.gens()[:9];h1,h2,h3=pp[0],pp[4],pp[8];aa=1/a;bb=1/b
E12=-(bb-1)*h1+(bb-aa)*h2;E13=(aa-1)*h1+(bb-aa)*h3;E23=-(aa-1)*h2+(bb-1)*h3
pi=matrix(R0,[[h1,-t*E12/(2*s*(bb-1)),E13/(2*s*(aa-1))],[-s*E12/(2*t*(aa-bb)),h2,-E23/(2*t*(aa-1))],[-s*E13/(2*(aa-bb)),t*E23/(2*(bb-1)),h3]])
sub=dict(zip(pp,pi.list()));n=vector(R0,[s,t,1]);eye=matrix.identity(R0,3);Lam=diagonal_matrix(K,[1,aa,bb])
checks={'source_zero_rotation':(pi.transpose()-pi.trace()*eye/2)*n==0,
 'target_zero_rotation':(pi.transpose()*Lam-(Lam*pi).trace()*eye/2)*n==0}
R=PolynomialRing(K,names=('H','V','Z'));H,V,Z=R.gens()
hom=R0.hom([H,0,0,0,V,0,0,0,Z]+[R(0)]*8,R)
Qs=[hom(q.subs(sub)) for q in d['residuals'] if q]
lift=vector(K,[1,0,-kap/la])
ell=kap*H+V-kap*Z/la
checks['tracezero_direct_source_lift']=sum(c*q for c,q in zip(lift,Qs))==ell
print('trace lift',[str(c.factor()) if c else '0' for c in lift],flush=True)
gr=[]
for row in d['scaled_pi_fields']:
 g=kap*row[0]+row[4]-kap*row[8]/la
 gr.append(hom(g.subs(sub)).subs({V:-kap*H+kap*Z/la}))
short=[H*(la**2*t*t*H+(kap-la*t*t)*Z),(-kap*H+kap*Z/la)*(Z-la**2*s*s*H),Z*(t*t*Z-la*(kap*s*s+t*t)*H)]
units=[]
for i,(g,sh) in enumerate(zip(gr,short)):
 mon=sh.monomials()[0];un=g.monomial_coefficient(mon)/sh.monomial_coefficient(mon);units.append(un)
 checks['trace_gradient_'+str(i)]=g==un*sh
 print('gradient unit',i,un.factor(),flush=True)
# Pointwise main-contact consequences: Z=la^2*s^2*H and t^2=kap*s^2/(la*s^2-1).
K0=FractionField(PolynomialRing(QQ,names=('a','b','s')));A,B,S=K0.gens();ka=A*(B-1)/(A-B);ll=ka+1
RT=PolynomialRing(K0,'t');tt=RT.gen();rel=tt**2-ka*S*S/(ll*S*S-1)
def conv(p):return RT(sum(K0(cf)*A**ex[0]*B**ex[1]*S**ex[2]*tt**ex[3] for ex,cf in p.dict().items()))
def red(c):return (conv(c.numerator())*(conv(c.denominator())%rel).inverse_mod(rel))%rel
qq=[q.subs({V:-kap*H+kap*Z/la,Z:la**2*s*s*H}) for q in Qs]
# Simultaneous substitution above leaves Z in V image on some backends; explicitly remove it.
qq=[q.subs({Z:la**2*s*s*H}) for q in qq]
z=ll*S*S;P=z*z-z+1
target=[(2*ka*ll*S*S, 2*ka*tt/S*(z*z-1)*P),
        (ll*tt,ka*ll*S*(2*z-1)/(z-1)*P),
        (2*ll*ll*S*S,2*ll*tt/S*(z*z-1)*P)]
coeff=[]
for i,q in enumerate(qq):
 cs=(red(q.monomial_coefficient(H)),red(q.monomial_coefficient(H*H)));coeff.append(cs)
 checks['whole_Q'+str(i)+'_linear']=cs[0]==target[i][0]
 checks['whole_Q'+str(i)+'_quadratic']=cs[1]==target[i][1]
det=(target[0][0]*target[1][1]-target[1][0]*target[0][1])%rel
positive=2*ka**2*ll*S*P**2/(z-1)
checks['protected_final_minor']=det==positive
checks['positive_material_lambda']=ll==B*(A-1)/(A-B)
checks['positive_quadratic_identity']=P==(z-QQ(1)/2)**2+QQ(3)/4
assert all(checks.values()),checks
rep={'status':'PASS_FULL_PARAMETER_PURE_D_D_MAIN_EXCLUSION','checks':checks,'count':len(checks),
 'trace_lift':[str(c.factor()) if c else '0' for c in lift],
 'gradient_units':[str(c.factor()) for c in units],
 'final_minor':'2*kappa^2*lambda*s*(z^2-z+1)^2/(z-1), z=lambda*s^2',
 'protected_region':'0<b<a<1; s*t != 0; h1*h2*h3 != 0 for the contradiction only',
 'seconds':time.time()-st}
save({'original_ring':R0,'pi':pi,'ring':R,'source_residuals':Qs,'trace':ell,'trace_lift':lift,
 'gradients':gr,'gradient_units':units,'short_gradients':short,'reduced_coefficients':coeff,'positive_minor':positive,'checks':checks},str(O/'dilation_pure_short_source_certificate.sobj'))
(O/'dilation_pure_short_source_reader.json').write_text(json.dumps(rep,indent=2),encoding='utf8');print(json.dumps(rep,indent=2),flush=True)
