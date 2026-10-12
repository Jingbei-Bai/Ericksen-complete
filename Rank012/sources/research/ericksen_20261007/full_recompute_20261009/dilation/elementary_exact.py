"""Fresh identities independent of archived state matrices or PASS records.

Every check is over a polynomial/rational function field. No material
specialization supplies a positive assertion.
"""
from pathlib import Path
import sympy as S
import json, hashlib, time
HERE=Path(__file__).absolute().parent
started=time.time(); checks=[]
def check(name, expression):
    print(name,flush=True)
    residual=S.factor(expression)
    ok=residual==0
    checks.append({'name':name,'passed':ok,'residual':str(residual)})
    if not ok: raise AssertionError((name,residual))

l1,l2,l3=S.symbols('l1 l2 l3', nonzero=True)
V=S.Matrix([[1,1,1],[l1,l2,l3],[1/l1,1/l2,1/l3]])
check('sim:projector_inverse_pressure_Vandermonde',V.det()+(l1-l2)*(l1-l3)*(l2-l3)/(l1*l2*l3))
sigma,w=S.symbols('sigma w',real=True)
Omega=S.Matrix([[0,-w,0],[w,0,0],[0,0,0]])
check('sim:unique_homothety_centre',(sigma*S.eye(3)+Omega).det()-sigma*(sigma*sigma+w*w))
# The diagonal entries are distinct: the three skew parameters vanish.
u,v,q=S.symbols('u v q')
A=S.Matrix([[0,-q,v],[q,0,-u],[-v,u,0]])
comm=S.diag(l1,l2,l3)*A-A*S.diag(l1,l2,l3)
for (i,j),ans in { (0,1):(l2-l1)*q,(0,2):(l1-l3)*v,(1,2):(l3-l2)*u}.items():
    check('sim:affine_simple_spectrum_centralizer_%d%d'%(i,j),comm[i,j]-ans)

p=S.Matrix(3,3,S.symbols('p0:9'))
pressure=S.Matrix([[p[1,2]-p[2,1],p[0,2],-p[0,1]],[-p[1,2],p[2,0]-p[0,2],p[1,0]],[p[2,1],-p[2,0],p[0,1]-p[1,0]]])
for j in range(3):check('dil:sum_projector_pressures_'+str(j),sum(pressure[:,j]))
def eps(i,j,k):
    return 0 if len({i,j,k})<3 else (1 if (i,j,k) in [(0,1,2),(1,2,0),(2,0,1)] else -1)
for i in range(3):
    for j,k in [(0,1),(0,2),(1,2)]:
        # d(exp(phi_i)theta_i)/exp(phi_i), in the theta basis.
        coefficient=sum(eps(j,k,m)*p[i,m] for m in range(3))
        coefficient += pressure[i,j]*int(k==i)-pressure[i,k]*int(j==i)
        check('dil:contact_structure_%d_%d%d'%(i,j,k),coefficient-eps(j,k,i)*p[i,i])

X,Y,Z=S.symbols('X Y Z')
F=(X+Y)**2-2*X*Z+2*Y*Z-3*Z**2
check('dil:nonaxial_XplusZ_face',F.subs(Z,-X)-Y**2)
check('dil:nonaxial_minusXplusYminusZ_face',F.subs(Z,-X+Y)-4*X*Y)
check('dil:nonaxial_XplusYplusZ_face',F.subs(Z,-X-Y)-4*Y*(-X-Y))

a,b,s,t,H=S.symbols('a b s t H',nonzero=True)
check('dil:pure_material_relation',a*(b-1)/(a-b)+1-b*(a-1)/(a-b))
kap=S.symbols('kappa',nonzero=True);ell=kap+1;z=ell*s*s;P=z*z-z+1
A0=2*kap*ell*s*s;B0=2*kap*t/s*(z*z-1)*P
A1=ell*t;B1=kap*ell*s*(2*z-1)*P/(z-1)
det=S.expand(A0*B1-A1*B0)
check('dil:pure_final_minor',det.subs(t*t,kap*s*s/(ell*s*s-1))-2*kap**2*ell*s*P**2/(z-1))
check('dil:pure_positive_quadratic',P-(z-S.Rational(1,2))**2-S.Rational(3,4))
zeta=S.symbols('zeta')
hermite=1+S.Rational(5,3)*zeta-zeta*zeta/12-S.Rational(7,6)*zeta**3-S.Rational(5,12)*zeta**4
for pt,val in [(0,1),(1,1),(-1,0),(-2,0)]:check('dil:Hermite_value_'+str(pt),hermite.subs(zeta,pt)-val)
check('dil:Hermite_double_negative_block',S.diff(hermite,zeta).subs(zeta,-1))
ann=zeta*(zeta-1)*(zeta+2)*(zeta+1)**2
check('dil:Hermite_idempotence_mod_annihilator',S.rem(hermite**2-hermite,ann,zeta))
check('dil:Hermite_positive_semisimplicity_mod_annihilator',S.rem(zeta*(zeta-1)*hermite,ann,zeta))

r,theta,zz,alpha,k,beta,c,ws,wt=S.symbols('r theta zz alpha k beta c ws wt',real=True)
angle=k*theta+beta*S.log(r)
map_=S.Matrix([alpha*r*S.cos(angle),alpha*r*S.sin(angle),c*zz])
Zmap=r*map_.diff(r)+ws*map_.diff(theta)+zz*map_.diff(zz)
target=map_+wt*S.Matrix([-map_[1],map_[0],0])
for i in range(3):check('dil:F5_exact_intertwining_'+str(i),(Zmap-target)[i].subs(wt,beta+k*ws))
F5=S.Matrix([[alpha,0,0],[alpha*beta,alpha*k,0],[0,0,c]])
C5=F5.T*F5
check('dil:F5_signed_J',F5.det()-alpha**2*k*c)
check('dil:F5_transverse_trace',C5[:2,:2].trace()-alpha**2*(1+beta**2+k*k))
check('dil:F5_transverse_determinant',C5[:2,:2].det()-alpha**4*k*k)
uu,vv,ll,sk=S.symbols('u v ell s_k',nonzero=True)
twist=(ll-uu/vv)*(vv/uu-ll)
check('dil:F5_full_spectrum_trace',(uu*vv/ll)*(1+twist+ll*ll)-(uu*uu+vv*vv))
quadratic=(1+ws*ws)*ll*ll-(uu/vv+vv/uu+2*sk*ws*wt)*ll+(1+wt*wt)
check('dil:F5_rate_quadratic',S.rem(S.expand((wt-sk*ll*ws)**2-twist-quadratic),sk*sk-1,sk))

report={'status':'PASS_EXACT_SYMBOLIC_IDENTITIES','checks_count':len(checks),'checks':checks,
        'scope':'Independent rational identities; state-derived blocks are checked separately.',
        'numerical_sampling_used':False,'seconds':time.time()-started,
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(HERE/'elementary_exact.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps({k:report[k] for k in ['status','checks_count','seconds']},indent=2))
