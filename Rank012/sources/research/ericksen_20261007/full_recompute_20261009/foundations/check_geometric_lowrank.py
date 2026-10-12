"""Exact independent calculations for Sections 2--4 and 9.
No sample spectrum is used. Run with ordinary Python + SymPy.
"""
from pathlib import Path
import json,time
import sympy as s
OUT=Path(__file__).absolute().parent; start=time.monotonic(); checks=[]
def check(name,expr):
    if isinstance(expr,s.MatrixBase): residues=list(expr.applyfunc(s.simplify))
    elif isinstance(expr,(list,tuple)):residues=[s.simplify(e) for e in expr]
    else:residues=[s.simplify(expr)]
    ok=all(e==0 for e in residues)
    checks.append(dict(name=name,passed=ok,scalar_identities=len(residues),nonzero_residuals=[str(e) for e in residues if e!=0]))
    print(name,ok,flush=True)
    (OUT/'geometric_lowrank_checks.json').write_text(json.dumps(dict(status='RUNNING',checks=checks),indent=2))
    assert ok,(name,residues)
a,b=s.symbols('a b',nonzero=True)
V=s.Matrix([[1,1,1],[1,a,b],[1,1/a,1/b]])
check('principal-projector Vandermonde',V.det()+(a-1)*(a-b)*(b-1)/(a*b))
p=s.Matrix(3,3,lambda i,j:s.Symbol('p%d%d'%(i+1,j+1)))
w=s.symbols('w1:4',nonzero=True)
def C(k,i,j):return -sum(s.LeviCivita(i,j,l)*p[k,l] for l in range(3))
def G(k,i,j,weights):return (weights[k]*C(k,i,j)-weights[i]*C(i,j,k)+weights[j]*C(j,k,i))/(2*weights[k])
check('Koszul torsion and metric compatibility',[G(k,i,j,w)-G(k,j,i,w)-C(k,i,j) for i in range(3) for j in range(3) for k in range(3)]+[w[k]*G(k,i,j,w)+w[j]*G(j,i,k,w) for i in range(3) for j in range(3) for k in range(3)])
sym={p[1,0]:p[0,1],p[2,0]:p[0,2],p[2,1]:p[1,2]};N=p.subs(sym)
def R(l,k,i,j):return sum(G(m,j,k,[1]*3)*G(l,i,m,[1]*3)-G(m,i,k,[1]*3)*G(l,j,m,[1]*3)-C(m,i,j)*G(l,m,k,[1]*3) for m in range(3))
Ric=s.Matrix(3,3,lambda k,j:sum(R(i,k,i,j) for i in range(3))).subs(sym)
expected=2*N*N-s.trace(N)*N+(s.trace(N)**2/2-s.trace(N*N))*s.eye(3)
check('constant symmetric frame Ricci tensor',Ric-expected)
r,z=s.symbols('r z',nonzero=True);betas=[1,a,b]
al=lambda i,j:1-betas[i]/betas[j]
R1=-al(1,0)/(al(1,2)*z);R2=-al(2,1)*r/al(2,0);R3=-al(0,2)*z/(al(0,1)*r)
check('opposite-pair product',R1*R2*R3-1)
check('cyclic alpha quotient',al(0,1)*al(1,2)*al(2,0)/(al(1,0)*al(2,1)*al(0,2))+1)
h1,h2,h3=s.symbols('h1 h2 h3');AA=a*(b-1);CC=b*(a-1);DD=a-b
q=lambda aa,bb:bb*(aa-1)/(aa*(bb-1))
check('second mu wall cyclic normalization',h1-q(b/a,1/a)*h3-(h1-DD*h3/CC))
check('third mu wall cyclic normalization',h2-q(1/b,a/b)*h1+AA/DD*(-h1-DD*h2/AA))
M=s.Matrix([[R1*(R1+1),-al(2,0)/al(1,0)*R2,-al(0,2)/al(1,2)*R1*R3],[-al(1,0)/al(2,0)*R1*R2,R2*(R2+1),-al(0,1)/al(2,1)*R3],[-al(1,2)/al(0,2)*R1,-al(2,1)/al(0,1)*R2*R3,R3*(R3+1)]])
M=M.applyfunc(s.cancel)
detM=sum(s.cancel(sign*M[0,i]*M[1,j]*M[2,k]) for i,j,k,sign in ((0,1,2,1),(1,2,0,1),(2,0,1,1),(2,1,0,-1),(1,0,2,-1),(0,2,1,-1)))
check('constant-ratio square matrix determinant',s.cancel(detM)-4)
sigma,u=s.symbols('sigma u',nonzero=True)
tr=sigma/(2*u);tr2=sigma**2/(2*u**2);prod=(tr**2-tr2)/2
check('eikonal two-eigenvalue product',prod+sigma**2/(8*u**2))
check('eikonal incompatible cube trace',tr**3-3*prod*tr-5*sigma**3/(16*u**3))
theta=s.symbols('theta',real=True);lam,m,mp,tp=s.symbols('lambda m mprime thetaprime',real=True)
rot=s.Matrix([[s.cos(theta),-s.sin(theta)],[s.sin(theta),s.cos(theta)]])
Gam=rot*s.diag(lam,m)*rot.T
Gamp=Gam.diff(theta)*tp+Gam.diff(m)*mp
check('translation orbit Gauss determinant',s.trigsimp(Gamp.det())+(lam-m)**2*tp**2)
rr=s.symbols('rho',positive=True);beta,gamma,alpha,delta=s.symbols('beta gamma alpha delta',real=True)
uu=beta+gamma/rr**2;vv=alpha*rr+delta/rr
check('cylindrical pressure integration',[s.diff(rr*(s.diff(uu,rr)+2*uu/rr),rr),s.diff(s.diff(vv,rr)+vv/rr,rr)])
a0,n11,n12,n22=s.symbols('a0 n11 n12 n22',real=True)
Gam=s.Matrix([[a0*rr**2+n11,n12],[n12,n22]])
check('cylindrical orbit determinant', (Gam-lam*s.diag(rr**2,1)).det()-((a0-lam)*(n22-lam)*rr**2+n11*(n22-lam)-n12**2))
lj,lk,kap=s.symbols('lambda_j lambda_k kappa',real=True)
res=(lj-lk)**2*(uu**2+(kap*rr*uu)**2)-(lj-n22-n11/rr**2)*(n22+n11/rr**2-lk)
check('cylindrical leading positive coefficient',s.expand(res).coeff(rr,2)-(lj-lk)**2*kap**2*beta**2)
check('cylindrical negative-four coefficient',s.expand(res.subs(beta,0)).coeff(rr,-4)-((lj-lk)**2*gamma**2+n11**2))
aa,bb,cc,dd,zz1,zz2,k=s.symbols('A11 A12 A21 A22 zeta1 zeta2 k',nonzero=True)
AA=s.Matrix([[aa,bb],[cc,dd]]);J=AA.det();zg=s.Matrix([[zz1,zz2]])
FF=AA.row_join(s.zeros(2,1)).col_join(s.Matrix([[zz1,zz2,k]]))
etanorm=((zg*AA.inv())*(zg*AA.inv()).T)[0]
check('nonabelian forward trace',s.trace(FF.T*FF)-(sum(e**2 for e in AA)+zz1**2+zz2**2+k**2))
check('nonabelian inverse trace',s.trace(FF.inv().T*FF.inv())-(sum(e**2 for e in AA)/J**2+(1+etanorm)/k**2))
ss,c0,d0,E,JJ=s.symbols('s c0 d0 E J',nonzero=True);av,bv=s.symbols('a0 b0',real=True)
tt=c0*ss+d0;bv0=JJ*tt/ss
PP=c0*E-d0-1-JJ**2*c0**2-c0*ss+(d0*E-2*JJ**2*c0*d0)/ss-JJ**2*d0**2/ss**2
check('harmonic-coordinate trace elimination',ss/tt*(1+PP+bv0**2)+ss-E)
check('exponential obstruction leading coefficient',s.expand(ss**4*s.diff(PP,ss)**2).coeff(ss,4)-c0**2)
x,y,nu,omega=s.symbols('x y nu omega',real=True);u0,bc,C0,n=s.symbols('u0 b c n',nonzero=True)
L=lambda f,var:(1-var**2)*s.diff(f,var,2)-2*var*s.diff(f,var)+6*f
P2=(3*x*x-1)/2;Q2=P2*s.log((1+x)/(1-x))/2-3*x/2
check('full Legendre basis',[L(P2,x),s.simplify(L(Q2,x))])
S=u0+bc*x*x;yy=C0*x/(u0*s.sqrt(S))+n
check('quadratic Legendre volume primitive',s.diff(yy,x)-C0/S**s.Rational(3,2))
Ty=1/u0-u0*bc*(y-n)**2/C0**2
check('reciprocal quadratic identity',s.simplify(Ty.subs(y,yy)-1/S))
check('inverse Legendre linear obstruction',s.expand(L(Ty,y)).coeff(y,1)-8*u0*bc*n/C0**2)
# Recurrence at epsilon=+/-1 checked without assumptions about a spatial axis.
for e in (-1,1):
    for n0 in range(5):
        coeff=[s.Symbol('a%d'%j) for j in range(7)]
        f=sum(coeff[j]*(y-e)**j for j in range(7))
        at=s.expand(L(f,y).subs(y,s.Symbol('t')+e)).coeff(s.Symbol('t'),n0)
        check('Legendre endpoint %s Taylor coefficient %s'%(e,n0),at-((6-n0*(n0+1))*coeff[n0]-2*e*(n0+1)**2*coeff[n0+1]))
U,V,Up,Vp,psip,st=s.symbols('U V Up Vp psi_prime sin_theta',nonzero=True)
Fdr=s.Matrix([[U,Up,0],[omega*U,U*psip,k*U/st],[V,Vp,0]])
check('dilation-rotation signed volume',Fdr.det()-k*U*(Up*V-U*Vp)/st)
t0=s.symbols('t',real=True);ma,mb,k0,c1=s.symbols('a m K0 c1',real=True)
Q=ma+mb*t0*t0;j=c1*s.sqrt(1+t0*t0)+k0*t0
check('azimuthal ODE full solution',(1+t0*t0)*s.diff(j,t0)-t0*j-k0)
F=s.Matrix([[ma/s.sqrt(Q),0,mb*t0/s.sqrt(Q)],[s.sqrt(Q)*(omega-t0*j/Q),k*s.sqrt(Q),j/s.sqrt(Q)],[0,0,nu]])
PP=ma**2+mb**2*t0**2+(k**2+omega**2)*Q**2+(1+t0*t0)*(c1**2*(1+t0*t0)+k0**2*t0*t0)-2*omega*k0*t0*t0*Q+(nu**2-E)*Q
check('phase trace identity',s.expand(Q*(s.trace(F.T*F)-E)-PP-2*c1*t0*(k0*(1+t0*t0)-omega*Q)*s.sqrt(1+t0*t0)))
check('phase sum-of-squares leading coefficient',s.expand(PP).coeff(t0,4)-(k*k*mb*mb+(k0-omega*mb)**2+c1*c1))
ell,offset,cc0,betaC=s.symbols('ell offset c0 betaC',real=True)
X0=1-x*x;W0=1-(offset+ell*x)**2
h=cc0+(betaC/(k*ma)-3*omega*(1-offset**2))*x+3*omega*offset*ell*x*x+omega*ell**2*x**3
for ep in (-1,1):
    check('constant norm endpoint %s'%ep,(omega**2*X0*W0**2+ell**2*X0**2+k**2*W0**2+h*h).subs(x,ep)-(k*k*W0.subs(x,ep)**2+h.subs(x,ep)**2))
check('constant norm surviving phase',s.cancel(h.subs({offset:0,ell:1,cc0:0,betaC:2*k*ma*omega})/(X0*(1-x*x)))+omega*x/(1-x*x))
cx,cy=s.symbols('cart_x cart_y',positive=True);ang=s.symbols('angle',real=True)
rho=s.sqrt(cx*cx+cy*cy);ca=s.cos(ang);sa=s.sin(ang)
ev=s.Matrix([(ca*cx-sa*cy)/rho,(ca*cy+sa*cx)/rho]);proj=ev*ev.T
div=s.Matrix([sum(s.diff(proj[i,j],(cx,cy)[j]) for j in range(2)) for i in range(2)])
pressure=s.cos(2*ang)*s.log(rho)+s.sin(2*ang)*s.atan(cy/cx)
check('log-cylinder projector local pressure from Cartesian derivatives',div-s.Matrix([s.diff(pressure,cx),s.diff(pressure,cy)]))
lamj,lamk,scale=s.symbols('lambda_j lambda_k scale',nonzero=True)
k2=lamj*lamk/scale**2;om2=(lamk-scale)*(scale-lamj)/scale**2
check('all-spectrum cylinder trace',scale*(1+k2+om2)-lamj-lamk)
check('all-spectrum cylinder determinant',scale**2*k2-lamj*lamk)
# Highest material-jet coefficients are obtained by the product rule.
uvec=s.Matrix(s.symbols('u1:4'));vvec=s.Matrix(s.symbols('v1:4'))
B=s.Matrix(3,3,lambda i,j:s.Symbol('B%d%d'%(min(i,j),max(i,j))));c=s.Matrix(3,3,lambda i,j:s.Symbol('c%d%d'%(min(i,j),max(i,j))))
q111,q112,q122,q222=s.symbols('W111 W112 W122 W222')
curlhigh=(q111*uvec+q112*vvec).cross(B*uvec)+(q112*uvec+q122*vvec).cross(B*vvec-c*uvec)-(q122*uvec+q222*vvec).cross(c*vvec)
expected=[uvec.cross(B*uvec),vvec.cross(B*uvec)+uvec.cross(B*vvec)-uvec.cross(c*uvec),vvec.cross(B*vvec)-vvec.cross(c*uvec)-uvec.cross(c*vvec),-vvec.cross(c*vvec)]
check('four independent highest material-jet wedge equations',[s.expand(curlhigh[i]).coeff(q)-expected[j][i] for j,q in enumerate((q111,q112,q122,q222)) for i in range(3)])
alpha=s.symbols('alpha',positive=True);q=s.symbols('q')
check('repeated spectrum characteristic polynomial',s.expand((q-alpha)**2*(q-alpha**-2))-(q**3-(2*alpha+alpha**-2)*q**2+(alpha**2+2/alpha)*q-1))
report=dict(status='PASS',checks=checks,scalar_identities=sum(c['scalar_identities'] for c in checks),elapsed_seconds=time.monotonic()-start,arithmetic='exact symbolic',sample_spectra_used=0)
(OUT/'geometric_lowrank_checks.json').write_text(json.dumps(report,indent=2))
print('PASS',len(checks),'groups',report['scalar_identities'],'scalar identities',flush=True)
