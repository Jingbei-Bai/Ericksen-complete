"""Independent exact recomputation of Section 5's explicit algebra.

No saved PASS record or Sage object is an input. All mathematical inputs below
are literal expressions transcribed from the maintained manuscript. A verified
consequence of an asserted source identity does not verify source provenance.
"""
from pathlib import Path
import sympy as s
import json, time, itertools, hashlib
OUT=Path(__file__).absolute().parent
START=time.monotonic(); checks=[]
def eq(name, expression, section, scope='identity'):
    answer=s.cancel(expression)==0
    checks.append(dict(name=name,passed=bool(answer),section=section,scope=scope))
    if not answer: raise AssertionError((name,s.factor(expression)))
def yes(name,condition,section,scope='identity'):
    checks.append(dict(name=name,passed=bool(condition),section=section,scope=scope))
    if not condition: raise AssertionError(name)

a,b,h,z,H,Z,k,r,q,m=s.symbols('a b h z H Z k r q m')
ra=(a-b)/(a*(b-1));qa=1+ra;ka=a*ra;ma=1+ka
eq('material q',qa-b*(a-1)/(a*(b-1)),'r3:algebra')
eq('material m',ma-(a-1)/(b-1),'r3:algebra')
eq('positive lambda',ma*qa-b/a*ma**2,'r3:algebra')
eq('positive N',ka*ra-a*ra**2,'r3:algebra')
for name,actual,expected in [('r',(b-a)/(b*(a-1)),-ra/qa),('q',1+(b-a)/(b*(a-1)),1/qa),('k',b*(b-a)/(b*(a-1)),-ka/ma),('m',1+(b-a)/(a-1),1/ma)]:
    eq('swap23 '+name,actual-expected,'r3:boundarysources')
Aa=b/a;Bb=1/a;rr=(Aa-Bb)/(Aa*(Bb-1));qq=1+rr
eq('cyclic r',rr+1/qa,'r3:curves')
eq('cyclic q',qq-ra/qa,'r3:curves')
eq('cyclic H',rr*z/h+(ra*z/qa)/(ra*h),'r3:curves')
eq('cyclic Z',rr/(qq*h)+1/(ra*h),'r3:curves')

G0=H*Z-H+Z;V0=(H+2)*(Z-2)*(2*H+1)*(Z-2*H)*(2*Z-1)*(2*Z-H)
eq('primary G0 center',G0.subs({H:1,Z:-1})+3,'r3:boundarysources')
eq('primary V0 center',V0.subs({H:1,Z:-1})-729,'r3:boundarysources')
x,y,w,t,v,E,X,Y,eta,ell,rho=s.symbols('x y w t v E X Y eta ell rho')
fA=(-8*v+152*q)*w+t*(-232*v-392*q)
fB=(392*v+232*q)*w+t*(-152*v+8*q)
coef=s.Matrix([[s.diff(fA,w),s.diff(fA,t)],[s.diff(fB,w),s.diff(fB,t)]])
eq('missing row determinant',coef.det()-92160*(v+q)**2,'r3:boundarysources')
eq('missing row w=t',fA.subs(v,-q)-160*q*(w-t),'r3:boundarysources')
normal=s.Matrix([[0,-64*w*(x+z)**2],[r*z*w**2*(x+z)/2,X]])
eq('single RR normal determinant',normal.det()-32*r*z*w**3*(x+z)**3,'r3:boundarysources', 'consequence of stated normal matrix')
J=s.groebner([w*w-k,x*z-10*(1+k),3*(x*x+z*z)-x*z],x,z,w,domain=s.QQ.frac_field(k,r))
eq('eight candidate sum',J.reduce((x+z)**2-s.Rational(70,3)*(1+k))[1],'r3:boundarysources')
det=32*r*z*w**3*(x+z)**3
dinv=x*w*(x+z)**3/(32*r*10*(1+k)*k**2*(s.Rational(70,3)*(1+k))**3)
eq('eight candidate determinant inverse',J.reduce(s.cancel(det*dinv-1))[1],'r3:boundarysources')

theta=(Z-H-1)*v-2*H*Z
base=s.Matrix([[ell*H*(v-2*H)/(H+1),ell*H*(v+2)/(H+1)],[-v*H,-H*Z*(v+2)/(H+1)]])
eq('ordinary active base determinant',base.det()+ell*H**2*(v+2)*theta/(H+1)**2,'r3:ordinaryincoming')
Nb,Lb=s.symbols('N_b Lambda_b')
eq('ordinary complementary minors',5*(5*Lb-2*Nb)-2*(11*Lb-5*Nb)-3*Lb,'r3:ordinaryincoming')
eq('ordinary resonance equivalence',(ma*qa-ka*ra)-(ka+qa),'r3:ordinaryincoming')

q1=eta**2-4*y*eta+64*y*y;hh=4*y*y-y*eta
q2=E*E-64*k*y*y-4*k*y*eta-k*eta*eta
for name,e in [('q1-4h',q1-4*hh-eta*eta-48*y*y),('negative h',-hh-y*eta+4*y*y),('E square',q2+k*q1-8*k*hh-E*E+32*k*y*y),('cubic nilpotence',y*q1+eta*hh-64*y**3)]: eq(name,e,'r3short:pure-axis')
GB=s.groebner([q1,hh,q2],E,eta,y,domain=s.QQ.poly_ring(k))
for i in range(5):
 for j in range(5-i):
  l=4-i-j;eq('degree4 ideal '+str((i,j,l)),GB.reduce(E**i*eta**j*y**l)[1],'r3short:pure-axis')
yes('sharp cubic survives',GB.reduce(E*y*y)[1]!=0,'r3short:pure-axis')
D9=29*rho*rho+49*q;D15=q-19*rho*rho
eq('inactive alternate pivot',D15.subs(rho*rho,-49*q/29)-960*q/29,'r3:boundarysources')
eq('inactive alternate quotient',(rho*rho+16*q).subs(rho*rho,-49*q/29)-415*q/29,'r3:boundarysources')
eq('paired wall complementary source',5*(5*m*ell**2-4*k)-4*(4*m*ell**2-5*k)-9*m*ell**2,'r3:boundarysources')
eq('resonance positive expression',r*r+r+s.Rational(13,4)-(r+s.Rational(1,2))**2-3,'r3:boundarysources')
eq('resonance equation', (1+k)*(1+r)-k*r-(1+k+r),'r3:boundarysources')
eq('resonance bound',r*(-1-r)+s.Rational(1,4)-(s.Rational(1,4)-r*(1+r)),'r3:boundarysources')
for name,A,B,identity in [('pair1',Z-3*H,H-3*Z,lambda A,B:-(3*A+B)/(8*H)-1),('pair2',Z-3,3*Z-1,lambda A,B:(B-3*A)/8-1),('pair3',H+3,3*H+1,lambda A,B:(3*A-B)/8-1)]: eq('supercritical '+name,identity(A,B),'r3:boundarysources')

# Balanced-line elimination is rederived from the printed rational vector field.
Ab=a*(b-1);Cb=b*(a-1);Db=a-b
NY=a*b*(a*b-1)-2*a*b*(b-1)*h-b*Db*h*h+2*a*b*(a-1)*z-3*a*Db*z*z
YY=-NY/(4*h*b*Db)
N3=a*b*(a*b-1)+2*a*b*(b-1)*h+3*b*Db*h*h-2*a*b*(a-1)*z+a*Db*z*z
mu=1-Db*z/Cb
ff=2*(1-2*Db**2*h*z/(Ab*Cb))*YY+N3/(2*b*Ab)-2*h*mu/Ab*((Ab+Db*h)*s.diff(YY,h)+Db*z*s.diff(YY,z))
Nf=s.cancel(8*(b-a)*h*(b-1)*(a-1)*a*b*b*ff)
yes('balanced numerator polynomial',s.denom(Nf)==1,'r3:balancedline')
restricted=s.cancel(Nf.subs(z,Cb*(Ab+Db*h)/(Ab*Db)))
eq('balanced constant coefficient',restricted.subs(h,0)-4*a**3*(b-1)**3*(a-1)*b*b/(b-a),'r3:balancedline')
U,V,W=s.symbols('U V W')
BQ=(V-W)**2-U*(V+W);BQ2=(U+W)**2-V*(U-W)
QM=s.hessian(BQ,(U,V,W))/2
eq('balanced quadric determinant',QM.det()+1,'r3:balancedline')
yes('balanced quadrics coprime',s.Poly(s.gcd(BQ,BQ2),U,V,W).total_degree()==0,'r3:balancedline')
u,vv,om=s.symbols('u vv omega')
pb=BQ.subs({U:Db*(b-1)*u*u,V:Db*(a-1)*vv*vv,W:(b-1)*(a-1)*om*om})
pc=BQ2.subs({U:Db*(b-1)*u*u,V:Db*(a-1)*vv*vv,W:(b-1)*(a-1)*om*om})
yes('balanced pulled quartics generic coprime',s.Poly(s.gcd(pb,pc),u,vv,om).total_degree()==0,'r3:balancedline','generic algebra; all-material follows finite morphism argument')

S,T=s.symbols('S w');L=H*S-Z;K=1+H-Z;c=L/(2*H*Z)
Qstar=2*Z*T*(2*H*Z*K*T-(L*L+H*S*S-Z))
Q=(2*H*Z*T+L)*(2*Z*K*T-(S-1)*(L+S))
eq('moment scalar translation',Qstar.subs(T,T+c)-Q,'r3:curves')
aq=s.expand(Q).coeff(T,2);bq=s.expand(Q).coeff(T);cq=Q.subs(T,0)
eq('moment leading Q',aq-4*H*Z*Z*K,'r3:Qcoefficients')
P=s.Matrix([[1,H*(T+bq/aq),-T,H/aq],[0,1,0,0],[0,-H,1,0],[0,-aq*T,0,1]])
eq('column transform det',P.det()-1,'r3:columntransform')
ys=s.symbols('y1:5');xs=P*s.Matrix(ys);y1,y2,y3,y4=ys
qr=T*(aq*y1*y2+bq*y2*y3+y3*y4)+cq*y2*(y3-H*y2)-y1*y4-H/aq*y4*(bq*y2+y4)
eq('transformed quadric',Q*xs[1]*xs[2]-xs[0]*xs[3]-qr,'r3:transformedmomentquadric')
weights=[1,0,1,0]
for i in range(4):
 for j in range(i,4):
  cc=s.expand(qr).coeff(ys[i],2) if i==j else s.expand(qr).coeff(ys[i]).coeff(ys[j])
  yes('quadric degree '+str((i,j)),cc==0 or s.degree(cc,T)<=weights[i]+weights[j],'r3:transformedmomentquadric')

A=s.Matrix(2,2,s.symbols('a0:4'));B=s.Matrix(2,2,s.symbols('b0:4'));C=s.Matrix(2,2,s.symbols('c0:4'));D=s.Matrix(2,2,s.symbols('d0:4'))
dd=A.det();AA=A.adjugate();RR=dd*D-C*AA*B
NN=A.row_join(B).col_join(C.row_join(D));LL=(-AA*B).col_join(dd*s.eye(2));JJ=RR.adjugate()
eq('fraction free determinant',s.expand(RR.det()-dd*NN.det()),'r3schur:basic')
for i,e in enumerate(NN*LL-s.zeros(2,2).col_join(RR)):eq('fraction free lift '+str(i),s.expand(e),'r3schur:basic')
for i,e in enumerate(LL*JJ-dd*NN.adjugate()[:,2:]):eq('kernel lift '+str(i),s.expand(e),'r3schur:columns')
block=(RR.det()*AA+AA*B*JJ*C*AA).row_join(-dd*AA*B*JJ).col_join((-dd*JJ*C*AA).row_join(dd*dd*JJ))
for i,e in enumerate(block-dd*dd*NN.adjugate()):eq('full adjugate '+str(i),s.expand(e),'r3schur:fulladj')
d0,d1,d2=s.symbols('d0 d1 d2');fs=s.symbols('f0:5');dp=d2*T*T+d1*T+d0;fp=sum(fs[i]*T**i for i in range(5))
af=d2**3*fs[1]-d1*d2*d2*fs[2]+d2*(d1*d1-d0*d2)*fs[3]+(-d1**3+2*d0*d1*d2)*fs[4]
bf=d2**3*fs[0]-d0*d2*d2*fs[2]+d0*d1*d2*fs[3]+(d0*d0*d2-d0*d1*d1)*fs[4]
eq('quartic pseudoremainder',s.rem(d2**3*fp-af*T-bf,dp,T),'r3schur:norm')
result=s.resultant(dp,fp,T);numnorm=d2*bf*bf-d1*af*bf+d0*af*af
eq('norm identity',s.expand(numnorm-d2**3*result),'r3schur:normidentity')
eq('norm leading drop',s.expand((numnorm-d2**3*result).subs(d2,0)),'r3schur:normidentity')
eq('double-root resultant',result.subs({d0:1,d1:-2,d2:1})-sum(fs)**2,'r3schur:normidentity')
yes('D0 degree',sum([12,20,20,25,25,26,6])==134,'r3:contact')
yes('projection degree arithmetic',481-134==347,'r3:contact')
yes('infinity degree arithmetic',487-6==481,'r3:contact')
out=dict(status='PASS',checks=checks,check_count=len(checks),elapsed_seconds=time.monotonic()-START,no_material_specialization=True,no_historical_certificate_input=True,
 scope='Exact printed identities and consequences. Complete-source claims, F/G construction, source support coverage and modular computation have separate statuses; this report does not certify all Rank 3.')
(OUT/'fresh_printed_algebra.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps({k:out[k] for k in ['status','check_count','elapsed_seconds']},indent=2))
