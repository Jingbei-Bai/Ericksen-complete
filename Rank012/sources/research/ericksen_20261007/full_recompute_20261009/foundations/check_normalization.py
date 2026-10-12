"""Independent printed cubic and constant-ratio substitutions in fresh fields."""
from sage.all import *
from pathlib import Path
import json,time
OUT=Path(__file__).absolute().parent; START=time.monotonic(); checks=[]
src=load(str(OUT/'state17_fresh.sobj'));P=src['ring'];v=P.gens();K=P.base_ring();a,b=K.gens();F=src['F']
p11,p12,p13,p21,p22,p23,p31,p32,p33,t0,t1,t2,t3,t4,t5,u0,u1=v
def check(name,ok,**extra):
    checks.append(dict(name=name,passed=bool(ok),**extra));print(name,bool(ok),flush=True);assert ok,name
A=a*(b-1);C=b*(a-1);D=a-b
C1=p11**2*p23-p11**2*p32+D/A*(p21*t3+p22*t1-p21**2*p23-p22*t4-p21**2*p32-p31*t2/2+2*p21*p22*p31)+D/C*(p11*p32*p33+p33*t4-p21*p31*p33)
C2=p11*p12*p32-p11*t3+p11*t5-p12**2*p13-p12**2*p31-p12*t1+D/A*(p13*p22**2-p22**2*p31-p32*t2/2+p21*p22*p32)+D/C*(-p33*t5-p12*p32*p33+p22*p31*p33)
C3=p11*p12*p33-p11*p13*p32-p11*p23*p31-p11*p31*p32-p11*t0-p12*p13**2-p13*t1+2*p11*p13*p23+p11*t2*C/(2*A)-p12*p13*p31+D/C*(-p12*p33**2+p21*p33**2)+D/A*(p21*p23**2-p22*t0-p23*t3+p21*p23*p32-p21*p22*p33)-C/A*p11*p21*p22
printed=[C1,C2,C3]
if printed!=src['C3']:print('cubic differences',list(f-g for f,g in zip(printed,src['C3'])),flush=True)
check('all three printed cubics match full fresh source',printed==src['C3'])
mu=[D*(p33/C-p22/A),p11-D*p33/C,-p11-D*p22/A]
check('cubic reconstruction coefficient columns',matrix(P,[[f.derivative(t) for t in (t4,t5,t0)] for f in printed])==diagonal_matrix(P,mu))
weights=[1]*9+[2]*6+[3]*2
check('Euler weighted homogeneity of all state components',all(sum(vv*f.derivative(vv)*ww for vv,ww in zip(v,weights))==(weights[j]+1)*f for row in F for j,f in enumerate(row)))
# Constant ratios: differentiate undivided pair relations, then substitute.
R=PolynomialRing(QQ,names=('a','b','r','z','x','y','w','tau0','tau1','tau2','tau3','tau4','tau5','ups0','ups1'));a,b,r,z,x,y,w,*extra=R.gens();L=R.fraction_field()
beta=(L.one(),L(a),L(b));alpha=lambda i,j:1-beta[i]/beta[j]
rr=[-alpha(1,0)/(alpha(1,2)*z),-alpha(2,1)*r/alpha(2,0),-alpha(0,2)*z/(alpha(0,1)*r)]
point=[1,y,rr[0]*x,rr[1]*y,r,w,x,rr[2]*w,z]+list(extra)
def phi(f):return sum((L(c.numerator()(a,b))/L(c.denominator()(a,b))*prod(point[i]**e for i,e in enumerate(mon)) for mon,c in f.dict().items()),L.zero())
derived=[phi(F[1][2])-rr[0]*phi(F[1][6]),phi(F[2][3])-rr[1]*phi(F[2][1]),phi(F[0][7])-rr[2]*phi(F[0][5])]
# These directional pair derivatives contain no free jet after the ratios are fixed.
check('constant-ratio opposite-pair derivatives eliminate free jets',all(not f.derivative(e) for f in derived for e in extra))
expected=matrix(L,[[rr[0]*(rr[0]+1),-alpha(2,0)/alpha(1,0)*rr[1],-alpha(0,2)/alpha(1,2)*rr[0]*rr[2]],[-alpha(1,0)/alpha(2,0)*rr[0]*rr[1],rr[1]*(rr[1]+1),-alpha(0,1)/alpha(2,1)*rr[2]],[-alpha(1,2)/alpha(0,2)*rr[0],-alpha(2,1)/alpha(0,1)*rr[1]*rr[2],rr[2]*(rr[2]+1)]])
coeff=matrix(L,[[f.derivative(e,2)/2 for e in (x,y,w)] for f in derived])
scales=[]
for i in range(3):
    j=next(j for j in range(3) if expected[i,j]);scale=coeff[i,j]/expected[i,j];scale.reduce()
    check('pair derivative square row %d'%(i+1),scale==1 and all(coeff[i,k]==expected[i,k] for k in range(3)),row_multiplier=str(scale))
    scales.append(scale)
check('opposite-pair square row multipliers nonzero protected',all(scales))
report=dict(status='PASS',checks=checks,elapsed_seconds=time.monotonic()-START,sample_spectra_used=0)
(OUT/'normalization_checks.json').write_text(json.dumps(report,indent=2))
print('COMPLETE',report['elapsed_seconds'])
