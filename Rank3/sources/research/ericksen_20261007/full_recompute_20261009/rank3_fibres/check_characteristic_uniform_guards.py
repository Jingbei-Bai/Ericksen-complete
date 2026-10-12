"""Turn all coefficient-table guard assertions into exact executable checks."""
from sage.all import *
from pathlib import Path
import json
out=Path(__file__).absolute().parent;p=out/'rebuild/round10';checks={}
d=load(str(p/'pencil_base_content_and_klein.sobj'));C=d['base_ring'];a,q,T,H=C.gens();r=q-1
guard=a*q*r*(1+a*r)*(a-1)*T*(T-1)*H*(H+T)*(H+1)*(1-H-T)
for i,f in enumerate(d['base_contents']):checks['base content protected '+str(i)]=all(not guard%g for g,e in f.factor())
c=load(str(p/'pencil_u_primitive_coefficient.sobj'))['coefficient']
B0=(T-2)*H**2+(3*T-4)*H+T*(T-1)
B1=2*(T-1)*H**2+(2*T**2-4*T+1)*H+T*(1-T)
C0=H**3+(-2*T**2+4*T)*H**2+(-3*T**3+6*T**2-2*T)*H-T**4+2*T**3-T**2
quot,rem=c.quo_rem(B0*B1*C0);checks['complete leading coefficient split']=not rem and all(not guard%g for g,e in quot.factor())
checks['B0 nonzero every T']=B0.subs({H:0})==T*(T-1)
checks['B1 nonzero every T']=B1.derivative(H,2)/2==2*(T-1)
checks['C0 monic every T']=C0.derivative(H,3)/6==1
rdata=load(str(p/'pencil_resultant_initial_certificate.sobj'));co=rdata['coefficients'][2];O=co.parent();aa,qq,tt=O.gens()
Psi=(qq-1)*tt*(aa*qq+1)-4
expected=-1024*tt**6*qq**8*(tt-2)**13*(tt-1)**13*aa**15*(qq-1)**40*(aa*qq+1)**5*Psi**2
checks['full generic resultant coefficient exact']=co==expected
checks['only three residual walls']=Psi.subs({aa:-1/qq})==-4
assert all(checks.values()),checks
(out/'characteristic_uniform_guards.json').write_text(json.dumps(dict(status='PASS',checks=checks,count=len(checks)),indent=2));print(checks)
