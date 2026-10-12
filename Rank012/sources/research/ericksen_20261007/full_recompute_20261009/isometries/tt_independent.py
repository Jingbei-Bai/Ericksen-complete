"""Fresh TT moving-frame jet elimination, starting from the two-metric block.
No archived symbolic expression or material sample is loaded.
"""
from pathlib import Path
import sympy as s
import json,time
OUT=Path(__file__).absolute().parent
start=time.time()
t,K,d,E,d0,B,a,b,P,Q,c,e=s.symbols('t K d E d0 B a b P Q c e')
C0=K+d/t**2;A=E-t**2-C0;B2=s.expand(A*C0-d0)
checks={};expressions={}
def red(x):
    n,de=s.fraction(s.cancel(x))
    n=s.rem(n,B**2-B2,B);de=s.rem(de,B**2-B2,B)
    return s.factor(s.cancel(n/de))
def eq(name,x):
    z=red(x);checks[name]=z==0
    if z!=0:expressions[name+'_residual']=str(z)
    print(name,checks[name],round(time.time()-start,1),flush=True)
Bp=s.diff(B2,t)/(2*B)
def dt(x):return s.diff(x,t)+s.diff(x,B)*Bp
u=-2*d/(d0*t**2)
N=d*(E-2*K)*t**2+K*t**6-2*d**2
w=N/(B*d0*t**4)
eq('inverse_divergence_b_coefficient',w+t*Bp/d0)
inv=u*a+w*b+C0*P/d0-Q
eq('affine_hessian_line',inv*B*d0*t**4-(-2*B*d*t**2*a+N*b+B*K*P*t**4+B*P*d*t**2-B*Q*d0*t**4))
def dn(x):return red(dt(x)*a+s.diff(x,a)*(c+2*b*b/t)+s.diff(x,b)*(e+(P-2*a)*b/t))
def dd(x):return red(dt(x)*b+s.diff(x,a)*(e+2*b*(P-a)/t)+s.diff(x,b)*(-c+(P-2*a)*(P-a)/t))
eq('third_symbol_11',s.diff(dn(inv),c)-u)
eq('third_symbol_12',s.diff(dn(inv),e)-w)
eq('third_symbol_21',s.diff(dd(inv),c)+w)
eq('third_symbol_22',s.diff(dd(inv),e)-u)
det=red(u*u+w*w)
eq('positive_third_determinant',det-(4*B*B*d*d*t**4+N*N)/(B*B*d0*d0*t**8))
# Form the homogeneous Hessian coefficient BEFORE any inverse expansion.
# Lower affine terms are immaterial for this coefficient and remain arbitrary.
homsub={P:0,Q:0,a:N,b:2*B*d*t**2}
rn=red(dn(inv).subs({c:0,e:0}).subs(homsub,simultaneous=True));rt=red(dd(inv).subs({c:0,e:0}).subs(homsub,simultaneous=True))
Vn=dt(A)*a+Bp*b+(A-C0)*(P-a)/t-2*B*b/t
Vt=Bp*a+dt(C0)*b+(A-C0)*b/t+2*B*(P-a)/t
curl=red((dn(Vt)-dd(Vn)+Vn*b/t+Vt*(P-a)/t).subs(homsub,simultaneous=True))
cc=red(s.diff(curl,c));ce=red(s.diff(curl,e));c0=red(curl.subs({c:0,e:0}))
leading=red(det*c0+cc*(-u*rn+w*rt)+ce*(-w*rn-u*rt))
expressions['adjugate_residual_line_direction_rho_squared']=str(leading)
expressions['literal_curl_line_direction_rho_squared']=str(red(leading/det))
D=red(-t**4*B**2)
H=red(-leading*d0**2*D**3/(B*t**6))
expressions['H_claim_reconstructed']=str(H)
checks['H_is_polynomial_in_t']=s.denom(H)==1
if checks['H_is_polynomial_in_t']:
    hp=s.Poly(H,t)
    checks['H_degree_24']=hp.degree()==24
    checks['H_leading_minus_K6']=s.factor(hp.LC()+K**6)==0
else:
    checks['H_degree_24']=False;checks['H_leading_minus_K6']=False
report={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'expressions':expressions,'seconds':time.time()-start,'scope':'generic parameters, neither pressures nor material spectrum specialized; P=Q=0 only extracts quadratic homogeneous Hessian coefficient'}
(OUT/'tt_independent.json').write_text(json.dumps(report,indent=2))
print(json.dumps({'status':report['status'],'checks':checks,'seconds':report['seconds']},indent=2),flush=True)
