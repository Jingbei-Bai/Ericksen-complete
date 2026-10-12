"""Exact all-parameter checks for the actual-open psi_q=0 branch."""
from pathlib import Path
import sympy as S,json
out=Path(__file__).resolve().parent
q,v,w,k,p,j=S.symbols('q v w k p j', nonzero=True)
I1,I2,KB,KC=S.symbols('I1 I2 KB KC')
P=v**3-I1*v**2+I2*v-j
A=j/v-P/(v**2*w**2)
FT=(v*I2-j*(1+v*w**2))*A-v*A*A-j*v*v
Tgraph=-v*(KB*A+KC)/(2*k*w)
res=S.factor(FT-Tgraph**2)
num,den=S.fraction(res)
pol=S.Poly(num,w)
assert pol.degree()==8
assert S.factor(pol.LC())==-4*j*j*k*k*v**3
checks={'square_residual_degree_w_8':pol.degree()==8,'leading_all_material_nonzero':S.factor(pol.LC())==-4*j*j*k*k*v**3}
# Formula for G on a_q=0; independently substitute its full spectral expression.
aq,az,T=S.symbols('aq az T')
K1=1/(2*q);K2=2*q/v; Ki1=2*q;Ki2=v/(2*q)
B=v*Ki2*az**2;AA=j/v-(v**3-I1*v*v+I2*v-j)/(v*B)
ss=I2/j-(1+B)/v
aa=S.Matrix([0,az]);rr=S.Matrix([-Ki2*az,0]);KK=S.diag(K1,K2)
G=AA/j*KK+v/B*(ss-2*AA/j)*(aa*aa.T)-T/(j*B)*(aa*rr.T+rr*aa.T)
checks['Gqq_A_over_2jq']=S.cancel(G[0,0]-AA/(2*j*q))==0
checks['Gqz_T_over_jv']=S.cancel(G[0,1]-T/(j*v))==0
beta,Tq,Aw,A0,T0=S.symbols('beta Tq Aw A0 T0')
forward=-w*Tq/k+k*beta*(A0+w*Aw)
hcompat=T0/q+k*k*Aw*beta
checks['vertical_harmonic_cancellation']=S.expand(forward.subs(Tq,hcompat)-(k*beta*A0-w*T0/(k*q)))==0
checks['inverse_elliptic_beta'] = S.cancel((v/(2*q)).subs(v,2*k*k*q)*(-KB/k)+k*KB)==0
assert all(checks.values())
data={'status':'PASS','checks':checks,'check_count':len(checks),'square_residual_numerator':str(num),
 'square_residual_denominator':str(den),'highest_w_coefficient':str(S.factor(pol.LC())),
 'scope':'symbolic k,j,I1,I2,KB,KC,v retained; p!=0 reduction is differentiation of inverse phase equation; no chosen spectrum'}
(out/'trs_zero_radial_phase_reader.json').write_text(json.dumps(data,indent=2),encoding='utf-8')
print(json.dumps({k:data[k] for k in ['status','check_count','highest_w_coefficient']},indent=2))
