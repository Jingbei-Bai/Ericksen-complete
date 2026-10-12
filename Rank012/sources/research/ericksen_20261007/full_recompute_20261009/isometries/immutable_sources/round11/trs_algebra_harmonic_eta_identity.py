"""Exact full tensor source identity for the whole T->R/S family.

Arbitrary symmetric inverse metric H and arbitrary one-form a are retained
with all first and second derivatives. No material or pressure is specialized.
No third-derivative solution or finite-jet truncation supplies this identity.
"""
from pathlib import Path
import sympy as S
import pickle,json,time
OUT=Path(__file__).resolve().parent;START=time.time()
q,k,p,I1,KC=S.symbols('q k p I1 KC', nonzero=False)
# Symbols above carry no nonvanishing assertions; guards are stated separately.
fields=('H11','H12','H22','a1','a2')
orders=((0,0),(1,0),(0,1),(2,0),(1,1),(0,2))
jets={(name,i,j):S.Symbol(f'{name}_{i}{j}') for name in fields for i,j in orders}
def D(expr,axis):
    ans=S.diff(expr,q) if axis==0 else S.S.Zero
    for (name,i,j),symbol in jets.items():
        if i+j<2:
            ans+=S.diff(expr,symbol)*jets[name,i+(axis==0),j+(axis==1)]
    return S.expand(ans)
H=[[jets['H11',0,0],jets['H12',0,0]],[jets['H12',0,0],jets['H22',0,0]]]
a=[jets['a1',0,0],jets['a2',0,0]]
v=2*k*k*q+p*p;vq=2*k*k
K=[[1/(2*q),0],[0,2*q/v]]
C=[[K[i][j]+v*a[i]*a[j] for j in range(2)] for i in range(2)]
L=[[v*K[i][j] for j in range(2)] for i in range(2)]
def divergence(T):
    return [sum(D(H[i][j]*T[j][b],i)+T[i][j]*D(H[i][j],b)/2
                for i in range(2) for j in range(2)) for b in range(2)]
def curl(cov):return D(cov[1],0)-D(cov[0],1)
def clean(expr):return S.factor(S.cancel(S.together(expr)))
dc=divergence(C);dl=divergence(L)
EC=curl(dc);EL=curl(dl)
flux=[v*sum(H[i][j]*a[j] for j in range(2)) for i in range(2)]
V=sum(D(flux[i],i) for i in range(2))
trace=sum(H[i][j]*C[i][j] for i in range(2) for j in range(2))+v-I1
eta_cov=[p*a[0],p*a[1]+k*K[1][1]]
eta_flux=[sum(H[i][j]*eta_cov[j] for j in range(2)) for i in range(2)]
lap_eta=sum(D(eta_flux[i],i) for i in range(2))
f=-2*p*k/v**2
F=curl(a)-f
rhs=a[1]*D(V,0)-a[0]*D(V,1)+(k*k/v)*D(trace,1)+2*F*V+D(F,0)*flux[0]+D(F,1)*flux[1]
identity=clean(EC-EL/v+4*k/v*lap_eta-rhs)
assert identity==0
checks={'full_offshell_source_identity':True,
        'phase_curl_primitive':clean(p*f+k*S.diff(K[1][1],q))==0,
        'v_derivative':S.diff(v,q)==vq,
        'pitch_flux_coefficient':clean(S.diff(f,q)+2*vq*f/v)==0}
# Verify that the actual phase parametrization makes all curl defects vanish.
u,w,h11,h12,h22,c111,c112,c122,c222=S.symbols('u w h11 h12 h22 c111 c112 c122 c222')
phase={
 jets['a1',0,0]:u/k,jets['a2',0,0]:w/k+p/(k*v),
 jets['a1',1,0]:h11/k,jets['a1',0,1]:h12/k,
 jets['a2',1,0]:h12/k+S.diff(p/(k*v),q),jets['a2',0,1]:h22/k,
 jets['a1',2,0]:c111/k,jets['a1',1,1]:c112/k,jets['a1',0,2]:c122/k,
 jets['a2',2,0]:c112/k+S.diff(p/(k*v),q,2),jets['a2',1,1]:c122/k,jets['a2',0,2]:c222/k,
}
for label,e in [('curl',F),('curl_q',D(F,0)),('curl_zeta',D(F,1))]:
    checks['actual_phase_'+label]=clean(e.subs(phase))==0
for label,e,want in [('eta_q',eta_cov[0],p*u/k),('eta_zeta',eta_cov[1],(1+p*w)/k)]:
    checks[label]=clean(e.subs(phase)-want)==0
assert all(checks.values())
payload={'variables':(q,k,p,I1,KC),'jets':jets,'H':H,'a':a,'K':K,'C_horizontal':C,'L':L,
         'E_C':EC,'E_L':EL,'vertical_flux':flux,'vertical_divergence':V,
         'constant_trace_source':trace,'eta_covector':eta_cov,'lap_eta':lap_eta,
         'curl_defect':F,'right_hand_source_combination':rhs,'identity':identity,
         'phase_substitution':phase,'checks':checks}
(OUT/'trs_algebra_harmonic_eta_identity.pkl').write_bytes(pickle.dumps(payload))
summary={'status':'PASS','checks':checks,'check_count':len(checks),'seconds':time.time()-START,
 'identity':'EC-EL/v+(4k/v)Delta_G eta = a2*dq(V)-a1*dzeta(V)+(k^2/v)*dzeta(trace)+2F*V+dq(F)*fluxq+dzeta(F)*fluxzeta',
 'definitions':'V=div(vHa), trace=tr(H(K+vaa))+v-I1, F=da_qzeta+2pk/v^2, d eta=p a+k K22 d zeta',
 'consequence':'On the full constant-trace, constant-forward-vertical source system EC=EL=0 implies Delta_G eta=0, for k!=0 and v>0.',
 'scope':'Every pitch including zero, every pressure constant and positive constant spectrum; no flatness of G needed for this identity.'}
(OUT/'trs_algebra_harmonic_eta_identity.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2),flush=True)
