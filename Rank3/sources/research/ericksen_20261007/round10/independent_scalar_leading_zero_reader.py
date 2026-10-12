"""Independent short audit of the scalar N-leading-zero actual-open branch."""
from pathlib import Path
import time,json
from sage.all import load,save,QQ,PolynomialRing,matrix

out=Path(__file__).resolve().parent;start=time.monotonic();checks={}
d=load(str(out/'general_curve_scalar_contact_source.sobj'))
saved=load(str(out/'general_curve_scalar_leading_zero.sobj'))
R=d['ring'];H,Z,S,N,W=R.gens();F=d['rows'][1]
A=F.derivative(N);B=F.subs({N:0,W:0})
L1=H*(H-1)*S+2*H+Z;L2=(H+2*Z)*S-Z*(Z+1)
checks['source_linear_N']=F==A*N+B and A.degree(N)==0 and F.degree(W)==0
checks['leading_factorization']=A==-2*H*Z*(Z-H)*(1+H-Z)*L1*L2
P=PolynomialRing(QQ,names=('H','Z'));HH,ZZ=P.gens();K=P.fraction_field()
U=PolynomialRing(K,'S');SS=U.gen()
constant=U.zero()
for mon,c in B.dict().items():
    assert mon[3]==0 and mon[4]==0
    constant+=c*HH**mon[0]*ZZ**mon[1]*SS**mon[2]
checks['constant_source']=constant==saved['source_constant']
rows=[]
for idx,(p,q) in enumerate([(HH*(HH-1),-(2*HH+ZZ)),(HH+2*ZZ,ZZ*(ZZ+1))]):
    # For L=p*S-q, Res_S(L,B)=p^degree(B)*B(q/p).
    # This direct identity avoids reusing the generator's resultant routine.
    result=K.zero()
    n=constant.degree()
    for j,c in enumerate(constant.list()):result+=c*q**j*p**(n-j)
    checks['resultant_direct_'+str(idx)]=result==saved['rows'][idx]['resultant']
    components=[]
    for j,(f,e,df,oldgcd) in enumerate(saved['rows'][idx]['factor_checks']):
        derivative=p*f.derivative(HH)+q*f.derivative(ZZ)
        gg=f.gcd(derivative)
        checks['factor_contact_'+str(idx)+'_'+str(j)]=derivative==df and gg==oldgcd
        protected=(f==HH if idx==0 else f==ZZ)
        checks['nonprotected_coprime_'+str(idx)+'_'+str(j)]=protected or gg==1
        components.append((f,e,derivative,gg,protected))
    reconstructed=K.one()
    for f,e,_,_,_ in components:reconstructed*=f**e
    checks['factor_product_'+str(idx)]=bool(result/reconstructed) and (result/reconstructed).derivative(HH)==0 and (result/reconstructed).derivative(ZZ)==0
    rows.append({'vector':(p,q),'direct_resultant':result,'components':components})

# Pivot strata: H=0 is excluded; H=1 cannot contain an interval on which
# H itself is the independent variable. On H+2Z=0, L2=0 implies Z(Z+1)=0.
checks['pivot1_values']=L1.subs({H:1})==Z+2
checks['pivot2_values']=L2.subs({H:-2*Z})==-Z*(Z+1)
assert all(checks.values()),[k for k,v in checks.items() if not v]
save({'checks':checks,'source':F,'leading':A,'constant':B,'walls':(L1,L2),'rows':rows,
      'scope':'No actual H-parametrized scalar-monodromy open interval can have A identically zero. No assertion that A is nonzero at every point, and no exclusion of the remaining A!=0 scalar branch.'},
     str(out/'independent_scalar_leading_zero_reader.sobj'))
report={'status':'EXACT_SOURCE_RESULTANT_CONTACT_PASS','checks':checks,'elapsed_seconds':time.monotonic()-start}
(out/'independent_scalar_leading_zero_reader.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report),flush=True)
