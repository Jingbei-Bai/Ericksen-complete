"""Independent direct-resultant/contact audit of the lambda D00-zero wall."""
from pathlib import Path
import json,time
from sage.all import load,save,QQ,PolynomialRing

out=Path(__file__).resolve().parent;start=time.monotonic();checks={}
d=load(str(out/'independent_scalar_lambda_lift.sobj'))
other=load(str(out/'general_curve_scalar_D00_zero.sobj'))
O=d['ring'];h,z,s,n=O.gens();K0=1+h-z
Pold=O(d['D'][0,0]/(-8*K0*z*h**2))
Rold=h**2*(2*z+1)*s+z**2*(1-2*h)
checks['D00_exact_source']=d['D'][0,0]==-8*K0*z*h**2*Pold
checks['B12_exact_source']=d['blocks'][1,2][0,0]==-4*K0*z*h*Rold
R=other['ring'];H,Z=R.gens();K=R.fraction_field();U=other['S_ring'];S=U.gen()
def cv(p):
    assert p.degree(n)<=0
    return sum((c*H**mm[0]*Z**mm[1]*S**mm[2] for mm,c in p.dict().items()),U.zero())
P=cv(Pold);Q=cv(Rold)
checks['source_polynomials']=P==other['P'] and Q==other['R_equation']
p=H**2*(2*Z+1);q=Z**2*(2*H-1)
checks['linear_slope']=Q==p*S-q
result=(-1)**P.degree()*sum((c*q**j*p**(P.degree()-j) for j,c in enumerate(P.list())),K.zero())
checks['direct_resultant']=result==other['resultant']
factorproduct=K.one()
for idx,(f,e,df,g) in enumerate(other['factor_checks']):
    contact=p*f.derivative(H)+q*f.derivative(Z);gg=f.gcd(contact)
    checks['contact_'+str(idx)]=contact==df and gg==g
    checks['coprime_or_old_guard_'+str(idx)]=gg==1 or f in [H,Z]
    factorproduct*=f**e
ratio=result/factorproduct
checks['factor_product']=bool(ratio) and ratio.derivative(H)==0 and ratio.derivative(Z)==0
Qhalf=sum((K(c.numerator()(H,-QQ(1)/2))/K(c.denominator()(H,-QQ(1)/2))*S**j for j,c in enumerate(Q.list())),U.zero())
checks['pivot_Zhalf']=Qhalf==QQ(1)/4-H/2
assert all(checks.values()),[key for key,val in checks.items() if not val]
save({'checks':checks,'P':P,'R':Q,'direct_resultant':result,'vector':(p,q),
      'scope':'Excludes D00 identically zero only on an actual scalar-monodromy open H interval with B01,B02 invertible. Isolated zeros and singular block charts are not excluded.'},
     str(out/'independent_scalar_D00_zero_reader.sobj'))
report={'status':'EXACT_SOURCE_CONTACT_PASS','checks':checks,'elapsed_seconds':time.monotonic()-start}
(out/'independent_scalar_D00_zero_reader.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report),flush=True)
