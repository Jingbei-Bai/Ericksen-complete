"""Independent exact review of all three one-pair projective P1 boundaries."""
from pathlib import Path
import json,time
from sage.all import load,save,prod

out=Path(__file__).resolve().parent;start=time.monotonic();checks={}
extra=load(str(out/'general_curve_extra_quartic_cubics.sobj'))
cert=load(str(out/'general_curve_single_pair_certificate.sobj'))
R=extra['ring'];K=R.base_ring();a,r,H,Z,S,N=K.gens();q=1+r;xs=R.gens()
allowed=[a,r,q,1+a*r,H,Z,H-Z,Z-1,H+1,1+H-Z,S-1,H*S-Z+1,(H+1)*S-Z]
def unit(c):
    c=K(c)
    if not c:return False
    return all(f.is_constant() or any(all((K(f)/g).derivative(v)==0 for v in K.gens()) for g in allowed) for p in [c.numerator(),c.denominator()] for f,e in p.factor())
forms=[]
for pair in range(3):
    pairforms=[]
    for j,qi in enumerate(cert['indices'][pair]):
        idx=next(i for i,spec in enumerate(extra['row_specs']) if spec[0]==qi)
        p=extra['leading_rows'][idx]
        restricted=R({tuple(mm):c for mm,c in p.dict().items() if all(mm[k]==0 for k in range(6) if k not in [2*pair,2*pair+1])})
        c=restricted.monomial_coefficient(xs[2*pair]**3)
        checks['source_coefficient_'+str(qi)]=c==cert['units'][pair][j]*cert['polynomials'][pair][j]
        checks['unit_protected_'+str(qi)]=unit(cert['units'][pair][j])
        checks['homogeneous_binary_'+str(qi)]=all(sum(mm)==3 for mm in restricted.dict())
        pairforms.append(restricted)
    combination=cert['polynomials'][pair][1]+(1 if pair==2 else -1)*cert['polynomials'][pair][0]
    expected=[(Z-H)*(S-1),(Z-1)*(H*S-Z+1),(H+1)*((H+1)*S-Z)][pair]
    checks['combination_identity_'+str(pair)]=combination==cert['targets'][pair]==expected
    checks['combination_protected_'+str(pair)]=unit(expected)
    forms.append(pairforms)
assert all(checks.values()),[key for key,val in checks.items() if not val]
save({'checks':checks,'binary_forms':forms,'scope':'For every allowed full-parameter base, each one-pair P1 has at most three common zeros of these necessary leaders. No assumption of analyticity, genericity, N nonzero, or any new pivot. This does not exclude all one-pair directions or establish the full fiber bound.'},str(out/'independent_single_pair_reader.sobj'))
report={'status':'EXACT_UNIFORM_THREE_P1_AUDIT_PASS','check_count':len(checks),'checks':checks,'elapsed_seconds':time.monotonic()-start}
(out/'independent_single_pair_reader.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report),flush=True)
