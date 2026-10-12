"""Independent source/saturation audit with attempted explicit guard identities."""
from pathlib import Path
import json,time
from sage.all import load,save,QQ,PolynomialRing
from cysignals.alarm import alarm,cancel_alarm,AlarmInterrupt

out=Path(__file__).resolve().parent;start=time.monotonic()
source=load(str(out/'general_curve_dimensionless_leading.sobj'))
old=load(str(out/'general_curve_zero_bilinear_blocks.sobj'))
R=PolynomialRing(QQ,names=('H','Z','S'));H,Z,S=R.gens()
guards=[H,Z,H-Z,Z-1,H+1,H-Z+1,S-1,H*S-Z+1,H*S-Z+S]
guard=R.one()
for g in guards:guard*=g
def cv(c):
    result=R.zero()
    for mm,cc in c.dict().items():
        assert mm[0]==mm[1]==mm[5]==0
        result+=cc*H**mm[2]*Z**mm[3]*S**mm[4]
    return result
reports=[];data=[]
for saved in old['data']:
    pair=saved['pair'];B=source['blocks'][pair]
    raw=[cv(B[0,0]),cv(B[0,1]),cv(B[1,0])]
    polys=[R(p) for p in saved['polynomials']];removed=[R(p) for p in saved['removed']]
    checks={'raw_source':raw==[R(p) for p in saved['raw']],
            'guard_exact':guard==R(saved['guard']),
            'source_divisions':all(p*c==x for p,c,x in zip(polys,removed,raw))}
    for i,c in enumerate(removed):
        rest=c
        for g in guards:
            while rest:
                q,rem=rest.quo_rem(g)
                if rem:break
                rest=q
        checks['removed_protected_'+str(i)]=rest.is_constant() and rest!=0
    I=R.ideal(polys)
    sat,power=I.saturation(R.ideal(guard))
    checks['independent_unit_saturation']=list(sat.groebner_basis())==[R.one()]
    target=guard**power
    checks['guard_power_membership']=I.reduce(target)==0
    item={'pair':pair,'ring':R,'raw':raw,'polynomials':polys,'removed':removed,
          'guard':guard,'power':power,'target':target,'checks':checks}
    report={'pair':list(pair),'checks':checks,'power':int(power),'lift_status':'NOT_ATTEMPTED'}
    alarm(20)
    try:
        multipliers=list(target.lift(I))
        checks['explicit_identity']=sum((c*p for c,p in zip(multipliers,polys)),R.zero())==target
        assert checks['explicit_identity']
        item['multipliers']=multipliers
        report.update(lift_status='EXACT_IDENTITY_PASS',multiplier_terms=[len(p.dict()) for p in multipliers],
                      multiplier_degrees=[int(p.total_degree()) for p in multipliers])
    except AlarmInterrupt:
        report['lift_status']='UNKNOWN_TIMEOUT_20_SECONDS'
    finally:cancel_alarm()
    assert all(checks.values()),(pair,checks)
    data.append(item);reports.append(report)
    save({'data':data,'reports':reports,'scope':'Exact all-spectrum pointwise exclusion of a zero bilinear block on the nine guarded dimensionless chart. Rank-one blocks are retained.'},str(out/'independent_zero_bilinear_blocks.sobj'))
    (out/'independent_zero_bilinear_blocks.json').write_text(json.dumps({'reports':reports,'elapsed_seconds':time.monotonic()-start},indent=2)+'\n')
    print(json.dumps(report),flush=True)
print('TOTAL',time.monotonic()-start,flush=True)
