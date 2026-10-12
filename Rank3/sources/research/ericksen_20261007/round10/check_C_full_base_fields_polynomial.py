"""Check both complete lifted fields and their general-ratio base determinant.

Run with SageMath's Python. All inputs and the resulting JSON are in the
adjacent research/ericksen_20261007/round10 directory of the bundle.
"""
from sage.all import *
from pathlib import Path
import time,json
started=time.monotonic()
# Locate the extracted bundle (or the source workspace) from this file.
root = next(parent for parent in Path(__file__).resolve().parents
            if (parent / 'research/ericksen_20261007/round10').is_dir())
p = root / 'research/ericksen_20261007/round10'
d=load(str(p/'independent_partial_D_Newton_gauges_reader.sobj'))
K0=d['ring'].base_ring()
K=PolynomialRing(QQ,names=('a','r','H','Z')).fraction_field()
a,r,H,Z=K.gens(); q=1+r; m=1+a*r
coefficient_map=K0.ring().hom([a,a*q/m,H/r,q*Z/r],K)
cache={}
def cf(value):
    value=K0(value)
    if value not in cache:
        cache[value]=coefficient_map(value.numerator())/coefficient_map(value.denominator())
    return cache[value]
Q=PolynomialRing(K,names=('l','rho','kap','s','w','E0','e1','e2','delta'))
l,rho,kap,s,w,E0,e1,e2,delta=Q.gens()
alpha=(Z-H)/(Z*(H+1)); beta=H*(Z-1)/(Z*(H+1)); shift=-2*H*(Z+1)/(Z*(H+1))
v=alpha*l+kap*(E0+shift*l)
n1=(alpha/beta)*s*rho+kap*e1; n2=beta*w+kap*e2
pnums=[Q(1),n1,n2,s,w,l]
Rs=[v,rho,Q(1)]
normal_zero={s:0,w:0,E0:0,e1:0,e2:0,delta:0,kap:0}
records=[]
denominator_factors=set()
for direction,pindex in [(0,0),(2,5)]:
    record=next(it for it in d['items'] if it['direction']==direction and it['p_index']==pindex)
    degree=max(max(int(ex[1]),int(ex[2])) for f in list(record['field_numerators'])+list(record['base_numerators']) for ex in f.dict())
    T=(v*rho)**degree*(l if pindex==5 else 1)
    def pull(poly,kind):
        result=Q(0)
        for powers,coef in poly.dict().items():
            powers=tuple(map(int,powers)); pd,rd=sum(powers[:6]),sum(powers[6:])
            kd,dd={'p':(7-2*pd-3*rd,3-pd-rd),'R':(8-2*pd-3*rd,3-pd-rd),'b':(5-2*pd-3*rd,2-pd-rd)}[kind]
            assert kd>=0 and dd>=0
            result+=cf(coef)*prod(pnums[j]**powers[j] for j in range(6))*prod(Rs[j]**powers[6+j] for j in range(3))*kap**kd*delta**dd*v**(degree-powers[1])*rho**(degree-powers[2])
        return result
    fp=[pull(record['field_numerators'][j],'p') for j in range(6)]
    fr=[pull(record['field_numerators'][j],'R') for j in range(6,9)]
    bh=r*pull(record['base_numerators'][0],'b'); bz=(r/q)*pull(record['base_numerators'][1],'b')
    A,B=fp[0],fr[2]
    dl=fp[5]-l*A; drho=fr[1]-rho*B; dv=fr[0]-v*B
    ds=fp[3]-s*A; dw=fp[4]-w*A
    derivative=lambda c:c.derivative(H)*bh+c.derivative(Z)*bz
    n0=dv-alpha*dl-l*derivative(alpha)
    n_1=dv*n1+v*(v*fp[1]-n1*A)-v*(alpha/beta)*(drho*s+rho*ds)-v*rho*s*derivative(alpha/beta)
    n_2=drho*n2+rho*(rho*fp[2]-n2*A)-rho*beta*dw-rho*w*derivative(beta)
    print('direction',direction,'constructed',time.monotonic()-started,flush=True)
    divided=[]
    for name,numerator in [('E0',n0),('e1',n_1),('e2',n_2)]:
        assert numerator.subs({kap:0})==0,(direction,name,'not divisible')
        quot=Q({tuple(ex[j]-(1 if j==2 else 0) for j in range(9)):co for ex,co in numerator.dict().items()})
        assert numerator==kap*quot
        divided.append(quot)
    ee0=divided[0]-(E0+shift*l)*(A-B)-shift*dl-l*derivative(shift)
    ee1=divided[1]-v*e1*(A-B)
    ee2=divided[2]-rho*e2*(A-B)
    nums=[dl,drho,kap*(A-B),ds,dw,ee0,ee1,ee2,delta*(2*B-3*A)]
    dens=[T,T,T,T,T,T,v*T,rho*T,T]
    assert all(num.subs(normal_zero)==0 for num in nums[3:])
    # All chart denominators are units at the cyclic centre. Check every
    # coefficient denominator there before using the formal local ring.
    special=K.ring().hom([-q/r,r,1,-1],K)
    coeffs={co for num in nums+[bh,bz,A] for co in num.coefficients()}
    bad=[str(co.denominator()) for co in coeffs if special(co.denominator())==0]
    assert not bad,(direction,bad)
    for co in coeffs:
        center_den=K(special(co.denominator()))
        for polynomial in [center_den.numerator(),center_den.denominator()]:
            for factor,power in polynomial.factor():
                denominator_factors.add(str(factor))
    bnums=[bh.subs(normal_zero),bz.subs(normal_zero)]
    bt=T.subs(normal_zero)
    boundary=[Q.fraction_field()(f)/bt for f in bnums]
    print('direction',direction,'regular full base; all coefficient denominators units; base',list(map(str,boundary)),flush=True)
    records.append({'direction':direction,'base':boundary,'degree':degree,'regular_components':len(nums)+3,'denominators_unit':True})
F=Q.fraction_field()
det=matrix(F,[rec['base'] for rec in records]).det()
cyclic=K.ring().hom([a,r,1,-1],K)
def ccy(co):return cyclic(co.numerator())/cyclic(co.denominator())
def pcy(poly):return Q({ex:ccy(co) for ex,co in poly.dict().items()})
detcy=F(pcy(det.numerator()))/pcy(det.denominator())
assert detcy==3*l/4,detcy
report={'status':'PASS','full_general_ratio_base':True,'base_fields':[list(map(str,r['base'])) for r in records], 'base_determinant':str(det),'cyclic_base_determinant':str(detcy),'coefficient_denominator_factors_at_resonance':sorted(denominator_factors),'log_tau':'-A','seconds':time.monotonic()-started}
(p/'check_C_full_base_fields_polynomial.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2),flush=True)
