"""Principal frame and ratio equations on the full incidence boundary."""
from sage.all import *
from pathlib import Path
import json

p=Path(__file__).resolve().parent
g=load(str(p/'independent_ratio_p_boundary_second_gauge.sobj'))
old=g['coefficient_field']
K=PolynomialRing(QQ,names=('a','r','eta','zeta','s','z','b0','c0')).fraction_field()
a,r,H,Z,s,z,b0,c0=K.gens();q=1+r;m=1+a*r
coefficient_map=old.ring().hom([a,a*q/m,H/r,q*Z/r],K)
alpha=(Z-H)/(Z*(H+1));beta=H*(Z-1)/(Z*(H+1))
ps=[K.one(),b0*z,s,c0*z,s*b0,c0]
rs=[alpha*c0,beta*b0,K.one()]
images=ps+rs
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def face_value(f,degree):
    out=K.zero()
    for exponent,coefficient in f.dict().items():
        if (sum(exponent[:6]),sum(exponent[6:]))!=degree:continue
        coefficient=old(coefficient)
        out+=(coefficient_map(coefficient.numerator())/coefficient_map(coefficient.denominator()))*prod(images[j]**int(exponent[j]) for j in range(9))
    return out
records=[]
expected_base=[(H*c0/(H+1),H*c0/(H+1)),(-H,-H*Z/(H+1))]
for index,direction in enumerate([0,2]):
    dp=[face_value(g['first_gauged_fields'][direction,j],(1,1)) for j in range(6)]
    normal=[dp[j]-ps[j]*dp[0] for j in range(6)]
    for j,f in enumerate(normal):check('normalized frame component %d %d'%(direction,j),f==0)
    base=[r*face_value(g['base_fields'][direction,0],(0,1)),r/q*face_value(g['base_fields'][direction,1],(0,1))]
    for j,f in enumerate(base):check('ratio component %d %d'%(direction,j),f==expected_base[index][j])
    records.append({'direction':direction+1,'normal_frame_components':list(map(str,normal)),
                    'ratio_components':list(map(str,base))})
det=matrix(K,expected_base).det()
check('horizontal determinant',det==H**2*(H-Z+1)*c0/(H+1)**2)
check('cyclic horizontal determinant',det.subs({H:1,Z:-1})==3*c0/4)
report={'status':'PASS','checks':checks,'source':'independent_ratio_p_boundary_second_gauge.sobj',
        'coefficient_ring':str(K),'incidence_frame':list(map(str,ps)),
        'incidence_direction':list(map(str,rs)),'records':records,
        'horizontal_determinant':str(det),'cyclic_determinant':str(3*c0/4)}
stem=Path(__file__).stem
(p/(stem+'.json')).write_text(json.dumps(report,indent=2)+'\n')
save(report,str(p/(stem+'.sobj')),compress=False)
print(json.dumps(report,indent=2))
