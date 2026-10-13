"""Original C-direction determinant and its protected base coefficient."""
from sage.all import *
from pathlib import Path
import json,time
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_delta_subcritical_Q0_direction.sobj'));e=load(str(p/'independent_delta_subcritical_even_directions.sobj'))['items'][1]
K=e['coefficient_field'];a,r,H,Z,t=K.gens();old=d['field'];ph=old.ring().hom([a,r,H,Z,-(t+1)/t,t],K);checks=[]
def cv(f):f=old(f);return ph(f.numerator())/ph(f.denominator())
def ck(n,v):
    if not bool(v):raise AssertionError(n)
    checks.append(n)
C=matrix(K,3,3,[cv(f) for f in d['C'].list()]);det=C.det()
ck('original C determinant source',det==cv(d['C_determinant']))
ck('even branch same C determinant',det==e['detC_direction'])
G=H*Z-H+Z
P=H**2*t**3+2*H*Z*t**2+2*H*Z*t-2*H*t**2+2*Z*t**2-t**3+Z**2-2*H*t+2*Z*t-3*t**2-3*t-1
unit=3*r*H**2*Z**3/(2*(Z-H)**3*(Z-1)**2*t*(H+1)*G)
ck('original determinant exact factor',det==unit*P)
poly=K.ring()(P);co={}
for ee,c in poly.dict().items():
    ex=tuple(map(int,ee));ck('no material in P '+str(ex),ex[0]==0 and ex[1]==0)
    if ex[2]==2:co[(ex[3],ex[4])]=c
ck('entire H2 coefficient is protected t3',co=={(0,3):QQ.one()})
guard=K.ring()(r*H*Z*(Z-H)*(Z-1)*t*(H+1)*G)
for f in (unit.numerator(),unit.denominator()):
    for fac,j in f.factor():
        if fac.total_degree():ck('unit factor '+str(fac),guard.quo_rem(fac)[1]==0)
record={'checks':len(checks),'seconds':time.monotonic()-st,'H2_coefficient':'t^3','scope':'Only a genuine invariant Newton component and the established projective-p descent rule can use this coefficient to exclude persistent Csing. Not pointwise detC nonvanishing on the line.'}
save({'checks':checks,'record':record,'C':C,'determinant':det,'factor':P,'unit':unit},str(p/'independent_subcritical_Csing_descent_reader.sobj'),compress=False)
(p/'independent_subcritical_Csing_descent_reader.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2),flush=True)
