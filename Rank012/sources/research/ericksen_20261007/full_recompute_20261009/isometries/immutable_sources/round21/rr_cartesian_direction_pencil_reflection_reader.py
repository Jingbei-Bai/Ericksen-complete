"""Independent reconstruction of both direction Jacobians and target reflection."""
from sage.all import *
from pathlib import Path
import time,signal,json
O=Path(__file__).resolve().parent;start=time.time();signal.alarm(100)
dp=load(str(O/'rr_cartesian_direction_paired_pencil.sobj'));R=dp['ring'];K=R.base_ring();B=K.ring();s,t,p,q,c=K.gens();aa,bb,xx=R.gens()
dm=load(str(O/'rr_null_physical_protected_graph_minus.sobj'));F=dm['F'];H=dm['H'];C=F.T*F;As=dm['axes'];At=dm['axet'];u=vector(R,[0,0,1])
g1=[];g2=[]
for j in range(3):
    ej=vector(R,[int(i==j) for i in range(3)]);dZ=As.cross_product(ej);du=dZ-u*(u*dZ)
    dF=matrix(R,3,3,lambda i,k:H[tuple(sorted((j,k)))][i]);dC=dF.T*F+F.T*dF
    g1.append(u*dC*u+2*du*C*u);g2.append(u*(dC*C+C*dC)*u+2*du*C*C*u)
D=g1[0]*g2[1]-g1[1]*g2[0]
Fi=F.inverse();Ci=Fi*Fi.T;qb=vector(R,[0,0,0])
for j in range(3):
    dF=matrix(R,3,3,lambda i,k:H[tuple(sorted((j,k)))][i]);dC=dF.T*F+F.T*dF;dCi=-Ci*dC*Ci
    for i in range(3):qb[i]+=dCi[j,i]
P=At*sum((H[j,j] for j in range(3)),vector(R,[0,0,0]))+F.det()*(As*qb)
cf=B.hom([s,-t,p,q,c],K)
def refl(f):
    out=R.zero()
    for ex,co in f.dict().items():
        cc=cf(co.numerator())/cf(co.denominator())
        out+=cc*aa**ex[0]*(-bb)**ex[1]*xx**ex[2]
    return out
J=D.derivative(xx)*P.subs({xx:0})-P.derivative(xx)*D.subs({xx:0})
checks={'source_orbit_annihilates_moments':g1[2]==0 and g2[2]==0,
        'D_has_no_chi_square':D.monomial_coefficient(xx*xx)==0,
        'D_full_target_reflection':D==refl(dp['D']),
        'paired_full_target_reflection':P==-refl(dp['paired']),
        'cubic_full_target_reflection':J==-refl(dp['cubic'])}
assert all(checks.values()),checks
report={'status':'PASS_DIRECTION_PENCIL_TWO_SHEET_REFLECTION','checks':checks,'seconds':time.time()-start}
(O/'rr_cartesian_direction_pencil_reflection_reader.json').write_text(json.dumps(report,indent=2),encoding='utf8');print(json.dumps(report,indent=2),flush=True)
