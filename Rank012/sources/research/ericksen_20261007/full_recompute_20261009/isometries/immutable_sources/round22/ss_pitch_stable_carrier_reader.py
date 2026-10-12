"""Whole-material first-order SS carrier at a pitch-matched RR cylinder.

Exact forward/inverse linear universal sources and fixed-pressure removal
are replayed on both twist charts. Completeness and the extra pitch direction
are justified by the accompanying physical-equivariance/dual-number proof.
"""
from sage.all import *
from pathlib import Path
import json,time,signal
O=Path(__file__).resolve().parent; start=time.time(); signal.alarm(90)
K=FractionField(PolynomialRing(QQ,names=('E','k','omega','c')))
E,k,om,c=K.gens()
R=LaurentPolynomialRing(K,names=('r','z','lr','A0','A1','A2','b','v','h','Gamma','B0','W0'))
r,z,lr,A0,A1,A2,b,v,h,Ga,B0,W0=R.gens()
dr=lambda f:R(f).derivative(r)+R(f).derivative(lr)/r
dz=lambda f:R(f).derivative(z)
def values(H):
    vr=dr(H[0,0])+dz(H[0,2])+(H[0,0]-H[1,1])/r
    vt=dr(H[0,1])+dz(H[1,2])+2*H[0,1]/r
    vz=dr(H[0,2])+dz(H[2,2])+H[0,2]/r
    return vector(R,[dr(r*vt),dz(r*vt),dz(vr)-dr(vz)]),r*vt
def first_tensor(x,y,w,tw):
    return matrix(R,[[2*((1+tw**2)*x+r*dr(x)+tw*r*dr(y)),2*tw*k*x+k*r*dr(y),r*(dz(x)+tw*dz(y))+E*dr(w)],
      [2*tw*k*x+k*r*dr(y),2*k*k*x,k*r*dz(y)],
      [r*(dz(x)+tw*dz(y))+E*dr(w),k*r*dz(y),2*E*dz(w)]])
checks={}; A=A0+A1*z+A2*z*z; beta=(1-om*om-k*k)/om; gamma=(1+om*om-k*k)/om
x=A/r**2; y=-beta*A/(2*r*r)+B0+b*z; w=W0+v*r*r+h*lr
C0=matrix(R,[[1+om*om,om*k,0],[om*k,k*k,0],[0,0,E]])
Ci=matrix(R,[[1,-om/k,0],[-om/k,(1+om*om)/(k*k),0],[0,0,1/E]])
H=first_tensor(x,y,w,om); HI=-Ci*H*Ci
checks['twist_first_trace']=H.trace()==0
checks['twist_first_trace_square']=(C0*H).trace()==0
checks['twist_first_determinant']=(Ci*H).trace()==0
rc,kc=values(H); ri,kb=values(HI)
checks['twist_all_source_curls']=rc==0
checks['twist_all_inverse_curls']=ri==0
checks['twist_pressure_C']=kc==-k*beta*A2
checks['twist_pressure_B']=kb==gamma*A2/(E*k)
checks['twist_protected_pressure_pair']=gamma-beta==2*om
yz=Ga*lr+B0+b*z+A/r**2
Cz=diagonal_matrix(R,[1,k*k,E]); Czi=diagonal_matrix(R,[1,1/(k*k),1/E])
Hiz=first_tensor(R(0),yz,w,R(0)); HIz=-Czi*Hiz*Czi
checks['zero_twist_first_trace']=Hiz.trace()==0
checks['zero_twist_first_trace_square']=(Cz*Hiz).trace()==0
checks['zero_twist_first_determinant']=(Czi*Hiz).trace()==0
rcz,kcz=values(Hiz); riz,kbz=values(HIz)
checks['zero_twist_all_source_curls']=rcz==0
checks['zero_twist_all_inverse_curls']=riz==0
checks['zero_twist_pressure_C']=kcz==2*k*(Ga+A2)
checks['zero_twist_pressure_B']=kbz==-2*(Ga+A2/E)/k
pmat=matrix(K,[[2*k,2*k],[-2/k,-2/(k*E)]])
checks['zero_twist_protected_pressure_determinant']=pmat.det()==4*(E-1)/E
F0=matrix(K,[[1,0,0],[om,k,0],[0,0,c]])
ez=vector(K,[0,0,1])
checks['matched_pitch_equivariance_forcing']=c*ez-F0*ez==0
# Five canonical fixed-pressure germs. Pure target constants are removed.
def remove_base(f):
    out={}
    for ex,co in R(f).dict().items():
        if ex[1] or ex[2]:continue
        ee=list(ex);ee[0]=0;ee=tuple(ee)
        out[ee]=out.get(ee,K.zero())+co
    return R(f)-R(out)
def radial(f):return r*dr(f)+z*dz(f)
xf=(A0+A1*z)/r**2
yf=remove_base(-beta*(A0+A1*z)/(2*r*r)+b*z)
wf=v*(r*r-1)+h*lr
am=[A0,A1,h,b,v]; weights=[-2,-1,-1,1,1]
def weight(f):return sum(wt*aa*R(f).derivative(aa) for aa,wt in zip(am,weights))
checks['twist_radial_x_weights']=radial(xf)==weight(xf)
checks['twist_radial_phase_weights']=remove_base(radial(yf))==weight(yf)
checks['twist_radial_axial_weights']=remove_base(radial(wf)-wf)==weight(wf)
yzf=remove_base((A0+A1*z)/r**2+b*z)
checks['zero_twist_radial_phase_weights']=remove_base(radial(yzf))==weight(yzf)
# Axial action on germs; constants removed. Pitch has no physical map variation.
checks['twist_axial_only_A1_to_A0_x']=dz(xf)==A1/r**2
checks['twist_axial_only_A1_to_A0_phase']=remove_base(dz(yf))==remove_base(-beta*A1/(2*r*r))
checks['zero_twist_axial_only_A1_to_A0_phase']=remove_base(dz(yzf))==remove_base(A1/r**2)
checks['axial_h_v_action_zero']=dz(wf)==0
assert all(checks.values()),checks
save({'field':K,'ring':R,'twist_H':H,'twist_Hinverse':HI,'zero_H':Hiz,'zero_Hinverse':HIz,
 'pressure_matrix_zero':pmat,'five_weights':weights,'extended_pitch_weight':-1,
 'twist_germs':[xf,yf,wf],'zero_germs':[R(0),yzf,wf],'checks':checks},str(O/'ss_pitch_stable_carrier_certificate.sobj'))
rep={'status':'PASS_FULL_SOURCE_MATCHED_PITCH_LINEAR_CARRIER','checks':checks,'count':len(checks),
 'extended_radial_weights':[-2,-1,-1,1,1,-1],
 'scope':'All positive simple materials; both omega!=0 and omega=0. Completeness/IFT from accompanying proof, not inferred from testing these germs.',
 'seconds':time.time()-start}
(O/'ss_pitch_stable_carrier_reader.json').write_text(json.dumps(rep,indent=2),encoding='utf8')
print(json.dumps(rep,indent=2),flush=True)
