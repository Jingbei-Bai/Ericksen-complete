"""Original-source transport for both Cartesian angular vectors and pressures."""
from sage.all import *
from pathlib import Path
import json
OUT=Path(__file__).absolute().parent;d=load(str(OUT/'homothety_first_jet_reduction.sobj'));R=d['ring'];K=R.base_ring();a,b,s,t=K.gens();P=R.gens()[:9];pi=d['pi'];n=vector(R,d['direction']);fp=d['scaled_pi_fields'];dn=d['direction_fields'];res=matrix(R,3,3,d['residuals']);A=pi-pi.trace()/2*matrix.identity(R,3);rr=res-res.trace()/2*matrix.identity(R,3)
def E(i,f):
 return sum(f.derivative(p)*fp[i][j] for j,p in enumerate(P))+R({m:c.derivative(s) for m,c in f.dict().items()})*dn[i][0]+R({m:c.derivative(t) for m,c in f.dict().items()})*dn[i][1]
omega=d['source_rotation_vector'];pressure=d['longitudinal_pressures']
lam=diagonal_matrix(K,[1,1/a,1/b]);target=(pi.transpose()*lam-(lam*pi).trace()/2*matrix.identity(R,3))*n
def pm(x):return matrix(R,[[x[1,2]-x[2,1],x[0,2],-x[0,1]],[-x[1,2],x[2,0]-x[0,2],x[1,0]],[x[2,1],-x[2,0],x[0,1]-x[1,0]]])
checks={};psres=pm(res)
for i in range(3):
 ww=vector(R,[E(i,z) for z in omega])+A.row(i).cross_product(omega)
 for j in range(3):
  checks[f'pressure_{i}_{j}']=E(i,pressure[j])==psres[j,i]
  checks[f'angular_{i}_{j}']=ww[j]==rr[i,j]
 cvec=lam[i,i]*pi.row(i)-vector(R,[int(i==j)*(lam*pi).trace()/2 for j in range(3)])
 tt=lam*vector(R,[E(i,z) for z in target])+cvec.cross_product(target)
 rt=lam*(lam[i,i]*res.row(i)-vector(R,[int(i==j)*(lam*res).trace()/2 for j in range(3)]))
 for j in range(3):checks[f'target_angular_{i}_{j}']=tt[j]==rt[j]
print(checks,flush=True)
assert all(checks.values()),checks
(OUT/'homothety_observation_transport_reader.json').write_text(json.dumps({'status':'PASS','check_count':len(checks),'checks':checks,'scope':'Whole-source defects of physical angular/pressure transport are exactly linear maps of the retained three first-jet quadratic residuals.'},indent=2))
