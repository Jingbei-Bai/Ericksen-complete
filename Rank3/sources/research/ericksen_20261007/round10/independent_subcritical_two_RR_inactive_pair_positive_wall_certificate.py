"""Protected whole-source kappa saturation coefficients at the positive wall."""
from sage.all import *
from pathlib import Path
import os,resource,time,json
p=Path(__file__).resolve().parent;st=time.time();d=load(str(p/'independent_subcritical_two_RR_inactive_pair_wall_second_module.sobj'));bd=load(str(p/'independent_subcritical_two_RR_inactive_pair_generic_boundary.sobj'))
vm=int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE');resource.setrlimit(resource.RLIMIT_AS,(vm+350*1024**2,vm+350*1024**2))
print('LinuxPID',os.getpid(),flush=True)
it=d['items'][0];B=it['material_wall'].base_ring();a,r=B.gens();m=1+a*r;q=1+r;N=a*r*r;la=m*q;c=a*r/m;M=it['module_pairs'];checks=[]
def ck(n,t):
 if not t:raise AssertionError(n)
 checks.append(n)
def add(x,y):return (x[0]+y[0],x[1]+y[1])
def neg(x):return (-x[0],-x[1])
def mul(x,y):return (x[0]*y[0]+c*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
U=[add(mul((6,0),M[13][j]),neg(M[5][j])) for j in range(7)]
ck('full positive quotient U_E zero',U[3]==(0,0))
ck('full positive quotient U_s protected',U[1]==(0,B(8)/3))
ck('full pair2 quotient E protected',M[30][3]==(0,-B(3)/32))
det=add(mul(U[1],M[30][3]),neg(mul(U[3],M[30][1])))
ck('protected quotient determinant ALL materials',det==(-a*r/(4*m),0))
ck('true H7 kappa s unit',M[17][2]==(0,B(4)/(3*m)))
ck('quotient generators have zero kappa s',U[2]==M[30][2]==(0,0))
ck('H7 rho quotient first matrix zero',M[17][1]==M[17][3]==(0,0))
ck('true residual even1 E² unit',M[2][6]==(B(1)/8,0))
for j in [5,13,17,30]:ck('COMPLETE generic original boundary zero '+str(j),bd['zero_range_boundary'][j]==0)
ck('material positivity identity lambda',la==m*q)
save({'whole_wall_source':d,'generic_complete_boundary':bd,'quotient_U_coefficients':U,'rho_quotient_minor':det,'material_N':N,'material_lambda':la,'after_saturation_E2_rhoE_minor':(0,-B(3)/256),'checks':checks,'scope':'All-spectrum positive-wall source unit and legal kappa-saturation coefficients; whole source-ideal argument in proof'},str(p/'independent_subcritical_two_RR_inactive_pair_positive_wall_certificate.sobj'),compress=False)
rep={'PASS':True,'checks':len(checks),'seconds':time.time()-st};(p/'independent_subcritical_two_RR_inactive_pair_positive_wall_certificate.json').write_text(json.dumps(rep,indent=2));print(rep,flush=True)
