"""Exact endpoint-Newton certificate for ALL 0<ord(L)/ord(tau)<1.

The same surviving monomial realizes the minimum at both endpoints.
Every other term is above at both endpoints, hence above for every
intermediate rational/real weight. No sampled slopes stand in for this.
"""
from sage.all import *
from pathlib import Path
import time,resource,json,sys
OUT=Path(__file__).resolve().parent;START=time.time();LWEIGHT=0
source=(OUT/'tr_pure_rotation_B_infinity_leaders.py').read_text()
exec(source[source.index("NF=QuadraticField"):source.index("vm=int")],globals())
vm=int(Path('/proc/self/statm').read_text().split()[0])*4096
resource.setrlimit(resource.RLIMIT_AS,(vm+650*1024**2,vm+650*1024**2))
data=load(str(OUT/'tr_pure_rotation_L_frobenius.sobj'))['coefficients']
def endpointcoef(z,side):
 n=z.numerator();den=z.denominator();ix=R.ngens()-1
 pick=min if side==0 else max
 dn=pick(mon[ix] for mon in n.dict());dd=pick(mon[ix] for mon in den.dict())
 nn=sum(K(c)*prod(K(R.gen(i))**mon[i] for i in range(ix)) for mon,c in n.dict().items() if mon[ix]==dn)
 de=sum(K(c)*prod(K(R.gen(i))**mon[i] for i in range(ix)) for mon,c in den.dict().items() if mon[ix]==dd)
 return dn-dd,nn/de
checks={};out={}
for label,f in data.items():
 results=[]
 for which,poly in [('numerator',f.numerator()),('denominator',f.denominator())]:
  LWEIGHT=0;o0,c0=leader(poly);LWEIGHT=1;o1,c1=leader(poly)
  l0,b0=endpointcoef(c0,0);l1,b1=endpointcoef(c1,1)
  checks[label+'_'+which+'_same_extremal_L_power']=l0==l1
  checks[label+'_'+which+'_same_tau_exponent']=o0==o1-l1
  checks[label+'_'+which+'_same_exact_coefficient']=b0==b1
  results.append((o0,l0,b0))
 (on,ln,cn),(od,ld,cd)=results
 extra=1 if label in ['EC','ED'] else 0
 out[label]={'tau_power_with_E':on-od+extra,'L_power':ln-ld,'coefficient_without_E_unit':cn/cd}
 print(label,out[label]['tau_power_with_E'],out[label]['L_power'],(cn/cd).factor(),flush=True)
I1symbol=K(P.parent().gen(0)) if False else None
checks['constant_leader']=out['1']['coefficient_without_E_unit']==6*(v**3+v*I2-j-P)/(v*j)
checks['EC_leader']=out['EC']['coefficient_without_E_unit']==-7*ii*w/j
checks['constant_weight_is_2']=out['1']['tau_power_with_E']==2 and out['1']['L_power']==0
checks['EC_weight_is_2_minus_alpha']=out['EC']['tau_power_with_E']==2 and out['EC']['L_power']==-1
checks['CC_weight_is_4_minus_2alpha']=out['CC']['tau_power_with_E']==4 and out['CC']['L_power']==-2
checks['CD_weight_is_3_minus_alpha']=out['CD']['tau_power_with_E']==3 and out['CD']['L_power']==-1
checks['ED_weight_is_3']=out['ED']['tau_power_with_E']==3 and out['ED']['L_power']==0
checks['DD_weight_is_5_minus_alpha']=out['DD']['tau_power_with_E']==5 and out['DD']['L_power']==-1
assert all(checks.values()),checks
save({'field':K,'leaders':out,'checks':checks},str(OUT/'tr_pure_rotation_B_infinity_L0_Newton_reader.sobj'))
summary={'status':'PASS','check_count':len(checks),'checks':checks,'scope':'All rational ramification weights 0<alpha<1 by identical extremal monomials at BOTH endpoints, not interpolation of sampled slopes.'}
(OUT/'tr_pure_rotation_B_infinity_L0_Newton_reader.json').write_text(json.dumps(summary,indent=2))
print('PASS',len(checks),round(time.time()-START,2),flush=True)
