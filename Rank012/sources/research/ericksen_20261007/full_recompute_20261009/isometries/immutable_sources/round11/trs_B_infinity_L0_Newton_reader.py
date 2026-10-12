"""Exact all-pitch Newton support for every 0<ord L/ord tau<1.
Each numerator/denominator is checked at both endpoint Newton faces.
"""
from pathlib import Path
import json
OUT=Path(__file__).resolve().parent
source=(OUT/'trs_B_infinity_L0_leaders.py').read_text()
exec(source[:source.index('results=[]')],globals())
weight=0;wall=0;lrpoly=Q(1);Lunit=ell
def endpointcoef(z,side):
 n=z.numerator();de=z.denominator();ix=R.ngens()-1;pick=min if side==0 else max
 dn=pick(m[ix] for m in n.dict());dd=pick(m[ix] for m in de.dict())
 nn=sum(K(co)*prod(K(R.gen(i))**m[i] for i in range(ix)) for m,co in n.dict().items() if m[ix]==dn)
 den=sum(K(co)*prod(K(R.gen(i))**m[i] for i in range(ix)) for m,co in de.dict().items() if m[ix]==dd)
 return dn-dd,nk(nn/den)
checks={};out=[]
for odd,pol in enumerate(data['frobenius'][1]):
 for press,f in sorted(pol.items()):
  if not f:continue
  key=str((odd,press));ends=[]
  for which,poly in [('numerator',f.numerator()),('denominator',f.denominator())]:
   weight=0;o0,c0=leader(poly);weight=1;o1,c1=leader(poly)
   l0,b0=endpointcoef(c0,0);l1,b1=endpointcoef(c1,1)
   checks[key+which+'_same_L_power']=l0==l1
   checks[key+which+'_same_tau_power']=o0==o1-l1
   checks[key+which+'_same_coefficient']=b0==b1
   ends.append((o0,l0,b0))
  (on,ln,cn),(od,ld,cd)=ends
  item={'U_parity':odd,'pressure':press,'tau_power_with_U':on-od-odd,'L_power':ln-ld,'coefficient_without_U_unit':nk(cn/cd)}
  out.append(item)
  print(odd,press,'order',on-od-odd,'L',ln-ld,'coefficient',item['coefficient_without_U_unit'].factor(),flush=True)
assert all(checks.values()),checks
save({'field':K,'leaders':out,'checks':checks},str(OUT/'trs_B_infinity_L0_Newton_reader.sobj'))
(OUT/'trs_B_infinity_L0_Newton_reader.json').write_text(json.dumps({'status':'PASS','check_count':len(checks),'checks':checks},indent=2))
print('PASS',len(checks),flush=True)
