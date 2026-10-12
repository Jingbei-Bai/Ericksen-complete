"""Independent whole-corner regularity audit; no relative-rate assumptions."""
from sage.all import *
from pathlib import Path
import json,time
OUT=Path(__file__).resolve().parent;START2=time.time();checks=[]
def ck(label,yes):
 assert yes,label
 checks.append(label)
src=OUT/'trs_algebra_frobenius_infinity_leading.py';ns={'__file__':str(src),'__name__':'null_audit_library'}
exec(src.read_text().split('items=[]')[0],ns)
dat=ns['d'];R=ns['R'];F=ns['F'];lead=ns['poly_lead'];rlead=ns['rational_lead']
v,p,j,I2,P,L=map(F,R.gens());old=dat['source']['field'].ring();Dold=old.gen(6)**2+old.gen(0);Dnew=R.gen(5)**2+R.gen(0)
def multiplicity(poly,div):
 nn=poly;count=0
 while nn:
  qq,rr=nn.quo_rem(div)
  if rr:return count
  nn=qq;count+=1
 return 10**9
units=[R.gen(0),R.gen(2),R.gen(5),R.gen(1)**2-R.gen(0)]
records=[]
for odd,columns in enumerate(dat['frobenius'][1]):
 for mon,x in sorted(columns.items()):
  if not x:continue
  label=str((odd,mon));nn=x.numerator();dd=x.denominator();vn=multiplicity(nn,Dold);vd=multiplicity(dd,Dold)
  ck(label+' full denominator no D factor',vd==0)
  if odd:ck(label+' full odd numerator contains D',vn>=1)
  dn=lead(dd);co=dn[1]
  for name,pol in [('numerator',co.numerator()),('denominator',co.denominator())]:
   rem=pol.quo_rem(Dnew)[1];ck(label+' leading denominator coefficient '+name+' nonzero at null',rem!=0)
   ck(label+' leading denominator coefficient '+name+' only physical units',all(ff in units for ff,m in rem.factor()))
  power,c,parity=rlead(x);power+=odd
  ck(label+' entire normalized source has no negative tau order',power<=-2)
  if mon!=(0,0) or odd:
   ck(label+' term vanishes at complete corner',power<-2 or (power==-2 and vn>=1))
  records.append({'odd':odd,'pressure':list(mon),'twice_B_order':int(power),'D_order':int(vn)})
power,c,parity=rlead(dat['frobenius'][1][0][(0,0)])
i1=(v**3+v*I2-j-P)/v**2
expected=-6*i1/(j*v)*(p*p*L*L+2*v*p*p-v*v)
ck('full generic no-pressure leader',c==expected)
ck('no-pressure parity independent of T sign',parity==0)
num=(c-6*i1*(v-p*p)/j).numerator()
ck('corner unit 6 I1 Q/j exactly',num.quo_rem(Dnew)[1]==0)
report={'status':'PASS','checks':len(checks),'seconds':float(time.time()-START2),
 'scope':'All approach rates to B=infinity, L^2+v=0; both square-root signs, every real pressure and pitch.',
 'corner_unit':'6*I1*(v-p^2)/j',
 'independent_extra_check':'Both numerator and denominator of the B-leading denominator coefficient remain physical units at the null-angle corner.',
 'records':records,'labels':checks}
(OUT/'trs_null_angle_independent_reader.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k not in ['labels','records']},indent=2),flush=True)
