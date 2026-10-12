"""Exact symmetric-variable test for the extreme-axis resultant guard.

1<s<t corresponds to y=s*t>1 and 0<z=(s-t)^2<(y-1)^2.
"""
from sage.all import *
from pathlib import Path
import json,time,signal,sys,resource
O=Path(__file__).resolve().parent;st=time.time();signal.alarm(120)
resource.setrlimit(resource.RLIMIT_CORE,(0,0))
vm=int(Path('/proc/self/statm').read_text().split()[0])*4096
resource.setrlimit(resource.RLIMIT_AS,(vm+1400*1024**2,vm+1400*1024**2))
tag=sys.argv[1] if len(sys.argv)>1 else 'plus';d=load(str(O/f'rr_null_physical_axis_resultant_sign_{tag}.sobj'))
P=d['ring'];s,t=P.gens();re=d['parts']['real']['polynomial'];imfac=d['parts']['imag']['factorization'];im=max(imfac,key=lambda z:z[0].total_degree())[0]
Q=PolynomialRing(QQ,names=('z','y'));z,y=Q.gens()
power=[Q(2),z+2*y]
for j in range(2,20):power.append((z+2*y)*power[-1]-y*y*power[-2])
def symmetric(f):
 assert f(s=t,t=s)==f
 out=Q(0)
 for (i,j),co in f.dict().items():
  if i<j:continue
  if i==j:out+=co*y**j
  else:
   assert (i-j)%2==0
   out+=co*y**j*power[(i-j)//2]
 assert Q.hom([(s-t)**2,s*t],P)(out)==f
 return out
rr=symmetric(re);ss=symmetric(im)
Y=PolynomialRing(QQ,'y');yv=Y.gen();Z=PolynomialRing(Y,'z');zv=Z.gen()
def convert(f):
 return Z([Y({j:co for (i,j),co in f.dict().items() if i==k}) for k in range(f.degree(z)+1)])
fr=convert(rr);fs=convert(ss)
print('degrees',fr.degree(),fs.degree(),time.time()-st,flush=True)
res=fr.resultant(fs);fac=res.factor()
save({'ring':Q,'real':rr,'imag_core':ss,'resultant':res,'factorization':fac,'source':str(O/f'rr_null_physical_axis_resultant_sign_{tag}.sobj')},str(O/f'rr_null_physical_axis_symmetric_resultant_{tag}.sobj'))
report={'status':'EXACT_MATERIAL_GUARD_ELIMINATION_NOT_YET_WHOLE_POSITIVITY','tag':tag,'degrees':[fr.degree(),fs.degree()],
 'factors':[{'degree':int(f.degree()),'multiplicity':int(m),'polynomial':str(f),'positive_shift_negative_coeffs':int(sum(v<0 for v in f(yv+1).coefficients()))} for f,m in fac],
 'seconds':time.time()-st}
print('factor degrees',[(f.degree(),m) for f,m in fac],time.time()-st,flush=True)
(O/f'rr_null_physical_axis_symmetric_resultant_{tag}.json').write_text(json.dumps(report,indent=2,default=str),encoding='utf8')
