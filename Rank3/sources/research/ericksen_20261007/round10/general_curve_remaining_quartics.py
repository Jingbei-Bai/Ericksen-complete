"""Three additional original quartics on the actual noncharacteristic graph.

Use true H-velocities so the nine-variable substitution is short. No
fixed spectrum or curve profile is chosen. Sources Q1,Q5,Q6 are retained.
"""
from pathlib import Path
import time,json
from sage.all import QQ,PolynomialRing,load,save,lcm,gcd
out=Path(__file__).resolve().parent;start=time.monotonic()
d=load(str(out.parent/'round8/rank3_joint_fiber_leading.sobj'))
v=load(str(out/'joint_fields_direct.sobj'))
C=PolynomialRing(QQ,names=('a','r','H','Z','S','N'));K=C.fraction_field()
a,r,H,Z,S,N=K.gens();q=1+r;b=a*q/(1+a*r)
R=PolynomialRing(K,names=('w1','m1','w2','m2','w3','m3'));w1,m1,w2,m2,w3,m3=R.gens()
images=[w3,w2,((1+H)*S-Z)*m3/(2*Z)+H*w3,w1,
 -(1+H*S-Z)*m2/(2*H)-Z*w2,-((1-S)*m1/2+Z*w1)/H,
 (Z-H*S)*m1/(H*Z),S*m2/Z,-m3/H]
O=PolynomialRing(K,names=d['ring'].variable_names());ev=O.hom(images,R)
def cc(c):return K(c.numerator()(a,b,H/r,q*Z/r))/K(c.denominator()(a,b,H/r,q*Z/r))
def pull(p):return ev(O({tuple(m):cc(c) for m,c in p.dict().items()}))
checks={}
for i in range(3):
 checks['h_velocity_'+str(i)]=pull(v['normalized_pi_fields'][i,3])==R.gen(2*i+1)/r
 checks['z_velocity_'+str(i)]=pull(v['normalized_pi_fields'][i,7])==q*S*R.gen(2*i+1)/r
assert all(checks.values())
full=load(str(out/'general_curve_dimensionless_full.sobj'))
checks['Q0_independent_replay']=pull(d['quartic_source_rows'][0])==R(full['rows'][3])
assert checks['Q0_independent_replay']
rows=[];leading=[];records=[]
for idx in [1,5,6]:
 f=pull(d['quartic_source_rows'][idx]);rows.append(f)
 top=R({m:c for m,c in f.dict().items() if sum(m)==f.total_degree()});leading.append(top)
 hist={};patterns=set()
 for m in f.dict():hist[sum(m)]=hist.get(sum(m),0)+1
 for m in top.dict():patterns.add(tuple(m[2*j]+m[2*j+1] for j in range(3)))
 rec={'source':idx,'terms':len(f.dict()),'histogram':hist,'highest_pair_degrees':sorted(patterns),'seconds':time.monotonic()-start}
 records.append(rec);print(json.dumps(rec),flush=True)
 payload={'ring':R,'parameters':K.gens(),'images':images,'rows':rows,'leading':leading,'records':records,'checks':checks,
  'source':'round8/rank3_joint_fiber_leading.sobj','scope':'Additional full necessary quartics; no dimension or rank conclusion.'}
 tmp=out/'general_curve_remaining_quartics_next.sobj';save(payload,str(tmp));tmp.replace(out/'general_curve_remaining_quartics.sobj')
(out/'general_curve_remaining_quartics.json').write_text(json.dumps({'records':records,'checks':checks,'elapsed':time.monotonic()-start},indent=2)+'\n')
