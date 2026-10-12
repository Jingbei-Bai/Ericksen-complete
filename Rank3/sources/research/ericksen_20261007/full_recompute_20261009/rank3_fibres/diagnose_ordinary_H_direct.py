"""Independent complete product-rule coefficients on the incoming family."""
from sage.all import *
from pathlib import Path
import json,time,hashlib
O=Path(__file__).absolute().parent;N=O.parent/'rank3_ordinary/rebuild/round10';OLD=O.parent.parent/'round10';st=time.monotonic()
receipt=O.parent/'rank3_ordinary/job_reports/independent_ratio_incoming12_generic_normal_two.json'
rec=json.loads(receipt.read_text());assert rec['status']=='PASS'
inputs=[it for it in rec['imports'] if Path(it['path']).name=='independent_ratio_global_R_derivative_degreefall.sobj']
assert inputs and len({it['path'] for it in inputs})==1
hp=Path(inputs[0]['path']).resolve();hh=hashlib.sha256(hp.read_bytes()).hexdigest()
assert all(it['sha256']==hh for it in inputs)
hd=load(str(hp))
np=N/'independent_ratio_incoming12_generic_nullwall_all.sobj';new=load(str(np))
U=new['ring'];K=U.base_ring();a,b,h,z=K.gens();l=U.gen();mod=new['relation'];kap=new['kappa'];cache={};checks=[]; entries=[]
# Actual monomial derivative evaluation at p=(1,0,0,0,0,l), RR=(-2l,0,-2).
# No normal-jet arrays enter this reconstruction.
orders=[(0,0,0),(1,0,0),(0,1,0),(0,0,1),(1,0,1),(0,1,1)]
def parts(ref):
 if ref in cache:return cache[ref]
 out={}
 for ee,c in hd['polynomial_pool'][ref].items():
  e=tuple(map(int,ee));d=sum(e)
  if e[2] or e[3] or e[7]:continue
  if d not in out:out[d]=[U.zero() for _ in orders]
  for j,(aa,bb,cc) in enumerate(orders):
   if e[1]!=aa or e[4]!=bb or e[6]<cc:continue
   out[d][j]+=K(c)*factorial(aa)*factorial(bb)*prod(e[6]-s for s in range(cc))*(-2*l)**(e[6]-cc)*l**e[5]*(-2)**e[8]
 cache[ref]=out;return out
for row in [11,13,23,25]:
 result=[U.zero(),U.zero()]
 for sgn,le,ri in hd['records'][row-10]['full_signed_product_DAG']:
  ff=parts(le);gg=parts(ri)
  for d,f in ff.items():
   for D,g in gg.items():
    if d+D==5:
     kp=(f[5]-f[4])*g[0]+(f[2]-f[1])*g[3]+f[3]*(g[2]-g[1])+f[0]*(g[5]-g[4])
     ep=f[4]*g[0]+f[1]*g[3]+f[3]*g[1]+f[0]*g[4]
     result[0]+=sgn*kp;result[1]+=sgn*kap*l*ep
    elif d+D==4:result[1]+=sgn*(f[3]*g[0]+f[0]*g[3])
 for col,f in zip([1,4],result):
  f=f%mod;nv=new['items'][row]['rhs'][col]
  assert f==nv,(row,col,'complete product-rule coefficient')
  checks.append({'row':row,'column':col,'equals_fresh':True})
  entries.append({'row':row,'column':col,'source_value':f})
assert len(entries)==8
payload={'ring':U,'relation':mod,'entries':entries,'H_input':{'path':str(hp),'sha256':hh},'nullwall_sha256':hashlib.sha256(np.read_bytes()).hexdigest()}
save(payload,str(O/'ordinary_H_direct_coefficients.sobj'),compress=False)
report={'status':'PASS','inputs':[payload['H_input'],{'path':str(np),'sha256':payload['nullwall_sha256']}],
        'checked_coefficients':8,'checks':checks,'seconds':time.monotonic()-st}
(O/'ordinary_H_direct_diagnostic.json').write_text(json.dumps(report,indent=2)+'\n')
print('PASS eight complete H coefficients',flush=True)
