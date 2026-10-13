"""Short complete source certificate at the critical direction (1,1)."""
from sage.all import *
from pathlib import Path
import json,time
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_delta_critical_point11_schur.sobj'));K=d['field'];a,r,H,Z,B,C=K.gens();q=1+r;items=d['items'];witnesses=[]
targets=[B*B+q,r*C*C-B*B]
for (i,j),target in zip(((0,2),(12,14)),targets):
 f=items[i]['delta']*items[j]['constant']-items[j]['delta']*items[i]['constant']
 quotient=f/target;assert f==target*quotient
 # A valid witness quotient must be a nonzero base polynomial times old
 # units and powers of nonzero B,C, with no remaining amplitude relation.
 fa=quotient.factor();print('WITNESS',i,j,'factor',fa,flush=True)
 witnesses.append({'source_rows':(i,j),'target':target,'minor':f,'quotient':quotient,'factorization':fa})
def canonical(f):
 f=K(f)
 def rr(poly):
  out=K.zero()
  for ee,co in poly.dict().items():
   e=tuple(map(int,ee));out+=co*a**e[0]*r**e[1]*H**e[2]*Z**e[3]*(-q)**(e[4]//2)*(-q/r)**(e[5]//2)*B**(e[4]%2)*C**(e[5]%2)
  return out
 return rr(f.numerator())/rr(f.denominator())
hdelta=canonical(items[0]['delta']);hconstant=canonical(items[0]['constant']);assert hdelta==0 and hconstant!=0
print('CANONICAL_CONSTANT',hconstant.factor(),flush=True)
save({'source_file':'independent_delta_critical_point11_schur.sobj','field':K,'witnesses':witnesses,'canonical_relations':targets,'canonical_H0_delta':hdelta,'canonical_H0_constant':hconstant,'scope':__doc__},str(p/'independent_delta_critical_point11_certificate.sobj'),compress=False)
(p/'independent_delta_critical_point11_certificate.json').write_text(json.dumps({'witness_pairs':[[0,2],[12,14]],'canonical_H0_constant':str(hconstant.factor()),'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n')
