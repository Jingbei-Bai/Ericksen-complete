"""Direct original-R source and all missing-p chart replay for the full-RR square forcing."""
from sage.all import *
from pathlib import Path
import json,time
p=Path(__file__).resolve().parent;st=time.monotonic();src=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'));d=load(str(p/'independent_supercritical_F12_module_identity.sobj'));B=d['ring'];BF=B.fraction_field();A,BB,C,b,c=B.gens();K=B.base_ring();a,r,H,Z=K.gens();q=1+r;alpha=(Z-H)/(Z*(H+1));beta=H*(Z-1)/(Z*(H+1));images=[A,BB,C,c*BB/b,b*C,c*A,alpha*c,beta*b,1];rows=[];checks=[]
def ck(s,z):assert z,s;checks.append(s)
for f in src['full_R_sources']:
 g=BF.zero()
 for ex,co in f.dict().items():
  ee=tuple(map(int,ex))
  if sum(ee[:6])==1 and sum(ee[6:])==2:g+=K(co)*prod(BF(images[j])**ee[j] for j in range(9))
 rows.append(g)
for block in d['modules']:
 mat=block['cofactors'];rr=block['source_rows']
 for j in range(2):ck('original_full_R_module',sum((mat[j,k]*rows[rr[k]] for k in range(2)),BF.zero())==block['targets'][j])
 # The module coefficients only use b and c among direction denominators.
 for z in mat.list():
  for ex in z.denominator().dict():ck('direction_denominator_only_RR_torus',all(ex[j]==0 for j in range(3)))
mis=load(str(p/'independent_RR_torus_missing_p_sources.sobj'));Mfield=mis['field'];ma,mr,mH,mZ,mx,mb,mc=Mfield.gens();km=K.ring().hom([ma,mr,mH,mZ],Mfield)
def pull(f,amps):
 vals=list(amps)+[mb,mc]
 def pp(pol):return sum((km(K(co).numerator())/km(K(co).denominator())*prod(vals[j]**int(ex[j]) for j in range(5)) for ex,co in pol.dict().items()),Mfield.zero())
 return pp(f.numerator())/pp(f.denominator())
for it,amps in zip(mis['items'],[(0,1,mx),(1,0,mx),(1,mx,0)]):
 for j,f in enumerate(rows):ck('original_missing_p_delta_column',pull(f,amps)==it['matrix'][j,3])
from itertools import product
for mask in product([0,1],repeat=3):
 ck('rankD3_support_equivalence',bool((mask[0] or mask[1]) and (mask[1] or mask[2]) and (mask[0] or mask[2]))==(sum(mask)>=2))
save({'checks':checks,'scope':__doc__,'source_files':['independent_delta_subcritical_Q0_full_source.sobj','independent_supercritical_F12_module_identity.sobj','independent_RR_torus_missing_p_sources.sobj']},str(p/'independent_supercritical_F12_source_reader.sobj'),compress=False)
(p/'independent_supercritical_F12_source_reader.json').write_text(json.dumps({'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n');print('PASS',len(checks),'seconds',time.monotonic()-st)
