"""Exact fresh-to-archive source provenance bridge, without using PASS files."""
from sage.all import *
from pathlib import Path
import json,time,hashlib
out=Path(__file__).absolute().parent;research=out.parent.parent;st=time.monotonic();checks=[]
newfile=out/'fresh_global_source_DAG.sobj';oldfile=research/'round10/independent_ratio_global_R_derivative_degreefall.sobj'
n=load(str(newfile));o=load(str(oldfile));K=n['coefficient_field'];R=PolynomialRing(K,names=n['variable_names']);oldK=o['coefficient_field'];mp=oldK.ring().hom(list(K.gens()),K)
def ck(n,v):
 if not v:raise AssertionError(n)
 checks.append(n)
def cv(c):c=oldK(c);return mp(c.numerator())/mp(c.denominator())
oldpool=[R({tuple(map(int,e)):cv(c) for e,c in d.items()}) for d in o['polynomial_pool']];newpool=[R(d) for d in n['polynomial_pool']]
for key in ['full_source_refs','units','sigma','source_order']:
 if key=='full_source_refs':ck(key,[oldpool[j] for j in o[key]]==[newpool[j] for j in n[key]])
 elif key=='units':ck(key,list(map(cv,o[key]))==list(n[key]))
 else:ck(key,tuple(o[key])==tuple(n[key]))
for key in ['full_first_gauged_field_refs','full_base_field_refs','D_refs','C_refs','L_refs']:
 for i in range(len(n[key])):
  for j in range(len(n[key][i])):ck(key+str((i,j)),oldpool[o[key][i][j]]==newpool[n[key][i][j]])
for i,(oldrec,newrec) in enumerate(zip(o['records'],n['records'])):
 od=oldrec['full_signed_product_DAG'];nd=newrec['full_signed_product_DAG']
 ck('DAG term count '+str(i),len(od)==len(nd))
 for j,((so,lo,ro),(sn,ln,rn)) in enumerate(zip(od,nd)):
  ck('complete DAG summand '+str((i,j)),so==sn and oldpool[lo]==newpool[ln] and oldpool[ro]==newpool[rn])
report=dict(status='PASS',checks=checks,count=len(checks),seconds=time.monotonic()-st,fresh_sha256=hashlib.sha256(newfile.read_bytes()).hexdigest(),historical_sha256=hashlib.sha256(oldfile.read_bytes()).hexdigest(),
 scope='Every coefficient of all six actual R sources, all original/gauged field components, and every signed factor in all18 complete derivative DAGs agrees with the fresh Section2 reconstruction. This transfers raw-source provenance to subsequent archived witness replays; it does not validate those witnesses without their arithmetic checks.')
(out/'fresh_boundary_source_bridge.json').write_text(json.dumps(report,indent=2));print(report['status'],len(checks),report['seconds'],flush=True)
