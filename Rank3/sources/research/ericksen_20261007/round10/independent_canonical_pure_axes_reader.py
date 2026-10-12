from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_ratio_canonical_partial_zero_probe.sobj'));old=load(str(p/'independent_ratio_canonical_amplitude_normal.sobj'));K=old['coefficient_field'];a,r=old['parameters'];R=PolynomialRing(K,names=('axis_amplitude','axis_unused'));amp=R.gen(0);checks=[];records=[]
def ck(n,v):
 if not bool(v):raise AssertionError(n)
 checks.append(n)
for it in d['items']:
 ax=int(it['axis']);rs=list(map(int,it['pivot_rows']));cs=list(map(int,it['pivot_columns']));ck('old_sources_only_'+str(ax),len(rs)==6 and max(rs)<28 and cs==list(range(6)));N=matrix(R,6,6)
 for ii,i in enumerate(rs):
  for j in range(6):
   f=old['matrix'][i,j];val=R.zero()
   for e,c in f.dict().items():
    e=tuple(map(int,e))
    if all(e[k]==0 for k in range(3) if k!=ax):val+=K(c)*amp**e[ax]
   N[ii,j]=val;oldval=R({tuple(map(int,e)):K(c) for e,c in it['matrix'][i,j].dict().items()});ck('axis_original_%d_%d_%d'%(ax,i,j),val==oldval)
 minor=N.det();saved=R({tuple(map(int,e)):K(c) for e,c in it['minor'].dict().items()});ck('axis_minor_'+str(ax),minor==saved);unit=K(minor.monomial_coefficient(amp**2));ck('pure_amplitude_square_'+str(ax),minor==unit*amp**2);ck('axis_material_unit_'+str(ax),unit!=0 and all(any(K(f)/u in QQ for u in (r,r+1)) for f,k in unit.numerator().factor()) and all(any(K(f)/u in QQ for u in (r,r+1)) for f,k in unit.denominator().factor()));records.append({'axis':ax,'rows':rs,'minor':str(minor.factor())})
report={'status':'PASS','count':len(checks),'seconds':time.monotonic()-st,'records':records,'scope':'Three nonzero pure-amplitude axes have full6 normal rank from original old28 rows only, all material spectra with r,q protected. Cyclic base; source/cone bridge needed for generic-base geometric use.'};save({'checks':checks,'records':records},str(p/'independent_canonical_pure_axes_reader.sobj'));(p/'independent_canonical_pure_axes_reader.json').write_text(json.dumps(report,indent=2)+'\n');print(report,flush=True)
