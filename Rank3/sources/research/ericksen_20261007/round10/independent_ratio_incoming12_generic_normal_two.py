"""Complete true second-normal jets on the whole shifted incoming12 family.

Six normals, tangent (l,u,v), all parameters. Full Q/R/H/E sources are used;
range corrections are NOT fixed here. Resume only uncomputed source rows.
"""
from sage.all import *
from pathlib import Path
from itertools import combinations_with_replacement
import os,resource,time,json
p=Path(__file__).resolve().parent;st=time.monotonic()
g=load(str(p/'independent_ratio_p_boundary_second_gauge.sobj'));hd=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'));ed=load(str(p/'independent_ratio_even_derivative_cyclic_direct.sobj'));oq=load(str(p.parent/'round8/rank3_joint_fiber_leading.sobj'));first=load(str(p/'independent_ratio_incoming12_generic_normal.sobj'))
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+900*1024**2,)*2)
K=g['coefficient_field'];a,b,h,z=K.gens();U=first['ring'];l,u,v=U.gens();lam=first['lambda'];A=lam*l
O=PolynomialRing(K,names=tuple('in12n2_p%d'%j for j in range(6))+tuple('in12n2_R%d'%j for j in range(3)));ov=O.gens()
def cp(f):return O({tuple(map(int,e)):K(c) for e,c in f.dict().items()})
exps=[(0,)*6]
for j in range(6):
 e=[0]*6;e[j]=1;exps.append(tuple(e))
for i,j in combinations_with_replacement(range(6),2):
 e=[0]*6;e[i]+=1;e[j]+=1;exps.append(tuple(e))
slots={e:i for i,e in enumerate(exps)};dcache={};factor_cache={}
def factor(dd,der=None):
 fkey=(str(der),tuple(sorted((tuple(map(int,e)),c) for e,c in dd.items())))
 if fkey in factor_cache:return factor_cache[fkey]
 out={}
 for ee,c in dd.items():
  e=tuple(map(int,ee));ne=e[1:5];sz=sum(ne)
  if sz>2:continue
  if der is not None:
   key=(c,str(der))
   if key not in dcache:dcache[key]=c.derivative(der)
   c=dcache[key]
  if not c:continue
  for j in range(min(e[6],2-sz)+1):
   key=(sum(e),ne+(j,));val=K(c)*binomial(e[6],j)*A**(e[6]-j)*l**e[5]*u**e[7]*v**e[8]
   out[key]=out.get(key,U.zero())+val
 result={e:f for e,f in out.items() if f};factor_cache[fkey]=result
 return result
def slot(d,ne,degree):
 tp=degree-d
 if tp<0 or sum(ne)+tp>2:return None
 return slots[ne+(tp,)]
def one(ff,degree):
 ans=[U.zero() for _ in exps]
 for (d,ne),f in ff.items():
  j=slot(d,ne,degree)
  if j is not None:ans[j]+=f
 return vector(U,ans)
def prod2(ff,gg,degree=5):
 ans=[U.zero() for _ in exps]
 for (dl,nl),fl in ff.items():
  for (dr,nr),fr in gg.items():
   ne=tuple(nl[j]+nr[j] for j in range(5));j=slot(dl+dr,ne,degree)
   if j is not None:ans[j]+=fl*fr
 return vector(U,ans)
outpath=p/'independent_ratio_incoming12_generic_normal_two.sobj'
rows=[];records=[];hc={}
if outpath.exists() and outpath.stat().st_size:
 previous=load(str(outpath));records=previous['source_records'];rows=[list(row) for row in previous['coefficient_matrix'].rows()]
while True:
 rowpath=p/('independent_ratio_incoming12_generic_normal_two.row%d.sobj'%len(rows))
 if not rowpath.exists() or not rowpath.stat().st_size:break
 savedrow=load(str(rowpath));assert savedrow['position']==len(rows)
 records.append(savedrow['record']);rows.append(list(savedrow['coefficients']))
done=len(rows);print('RESUME',done,flush=True)
def hf(ref):
 if ref not in hc:hc[ref]=factor(hd['polynomial_pool'][ref])
 return hc[ref]
def add(rec,val):
 assert val[0]==0,('family NOT on source',rec)
 assert list(val[1:7])==list(first['matrix'].row(len(rows)))
 rows.append(list(val[1:]));records.append(rec)
 print(rec,'nonzero',sum(bool(f) for f in val),'sec',round(time.monotonic()-st,1),flush=True)
 rowpath=p/('independent_ratio_incoming12_generic_normal_two.row%d.sobj'%(len(rows)-1))
 temp=p/('independent_ratio_incoming12_generic_normal_two.row%d.partial.sobj'%(len(rows)-1))
 save({'ring':U,'position':len(rows)-1,'record':rec,'coefficients':rows[-1]},str(temp),compress=False);os.replace(temp,rowpath)
Q={qi:cp(oq['quartic_source_rows'][qi]) for qi in (0,1,5,6)};src=[cp(f) for f in g['six_full_sources']]
for pos,(qi,f) in enumerate(Q.items()):
 if pos>=done:add({'kind':'Q','index':qi},one(factor(f.dict()),int(f.total_degree())))
for j,f in enumerate(src):
 if 4+j>=done:add({'kind':'R','index':j},one(factor(f.dict()),4))
for pos,item in enumerate(hd['records']):
 if 10+pos<done:continue
 val=sum((sgn*prod2(hf(le),hf(ri)) for sgn,le,ri in item['full_signed_product_DAG']),vector(U,[0]*len(exps)))
 add({'kind':'H','direction':item['direction'],'index':item['source_position']},val)
hat=[[cp(g['first_gauged_fields'][i,k]) for k in range(9)] for i in range(3)];base=[[cp(g['base_fields'][i,k]) for k in range(2)] for i in range(3)]
for pos,item in enumerate(ed['items']):
 if 28+pos<done:continue
 i=item['direction'];qi=item['Q_source_index'];f=Q[qi];val=vector(U,[0]*len(exps))
 for k in range(9):val+=prod2(factor(hat[i][k].dict()),factor(f.derivative(ov[k]).dict()))
 for k,bp in enumerate((h,z)):val+=prod2(factor(base[i][k].dict()),factor(f.dict(),bp))
 for co,(_,rpos,sgn) in zip(item['global_support_coefficient_dicts'],ed['supports'][i]):
  cpoly=O({tuple(map(int,e)):K(c) for e,c in co.items()})*g['source_unit_coefficients'][i]
  val-=sgn*prod2(factor(cpoly.dict()),factor(src[rpos].dict()))
 add({'kind':'E','direction':i,'index':qi},val)
payload={'coefficient_field':K,'ring':U,'lambda':lam,'normal_positions':first['normal_positions'],'normal_exponents':exps[1:],'source_records':records,'coefficient_matrix':matrix(U,rows),'scope':__doc__,'qualification':'All tangent values preserved. Need true cokernel/range compensation before any no-lift inference.'}
temp=p/'independent_ratio_incoming12_generic_normal_two.partial.sobj';save(payload,str(temp),compress=False);os.replace(temp,outpath)
(p/'independent_ratio_incoming12_generic_normal_two.json').write_text(json.dumps({'sources':len(rows),'first_normal_replayed':True,'seconds_this_run':time.monotonic()-st,'scope':__doc__},indent=2)+'\n')
print('DONE',time.monotonic()-st,flush=True)
