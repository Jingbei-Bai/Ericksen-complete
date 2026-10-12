from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();checks={}
def ck(n,v):
 if not v:raise AssertionError(n)
 checks[n]=True
d=load(str(p/'joint_fields_direct.sobj'));f=load(str(p/'independent_ratio_weighted_degreefall.sobj'));red=load(str(p/'independent_ratio_diagonal_Q0_reduction.sobj'));g=load(str(p/'independent_ratio_generic_canonical_normal_28.sobj'));saved=load(str(p/'independent_ratio_generic_canonical_whole_cone.sobj'));K=g['coefficient_field'];a,r,H,Z=g['parameters'];q=1+r;ok=d['ring'].base_ring();h,z=d['parameters'][2:];mp=ok.ring().hom(g['parameter_map'],K);sq=g['RR_squares'];ker=g['kernel_rows'];cols=g['column_RR_exponents'];wt=(1,)*6+(2,)*3;T=PolynomialRing(K,names=('indWholeK0','indWholeK1','indWholeK2'));ks=T.gens();cache={};dcache={}
def cv(c):
 c=ok(c)
 if c not in cache:cache[c]=mp(c.numerator())/mp(c.denominator())
 return cache[c]
def dc(c,j):
 key=(ok(c),j)
 if key not in dcache:dcache[key]=cv(key[0].derivative((h,z)[j]))
 return dcache[key]
def clean(row):return {tuple(map(int,e)):c for e,c in row.items() if c}
def original(poly):return {tuple(map(int,e)):ok(c) for e,c in poly.dict().items() if c}
def weight(e):return sum(a*b for a,b in zip(e,wt))
pmap=[]
for j in range(6):
 nz=[i for i in range(3) if ker[i,j]];ck('single_kernel_support_'+str(j),len(nz)==1);pmap.append((nz[0],ker[nz[0],j]))
def evaluate(row,base_derivative=None):
 buckets=[{} for _ in range(4)]
 for e,c in row.items():
  re=list(e[6:]);ke=[0,0,0];value=cv(c) if base_derivative is None else dc(c,base_derivative)
  for j in range(6):
   i,fac=pmap[j];ke[i]+=e[j];value*=fac**e[j]
   for l in range(3):re[l]+=e[j]*cols[j][l]
  value*=prod(sq[j]**(re[j]//2) for j in range(3));code=(re[0]%2)+2*(re[1]%2);ke=tuple(ke);buckets[code][ke]=buckets[code].get(ke,K.zero())+value
 return tuple(T(clean(v)) for v in buckets)
def mul(u,v):
 ans=[T.zero() for _ in range(4)]
 for i,fi in enumerate(u):
  if not fi:continue
  for j,gj in enumerate(v):
   if not gj:continue
   factor=(sq[0] if (i&j)&1 else 1)*(sq[1] if (i&j)&2 else 1);ans[i^j]+=factor*fi*gj
 return tuple(ans)
def plus(u,v):return tuple(x+y for x,y in zip(u,v))
def compare_stored(element,sd):
 got={}
 for code,poly in enumerate(element):
  for e,c in poly.dict().items():got[tuple(map(int,e))+(code&1,(code>>1)&1)]=K(c)
 return got=={tuple(map(int,e)):K(c) for e,c in sd.items() if c}
ss=[original(f['Q0'])]+[original(it['full']) for it in f['items']];tops=[]
for j,row in enumerate(ss):
 ck('complete_source_weight_'+str(j),all(weight(e)<=5 for e in row));top={e:c for e,c in row.items() if weight(e)==5};tops.append(top);ck('raw_top_source_'+str(j),top==clean(saved['raw_original_top_sources'][j]))
hat=[];base=[];hatv=[];basev=[]
for i in range(3):
 row=[]
 for j in range(9):
  vals=original(d['fields'][i,j])
  if j==6+i:
   for e,c in ss[0].items():vals[e]=vals.get(e,ok.zero())-red['items'][i]['kappa']*c
   vals=clean(vals)
  ck('complete_hat_weight_'+str(i)+'_'+str(j),all(weight(e)<=wt[j]+2 for e in vals));top={e:c for e,c in vals.items() if weight(e)==wt[j]+2};ck('raw_hat_source_'+str(i)+'_'+str(j),top==clean(saved['raw_original_principal_hat_fields'][i][j]));row.append(top)
 hat.append(row);hatv.append([evaluate(v) for v in row]);br=[]
 for j in range(2):
  vals=original(d['base_fields'][i,j]);ck('base_weight_'+str(i)+'_'+str(j),all(weight(e)<=2 for e in vals));top={e:c for e,c in vals.items() if weight(e)==2};ck('raw_base_source_'+str(i)+'_'+str(j),top==clean(saved['raw_original_principal_base_fields'][i][j]));br.append(top)
 base.append(br);basev.append([evaluate(v) for v in br]);print('fields',i,'seconds',round(time.monotonic()-st,2),flush=True)
 for j in range(9):ck('hat_quotient_source_'+str(i)+'_'+str(j),compare_stored(hatv[i][j],saved['principal_hat_field_restrictions'][i][j]))
 for j in range(2):ck('base_quotient_source_'+str(i)+'_'+str(j),compare_stored(basev[i][j],saved['principal_base_field_restrictions'][i][j]))
values=[]
for j,row in enumerate(tops):
 val=evaluate(row);ck('whole_old_cone_'+str(j),all(f==0 for f in val));values.append(val)
for j,row in enumerate(tops):
 derivatives=[]
 for l in range(9):
  dr={}
  for e,c in row.items():
   if not e[l]:continue
   ne=list(e);ne[l]-=1;dr[tuple(ne)]=e[l]*c
  val=evaluate(dr);derivatives.append(val);ck('fiber_derivative_'+str(j)+'_'+str(l),compare_stored(val,saved['derivative_evaluations'][j]['fiber'][l]))
 bd=[evaluate(row,l) for l in range(2)]
 for l in range(2):ck('base_derivative_'+str(j)+'_'+str(l),compare_stored(bd[l],saved['derivative_evaluations'][j]['base'][l]))
 for i in range(3):
  val=tuple(T.zero() for _ in range(4))
  for l in range(9):val=plus(val,mul(derivatives[l],hatv[i][l]))
  for l in range(2):val=plus(val,mul(bd[l],basev[i][l]))
  ck('whole_genuine_derivative_'+str(j)+'_'+str(i),all(f==0 for f in val));values.append(val)
 print('whole derivatives',j,'seconds',round(time.monotonic()-st,2),flush=True)
ck('all28',len(values)==28);report={'status':'PASS','count':len(checks),'checks':checks,'seconds':time.monotonic()-st,'scope':'All spectral and base parameters retained. Independently rebuilt original tops and principal full-Q0-gauged fields, then contracted exact derivatives in a four-basis quadratic extension algebra, including all original h,z coefficient derivatives. All28 complete tau-zero cone restrictions vanish. No J=1, no fixed spectrum, no truncation at cubic normal order.'}
(p/'independent_generic_RR_whole_cone_reader.json').write_text(json.dumps(report,indent=2)+'\n');save({'checks':checks},str(p/'independent_generic_RR_whole_cone_reader.sobj'));print({k:v for k,v in report.items() if k!='checks'},flush=True)
