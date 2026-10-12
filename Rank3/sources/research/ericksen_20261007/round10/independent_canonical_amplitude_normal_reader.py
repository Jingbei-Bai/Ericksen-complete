from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();checks={}
def ck(n,v):
 if not v:raise AssertionError(n)
 checks[n]=True
saved=load(str(p/'independent_ratio_canonical_amplitude_normal.sobj'));whole=load(str(p/'independent_ratio_RR_canonical_whole_cone.sobj'));d=load(str(p/'joint_fields_direct.sobj'));fd=load(str(p/'independent_ratio_weighted_degreefall.sobj'));red=load(str(p/'independent_ratio_diagonal_Q0_reduction.sobj'));sec=load(str(p/'independent_ratio_RR_cyclic_full_matrix.sobj'));tan=load(str(p/'independent_ratio_RR_canonical_tangent_obstruction.sobj'));K=sec['coefficient_field'];a,r=sec['parameters'];q=1+r;A=-q/r;B=-q;ok=d['ring'].base_ring();h,z=d['parameters'][2:];mp=ok.ring().hom(sec['old_parameter_map'],K);R=saved['amplitude_ring'];tt=R.gens();wt=(1,)*6+(2,)*3;cache={};dcache={}
def cv(c):
 c=ok(c)
 if c not in cache:cache[c]=mp(c.numerator())/mp(c.denominator())
 return cache[c]
def dc(c,j):
 key=(ok(c),j)
 if key not in dcache:dcache[key]=cv(key[0].derivative((h,z)[j]))
 return dcache[key]
def raw(f):return {tuple(map(int,e)):ok(c) for e,c in f.dict().items() if c}
def clean(d):return {e:c for e,c in d.items() if c}
def weight(e):return sum(x*y for x,y in zip(e,wt))
def add(out,e,c):out[e]=out.get(e,K.zero())+c
ss=[raw(fd['Q0'])]+[raw(it['full']) for it in fd['items']];hat=[];base=[]
for i in range(3):
 row=[raw(d['fields'][i,j]) for j in range(9)]
 for e,c in ss[0].items():row[6+i][e]=row[6+i].get(e,ok.zero())-red['items'][i]['kappa']*c
 row[6+i]=clean(row[6+i]);hat.append(row);base.append([raw(d['base_fields'][i,j]) for j in range(2)])
lower=[clean({e:cv(c) for e,c in row.items() if weight(e)==4}) for row in ss]
for j,row in enumerate(ss):
 for i in range(3):
  out={}
  for e,c in row.items():
   for l,field in enumerate(hat[i]):
    if e[l]==0:continue
    for u,v in field.items():
     ne=tuple(e[k]+u[k]-(1 if k==l else 0) for k in range(9))
     if weight(ne)==6:add(out,ne,e[l]*cv(c)*cv(v))
   for l,field in enumerate(base[i]):
    for u,v in field.items():
     ne=tuple(e[k]+u[k] for k in range(9))
     if weight(ne)==6:add(out,ne,dc(c,l)*cv(v))
  lower.append(clean(out))
 print('subtop source',j,'seconds',round(time.monotonic()-st,2),flush=True)
rows=[]
for idx,(top,low) in enumerate(zip(whole['complete_weighted_top_sources'],lower)):
 ck('all_true_tau_source_'+str(idx),low==clean({tuple(map(int,e)):K(c) for e,c in saved['complete_tau_coefficient_sources'][idx].items()}));pa=tan['row_parities'][idx];unit=sec['row_scales'][idx];out=[R.zero() for _ in range(6)]
 # Alternative reduction: c0^(1-n)=(RR1*RR2)^((1-n)/2).
 # Absorb this as integer powers of each root before squaring, rather than
 # using creator's parity/J formula. c0 is fixed when differentiating X,Y.
 for tp,source in ((0,top),(1,low)):
  for ee,c in source.items():
   e=tuple(map(int,ee));n=sum(e[:6])+tp;ck('normal_odd_'+str(idx)+'_'+str(tp)+'_'+str(e),n%2==1);rootexp=(e[6]+e[3]+e[5]+tp-pa[0],e[7]+e[1]+e[4]+tp-pa[1]);adjust=(n-1)//2;ex=(rootexp[0]-adjust,rootexp[1]-adjust);ck('fourth_root_reduction_even_'+str(idx)+'_'+str(tp)+'_'+str(e),ex[0]%2==ex[1]%2==0);value=K(c)/unit*A**(ex[0]//2)*B**(ex[1]//2);ke=(e[0]+e[5],e[1]+e[3],e[2]+e[4]);mon=prod(tt[j]**ke[j] for j in range(3))
   if tp==1:out[3]+=value*mon;continue
   for j in range(3):
    if e[j]:out[j]+=e[j]*value*mon/tt[j]
   out[4]+=-QQ(rootexp[0])/2*value*mon;out[5]+=-QQ(rootexp[1])/2*value*mon
 rows.append([R(v) for v in out])
M=matrix(R,rows);ck('all_normal_entries',all(M[i,j]==saved['matrix'][i,j] for i in range(28) for j in range(6)));at0=matrix(K,28,6,[f(0,0,0) for f in M.list()]);ck('origin_matches_four_columns',at0.matrix_from_columns([0,1,2,3])==tan['constant_matrix'].matrix_from_columns([0,1,2,6]));ck('origin_RR_columns_zero',all(v==0 for v in at0.column(4)) and all(v==0 for v in at0.column(5)))
report={'status':'PASS','count':len(checks),'checks':checks,'seconds':time.monotonic()-st,'scope':'Complete cyclic-base whole-amplitude normal Jacobian from original f4/E6 tau coefficients and full already-audited f5/E7 tops; independent fourth-root exponent normalization, retaining derivative of the RR product and all h,z source coefficients. All amplitudes and material spectra symbolic. No rank/no-lift assertion away from kernel origin, and no generic-base spreading assertion.'}
(p/'independent_canonical_amplitude_normal_reader.json').write_text(json.dumps(report,indent=2)+'\n');save({'checks':checks},str(p/'independent_canonical_amplitude_normal_reader.sobj'));print({k:v for k,v in report.items() if k!='checks'},flush=True)
