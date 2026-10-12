from sage.all import *
from pathlib import Path
import json,time
p=Path(__file__).resolve().parent;st=time.monotonic();checks={}
def ck(n,v):
 if not v:raise AssertionError(n)
 checks[n]=True
d=load(str(p/'joint_fields_direct.sobj'));f=load(str(p/'independent_ratio_weighted_degreefall.sobj'));red=load(str(p/'independent_ratio_diagonal_Q0_reduction.sobj'));sec=load(str(p/'independent_ratio_RR_cyclic_full_matrix.sobj'));saved=load(str(p/'independent_ratio_RR_canonical_whole_cone.sobj'));cubic=load(str(p/'independent_ratio_RR_canonical_cubic_sources.sobj'))
K0=d['ring'].base_ring();aa,bb,hh,zz=d['parameters'];K=sec['coefficient_field'];a,r=sec['parameters'];q=1+r;mp=K0.ring().hom(sec['old_parameter_map'],K);wt=(1,)*6+(2,)*3
cache={};dcache={}
def ev(v):
 v=K0(v)
 if v not in cache:cache[v]=mp(v.numerator())/mp(v.denominator())
 return cache[v]
def deriv(v,i):
 key=(K0(v),i)
 if key not in dcache:dcache[key]=ev(key[0].derivative((hh,zz)[i]))
 return dcache[key]
def pdict(x):return {tuple(map(int,e)):K0(c) for e,c in x.dict().items() if c}
def weight(e):return sum(x*y for x,y in zip(e,wt))
def compact(x):return {e:c for e,c in x.items() if c}
def add(out,e,c):out[e]=out.get(e,K.zero())+c
ss=[pdict(f['Q0'])]+[pdict(i['full']) for i in f['items']];qzero=ss[0];fields=[];base=[]
for i in range(3):
 fr=[pdict(d['fields'][i,j]) for j in range(9)];kap=red['items'][i]['kappa']
 for e,c in qzero.items():fr[6+i][e]=fr[6+i].get(e,K0.zero())-kap*c
 fr[6+i]=compact(fr[6+i]);fields.append(fr);base.append([pdict(d['base_fields'][i,j]) for j in range(2)])
 for j,row in enumerate(fr):ck('field_weight_'+str(i)+'_'+str(j),all(weight(e)<=wt[j]+2 for e in row))
 ck('base_weight_'+str(i),all(weight(e)<=2 for row in base[-1] for e in row))
raw=[compact({e:ev(c) for e,c in row.items() if weight(e)==5}) for row in ss]
for j,row in enumerate(ss):
 ck('source_weight_'+str(j),all(weight(e)<=5 for e in row))
 for i in range(3):
  out={}
  for e,c in row.items():
   for k,fr in enumerate(fields[i]):
    if e[k]==0:continue
    for u,v in fr.items():
     ne=tuple(e[t]+u[t]-(1 if t==k else 0) for t in range(9))
     if weight(ne)==7:add(out,ne,e[k]*ev(c)*ev(v))
   for k,fr in enumerate(base[i]):
    for u,v in fr.items():
     ne=tuple(e[t]+u[t] for t in range(9))
     if weight(ne)==7:add(out,ne,deriv(c,k)*ev(v))
  raw.append(compact(out))
 print('source rows',7+3*(j+1),'seconds',round(time.monotonic()-st,2),flush=True)
ck('source_count',len(raw)==28)
# Independent evaluation in K[R1,R2]/(R1^2+q/r,R2^2+q), R3=1.
# This keeps all four complex embeddings and never replaces J=R1*R2 by 1.
records=[]
for i,row in enumerate(raw):
 old={tuple(map(int,e)):K(c) for e,c in saved['complete_weighted_top_sources'][i].items() if c};ck('whole_source_'+str(i),row==old)
 cpart={e:c for e,c in row.items() if sum(e[:6])==3};olds={tuple(map(int,e)):K(c) for e,c in cubic['raw_cubic_terms'][i].items() if c};ck('cubic_source_'+str(i),cpart==olds)
 out={};hist={}
 for e,c in row.items():
  n=sum(e[:6]);hist[n]=hist.get(n,0)+1
  R1=e[6]+e[3]+e[5];R2=e[7]+e[1]+e[4]
  key=(e[0]+e[5],e[1]+e[3],e[2]+e[4],R1%2,R2%2)
  add(out,key,c*(-q/r)**(R1//2)*(-q)**(R2//2))
 ck('entire_quadratic_algebra_zero_'+str(i),not compact(out));ck('normal_parity_'+str(i),all(n>0 and n%2==1 for n in hist))
 records.append({'row':i,'terms':len(row),'normal_histogram':hist})
# Recheck exact spectrum protection for full kernel and Schur row choices.
tan=load(str(p/'independent_ratio_RR_canonical_tangent_obstruction.sobj'));N0=tan['constant_matrix'];B=tan['obstruction_matrix'];ids=list(tan['pivot_rows']);cols=list(tan['pivot_columns']);bs=list(tan['rank_minor_rows']);C=N0.matrix_from_columns(cols);Kmat=tan['kernel_rows'].transpose()
ck('kernel_exact',N0*Kmat==0);dp=C.matrix_from_rows(ids).det();db=B.matrix_from_rows(bs).det();ck('pivot_formula',dp==-q*q/(r*r));ck('six_formula',db==QQ(2048)/729*q**6/r**5)
ck('ten_distinct_rows',len(set(ids+bs))==10);ck('ten_product_formula',dp*db==-QQ(2048)/729*q**8/r**7)
report={'status':'PASS','count':len(checks),'checks':checks,'records':records,'seconds':time.monotonic()-st,'scope':'Full28 weighted-top sources independently rebuilt from original full fields including all h,z coefficient derivatives and full Q0 field gauge; complete cone identities checked in the two quadratic RR extension algebra, all four embeddings. Fixed cyclic base only. Combined with protected rank10 generator matrix proves local tau membership, not generic-base exclusion.'}
(p/'independent_RR_whole_cone_reader.json').write_text(json.dumps(report,indent=2)+'\n');save({'checks':checks,'records':records,'combined_minor':dp*db},str(p/'independent_RR_whole_cone_reader.sobj'));print({k:v for k,v in report.items() if k not in ('checks','records')},flush=True)

