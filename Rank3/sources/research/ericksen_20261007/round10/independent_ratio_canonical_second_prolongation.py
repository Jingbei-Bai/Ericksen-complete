"""Genuine second source prolongations, two weighted layers, exact base jets.

No repeated old-normal template: compute hatFi(hatFj S), including two
h,z coefficient derivatives, then its six-column canonical normal form.
"""
from sage.all import *
from pathlib import Path
import time,json,sys
p=Path(__file__).resolve().parent;start=time.monotonic();d=load(str(p/'joint_fields_direct.sobj'));fd=load(str(p/'independent_ratio_weighted_degreefall.sobj'));kd=load(str(p/'independent_ratio_diagonal_Q0_reduction.sobj'));sec=load(str(p/'independent_ratio_RR_cyclic_full_matrix.sobj'));oldnormal=load(str(p/'independent_ratio_canonical_amplitude_normal.sobj'));oldtop=load(str(p/'independent_ratio_RR_canonical_whole_cone.sobj'));tangent=load(str(p/'independent_ratio_RR_canonical_tangent_obstruction.sobj'))
oldK=d['ring'].base_ring();a,b,h,z=d['parameters'];K=sec['coefficient_field'];aa,r=sec['parameters'];q=1+r;Jsquare=q*q/r;cvhom=oldK.ring().hom(sec['old_parameter_map'],K);w=(1,)*6+(2,)*3;R=PolynomialRing(K,names=('n2_t0','n2_t1','n2_t2'));tt=R.gens();cache={};checks=[]
def ck(n,v):assert v,n;checks.append(n)
def wt(e):return sum(e[i]*w[i] for i in range(9))
def ev(c):c=oldK(c);return cvhom(c.numerator())/cvhom(c.denominator())
def jet(c,order=1):
 c=oldK(c)
 if c not in cache:
  ch=c.derivative(h);cz=c.derivative(z);cache[c]=(ev(c),ev(ch),ev(cz),ev(ch.derivative(h)),ev(ch.derivative(z)),ev(cz.derivative(z)))
 return cache[c]
Q={tuple(map(int,e)):oldK(c) for e,c in fd['Q0'].dict().items()};kap=[it['kappa'] for it in kd['items']];fields=[];bases=[]
for i in range(3):
 row=[]
 for j in range(9):
  vals={tuple(map(int,e)):oldK(c) for e,c in d['fields'][i,j].dict().items()}
  if j==6+i:
   for e,c in Q.items():vals[e]=vals.get(e,oldK.zero())-kap[i]*c
  vals={e:c for e,c in vals.items() if c and wt(e)-w[j] in (1,2)};row.append({e:jet(c)[:3] for e,c in vals.items()})
 fields.append(row);bases.append([{tuple(map(int,e)):jet(c)[:3] for e,c in d['base_fields'][i,j].dict().items() if wt(e) in (1,2)} for j in range(2)])
print('FIELDS',len(cache),time.monotonic()-start,flush=True)
def first_action(vals,i):
 out={}
 def add(e,cs):
  if e not in out:out[e]=[K.zero(),K.zero(),K.zero()]
  for j in range(3):out[e][j]+=cs[j]
 for e,c in vals.items():
  cj=jet(c)
  for l in range(9):
   if not e[l]:continue
   ee=list(e);ee[l]-=1
   for f,vj in fields[i][l].items():
    nu=tuple(ee[k]+f[k] for k in range(9))
    if wt(nu) not in (6,7):continue
    add(nu,[e[l]*cj[0]*vj[0],e[l]*(cj[1]*vj[0]+cj[0]*vj[1]),e[l]*(cj[2]*vj[0]+cj[0]*vj[2])])
  for l in range(2):
   dj=(cj[1],cj[3],cj[4]) if l==0 else (cj[2],cj[4],cj[5])
   for f,vj in bases[i][l].items():
    nu=tuple(e[k]+f[k] for k in range(9))
    if wt(nu) not in (6,7):continue
    add(nu,[dj[0]*vj[0],dj[1]*vj[0]+dj[0]*vj[1],dj[2]*vj[0]+dj[0]*vj[2]])
 return {e:tuple(cs) for e,cs in out.items() if any(cs)}
def next_action(vals,i):
 out={}
 def add(e,c):out[e]=out.get(e,K.zero())+c
 for e,cj in vals.items():
  for l in range(9):
   if not e[l]:continue
   ee=list(e);ee[l]-=1
   for f,vj in fields[i][l].items():
    nu=tuple(ee[k]+f[k] for k in range(9))
    if wt(nu) in (8,9):add(nu,e[l]*cj[0]*vj[0])
  for l in range(2):
   for f,vj in bases[i][l].items():
    nu=tuple(e[k]+f[k] for k in range(9))
    if wt(nu) in (8,9):add(nu,cj[1+l]*vj[0])
 return {e:c for e,c in out.items() if c}
def normal(top,low,parity=None,scale=K.one()):
 if parity is None:
  e=next(iter(top));n=sum(e[:6]);jp=((n-1)//2)%2;rho=(e[6]+e[3]+e[5],e[7]+e[1]+e[4],e[8]+e[0]+e[2]);parity=tuple((rho[j]-jp)%2 for j in range(3))
 sq=(-q/r,-q,K.one());out=[{} for _ in range(6)];cone={}
 def add(dic,e,c):dic[e]=dic.get(e,K.zero())+c
 for e,c in top.items():
  n=sum(e[:6]);assert n%2==1;(jp:=((n-1)//2)%2);rho=(e[6]+e[3]+e[5],e[7]+e[1]+e[4],e[8]+e[0]+e[2]);dif=tuple(rho[j]-parity[j] for j in range(3));assert all(v%2==jp for v in dif)
  value=c/scale*Jsquare**(-((n-1)//4))
  for j in range(3):value*=sq[j]**((dif[j]-jp)//2)
  ke=(e[0]+e[5],e[1]+e[3],e[2]+e[4]);add(cone,ke,value)
  for j in range(3):
   if e[j]:ne=list(ke);ne[j]-=1;add(out[j],tuple(ne),e[j]*value)
  add(out[4],ke,-K(dif[0])/2*value);add(out[5],ke,-K(dif[1])/2*value)
 for e,c in low.items():
  n=sum(e[:6])+1;jp=((n-1)//2)%2;rho=(e[6]+e[3]+e[5]+1,e[7]+e[1]+e[4]+1,e[8]+e[0]+e[2]+1);dif=tuple(rho[j]-parity[j] for j in range(3));assert all(v%2==jp for v in dif)
  value=c/scale*Jsquare**(-((n-1)//4))
  for j in range(3):value*=sq[j]**((dif[j]-jp)//2)
  ke=(e[0]+e[5],e[1]+e[3],e[2]+e[4]);add(out[3],ke,value)
 ck('whole cone vanishes '+str(len(checks)),not any(cone.values()))
 return [R({e:c for e,c in dd.items() if c}) for dd in out],parity
for source_index in ([int(v) for v in sys.argv[1:]] or [0]):
 source_index=source_index;all_sources=[fd['Q0']]+[it['full'] for it in fd['items']];f=all_sources[source_index];vals={tuple(map(int,e)):oldK(c) for e,c in f.dict().items() if wt(e) in (4,5)};first=[];new=[]
 for j in range(3):
  v=first_action(vals,j);first.append(v);top={e:c[0] for e,c in v.items() if wt(e)==7 and c[0]};low={e:c[0] for e,c in v.items() if wt(e)==6 and c[0]};rowid=7+3*source_index+j;row,_=normal(top,low,tangent['row_parities'][rowid],sec['row_scales'][rowid]);oldrow=[R({tuple(map(int,e)):c for e,c in ff.dict().items()}) for ff in oldnormal['matrix'].row(rowid)];ck(f'old normal source {source_index}/{j}',row==oldrow);print('FIRST VERIFIED',source_index,j,len(v),time.monotonic()-start,flush=True)
  for i in range(3):
   v2=next_action(v,i);top={e:c for e,c in v2.items() if wt(e)==9};low={e:c for e,c in v2.items() if wt(e)==8};row,pa=normal(top,low);rec={'source_index':source_index,'inner':j,'outer':i,'top_terms':len(top),'tau_terms':len(low),'normal_terms':[len(f.dict()) for f in row],'seconds':time.monotonic()-start};print(rec,flush=True);new.append({'source_index':source_index,'inner':j,'outer':i,'top':top,'tau':low,'parity':pa,'normal_row':row,'record':rec})
  save({'source_index':source_index,'coefficient_field':K,'parameters':(aa,r),'amplitude_ring':R,'first_source_jets':first,'items':new,'checks':checks,'scope':__doc__},str(p/f'independent_ratio_canonical_second_prolongation_source{source_index}.sobj'),compress=False)
 (p/f'independent_ratio_canonical_second_prolongation_source{source_index}.json').write_text(json.dumps({'status':'EXACT_SECOND_PROLONGATION_NORMAL_ROWS','checks':len(checks),'records':[it['record'] for it in new],'seconds':time.monotonic()-start},indent=2)+'\n');print('DONE',source_index,len(checks),time.monotonic()-start,flush=True)
 
