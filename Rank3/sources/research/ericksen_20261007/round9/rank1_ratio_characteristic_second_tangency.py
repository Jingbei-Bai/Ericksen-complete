"""Actual second tangencies of the homogeneous pencil relation L=0.

Computed directly in the original 17-state polynomial source, before any
u elimination.  The E1/p23!=0 t graph is substituted with p23 denominator.
"""
from pathlib import Path
import pickle,time,json
from sage.all import load,save,PolynomialRing,matrix
out=Path(__file__).resolve().parent;base=out.parent;start=time.monotonic()
dd=load(str(out/'rank1_ratio_characteristic_quartics_partial.sobj'))
R=dd['ring'];K=R.base_ring();a,b,c,h=dd['parameters'];p12,p13,p21,s,p31=R.gens()
r=(a-b)/(a*(b-1));q=1+r;z=q*h+c
tn=dd['t_numerators_over_p23'];assert len(tn)==6
old=pickle.loads((base/'round8/rank3_source_minor_candidates.sobj').read_bytes())
O=PolynomialRing(K,names=old['ring'].variable_names());vv=O.gens();cache={}
def ck(co):
 if co not in cache:
  def ev(p):return sum((K(v)*a**m[0]*b**m[1] for m,v in p.dict().items()),K.zero())
  cache[co]=ev(co.numerator())/ev(co.denominator())
 return cache[co]
def lift(p):return O({tuple(m):ck(co) for m,co in p.dict().items()})
F=[[lift(f) for f in row] for row in old['F']]
DL=[F[i][8]-q*F[i][4]-c*F[i][0] for i in range(3)]
DD=[]
for i in range(3):
 for j in range(3):
  f=sum((DL[i].derivative(vv[k])*F[j][k] for k in range(15)),O.zero())
  DD.append(f)
print('source_DD_ready',round(time.monotonic()-start,2),flush=True)
pics=[R.one(),p12,p13,p21,R(h),s,p31,-z*s/(q*h),R(z)]
low=[];uc=[];records=[]
for ii,f in enumerate(DD):
 f0=R.zero();cu=[K.zero(),K.zero()]
 for m,co in f.dict().items():
  if sum(m[15:]):
   assert sum(m[15:])==1 and not any(m[:15])
   cu[0 if m[15] else 1]+=co
  else:
   term=R(co)
   for x,e in zip(pics,m[:9]):
    if e:term*=x**e
   dt=sum(m[9:15]);assert dt<=1
   term*=s**(1-dt)
   for n,e in zip(tn,m[9:15]):
    if e:term*=n**e
   f0+=term
 low.append(f0);uc.append(cu)
 records.append({'row':ii,'directions':[ii//3,ii%3],'terms':len(f0.dict()),'degree':int(f0.total_degree()) if f0 else -1})
 print('DD_pullback',records[-1],round(time.monotonic()-start,2),flush=True)
U=matrix(K,uc);minors=[]
for i in range(9):
 for j in range(i+1,9):
  m=U.matrix_from_rows([i,j]).det()
  if m:minors.append({'rows':[i,j],'factorization':str(m.factor())})
payload={'ring':R,'parameters':(a,b,c,h),'t_numerators_over_p23':tn,'scaled_lower_rows':low,
 'u_coefficient_matrix':U,'source_second_tangencies':DD,'original_extended_ring':O,
 'records':records,'nonzero_u_minors':minors,'scale':'Original DD L = low/p23 + U*(u0,u1).',
 'source':'rank3_source_minor_candidates.sobj original17 Fi; L=p33-q*p22-c*p11'}
save(payload,str(out/'rank1_ratio_characteristic_second_tangency.sobj'))
report={'status':'EXACT_ACTUAL_SECOND_TANGENCIES','records':records,'u_rank_generic':int(U.rank()),
 'nonzero_u_minors':minors,'elapsed_seconds':time.monotonic()-start,
 'scope':'Necessary actual-open relations, all spectra and c. No dimension assertion.'}
(out/'rank1_ratio_characteristic_second_tangency.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2),flush=True)
