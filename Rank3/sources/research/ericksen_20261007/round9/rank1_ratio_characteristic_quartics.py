"""Exact original-Q4 pullback to the actual E1 pencil's protected t graph.

All state denominators are powers of p23; no new spectral/base pivot is used.
This exports necessary equations and does not assert a fiber dimension bound.
"""
from pathlib import Path
import pickle,time,json
from sage.all import QQ,PolynomialRing,load,save,vector
out=Path(__file__).resolve().parent;base=out.parent;start=time.monotonic()
d=load(str(out/'rank1_ratio_characteristic_cramer.sobj'))
Rold=d['ring'];K=Rold.base_ring();a,b,c,h=d['parameters'];r=(a-b)/(a*(b-1));q=1+r;z=q*h+c
R=PolynomialRing(K,names=('p12','p13','p21','p23','p31'))
p12,p13,p21,p23,p31=R.gens()
hom=Rold.hom([*R.gens(),*[R.zero()]*6],R)
M=d['t_coefficient_matrix'];bb=d['constant_terms'];N=d['t_numerators'];det=d['t_denominator']
expected=4*q*c*c*(1+r*h)*(1-r*z/q)*(1-r*c/q)**2*Rold.gen(3)/h
checks={'protected_minor':det==expected}
for i in range(6):
 checks[f'cramer_identity_{i}']=not(sum((M[i,j]*N[j] for j in range(6)),Rold.zero())+bb[i]*det)
 print('cramer_audit',i,checks[f'cramer_identity_{i}'],round(time.monotonic()-start,2),flush=True)
assert all(checks.values()),checks
scalar=K(hom(det).monomial_coefficient(p23))
tn=[hom(n)/scalar for n in N] # t_j=tn_j/p23, scalar is protected by identity above.
assert all(x in R for x in tn)
tn=list(map(R,tn))
old=pickle.loads((base/'round8/rank3_source_minor_candidates.sobj').read_bytes())
normal=pickle.loads((base/'round8/rank3_k_quartic_diagonal_normal.pkl').read_bytes())
pics=[R.one(),p12,p13,p21,R(h),p23,p31,-z*p23/(q*h),R(z)]
cache={}
def ck(co):
 if co not in cache:
  def ev(pp):return sum((K(v)*a**m[0]*b**m[1] for m,v in pp.dict().items()),K.zero())
  cache[co]=ev(co.numerator())/ev(co.denominator())
 return cache[co]
def img_pi(m):
 ans=R.one()
 for p,e in zip(pics,m[:9]):
  if e:ans*=p**e
 return ans
low=[];uc=[];records=[]
for ii,poly in enumerate(old['Q4']):
 f=R.zero();cu=[R.zero(),R.zero()]
 for m,co in poly.dict().items():
  term=ck(co)*img_pi(m)
  if sum(m[15:]):
   assert sum(m[15:])==1 and sum(m[9:15])==0
   cu[0 if m[15] else 1]+=term
  else:
   nt=sum(m[9:15]);assert nt<=2
   term*=p23**(2-nt)
   for n,e in zip(tn,m[9:15]):
    if e:term*=n**e
   f+=term
 low.append(f);uc.append(cu)
 print('Q_lower',ii,'terms',len(f.dict()),'seconds',round(time.monotonic()-start,2),flush=True)
 save({'ring':R,'parameters':(a,b,c,h),'t_numerators_over_p23':tn,'Q_scaled_lower':low,'Q_u_coefficients':uc,'checks':checks},str(out/'rank1_ratio_characteristic_quartics_partial.sobj'))
L=normal['coefficient_transform'];indices=normal['source_quartic_indices']
def symco(e):return K(str(e))
un=[-sum((symco(L[i,j])*low[indices[j]] for j in range(6)),R.zero()) for i in range(2)]
for i in range(2):
 for k in range(2):
  checks[f'normal_u_{i}_{k}']=sum((symco(L[i,j])*uc[indices[j]][k] for j in range(6)),R.zero())==int(i==k)
assert all(checks.values()),checks
Q=[]
for ii in range(21):
 f=low[ii]+uc[ii][0]*un[0]+uc[ii][1]*un[1]
 Q.append(f)
 records.append({'row':ii,'degree':int(f.total_degree()) if f else -1,'terms':len(f.dict())})
 print('Q_final',ii,records[-1],round(time.monotonic()-start,2),flush=True)
payload={'ring':R,'parameters':(a,b,c,h),'q':q,'r':r,'z':z,'p32':pics[7],
 't_numerators_over_p23':tn,'u_numerators_over_p23_squared':un,'Q_scaled_lower':low,
 'Q_u_coefficients':uc,'quartic_source_rows':Q,'source_row_scale':'p23^2',
 'remaining_tangencies':[hom(f)/scalar for f in d['residual_numerators']],
 'normal_source_indices':indices,'normal_coefficient_transform':L,'checks':checks,'records':records}
save(payload,str(out/'rank1_ratio_characteristic_quartics.sobj'))
report={'status':'PASS_EXACT_ORIGINAL_Q4_PULLBACK','checks':checks,'records':records,
 'denominators':'Only original spectrum/base guards, K=1-r*c/q, and p23; t_j=tn_j/p23, u_j=un_j/p23^2.',
 'scope':'Necessary source equations on the E1 characteristic pencil, all spectra and all c. No fiber-dimension or exclusion conclusion.',
 'elapsed_seconds':time.monotonic()-start}
(out/'rank1_ratio_characteristic_quartics.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2),flush=True)
