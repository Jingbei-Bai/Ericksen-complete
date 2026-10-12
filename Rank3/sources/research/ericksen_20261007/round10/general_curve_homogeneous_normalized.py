"""Remove exact spectral row/column factors from homogeneous necessary sources.

For F3 the odd pair contributes one free direction factor; the remaining
pair-square contributes c_j. Over an algebraic closure one invertible
diagonal rescaling of both entries of each pair removes all these c_j.
No positive-spectrum point is chosen.
"""
from pathlib import Path
import time,json
from sage.all import *
out=Path(__file__).resolve().parent;start=time.monotonic()
d=load(str(out/'general_curve_extra_quartic_cubics.sobj'))
q4=load(str(out/'general_curve_remaining_quartics.sobj'))
R=d['ring'];K=R.base_ring();a,r,H,Z,S,N=K.gens();q=1+r
c=[1/(a*r),K.one(),1/(1+a*r)]
rho=[a*r,a*r,K.one(),K.one(),q*(1+a*r)/r,q*(1+a*r)/r]
C=PolynomialRing(QQ,names=('H','Z','S','N'));L=C.fraction_field()
T=PolynomialRing(L,names=R.variable_names());checks={};cubics=[];quartics=[];records=[];odds=[]
def remove_spectrum(v):
 assert not v.derivative(a) and not v.derivative(r)
 num=v.numerator();den=v.denominator()
 # Evaluation at a=r=0 is only the canonical restriction of an already
 # proved a,r-independent rational expression, never a spectrum sample.
 return L(num(0,0,*L.gens()))/L(den(0,0,*L.gens()))
for i,p in enumerate(d['leading_rows']):
 vals={};oddset=set()
 for m,co in p.dict().items():
  deg=[int(m[2*j]+m[2*j+1]) for j in range(3)]
  odd=[j for j,x in enumerate(deg) if x%2];assert len(odd)==1
  oddset.add(odd[0]);square=[j for j in range(3) if deg[j]-(j==odd[0])==2];assert len(square)==1
  v=co/(rho[i]*c[square[0]])
  checks[f'cubic_{i}_{m}']=not v.derivative(a) and not v.derivative(r)
  vals[tuple(m)]=remove_spectrum(v)
 assert len(oddset)==1;odds.append(oddset.pop());f=T(vals);cubics.append(f)
 records.append({'kind':'cubic','source':int(d['row_specs'][i][0]),'odd_pair':odds[-1],'terms':f.number_of_terms()})
 print(records[-1],flush=True)
quartic_scales=[]
for i,p in enumerate(q4['leading']):
 vals={}
 for m,co in p.dict().items():
  deg=[int(m[2*j]+m[2*j+1]) for j in range(3)];assert all(x%2==0 for x in deg)
  vals[tuple(m)]=co/prod(c[j]**(deg[j]//2) for j in range(3))
 first=next(iter(vals.values()))
 # A common spectral factor may remain. Extract it from one coefficient
 # by separating irreducible numerator and denominator factors in a,r.
 spec=K.one()
 for pol,sgn in ((first.numerator(),1),(first.denominator(),-1)):
  for fac,e in pol.factor():
   if fac.degree(K.ring().gen(0)) or fac.degree(K.ring().gen(1)):
    assert not any(fac.degree(K.ring().gen(j)) for j in (2,3,4,5)),str(fac)
    spec*=K(fac)**(sgn*e)
 pure={}
 for m,v in vals.items():
  v/=spec;checks[f'quartic_{i}_{m}']=not v.derivative(a) and not v.derivative(r)
  pure[m]=remove_spectrum(v)
 f=T(pure);quartics.append(f);quartic_scales.append(spec)
 records.append({'kind':'quartic','source':[1,5,6][i],'terms':f.number_of_terms(),'spectral_scale':str(spec.factor())})
 print(records[-1],flush=True)
assert all(checks.values())
payload={'ring':T,'cubics':cubics,'odd_pairs':odds,'cubic_sources':[int(x[0]) for x in d['row_specs']],
 'quartics':quartics,'quartic_sources':[1,5,6],'spectral_pair_factors':c,'spectral_row_factors':rho,'spectral_quartic_factors':quartic_scales,'checks':checks,'records':records,
 'scope':'Exact full-spectrum homogeneous source transformation by invertible pair scaling over algebraic closure. No dimension conclusion.'}
save(payload,str(out/'general_curve_homogeneous_normalized.sobj'))
(out/'general_curve_homogeneous_normalized.json').write_text(json.dumps({'records':records,'checks_passed':sum(checks.values()),'elapsed':time.monotonic()-start},indent=2)+'\n')
print('DONE',time.monotonic()-start,flush=True)
