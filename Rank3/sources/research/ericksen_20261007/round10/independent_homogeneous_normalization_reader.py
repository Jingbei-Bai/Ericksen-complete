"""Full original-Q pullback and coherent all-spectrum pair normalization audit.

Polynomial parameter substitution uses one Laurent clearing per coefficient,
independent of the original generator's repeated rational substitutions.
"""
from pathlib import Path
from functools import lru_cache
import time,json
from sage.all import load,save,prod

out=Path(__file__).resolve().parent;base=out.parent;start=time.monotonic();checks={}
normal=load(str(out/'general_curve_homogeneous_normalized.sobj'))
extra=load(str(out/'general_curve_extra_quartic_cubics.sobj'))
q4=load(str(out/'general_curve_remaining_quartics.sobj'))
original=load(str(base/'round8/rank3_joint_fiber_leading.sobj'))
joint=load(str(out/'joint_fields_direct.sobj'))
full=load(str(out/'general_curve_dimensionless_full.sobj'))
R=extra['ring'];K=R.base_ring();CR=K.ring();a,r,H,Z,S,N=K.gens()
aa,rr,hh,zz,ss,nn=CR.gens();qq=1+rr;mm=1+aa*rr;q=1+r
c=normal['spectral_pair_factors'];rho=normal['spectral_row_factors'];sq=normal['spectral_quartic_factors']
def spectral_unit(v):
    v=K(v)
    if not v:return False
    return all(CR(f).is_constant() or any(all((K(f)/g).derivative(x)==0 for x in K.gens()) for g in [a,r,1+r,1+a*r]) for p in [v.numerator(),v.denominator()] for f,e in p.factor())
checks['pair_factors_exact']=c==[1/(a*r),K.one(),1/(1+a*r)]
checks['all_spectral_units']=all(spectral_unit(v) for v in c+rho+sq)
def back(v):return K(v.numerator()(H,Z,S,N))/K(v.denominator()(H,Z,S,N))
for i,p in enumerate(extra['leading_rows']):
    f=normal['cubics'][i]
    checks['cubic_support_'+str(i)]=set(map(tuple,p.dict()))==set(map(tuple,f.dict()))
    for mon,co in p.dict().items():
        degrees=[sum(mon[2*j:2*j+2]) for j in range(3)]
        odd=[j for j,degree in enumerate(degrees) if degree%2]
        square=[j for j in range(3) if degrees[j]-(j in odd)==2]
        assert len(odd)==len(square)==1
        v=back(f.monomial_coefficient(f.parent().monomial(*mon)))
        key=str(i)+'_'+','.join(map(str,mon))
        checks['cubic_coefficient_'+key]=co==rho[i]*c[square[0]]*v
        checks['cubic_spectrum_free_'+key]=v.derivative(a)==0 and v.derivative(r)==0
        checks['coherent_odd_pair_'+str(i)]=normal['odd_pairs'][i]==odd[0]
for i,p in enumerate(q4['leading']):
    f=normal['quartics'][i]
    checks['quartic_support_'+str(i)]=set(map(tuple,p.dict()))==set(map(tuple,f.dict()))
    for mon,co in p.dict().items():
        degrees=[sum(mon[2*j:2*j+2]) for j in range(3)]
        assert all(degree%2==0 for degree in degrees)
        v=back(f.monomial_coefficient(f.parent().monomial(*mon)))
        key=str(i)+'_'+','.join(map(str,mon))
        checks['quartic_coefficient_'+key]=co==sq[i]*prod([c[j]**(degrees[j]//2) for j in range(3)])*v
        checks['quartic_spectrum_free_'+key]=v.derivative(a)==0 and v.derivative(r)==0
print('coherent normalization coefficients PASS',round(time.monotonic()-start,2),flush=True)

images=extra['images']
checks['same_all_pair_coordinates']=images==q4['images']
@lru_cache(None)
def evaluate_polynomial(poly):
    if not poly:return CR.zero(),CR.one()
    terms=poly.dict();bdeg=max(int(ex[1]) for ex in terms);rdeg=max(int(ex[2]+ex[3]) for ex in terms)
    answer=CR.zero()
    for ex,co in terms.items():
        ia,ib,ih,iz=map(int,ex)
        answer+=co*aa**(ia+ib)*qq**(ib+iz)*hh**ih*zz**iz*mm**(bdeg-ib)*rr**(rdeg-ih-iz)
    return answer,mm**bdeg*rr**rdeg
@lru_cache(None)
def coefficient(v):
    p,pden=evaluate_polynomial(v.numerator());q,qden=evaluate_polynomial(v.denominator())
    return K(p*qden)/K(q*pden)
@lru_cache(None)
def state_monomial(ex):return prod([images[j]**int(e) for j,e in enumerate(ex) if e])
def pull(poly):
    return sum((coefficient(v)*state_monomial(tuple(ex)) for ex,v in poly.dict().items()),R.zero())
for i in range(3):
    checks['original_h_velocity_'+str(i)]=pull(joint['normalized_pi_fields'][i,3])==R.gen(2*i+1)/r
    checks['original_z_velocity_'+str(i)]=pull(joint['normalized_pi_fields'][i,7])==q*S*R.gen(2*i+1)/r
checks['original_Q0']=pull(original['quartic_source_rows'][0])==full['rows'][3]
print('velocity and Q0 provenance PASS',round(time.monotonic()-start,2),flush=True)
for i,index in enumerate([1,5,6]):
    p=pull(original['quartic_source_rows'][index])
    checks['full_original_Q'+str(index)]=p==q4['rows'][i]
    top=R({tuple(ex):v for ex,v in p.dict().items() if sum(ex)==4})
    checks['original_highest_Q'+str(index)]=top==q4['leading'][i]
    print('full original Q',index,'PASS',round(time.monotonic()-start,2),flush=True)
checks['complete_three_quartics']=len(q4['rows'])==len(q4['leading'])==3 and normal['quartic_sources']==[1,5,6]
assert all(checks.values()),[key for key,val in checks.items() if not val]
payload={'checks':checks,'pair_factors':c,'cubic_row_factors':rho,'quartic_row_factors':sq,
 'scope':'Full original Q1,Q5,Q6 pullback and all homogeneous normalization identities, with the same invertible pair scaling x_old=gamma*x_new and gamma_j^2=1/c_j over C. Every spectral factor is an old guard. No assertion that real square roots exist, and no fiber dimension conclusion.'}
save(payload,str(out/'independent_homogeneous_normalization_reader.sobj'))
report={'status':'EXACT_FULL_SOURCE_AND_PAIR_NORMALIZATION_PASS','checks':checks,'check_count':len(checks),'elapsed_seconds':time.monotonic()-start}
(out/'independent_homogeneous_normalization_reader.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({key:val for key,val in report.items() if key!='checks'}),flush=True)
