"""Independent regular-lambda source -> projection -> actual-contact audit."""
from pathlib import Path
import time,json
from sage.all import load,save,QQ,PolynomialRing,prod

out=Path(__file__).resolve().parent;start=time.monotonic();checks={}
red=load(str(out/'independent_scalar_lambda_reduction.sobj'))
reg=load(str(out/'general_curve_scalar_lambda_regular.sobj'))
prim=load(str(out/'general_curve_scalar_lambda_primitive.sobj'))
lc=load(str(out/'general_curve_scalar_lambda_contact_L.sobj'))
proj=load(str(out/'general_curve_scalar_lambda_projection.sobj'))
fac=load(str(out/'general_curve_scalar_lambda_projection_factors.sobj'))
actual=load(str(out/'general_curve_scalar_projection_contact.sobj'))
R=reg['ring'];H,Z,S=R.gens();K=R.fraction_field()
lam,tau=reg['lambda'],reg['tau'];oldL,oldTau=red['ring'].gens()[3:]
aff=red['affine_sources']
checks['first_cramer_equation']=K(aff[0](H,Z,S,lam,tau))==0
checks['second_cramer_equation']=K(aff[1](H,Z,S,lam,tau))==0
def base(p):return R(p(H,Z,S,0,0))
checks['lambda_pivot']=reg['lambda_pivot']==base(aff[0].derivative(oldL))
checks['tau_pivot']=reg['tau_pivot']==base(aff[1].derivative(oldTau))
oldguards=[H,Z,H-Z,1-Z,1+H,1+H-Z,1-S,1+H*S-Z,(1+H)*S-Z]
allowed=oldguards+[reg['P'],reg['L1'],reg['L2']]
def unit(poly):
    poly=R(poly)
    if not poly:return False
    for p,e in poly.factor():
        if p.is_constant():continue
        if not any(all((K(p)/g).derivative(x)==0 for x in R.gens()) for g in allowed):return False
    return True
for name,f in [('lambda',lam),('tau',tau),('N',reg['N'])]:checks[name+'_denominator_guarded']=unit(f.denominator())
checks['both_pivots_guarded']=unit(reg['lambda_pivot']) and unit(reg['tau_pivot'])
checks['N_reconstruction']=reg['N']==K(red['N_numerator'](H,Z,S,lam,tau))/K(red['N_denominator'](H,Z,S,0,0))
for i,p in enumerate([aff[2],red['quadratic_source']]):
    val=K(p(H,Z,S,lam,tau))
    checks['regular_original_source_'+str(i)]=val.numerator()==reg['rows'][i]*reg['removed'][i]
    checks['regular_strip_guarded_'+str(i)]=unit(reg['removed'][i]) and unit(val.denominator())
    checks['primitive_source_'+str(i)]=reg['rows'][i]==prim['source_units'][i]*prim['rows'][i] and QQ(prim['source_units'][i])!=0

Lring=lc['ring'];LH,LZ,L=Lring.gens();LF=Lring.fraction_field()
for i,p in enumerate(prim['rows']):
    val=LF(p(LH,LZ,(LZ+L)/LH))
    checks['L_coordinate_source_'+str(i)]=val==lc['source_units'][i]*lc['rows'][i]/lc['source_denominators'][i]
    den=lc['source_denominators'][i]
    checks['L_den_only_H_'+str(i)]=all(f.is_constant() or (LF(f)/LH).derivative(LH)==0 and (LF(f)/LH).derivative(LZ)==0 and (LF(f)/LH).derivative(L)==0 for f,e in den.factor())
checks['projection_source_order']=proj['rows']==lc['rows']
print('source chain replayed',round(time.monotonic()-start,2),flush=True)
# Independently recompute the univariate elimination from the verified rows.
recomputed=lc['rows'][0].resultant(lc['rows'][1],L)
checks['projection_resultant']=recomputed==proj['resultant']==fac['resultant'] and bool(recomputed)
factorproduct=prod([f**e for f,e in fac['factors']])
ratio=LF(recomputed/factorproduct)
checks['complete_factor_product']=bool(ratio) and all(ratio.derivative(x)==0 for x in Lring.gens())
checks['all_factors_independent_L']=all(f.degree(L)==0 for f,e in fac['factors'])

B=PolynomialRing(QQ,names=('H','Z'));BH,BZ=B.gens()
def bpoly(f):return B(f(BH,BZ,0))
oldplane=[BH,BZ,BH-BZ,1-BZ,BH+1,1+BH-BZ]
def coefficients(p):
    ans=[B.zero()]*(p.degree(L)+1)
    for mm,c in p.dict().items():ans[mm[2]]+=c*BH**mm[0]*BZ**mm[1]
    return ans
coeffs=[coefficients(p) for p in lc['rows']]
def rem(p,g):return p.quo_rem(g)[1]
def independent_eval(cs,num,den,g):
    # Homogenized evaluation using independently built modular power tables,
    # rather than the generating script's Horner recurrence.
    n=len(cs)-1;np=[B.one()];dp=[B.one()]
    for _ in range(n):
        np.append(rem(np[-1]*num,g));dp.append(rem(dp[-1]*den,g))
    result=B.zero()
    for i,c in enumerate(cs):result=rem(result+rem(c*np[i],g)*dp[n-i],g)
    return result
records=[];seen=[]
for index,(ff,e) in enumerate(fac['factors']):
    g=bpoly(ff)
    if any(g.quo_rem(guard)[1]==0 for guard in oldplane):continue
    seen.append(index);den=g.derivative(BZ);num=-BH*g.derivative(BH)-BZ*den
    checks['R_RZ_'+str(index)]=g.gcd(den).is_constant()
    residues=[independent_eval(cs,num,den,g) for cs in coeffs]
    saved=next(item for item in actual['certificates'] if item['factor_index']==index)
    checks['contact_numerator_'+str(index)]=num==saved['num_L'] and den==saved['den_L']
    checks['source_residues_'+str(index)]=residues==saved['source_residues']
    gg=g
    for p in residues:gg=gg.gcd(p)
    checks['excludes_interval_'+str(index)]=gg.is_constant()
    records.append({'index':index,'factor':g,'residues':residues,'joint_gcd':gg})
    print('contact factor',index,'PASS',round(time.monotonic()-start,2),flush=True)
checks['all_eleven_factors']=len(seen)==11 and seen==[item['factor_index'] for item in actual['certificates']]
checks['complete_flag']=actual['complete'] is True
assert all(checks.values()),[key for key,val in checks.items() if not val]
save({'checks':checks,'records':records,'resultant':recomputed,'factor_unit':ratio,
 'scope':'Exact source/projection/contact certificate. Together with independently excluded A=0 and D00=0 open walls, excludes every actual scalar-monodromy open interval in the two-invertible-block chart. At least-two-rank-one-block charts remain outside this conclusion.'},
 str(out/'independent_scalar_regular_projection_reader.sobj'))
report={'status':'EXACT_COMPLETE_REGULAR_SCALAR_AUDIT_PASS','checks':checks,'check_count':len(checks),'elapsed_seconds':time.monotonic()-start}
(out/'independent_scalar_regular_projection_reader.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({key:val for key,val in report.items() if key!='checks'}),flush=True)
