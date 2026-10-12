"""Entire Mobius parameter curve: original seven directional sources to unit."""
from pathlib import Path
import json,time
from sage.all import load,save,prod

out=Path(__file__).resolve().parent;start=time.monotonic();checks={}
raw=load(str(out/'general_curve_dimensionless_full.sobj'))
unit=load(str(out/'general_curve_unit_gradient_directions.sobj'))
input=load(str(out/'general_curve_B_mobius_directions_input.sobj'))
cert=load(str(out/'general_curve_B_mobius_unit_certificate.sobj'))
R=input['ring'];F=R.base_ring();C=F.ring();H=F.gen();ts=R.gens();Z=H/(H+1);S=1/(H+1)**2;N=-2/(H+1)**3
oldK=raw['ring'].base_ring();a,r,h,z,s,n=oldK.gens()
checks['exact_entire_actual_curve']=tuple(input['curve'])==(Z,S,N)==tuple(cert['curve'])
checks['first_contact_derivative']=Z.derivative(H)==S
checks['second_contact_derivative']=S.derivative(H)==N
checks['Mobius_plane_equation']=H*Z-H+Z==0
checks['Mobius_B_equation']=H*Z*(H**2-1)*S**2+2*(Z**2-H**2)*S-H*Z*(Z**2-1)==0
images=[ts[0],R.one(),ts[1],R.one(),ts[2],R.one()]
sources=[];scales=[]
for i,f in enumerate(raw['rows']):
    degree=2 if i<3 else 3
    oldscale=(1+r)/r if i<3 else unit['old_cubic_scales'][i-3]
    original=R.zero()
    for m,c in f.dict().items():
        if sum(m)!=degree:continue
        co=c/oldscale
        assert not co.derivative(a) and not co.derivative(r)
        fc=F(co.numerator()(0,0,H,Z,S,N))/F(co.denominator()(0,0,H,Z,S,N))
        original+=fc*prod(images[j]**int(m[j]) for j in range(6))
    checks[f'original_full_source_leader_{i}']=original==input['sources'][i]
    scale=F(input['scales'][i]);scales.append(scale)
    checks[f'primitive_exact_source_scale_{i}']=original==scale*input['primitive'][i]
    checks[f'scalar_is_nonzero_QQ_H_{i}']=bool(scale) and scale.numerator().parent()==C and scale.denominator().parent()==C
    checks[f'certificate_original_primitive_{i}']=cert['sources'][i]==input['primitive'][i]
    sources.append(input['primitive'][i])
checks['all_seven_sources_present']=len(sources)==len(cert['cofactors'])==len(cert['integer_H_cofactors'])==7
identity=sum((c*f for c,f in zip(cert['cofactors'],sources)),R.zero())
checks['original_seven_source_UNIT_identity']=identity==1
den=C(cert['H_denominator'])
checks['H_denominator_nonzero_polynomial']=bool(den) and den.degree()==12
for i,c in enumerate(cert['cofactors']):
    checks[f'cofactor_denominator_cleared_{i}']=c*den==cert['integer_H_cofactors'][i]
    checks[f'cleared_cofactor_coefficients_are_polynomials_{i}']=all(v.denominator()==1 for v in cert['integer_H_cofactors'][i].coefficients())
checks['polynomial_H_identity']=sum((c*f for c,f in zip(cert['integer_H_cofactors'],sources)),R.zero())==R(den)
checks['cofactor_terms_record']=[len(c.dict()) for c in cert['cofactors']]==[1,2,2,1,0,1,1]
assert all(checks.values()),[k for k,v in checks.items() if not v]
scope='Entire persistent real Mobius parameter curve excluded for all-mi-nonzero fullsupport infinity directions by exact original-seven-source unit identity over QQ(H). Exceptional scalar H roots are finite and cannot carry an H interval. All original spectra retained. No full rank classification.'
save({'checks':checks,'original_replayed_primitives':sources,'source_scales':scales,'H_denominator':den,'scope':scope},str(out/'independent_B_mobius_unit_reader.sobj'))
summary={'checks_passed':len(checks),'H_denominator':str(den.factor()),'elapsed_seconds':time.monotonic()-start,'scope':scope}
(out/'independent_B_mobius_unit_reader.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)
