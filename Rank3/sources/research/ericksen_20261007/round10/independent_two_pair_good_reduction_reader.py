"""Exact integer/Gauss certificate audit; not a sampled-spectrum argument."""
from pathlib import Path
import json,time,sys
from sage.all import load,save,QQ,ZZ,GF,PolynomialRing,gcd
out=Path(__file__).resolve().parent;start=time.monotonic();checks={};cert=[]
reverse='reverse' in sys.argv[1:]
prefix='general_curve_two_pair_reverse' if reverse else 'general_curve_two_pair_endpoint'
output='independent_two_pair_reverse_good_reduction_reader' if reverse else 'independent_two_pair_good_reduction_reader'
d=load(str(out/(prefix+'_big_certificate.sobj')))
R=d['ring'];H,Z=R.gens();UH=PolynomialRing(QQ,'h');hh=UH.gen()
for item in d['items']:
    missing=item['missing'];tag=str(missing)
    pr=load(str(out/(prefix+'_projection_'+tag+'.sobj')))
    branch=pr['items'][2];source=branch['source_factor'];OS=source.parent().gen(2)
    F=item['F'];cs=item['source_coefficients'];degree=item['source_S_degree']
    rawlarge=[f for f,e in branch['factors'] if len(f.dict())>100]
    checks['one_big_factor_'+tag]=len(rawlarge)==1
    checks['original_factor_'+tag]=F*item['F_unit']==R(rawlarge[0](H,Z,0)) and bool(item['F_unit'])
    checks['F_integer_'+tag]=all(c in ZZ for c in F.coefficients())
    checks['source_integer_'+tag]=all(c in ZZ for p in cs for c in p.coefficients())
    rec=source.parent().zero()
    for j,cc in enumerate(cs):
        for mon,c in cc.dict().items():rec+=c*source.parent().gen(0)**mon[0]*source.parent().gen(1)**mon[1]*OS**j
    checks['integer_T_source_reconstruction_'+tag]=rec*item['source_unit']==source and bool(item['source_unit']) and degree==source.degree(OS) and len(cs)==degree+1
    byZ=[UH.zero()]*(F.degree(Z)+1)
    for mon,c in F.dict().items():byZ[mon[1]]+=c*hh**mon[0]
    content=gcd(byZ)
    # QQ[H] gcds are units up to a nonzero rational; normalization may depend
    # on coefficient order. The property needed is exactly constant nonzero.
    saved_content=item['content_Z_F']
    checks['no_H_only_factor_'+tag]=bool(content) and content.is_constant() and bool(saved_content) and saved_content.is_constant()
    gr=item['good_reduction'];p=ZZ(gr['prime']);h0=ZZ(gr['H_value'])
    checks['prime_'+tag]=p.is_prime()
    FF=GF(p);A=PolynomialRing(FF,'z');z=A.gen()
    def ev(poly):return sum((FF(c)*FF(h0)**mon[0]*z**mon[1] for mon,c in poly.dict().items()),A.zero())
    ff=ev(F);nn=ev(-F.derivative(H));dd=ev(F.derivative(Z))
    # A different exact evaluation from the exporter's Horner recurrence.
    tt=sum((ev(co)*nn**j*dd**(degree-j) for j,co in enumerate(cs)),A.zero())
    checks['F_reduction_'+tag]=ff==gr['F_reduction']
    checks['T_integer_homomorphism_'+tag]=tt==gr['T_reduction']
    checks['Z_degree_preserved_'+tag]=ff.degree()==F.degree(Z)==gr['retained_Z_degree']
    checks['derivative_homomorphism_'+tag]=dd==ff.derivative()
    aa,bb=gr['bezout'];da,db=gr['derivative_bezout']
    checks['contact_Bezout_one_'+tag]=A(aa)*ff+A(bb)*tt==1
    checks['derivative_Bezout_one_'+tag]=A(da)*ff+A(db)*dd==1
    cert.append({'missing':missing,'integer_F':F,'integer_T_source_coefficients':cs,'content_Z':content,'F_reduction':ff,'T_reduction':tt,'prime':p,'H_value':h0,'source_identity':rec*item['source_unit']==source})
checks['complete_three_cases']=d['complete'] and len(cert)==3 and [v['missing'] for v in cert]==[0,1,2]
assert all(checks.values()),[k for k,v in checks.items() if not v]
save({'checks':checks,'certificates':cert,'scope':'Characteristic-zero bivariate coprimality certified by primitive integer content and one degree-preserving exact good reduction per polynomial. Not a generic-spectrum inference. Upstream endpoint/resultant provenance is audited separately.'},str(out/(output+'.sobj')))
report={'status':'EXACT_INTEGER_GOOD_REDUCTION_COPRIMALITY_PASS','check_count':len(checks),'checks':checks,'elapsed_seconds':time.monotonic()-start}
(out/(output+'.json')).write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='checks'}),flush=True)
