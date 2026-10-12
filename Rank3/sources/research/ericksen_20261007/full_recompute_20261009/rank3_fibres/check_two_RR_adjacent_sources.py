"""Exact adjacent source-list/normalization links in the complete two-RR chain."""
from sage.all import *
from pathlib import Path
import json,hashlib
out=Path(__file__).absolute().parent;research=out.parent.parent;checks=[];imports=[]
def old(n):
 p=research/'round10'/(n+'.sobj');imports.append(dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()));return load(str(p))
def ck(n,v):assert bool(v),n;checks.append(n)
d=old('independent_subcritical_two_RR_delta_even_elimination');sc=old('independent_subcritical_two_RR_delta_material_scaling');P=d['ring'];t,w,s,e=P.gens()
raw=[]
for rec in d['necessary_eliminants']:
 for key,f in sorted(rec['base_coefficients'].items()):raw.append((('F_even',rec['even_row'],key),f))
ev=d['even_sources'];f=ev[0].derivative(e)*ev[1]-ev[1].derivative(e)*ev[0];raw.append((('old_even_compat',),P(f/f.factor().unit())))
ck('complete scaling source count',len(raw)==len(sc['items']))
for j,((lab,f),it) in enumerate(zip(raw,sc['items'])):
 ck('actual scaling label '+str(j),lab==it['label']);ck('actual scaling source '+str(j),f==it['original_source'])
ck('all primitive scaling outputs',[it['primitive'] for it in sc['items']]==sc['polynomials'])
curve=old('independent_subcritical_two_RR_D_curve_sources')
for j,it in enumerate(curve['items']):ck('actual D curve input '+str(j),it['original_source']==sc['polynomials'][it['source_index']])
ck('complete Dcurve generator list',list(dict.fromkeys(it['primitive'] for it in curve['items'] if it['primitive']))==curve['sources'])
fact=old('independent_subcritical_two_RR_D_curve_factor_ideals')
ck('compatibility is actual final Dcurve source',fact['compatibility']==curve['sources'][-1])
edata=old('independent_subcritical_two_RR_even_resultants');eq=old('independent_subcritical_two_RR_D_curve_even_resultants')
F=edata['flat_ring'];BF=PolynomialRing(QQ,names=('r','c','T','y'));r,c,T,y=BF.gens();BFF=BF.fraction_field();ph=F.hom([c*c/r,r,BFF.zero(),BFF.zero(),-y/(1+y),c*T*y,c*T],BFF)
def bc(f):
 def pc(v):return BF({tuple(map(int,e)):QQ(co) for e,co in v.dict().items()}) if hasattr(v,'dict') else BF(v)
 return BFF(pc(f.numerator()))/pc(f.denominator())
for j,it in enumerate(eq['items']):
 ck('original quadratic-even resultant chain '+str(j),it['original_resultant']==edata['items'][2]['items'][j]['flat_resultant'])
 val=ph(it['state_factor']);ck('actual Dcurve even pullback '+str(j),val==bc(it['pullback']))
 hp=eq['ring'].hom([T*T,y],BF);num=bc(it['content'])*T**it['T_power']*y**it['y_power']*(y+1)**it['yplus_power']*it['rational_content']*hp(it['primitive'])
 ck('complete even source primitive reconstruction '+str(j),val.numerator()==num and val.denominator()==bc(it['denominator']))
 for f,n in bc(it['content']).factor():
  ck('removed even material content protected '+str((j,f)),any(BFF(f)/g in QQ and BFF(f)/g for g in [r,c,r+1,c*c+1]))
 for f,n in val.denominator().factor():
  ck('even pullback denominator protected '+str((j,f)),any(BFF(f)/g in QQ and BFF(f)/g for g in [r,c,r+1,c*c+1,T,y,y+1]))
cert=old('independent_subcritical_two_RR_D_curve_exclusion_certificate');AP=old('independent_subcritical_two_RR_D_curve_A_projection');Q=fact['ring'];N,X,Y=Q.gens();h=eq['ring'].hom([X,Y],Q)
ck('R witness actual full source list',cert['R_sources']==fact['R_sources']+[fact['compatibility'],h(eq['items'][0]['primitive'])])
ck('A witness actual full source list',cert['A_sources']==AP['A_sources'])
projection=old('independent_subcritical_two_RR_D_curve_even_projection');uy=projection['ring'].gen();ay=cert['A_target'];actualAy=projection['ring'](ay(0,uy))
for j,it in enumerate(cert['A_even_Bezout']):ck('A even Bezout actual projected target '+str(j),it['A_y']==actualAy)
tr=old('independent_three_D_F_even_transport')
for ci,it in enumerate(tr['items']):
 coef=[]
 for rec in it['records']:coef.extend(f for key,f in sorted(rec['base_coefficients'].items()))
 ck('all actual transported coefficient source rows '+str(ci),coef==it['coefficient_sources'])
 ck('actual transported reference sources '+str(ci),tr['reference_sources']==sc['polynomials'][:10])
# The kappa elimination also uses an intermediate selected source list.
k=old('independent_subcritical_two_RR_kappa_compatibility');fs=k['original_primitive_rows']
ck('kappa actual three linear source selection',k['linear_rows']==[fs[0],fs[2],fs[3]])
report=dict(status='PASS',count=len(checks),checks=checks,imports=imports,scope='All adjacent coefficient/source arrays in representative Dcurve and other-support F-even transport; exact primitive reconstruction with only protected material/support contents; complete witness source lists and pure-kappa source selection.')
(out/'two_RR_adjacent_sources.json').write_text(json.dumps(report,indent=2));print('PASS',len(checks))
