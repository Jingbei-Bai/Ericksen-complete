from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();checks={}
def ck(n,v):
 if not v:raise AssertionError(n)
 checks[n]=True
bd=load(str(p/'independent_ratio_RR_boundary_leading.sobj'));cl=load(str(p/'independent_ratio_RR_boundary_classification.sobj'));pt=load(str(p/'independent_ratio_two_RR_principal_tangency.sobj'));hi=load(str(p/'independent_ratio_two_RR_genuine_higher.sobj'));ex=load(str(p/'independent_ratio_two_RR_mixed_exclusion.sobj'));fd=load(str(p/'independent_ratio_weighted_degreefall.sobj'));raw=load(str(p/'joint_fields_direct.sobj'));red=load(str(p/'independent_ratio_diagonal_Q0_reduction.sobj'));orig=load(str(p.parent/'round8/rank3_joint_fiber_leading.sobj'))
K=bd['coefficient_field'];a,r,H,Z=bd['parameters'];q=1+r;G=bd['physical_ring'];rr=G.gens()[:3];pp=G.gens()[3:];ok=raw['ring'].base_ring();h,z=raw['parameters'][2:];phi=ok.ring().hom(bd['parameter_map'],K);wt=(1,)*6+(2,)*3;cache={};dcache={}
def cv(c):
 c=ok(c)
 if c not in cache:cache[c]=phi(c.numerator())/phi(c.denominator())
 return cache[c]
def dc(c,j):
 key=(ok(c),j)
 if key not in dcache:dcache[key]=cv(key[0].derivative((h,z)[j]))
 return dcache[key]
def rd(f):return {tuple(map(int,e)):ok(c) for e,c in f.dict().items() if c}
def weight(e):return sum(a*b for a,b in zip(e,wt))
def clean(d):return {e:c for e,c in d.items() if c}
def top(f,n):return {e:c for e,c in rd(f).items() if weight(e)==n}
def mapped(di):return G({tuple(e[6:])+tuple(e[:6]):cv(c) for e,c in di.items()})
alpha=(Z-H)/(Z*(H+1));beta=H*(Z-1)/(Z*(H+1));Ls=(rr[0]*pp[0]-alpha*rr[2]*pp[5],rr[0]*pp[1]-alpha/beta*rr[1]*pp[3],rr[1]*pp[2]-beta*rr[2]*pp[4]);ck('global_linear_forms',tuple(bd['linear_forms'])==Ls)
for idx,j in bd['quadratic_linear_indices'].items():ck('original_global_square_'+str(idx),mapped(top(orig['quartic_source_rows'][idx],6))==bd['quadratic_units'][idx]*Ls[j]**2)
ck('original_Q0',mapped(top(orig['quartic_source_rows'][0],5))==bd['original_Q0_top'])
Rtops=[]
for j,it in enumerate(fd['items']):
 f=top(it['full'],5);Rtops.append(f);ck('six_R_source_dict_'+str(j),f=={tuple(map(int,e)):ok(c) for e,c in bd['original_R_top_dictionaries'][j][1].items()});ck('six_R_source_poly_'+str(j),mapped(f)==hi['original_mapped_R_tops'][j][1])
# All six support charts exhaust the three global linear constraints.
for bi,item in enumerate(bd['items']):
 ch=item['chart'];T0=bd['chart_ring'];free=ch['free'];ev=G.hom(list(ch['RR'])+list(ch['p']),T0)
 ck('boundary_L_zero_'+str(bi),all(ev(f)==0 for f in Ls));J=matrix(T0,[[ev(f.derivative(v)) for v in pp] for f in Ls]);B=matrix(T0,[[f.derivative(v) for v in free] for f in ch['p']]);ck('complete_chart_dimension_'+str(bi),B.rank()==len(free) and J.rank()+len(free)==6)
 for j,f in enumerate(Rtops):ck('boundary_R_pullback_'+str(bi)+'_'+str(j),ev(mapped(f))==item['sources'][j]['polynomial'])
 ci=cl['items'][bi];expected=ci['expected_equations'];units=ci['units'];FT=T0.fraction_field()
 for source in item['sources']:
  name=source['source'];value=source['polynomial'];ck('boundary_factored_'+str(bi)+'_'+name,value==0 if name not in expected else FT(value)==units[name]*expected[name])
  if name in units:
   for side in (units[name].numerator(),units[name].denominator()):
    di=side.dict();ck('boundary_unit_single_'+str(bi)+'_'+name+'_'+str(len(checks)),len(di)==1);ee,co=next(iter(di.items()));ck('boundary_unit_RR_only_'+str(bi)+'_'+name+'_'+str(len(checks)),all(int(e)==0 for e in ee[1:]))
    for po in (co.numerator(),co.denominator()):
     for fac,e in po.factor():ck('boundary_unit_guard_'+str(bi)+'_'+name+'_'+str(fac),bd['baseguard'].quo_rem(fac)[1]==0)
ck('two_support_determinants',tuple(cl['two_RR_minors'])==(3*(H*Z-H+Z),-3*(H*Z-H+Z),-3*(H*Z-H+Z)))
# Rebuild complete principal fields before any support substitution.
principal=[]
for i in range(3):
 fs=[]
 for j in range(9):
  f=rd(raw['fields'][i,j])
  if j==6+i:
   for e,c in rd(fd['Q0']).items():f[e]=f.get(e,ok.zero())-red['items'][i]['kappa']*c
  f=clean(f);ck('full_hat_bound_'+str(i)+'_'+str(j),all(weight(e)<=wt[j]+2 for e in f));f={e:c for e,c in f.items() if weight(e)==wt[j]+2};fs.append(f);ck('principal_field_source_'+str(i)+'_'+str(j),f==pt['full_principal_fields_original_coordinates'][i]['fields'][j])
 bs=[top(raw['base_fields'][i,j],2) for j in range(2)];principal.append((fs,bs))
 for j in range(2):ck('principal_base_source_'+str(i)+'_'+str(j),bs[j]==pt['full_principal_fields_original_coordinates'][i]['bases'][j])
print('global sources and six charts PASS',round(time.monotonic()-st,2),flush=True)
T=pt['ring'];v,x,y=T.gens();FT=T.fraction_field()
def reduce_v(f,sigma):
 out={}
 for ee,c in T(f).dict().items():
  e=tuple(map(int,ee));ne=(e[0]%2,e[1],e[2]);out[ne]=out.get(ne,K.zero())+c*sigma**(e[0]//2)
 return T(clean(out))
def pull(di,ch,derivative=None):
 vals=tuple(ch['p'])+tuple(ch['RR']);f=T.zero()
 for e,c in di.items():
  if any(n and value==0 for n,value in zip(e,vals)):continue
  term=T(cv(c) if derivative is None else dc(c,derivative))
  for n,value in zip(e,vals):term*=value**n
  f+=term
 return reduce_v(f,ch['v_squared'])
for m,it in enumerate(pt['items']):
 ch=it['chart'];ck('matching_chart_'+str(m),ch==hi['items'][m]['chart'] and ch==ex['items'][m]['chart']);ev=G.hom(list(ch['RR'])+list(ch['p']),T)
 plane=hi['items'][m]['plane_condition'];ck('necessary_quadratic_chart_'+str(m),reduce_v(ev(plane),ch['v_squared'])==0)
 pivotpair=ch['necessary_nonzero_factors'];ident=2*pivotpair[0]-pivotpair[1];ck('pivot_cannot_simultaneous_zero_'+str(m),ident==(-3*H if m==0 else (-3 if m==1 else 3)))
 # If either factor is zero then the nonzero two-RR support relation forces
 # both zero, contradicted by the preceding protected identity.
 for i,(fields,bases) in enumerate(principal):
  fe=[pull(f,ch) for f in fields];be=[pull(f,ch) for f in bases];ck('full_field_pullback_'+str(m)+'_'+str(i),all(f==g for f,g in zip(fe,it['fiber_fields_reduced'][i])));ck('base_chain_pullback_'+str(m)+'_'+str(i),r*be[0]==it['base_fields_reduced'][i][0] and r/q*be[1]==it['base_fields_reduced'][i][1])
  for j,f in enumerate(Rtops):
   ans=T.zero()
   for l in range(9):
    df={}
    for e,c in f.items():
     if e[l]:ne=list(e);ne[l]-=1;df[tuple(ne)]=e[l]*c
    ans+=pull(df,ch)*fe[l]
   for l in range(2):ans+=pull(f,ch,l)*be[l]
   ans=reduce_v(ans,ch['v_squared']);ck('genuine_higher_source_'+str(m)+'_'+str(j)+'_'+str(i),ans==hi['items'][m]['genuine_R_derivatives'][j]['values'][i])
 cert=ex['items'][m]
 for j,c in enumerate(cert['source_amplitude_certificates']):
  vals={it['source']:it['values'] for it in hi['items'][m]['genuine_R_derivatives']};fs=[vals[c['source_name']][i] for i in c['directions']];ck('amplitude_original_sources_'+str(m)+'_'+str(j),fs==c['original_reduced_higher_sources']);ck('amplitude_identity_'+str(m)+'_'+str(j),sum((co*f for co,f in zip(c['bezout_coefficients'],fs)),FT.zero())==c['amplitude'])
  for unit in c['source_units']:
   for po in (unit.numerator(),unit.denominator()):
    di=po.dict();ck('higher_unit_monomial_'+str(m)+'_'+str(j)+'_'+str(len(checks)),len(di)==1);ee,cc=next(iter(di.items()));ck('higher_unit_v_only_'+str(m)+'_'+str(j)+'_'+str(len(checks)),tuple(map(int,ee[1:]))==(0,0))
    for pol in (cc.numerator(),cc.denominator()):
     for fac,e in pol.factor():ck('higher_unit_protected_'+str(m)+'_'+str(j)+'_'+str(fac),cert['guard'].quo_rem(fac)[1]==0)
 print('entire two-RR chart',m,'PASS seconds',round(time.monotonic()-st,2),flush=True)
report={'status':'PASS','count':len(checks),'checks':checks,'seconds':time.monotonic()-st,'scope':'Original global square sources, complete six RR-support linear charts and old R factorizations, all original principal fields and all54 true derivative evaluations, six explicit amplitude identities independently replayed. For exactly two nonzero RR entries every nonzero-p mixed initial point is excluded on the stated generic base; p0 is the already treated RR plane. No support-invariance assumption, no real positivity, and no dropped quadratic-chart pivot.'}
(p/'independent_two_RR_mixed_exclusion_reader.json').write_text(json.dumps(report,indent=2)+'\n');save({'checks':checks},str(p/'independent_two_RR_mixed_exclusion_reader.sobj'));print({k:v for k,v in report.items() if k!='checks'},flush=True)
