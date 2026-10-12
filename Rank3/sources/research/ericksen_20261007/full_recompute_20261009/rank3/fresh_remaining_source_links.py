"""Close polynomial-to-matrix and critical-direction intermediate links."""
from sage.all import *
from pathlib import Path
import json,time
out=Path(__file__).absolute().parent;research=out.parent.parent;st=time.monotonic();checks=[]
def ck(n,v):
 if not v:raise AssertionError(n)
 checks.append(n)
def old(name):return load(str(research/'round10'/(name+'.sobj')))
def regenerated(script,name=None):return load(str(out/'source_replays'/script/'generated/round10'/((name or script)+'.sobj')))
def bind_generic_w_sources(fg,joint,bridge,pairs):
 R=fg['ring'];K=R.base_ring();a,b,h,z=K.gens();x=R.gens()
 r=(a-b)/(a*(b-1));q=1+r;H=r*h;Z=r*z/q
 u=2*(H-Z)/(H*Z)
 sources=[R(f) for f in fg['sources']]
 F=x[0]*vector(R,list(fg['fields'].row(0)))
 F[6]-=u*sources[3]
 B=x[0]*vector(R,list(fg['base_fields'].row(0)))
 ck('generic W complete first partial gauge dimensions',len(F)==9 and len(B)==2)
 full_pairs=[x[5]*sources[0]+x[4]*sources[1],x[3]*sources[2]+x[2]*sources[5],x[1]*sources[3]-x[0]*sources[4]]
 ck('complete three global pair source order',[int(it['pair']) for it in pairs['items']]==[0,1,2])
 for i,it in enumerate(pairs['items']):
  target=it['full'];S=target.parent();SK=S.base_ring();ph=K.ring().hom(SK.gens(),SK)
  mapped=S({tuple(map(int,e)):ph(c.numerator())/ph(c.denominator()) for e,c in full_pairs[i].dict().items()})
  ck('complete global pair identity from fresh six R '+str(i),mapped==target)
 polys=[R(joint['quartic_source_rows'][5]),full_pairs[1],full_pairs[2]]
 ck('generic W three original source names',[it['name'] for it in bridge['items']]==['GQ5','GP1','GP2'])
 for i,(f,it) in enumerate(zip(polys,bridge['items'])):
  target=it['original_derivative'];S=target.parent();SK=S.base_ring();k,rr,HH,ZZ=SK.gens()
  images=(k/rr,k*(1+rr)/(rr*(1+k)),HH/rr,(1+rr)*ZZ/rr)
  phi=K.ring().hom(images,SK)
  cache={}
  def coefficient(c):
   c=K(c)
   if c not in cache:cache[c]=phi(c.numerator())/phi(c.denominator())
   return cache[c]
  def mapped(poly):return S({tuple(map(int,e)):coefficient(c) for e,c in poly.dict().items()})
  fm=mapped(f);Fm=[mapped(value) for value in F];Bm=[mapped(value) for value in B]
  print('Derive complete fresh generic W source',it['name'],flush=True)
  actual=sum((Fm[j]*fm.derivative(S.gen(j)) for j in range(9)),S.zero())
  for j,(param,scale) in enumerate(((HH,rr),(ZZ,rr/(1+rr)))):
   df=S({tuple(map(int,e)):scale*c.derivative(param) for e,c in fm.dict().items() if c.derivative(param)})
   actual+=Bm[j]*df
  ck('generic W whole fresh derivative '+it['name'],actual==target)
  ck('generic W exact scale offsets '+it['name'],tuple(it['offset'])==((13,6) if i==0 else (14,6)))
  for e in actual.dict():
   dg,rd=sum(map(int,e)),sum(map(int,e[6:]))
   ck('generic W complete nonnegative delta exponent '+it['name']+' '+str((dg,rd)),6-dg>=0)
  print('Fresh generic W derivative identity',it['name'],'terms',len(actual.dict()),flush=True)
name='independent_ratio_RR_cyclic_full_matrix';a=regenerated(name);b=old(name)
for key in ['matrix','row_scales','old_parameter_map','guard']:
 ck('fresh cyclic full28 '+key,a[key]==b[key])
# An even polynomial and its coefficient matrix are separate source objects.
d=old('independent_critical_omega_even_matrix');P=d['amplitude_ring'];R=d['coefficient_ring'];aa,xx,TT=P.gens();ra,rx=R.gens()
for item in d['records']:
 rows=[]
 for poly in item['polynomials']:
  row=[R.zero() for _ in range(3)]
  for e,c in poly.dict().items():
   ck('omega full degree at most2',int(e[2])<=2)
   row[2-int(e[2])]+=R(c)*ra**int(e[0])*rx**int(e[1])
  rows.append(row)
 ck('omega literal coefficient matrix '+item['direction_type'],matrix(R,rows)==item['matrix'])
# Exact source coefficients were reconstructed before these creators were run.
for name,keys in [
 ('independent_critical_two_RR_R_coefficients',['items']),
 ('independent_critical_two_RR_R_scaled',['items']),
]:
 a=regenerated(name);b=old(name)
 for key in keys:ck('fresh twoRR source '+name+' '+key,a[key]==b[key])
name='independent_critical_two_RR_direction_certificate';a=regenerated(name);scaled=regenerated('independent_critical_two_RR_R_scaled')
for item,src in zip(a['items'],scaled['items']):
 ck('fresh QQ witness source equals coefficients',item['sources']==src['QQ_sources'])
 for w in item['items']:ck('fresh exact direction membership',sum(c*f for c,f in zip(w['cofactors'],item['sources']))==w['target'])
a=regenerated(name,'independent_critical_two_RR_direction_H_probe');b=old('independent_critical_two_RR_direction_H_probe')
for j,(v,w) in enumerate(zip(a['items'],b['items'])):
 ck('all25 critical direction rows case'+str(j),v['matrix']==w['matrix'])
 ck('critical direction defining square case'+str(j),v['square']==w['square'])
 for rec in v['first_nonzero']:
  ck('fresh exact source determinant',v['matrix'].matrix_from_rows(rec['rows']).det()==rec['determinant'])
  ck('fresh direction quadratic remainder',rec['determinant']%(v['ring'].gen()**2-v['square'])==rec['remainder'])
name='pure_x_super_reader_generic_h_kappa_bridge';a=regenerated(name);b=old(name)
for key in ['source_first_two','support','generic_combination','incidence_cofactor','true_h_kappa_zero_divided_source','linear','quadratic','selected_minor','selected_monomials']:
 ck('whole generic h source linkage '+key,a[key]==b[key])
bind_generic_w_sources(load(str(out/'fresh_global_polynomials.sobj')),load(str(out/'fresh_joint_fields.sobj')),old('pure_x_root_generic_w_bridge'),old('independent_ratio_p_boundary_global_pair_syzygies'))
fg=load(str(out/'fresh_global_polynomials.sobj'));raw=old('independent_supercritical_first_boundary_probe');R=raw['ring'];K=R.base_ring();aa,rr,HH,ZZ=K.gens();OK=fg['ring'].base_ring();phi=OK.ring().hom([aa,aa*(1+rr)/(1+aa*rr),HH/rr,(1+rr)*ZZ/rr],K);cache={}
def kc(c):
 c=OK(c)
 if c not in cache:cache[c]=phi(c.numerator())/phi(c.denominator())
 return cache[c]
def source_map(poly,bi=None):return R({tuple(map(int,e)):kc(c) for e,c in poly.dict().items() if bi is None or (sum(e[:6]),sum(e[6:]))==bi})
for j,p in enumerate(fg['sources']):ck('fresh firstboundary full F12 '+str(j),source_map(p,(1,2))==raw['F12'][j])
small=old('independent_supercritical_RR_small_support_sources')
for j in range(18):
 row=load(str(out/'critical_rows'/('row'+str(j)+'.sobj')));p=row['source_faces'][1,3]
 ck('fresh firstboundary full H13 '+str(j),source_map(p)==raw['H13'][j])
 ck('fresh RRsmall-support original H13 '+str(j),p==fg['ring'](small['H13_original'][j]))
F=PolynomialRing(QQ,names=('mat_a','mat_b')).fraction_field();a,b=F.gens();r=(a-b)/(a*(b-1));q=1+r;m=1+a*r
ck('omega material sign first v',a*r**2!=0)
ck('omega material sign second v',(-q)*(a*(1-q)-1)==m*q)
X=-q/r;ck('omega material sign third v',a*X*(a-X-1)==-a*m*q/r**2)
ck('omega protected mq positive-square identity',m*q==b*(a-1)**2/(a*(b-1)**2))
ck('omega c both-conjugate rational lower bound',232135214-27335384*QQ(13)/2==54455218)
record=dict(status='PASS',checks=checks,count=len(checks),seconds=time.monotonic()-st,scope='Full cyclic normal matrix regeneration, all omega polynomial-to-coefficient-matrix links, complete critical twoRR direction-coefficient/scaling/source-probe closure. Actual original-to-twoRR25x5 linkage is in fresh_two_RR_sources.json.')
(out/'fresh_remaining_source_links.json').write_text(json.dumps(record,indent=2));print('PASS',len(checks),record['seconds'],flush=True)
