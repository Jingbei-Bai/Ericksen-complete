"""Independent original-source reader for the subcritical RR2 missing-row atlas.

Only full original R/Q0/even polynomials and the complete signed H15 DAG
are used to construct the necessary scalar identities. No horizontal p
descent, numerical material specialization, or frozen-fiber no-lift is used.
"""
from sage.all import *
from pathlib import Path
import os, resource, time, json

p=Path(__file__).resolve().parent; st=time.time(); checks=[]
src=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'))
even=load(str(p/'independent_delta_subcritical_even_full_source.sobj'))
g=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'))
previous=load(str(p/'independent_subcritical_two_RR_inactive_rankone_classification.sobj'))
partial=load(str(p/'independent_partial_D_Newton_gauges_reader.sobj'))
vm=int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')
resource.setrlimit(resource.RLIMIT_AS,(vm+650*1024**2,vm+650*1024**2))
print('PID',os.getpid(),'VM_cap_MB',(vm//1024**2)+650,flush=True)
O=src['ring']; K=O.base_ring(); a,r,H,Z=K.gens(); q=1+r; m=1+a*r
al=(Z-H)/(Z*(H+1)); be=H*(Z-1)/(Z*(H+1))
P=PolynomialRing(K,names=('w','t')); w,t=P.gens(); F=P.fraction_field()
ev=O.hom([1,w,0,0,0,0,0,t,1],F)
flows=[{5:-1/al},{3:-be/(al*t)},{4:-1/be}]
ph=K.ring().hom([a,r,1,-1],K)
def ck(n,b):
 if not b: raise AssertionError(n)
 checks.append(n)
def cy(x):
 x=K(x); nn=ph(x.numerator()); dd=ph(x.denominator())
 ck('protected cyclic denominator '+str(len(checks)),dd!=0)
 return nn/dd
def face(f,i,j):
 return O({e:c for e,c in f.dict().items() if (sum(e[:6]),sum(e[6:]))==(i,j)})
def first(f):
 return [sum((F(c)*ev(f.derivative(O.gen(h))) for h,c in flow.items()),F.zero()) for flow in flows]

# Subcritical delta=0: all source constants are checked directly, not omitted.
rows={3:first(face(src['full_R_sources'][3],3,1)),6:first(face(src['full_Q0'],1,2))}
ck('full original R3 purep4 zero',ev(face(src['full_R_sources'][3],4,0))==0)
ck('full original Q0 top zero',ev(face(src['full_Q0'],2,1))==0)
old=g['coefficient_field']; pool=g['polynomial_pool']
mapping=old.ring().hom([a,a*q/m,H/r,q*Z/r],K); converted={}
def cv(c):
 if c not in converted: converted[c]=mapping(c.numerator())/mapping(c.denominator())
 return converted[c]
out={(3,2):{},(4,1):{}}
for sign,li,ri in g['records'][15]['full_signed_product_DAG']:
 for ex,c in pool[li].items():
  bx=(sum(ex[:6]),sum(ex[6:])); targets={(i-bx[0],j-bx[1]):(i,j) for i,j in out if i>=bx[0] and j>=bx[1]}
  if not targets: continue
  for ey,d in pool[ri].items():
   key=targets.get((sum(ey[:6]),sum(ey[6:])))
   if key is None: continue
   ee=tuple(int(x)+int(y) for x,y in zip(ex,ey)); dd=out[key]
   dd[ee]=dd.get(ee,old.zero())+sign*c*d
rawH={ij:O({e:cv(c) for e,c in dd.items() if c}) for ij,dd in out.items()}
rows[22]=first(rawH[(3,2)])
ck('full signed H15 purep4RR zero',ev(rawH[(4,1)])==0)
for j,row in rows.items(): ck('true source e2 column zero '+str(j),row[2]==0)
ck('Q0 e1 coefficient cyclic unit',cy(K(rows[6][1]))!=0)
minors=[]
for j in (3,22):
 det=matrix(F,[rows[6][:2],rows[j][:2]]).det()
 num=P(det.numerator()); den=P(det.denominator())
 ck('original two-row minor '+str(j),det==previous['range_minors'][0 if j==3 else 1]['minor'])
 ck('two-row numerator complete linear w,t '+str(j),set(tuple(map(int,e)) for e in num.dict())=={(1,0),(0,1)})
 ck('two-row denominator contains no state '+str(j),den.degree(w)==den.degree(t)==0)
 co={tuple(map(int,e)):c for e,c in num.dict().items()}
 minors.append({'source_rows':[6,j],'minor':det,'wco':co[(1,0)],'tco':co[(0,1)],'den':den})
Rdet=minors[0]['wco']*minors[1]['tco']-minors[1]['wco']*minors[0]['tco']
ck('range coefficient determinant exact',Rdet==previous['range_base_determinant'])
Rcy=cy(Rdet); ck('range coefficient determinant all materials',Rcy==K(25)/(4096*49))

# Reconstruct complete even equations; after e0=e1=0 and w unit use E=e2/w.
# The endpoint w=0 is retained before any division: all three square units
# force every finite e component to zero there.
EP=PolynomialRing(K,names=('w','t','e0','e1','e2')); ew,et,e0,e1,e2=EP.gens()
def evalp(f,values,ring):
 return sum((ring(c)*prod(ring(values[h])**int(ex[h]) for h in range(9)) for ex,c in f.dict().items()),ring.zero())
items=[]; complete=[]
for j,it in enumerate(even['items']):
 vals=[1,ew,0,0,0,0,0,et,1]
 sq=K(it['square_unit']); comp=int(it['square_component'])
 lin=[evalp(f,vals,EP) for f in it['A_p_quadratic']]
 pure=evalp(it['pure_p4'],vals,EP)
 ef=EP(sq)*[e0,e1,e2][comp]**2+sum(f*e for f,e in zip(lin,[e0,e1,e2]))+pure
 A=K(lin[2]/ew); B=K(pure/ew**2)
 ck('complete even support '+str(j),ef==sq*[e0,e1,e2][comp]**2+A*ew*e2+B*ew**2)
 ck('full even square protected '+str(j),cy(sq)!=0)
 ck('pure-axis endpoint e2 coefficient '+str(j),ef.subs({ew:0,e0:0,e1:0})==(sq*e2**2 if comp==2 else 0))
 items.append((sq,A,B)); complete.append(ef)
u0,A0,B0=items[0]; u1,A1,B1=items[1]
W1=A0/u0; W0=B0/u0; L1=-A1/u1; L0=-B1/u1
necessary=L0**2-W1*L0*L1+W0*L1**2
ck('independent complete even resultant exact',necessary==previous['even_rangezero_resultant'])
ck('base-only even resultant every material',cy(necessary)==1024*a**2*r**2)
E=PolynomialRing(K,'E').gen(); EE=E.parent()
ck('literal resultant multiplyback',EE(L1*E+L0).resultant(EE(E**2+W1*E+W0))==necessary)

# Exact inactive support enumeration from the original L binomials and D.
# These binomials have protected coefficients, so matching their two term
# supports is necessary even when their amplitudes are arbitrary complex.
L=src['L']; D=src['D']; support=[]
for rm in range(1,8):
 rr=[bool(rm&(1<<i)) for i in range(3)]
 if sum(rr)!=2: continue
 for pm in range(1,64):
  pp=[bool(pm&(1<<i)) for i in range(6)]
  def survives(f):
   return [e for e,c in f.dict().items() if c and all(pp[h] or not e[h] for h in range(6))]
  impossible=False
  for h in range(3):
   active=[]
   for j in range(3):
    if rr[j]: active += survives(L[h,j])
   if len(active)==1: impossible=True
  if impossible: continue
  dg=[any(survives(D[i,j]) for i in range(6)) for j in range(3)]
  if all(dg): continue
  missing=rr.index(False); nz=[h for h,b in enumerate(pp) if b]
  activegroups=[pp[0] or pp[1],pp[2] or pp[3],pp[4] or pp[5]]
  if sum(activegroups)==1:
   ck('rankone exactly missing row '+str((rm,pm)),activegroups[missing] and len(nz)<=2)
   typ='pureaxis' if len(nz)==1 else 'mixed_missing_row'
  else:
   ck('ranktwo exact paired support '+str((rm,pm)),sum(activegroups)==2 and len(nz)==2 and not activegroups[missing])
   typ='paired'
  olditem=next(x for x in partial['support_records'] if x['RR_mask']==rm and x['p_mask']==pm)
  if typ=='paired': ck('paired has two true available partial fields '+str((rm,pm)),len(olditem['available_horizontal_directions'])==2)
  support.append({'RR_mask':rm,'p_mask':pm,'kind':typ,'D_rank':sum(activegroups),'horizontal':olditem['available_horizontal_directions']})
ck('all three RR charts atlas counts',len(support)==12)
ck('paired ranktwo three charts',sum(x['kind']=='paired' for x in support)==3)
ck('pureaxis six charts',sum(x['kind']=='pureaxis' for x in support)==6)
ck('mixed missingrow three charts',sum(x['kind']=='mixed_missing_row' for x in support)==3)

save({'original_rows':rows,'raw_H15_layers':rawH,'minors':minors,'range_base_determinant':Rdet,'complete_even_equations':complete,'even_items':items,'even_resultant':necessary,'support_records':support,'checks':checks},str(p/'independent_subcritical_two_RR_inactive_rankone_raw_reader.sobj'),compress=False)
report={'PASS':True,'checks':len(checks),'seconds':time.time()-st,'range_base_determinant_cyclic':str(Rcy.factor()),'complete_even_resultant_cyclic':str(cy(necessary).factor()),'atlas_counts':{'paired':3,'pureaxis':6,'mixed':3},'no_horizontal_descent_on_rankone':True,'all_material_parameters_symbolic':True}
(p/'independent_subcritical_two_RR_inactive_rankone_raw_reader.json').write_text(json.dumps(report,indent=2)+'\n')
print('PASS',report,flush=True)
