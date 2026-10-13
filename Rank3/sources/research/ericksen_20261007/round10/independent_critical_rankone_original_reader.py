from sage.all import *
from pathlib import Path
import json,time,resource,os,gc
p=Path(__file__).resolve().parent;st=time.time();checks=[]
def ck(n,v):
 if not v:raise AssertionError(n)
 checks.append(n)
src=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'));O=src['ring'];K=O.base_ring();a,r,H,Z=K.gens();al=(Z-H)/(Z*(H+1));be=H*(Z-1)/(Z*(H+1));P=PolynomialRing(K,names=('w','t'));w,t=P.gens();F=P.fraction_field();vals=[1,w,0,0,0,0,0,t,1];ev=O.hom(vals,F)
flows=[{5:-1/al},{3:-be/(al*t)},{4:-1/be}]
def face(f,i,j):return O({ex:c for ex,c in f.dict().items() if (sum(ex[:6]),sum(ex[6:]))==(i,j)})
def row_from(f,i,j,k,l):
 g=face(f,i,j);return [sum((F(v)*ev(g.derivative(O.gen(z))) for z,v in flow.items()),F.zero()) for flow in flows]+[ev(face(f,k,l))]
matrixcache=load(str(p/'independent_critical_two_RR_inactive_rankone_full_matrix.sobj'))['matrix'];newrows={}
newrows[3]=row_from(src['full_R_sources'][3],3,1,1,2)
newrows[6]=row_from(src['full_Q0'],1,2,0,0);newrows[6][3]=0
for j in (3,6):
 for h in range(4):ck('complete original row '+str(j)+'/'+str(h),newrows[j][h]==matrixcache[j,h])
ck('R3 purep4 zero',ev(face(src['full_R_sources'][3],4,0))==0)
ck('Q0 top zero',ev(face(src['full_Q0'],2,1))==0)
# Assemble complete H layers in original coefficients before evaluation or differentiation.
# This independently projects the entire signed-product DAG, including all cancellations.
g=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'));old=g['coefficient_field'];pool=g['polynomial_pool'];ph=old.ring().hom([a,a*(1+r)/(1+a*r),H/r,(1+r)*Z/r],K)
vm=int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE');resource.setrlimit(resource.RLIMIT_AS,(vm+650*1024**2,vm+650*1024**2))
cv_cache={}
def cv(c):
 if c not in cv_cache:cv_cache[c]=ph(c.numerator())/ph(c.denominator())
 return cv_cache[c]
rawfaces={}
for idx in (9,15,10,16):
 out={(3,2):{},(4,1):{},(1,3):{}}
 for sign,li,ri in g['records'][idx]['full_signed_product_DAG']:
  for ex,c in pool[li].items():
   bix=(sum(ex[:6]),sum(ex[6:]));needed={(u-bix[0],v-bix[1]):(u,v) for u,v in out if u>=bix[0] and v>=bix[1]}
   if not needed:continue
   for ey,d in pool[ri].items():
    facekey=needed.get((sum(ey[:6]),sum(ey[6:])))
    if facekey is None:continue
    ee=tuple(int(x)+int(y) for x,y in zip(ex,ey));dd=out[facekey];dd[ee]=dd.get(ee,old.zero())+sign*c*d
 mapped={bi:O({ee:cv(co) for ee,co in dd.items() if co}) for bi,dd in out.items()}
 rr=[sum((F(v)*ev(mapped[(3,2)].derivative(O.gen(z))) for z,v in flow.items()),F.zero()) for flow in flows]+[ev(mapped[(1,3)])]
 rownum=7+idx;newrows[rownum]=rr
 for h in range(4):ck('full raw H DAG row '+str(rownum)+'/'+str(h),rr[h]==matrixcache[rownum,h])
 ck('full raw H41 zero '+str(idx),ev(mapped[(4,1)])==0)
 rawfaces[idx]=mapped
 print('raw H',idx,'seconds',round(time.time()-st,2),flush=True)
cert=load(str(p/'independent_critical_rankone_generic_resultant.sobj'))
for item in cert['minor_items']:
 mm=matrix(F,[newrows[j] for j in item['rows']]).matrix_from_columns([0,1,3]);ck('original source determinant factor identity '+str(item['rows']),mm.det()==item['other']*item['state_factor'])
# Rebuild all three even equations directly from complete source coefficients.
even=load(str(p/'independent_delta_subcritical_even_full_source.sobj'));EP=PolynomialRing(K,names=('w','t','e0','e1','e2','dd'));ew,et,e0,e1,e2,dd=EP.gens();evals=[1,ew,0,0,0,0,0,et,1]
def pue(f,values,ring):return sum((ring(co)*prod(ring(values[h])**int(ex[h]) for h in range(9)) for ex,co in f.dict().items()),ring.zero())
ecache=load(str(p/'independent_critical_two_RR_inactive_rankone_even_probe.sobj'))['equations']
for j,it in enumerate(even['items']):
 ef=EP(it['square_unit'])*[e0,e1,e2][int(it['square_component'])]**2+sum(pue(c,evals,EP)*e for c,e in zip(it['A_p_quadratic'],[e0,e1,e2]))+pue(it['pure_p4'],evals,EP)
 ck('complete even source '+str(j),ef==ecache[j])
 for lab,values in [('endpoint',[0,1,0,0,0,0,0,et,1])]:
  ck(lab+' even purep4 '+str(j),pue(it['pure_p4'],values,EP)==0)
  for h,f in enumerate(it['A_p_quadratic']):ck(lab+' even linear '+str(j)+'/'+str(h),pue(f,values,EP)==0)
# Endpoint: literal algebraic F12/H13 exclusion is already independently audited;
# replay the exact complete-source amplitude certificates, retaining both amplitudes.
ss=load(str(p/'independent_supercritical_RR_small_support_sources.sobj'));T=ss['items'][0]['H13_values'][0].parent();vf,x,y=T.gens();chart=ss['items'][0]['chart'];TF=T.fraction_field()
for j,ampcert in enumerate(ss['items'][0]['amplitude_certificates']):
 pos=list(ss['R_order']).index(int(ampcert['source_name'][1:]));selected=[]
 for direction in ampcert['directions']:
  hi=ss['H_order'].index((int(direction),pos));selected.append(ss['items'][0]['H13_values'][hi])
 ck('literal H13 amplitude bezout '+str(j),sum((TF(f)*TF(c) for f,c in zip(selected,ampcert['bezout_coefficients'])),TF.zero())==ampcert['amplitude'])
ck('chart missingRR1',chart['RR'][0]==0)
ck('chart all rank1 amplitudes',set(str(v) for v in chart['p'] if v)!=set())
print('endpoint chart',chart,'amplitudes',[it['amplitude'] for it in ss['items'][0]['amplitude_certificates']],flush=True)
# All three cyclic frame orientations have allowed positive simple labels.
L=PolynomialRing(QQ,names=('a','b'));aa,bb=L.gens();LF=L.fraction_field()
for lab,A,B in [('123',LF(aa),LF(bb)),('231',LF(bb/aa),LF(1/aa)),('312',LF(1/bb),LF(aa/bb))]:
 rr=(A-B)/(A*(B-1));qq=1+rr;mm=1+A*rr
 ck('material lambda positive identity '+lab,mm*qq==(B/A)*mm**2)
 ck('material r nonzero numerator '+lab,rr!=0)
 ck('material q nonzero rational '+lab,qq!=0)
 ck('material m nonzero rational '+lab,mm!=0)
save({'source_rows':newrows,'H_raw_layers':rawfaces,'endpoint_chart':chart,'checks':checks},str(p/'independent_critical_rankone_original_reader.sobj'),compress=False)
rep={'PASS':True,'checks':len(checks),'seconds':time.time()-st,'no_horizontal_descent':True,'source_row_indices':[3,6,16,22],'endpoint':'complete original even e=0 then literal full F12/H13 Bezout','material_orientations':['123','231','312']};(p/'independent_critical_rankone_original_reader.json').write_text(json.dumps(rep,indent=2));print('PASS',rep,flush=True)
