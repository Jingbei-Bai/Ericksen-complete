from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(n,v):
 if not bool(v):raise AssertionError(n)
 checks.append(n)
old=load(str(p/'independent_ratio_canonical_amplitude_normal.sobj'));new=load(str(p/'independent_second_prolongation_reader.sobj'));compact=load(str(p/'independent_ratio_canonical_sumzero_five_minor_compact.sobj'))
K=old['coefficient_field'];a,r=old['parameters'];q=1+r
R=PolynomialRing(K,names=('audit_x','audit_y'));x,y=R.gens()
def pull(f):
 out=R.zero()
 for e,c in f.dict().items():
  e=tuple(map(int,e));out+=K(c)*x**e[0]*y**e[1]*(-x-y)**e[2]
 return out
M=matrix(R,[[pull(f) for f in row] for mat in (old['matrix'],new['normal_matrix']) for row in mat.rows()]);ck('91x6_reconstructed',M.dimensions()==(91,6))
saved=load(str(p/'independent_ratio_canonical_sumzero_second_normal.sobj'))['matrix']
def reb2(f):return R({tuple(map(int,e)):K(c) for e,c in f.dict().items()})
for i in range(91):
 for j in range(6):ck('full91_original_sumzero_%d_%d'%(i,j),M[i,j]==reb2(saved[i,j]))
N=M.matrix_from_columns([0,1,2,4,5]);targets=[x*(x+y),y*(x+y),x*x,x*y];allminors={};records=[]
for n,it in enumerate(compact['items']):
 seed=tuple(map(int,it['seed']));A=N.matrix_from_rows(seed);w=vector(R,[(-1)**j*A.matrix_from_columns([k for k in range(5) if k!=j]).det() for j in range(5)])
 ck('signed_cofactor_wedge_kernel_'+str(n),all(v==0 for v in A*w));value=R.zero()
 for i,c in it['terms']:
  c=K(c);ck('old_sources_only_%d_%d'%(n,i),int(i)<28);ck('coefficient_guard_%d_%d'%(n,i),c.denominator().divides(K.ring()(r**4*q**4*12**4)))
  dm=sum(N[i,j]*w[j] for j in range(5));direct=N.matrix_from_rows(seed+(int(i),)).det();ck('signed_minor_direct_%d_%d'%(n,i),dm==direct);allminors[(seed,int(i))]=dm;value+=c*dm
 ck('five_minor_target_'+str(n),value==targets[n]);ck('stored_target_'+str(n),reb2(it['target'])==targets[n]);records.append({'seed':seed,'target':str(targets[n]),'source_rows':[int(i) for i,c in it['terms']]})
print('FOUR ORIGINAL FIVE MINORS PASS',len(checks),round(time.monotonic()-st,2),flush=True)
axis=load(str(p/'independent_ratio_canonical_sumzero_axis.sobj'));wall=load(str(p/'independent_ratio_canonical_sumzero_axis_spectral_wall.sobj'))
Y=PolynomialRing(K,names=('audit_y_axis','audit_unused_axis'));yy=Y.gen(0)
def toaxis(f):return sum((K(c)*yy**int(e[1]) for e,c in f.dict().items() if int(e[0])==0),Y.zero())
MA=matrix(Y,91,6,[toaxis(f) for f in M.list()])
def reby(f):return Y({(int(e[0]),0):K(c) for e,c in f.dict().items()})
for i in range(91):
 for j in range(6):ck('axis_original_M_%d_%d'%(i,j),MA[i,j]==reby(axis['normal_matrix'][i,j]))
it=compact['items'][1];v=sum((K(c)*MA.row(i) for i,c in it['terms']),vector(Y,[0]*6));B=matrix(Y,[MA.row(i) for i in it['seed']]+[v]);BNT=B.matrix_from_columns([0,1,2,4,5]);ck('axis_non_tau_pivot',BNT.det()==yy**2)
for i in range(5):
 for j in range(6):ck('axis_source_block_%d_%d'%(i,j),B[i,j]==reby(axis['five_source_block'][i,j]))
nums=[]
for j in range(5):
 C=matrix(Y,BNT);C.set_column(j,-B.column(3));nums.append(C.det())
ln=vector(Y,[nums[0],nums[1],nums[2],yy**2,nums[3],nums[4]]);ck('axis_cramer_kernel',all(v==0 for v in B*ln))
for j in range(6):ck('axis_cramer_original_'+str(j),ln[j]==reby(axis['lift_numerator'][j]))
needed=sorted(set(int(i) for i,c in axis['exact_identity_terms'])|set(int(i) for i,c in wall['exact_identity_terms']))
res={}
for i in needed:
 res[i]=sum((MA[i,j]*ln[j] for j in range(6)),Y.zero());ck('axis_residual_original_'+str(i),res[i]==reby(axis['residuals'][i]))
print('ORIGINAL CRAMER/SELECTED RESIDUALS PASS',len(checks),round(time.monotonic()-st,2),flush=True)
coeffden=[];val=Y.zero()
for i,c in axis['exact_identity_terms']:
 c=K(c);coeffden.append(str(c.denominator().factor()));val+=c*res[int(i)]
 for f,e in c.denominator().factor():ck('generic_axis_coefficient_guard_%d_%s'%(i,f),any(K(f)/v in QQ for v in (r,q,a,a-1,1+a*r,a*r+r+1)))
ck('generic_axis_residual_identity',val==yy**2)
print('GENERIC AXIS PASS',len(checks),coeffden,round(time.monotonic()-st,2),flush=True)
P=PolynomialRing(QQ,'wall_r');rr=P.gen();L=P.fraction_field();theta=K.ring().hom([-(1+rr)/rr,rr],L)
def matmap(c):c=K(c);ck('wall_original_source_den_nonzero',theta(c.denominator())!=0);return theta(c.numerator())/theta(c.denominator())
YW=PolynomialRing(L,'wall_y');yw=YW.gen()
def wallpoly(f):return YW({int(e[0]):matmap(c) for e,c in f.dict().items()})
MW=matrix(YW,[[wallpoly(MA[i,j]) for j in range(6)] for i in needed]);lw=vector(YW,[wallpoly(f) for f in ln]);rw={}
for ii,i in enumerate(needed):
 rw[i]=sum((MW[ii,j]*lw[j] for j in range(6)),YW.zero());ck('wall_residual_specialization_'+str(i),rw[i]==wallpoly(res[i]))
val=YW.zero();wden=[];Q=rr**4+2*rr**3-QQ(269)/32*rr**2-QQ(301)/32*rr-QQ(227)/64
for i,c in wall['exact_identity_terms']:
 c=L(c);wden.append(str(c.denominator().factor()));val+=c*rw[int(i)]
 for f,e in c.denominator().factor():ck('wall_identity_den_%d_%s'%(i,f),any(L(f)/v in QQ for v in (rr,rr+1,Q)))
ck('wall_residual_identity',val==yw**2)
u=(rr+QQ(1)/2)**2;ck('wall_positive_identity',-128*Q==145+1236*u+128*u*(QQ(1)/4-u));aa=-(1+rr)/rr;bb=aa*(1+rr)/(1+aa*rr);ck('wall_b_equals_a_squared',bb==aa**2)
for j,c in enumerate(wall['parameter_map']):ck('wall_parameter_map_'+str(j),L(c)==[aa,rr][j])
report={'status':'PASS','count':len(checks),'seconds':time.monotonic()-st,'five_minor_records':records,'axis_generic_terms':[(int(i),str(c)) for i,c in axis['exact_identity_terms']],'axis_generic_denominators':coeffden,'axis_wall_terms':[(int(i),str(c)) for i,c in wall['exact_identity_terms']],'axis_wall_denominators':wden,'wall_Q':str(Q),'scope':'All91 normal rows rebuilt from independently audited old28 and new63. Four exact5-minor source identities prove non-tau rank5 away from origin on sum(t)=0 for all protected spectra. Axis x=0,y!=0 full6-rank proved by original-source Cramer residual identity, separately on ar+r+1!=0 and ar+r+1=0. On positive-material wall -1<r<0, -128Q=145+1236u+128u(1/4-u)>0 for u=(r+1/2)^2. Cyclic base only; no assertion for x!=0 full6-rank or whole rank3.'}
(p/'independent_sumzero_short_certificate_reader.json').write_text(json.dumps(report,indent=2)+'\n');save({'checks':checks,'records':records,'normal_matrix':M,'axis_matrix':MA,'axis_lift':ln,'axis_residuals':res},str(p/'independent_sumzero_short_certificate_reader.sobj'));print(report,flush=True)
