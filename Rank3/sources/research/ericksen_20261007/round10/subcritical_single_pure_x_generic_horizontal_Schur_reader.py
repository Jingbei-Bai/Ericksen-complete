"""Full relative shifted conormal closure: all H,Z,t,l,kappa remain symbolic."""
from sage.all import *
from pathlib import Path
import resource,os,time,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'pure_x_root_generic_horizontal_normals.sobj'));hd=load(str(p/'pure_x_super_reader_generic_h_compact_certificate.sobj'))
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+550*1024**2,)*2);print('PID',os.getpid(),flush=True)
P=d['ring'];K=P.base_ring();k,r,H,Z=K.gens();q=1+r;m=1+k;al=(Z-H)/(Z*(H+1));be=H*(Z-1)/(Z*(H+1));shift=d['shift'];T=PolynomialRing(K,names=('t','l','kap'));t,l,kap=T.gens();oldv=(al+kap*shift)*l;olde0=shift*l;checks=[]
def ck(n,b):assert bool(b),n;checks.append(n)
src=[P(f) for f in d['sources']];ck('exact 29 complete source conormals',len(src)==29)
def normalrow(f):
 out=[T.zero() for unused in range(7)]
 for ee,c in f.dict().items():
  ck('exact entire horizontal family zero',sum(ee[:7])==1)
  jj=next(j for j in range(7) if ee[j]);out[jj]+=T({tuple(int(ee[j]) for j in (7,8,9)):c})
 return vector(T,out)
N=[normalrow(f) for f in src]
co=P(hd['incidence_cofactor']);cf=T.zero()
for ee,c in co.dict().items():
 if any(ee[j] for j in (0,1,2,4,5,6)):continue
 cf+=c*olde0**int(ee[3])*t**int(ee[7])*oldv**int(ee[8])*kap**int(ee[9])
vec=vector(K,hd['generic_combination']);multipliers=(t,T.one(),T.zero(),T.one(),T.zero(),oldv,t);sumN=vector(T,7)
for a,mu,(j,unused) in zip(vec,multipliers,hd['support']):sumN+=a*mu*N[int(j)]
incN=vector(T,(oldv,0,-al/be*t,0,-kap,0,0))
def kd(f):
 ck('whole conormal true kappa divisibility',all(ee[2]>=1 for ee in f.dict()))
 return T({tuple(int(a)-(1 if j==2 else 0) for j,a in enumerate(ee)):c for ee,c in f.dict().items()})
rawh=vector(T,[kd(f) for f in sumN-cf*incN]);correction=T.zero();hN=rawh
if all(ee[2]>=1 for ee in cf.dict()):
 correction=kd(cf);hN=rawh+correction*incN
 ck('actual necessary incidence correction identity',hN==vector(T,[kd(f) for f in sumN]))
print('H CORRECTION TERMS',len(correction.dict()),'CORRECTED FIRSTNORMAL TERMS',[len(f.dict()) for f in hN],flush=True)
N.append(hN);RR=(3,6,25,28);RC=(1,4,6,2);NC=(0,3,5);JR=matrix(T,[[n[j] for j in RC] for n in N]);JN=matrix(T,[[n[j] for j in NC] for n in N]);A=JR.matrix_from_rows(RR);B=JN.matrix_from_rows(RR)
print('RANGE ENTRY TERMS',[[len(f.dict()) for f in row] for row in A],'normal terms',[[len(f.dict()) for f in row] for row in B],flush=True)
det=A.det();print('RANGE DET TERMS',len(det.dict()),'sec',time.monotonic()-st,flush=True)
adj=A.adjugate();Schur=det*JN.matrix_from_rows((26,27,29))-JR.matrix_from_rows((26,27,29))*adj*B
for i in range(3):
 for j in range(3):ck('COMPLETE GENERIC horizontal Schur exactly zero '+str((i,j)),Schur[i,j]==0)
ck('complete raw corrected h normal recorded',hN==N[29])
CB=PolynomialRing(QQ,names=('k','r')).fraction_field();kk,rr=CB.gens();BT=PolynomialRing(CB,names=T.variable_names())
for hf in (1,2):
 ph=K.ring().hom([kk,rr,hf,-1],CB)
 def sp(f):return BT({tuple(map(int,ee)):ph(c.numerator())/ph(c.denominator()) for ee,c in f.dict().items()})
 expect=-32*(1+rr)*kk/(1701*rr*(1+kk)**2) if hf==2 else -3*(1+rr)*kk/(224*rr*(1+kk)**2)
 ck('relative range unit specializes to protected each-material H'+str(hf),sp(det).constant_coefficient()==expect)
 ck('incidence correction leaves actual H'+str(hf)+' h source unchanged',sp(correction)==0)
 for j,a in enumerate(vec):ck('generic h coefficient denominator is protected H'+str(hf)+' '+str(j),ph(a.denominator()) in QQ and ph(a.denominator())!=0)
report={'status':'PASS','checks':len(checks),'full_relative_H_Z_t_l_kappa':True,'true_generic_four_range_Schur_zero':True,'necessary_h_incidence_correction_terms':len(correction.dict()),'range_determinant_terms':len(det.dict()),'H1_H2_each_material_range_units':True,'H1_H2_h_source_unchanged':True,'seconds':time.monotonic()-st,'scope':__doc__}
save({'ring':P,'coefficient_ring':T,'shift':shift,'complete_source_conormals':N,'h_actual_incidence_correction':correction,'h_original_conormal':rawh,'h_corrected_conormal':hN,'range_rows':RR,'range_cols':RC,'normal_cols':NC,'range_jacobian':A,'range_normal_jacobian':B,'range_determinant':det,'complete_generic_Schur':Schur,'checks':checks,'report':report,'scope':__doc__},str(p/'subcritical_single_pure_x_generic_horizontal_Schur_reader.sobj'),compress=False)
(p/'subcritical_single_pure_x_generic_horizontal_Schur_reader.json').write_text(json.dumps(report,indent=2)+'\n');print(report,flush=True)
