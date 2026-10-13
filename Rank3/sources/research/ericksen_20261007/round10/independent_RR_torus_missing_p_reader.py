"""Independent source/minor replay for all three missing-p RR-torus supports."""
from sage.all import *
from pathlib import Path
import time,json,itertools
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(n,v):
 if not bool(v):raise AssertionError(n)
 checks.append(n)
src=load(str(p/'independent_RR_torus_missing_p_sources.sobj'));amp=load(str(p/'independent_RR_torus_missing_p_critical_amplitudes.sobj'));nor=load(str(p/'independent_RR_torus_missing_p_super_normal.sobj'));cl=load(str(p/'independent_RR_torus_missing_p_closure.sobj'));old=src['field'];oa,orr,oH,oZ,ox,ob,oc=old.gens();K=nor['field'];a,r,H,Z,y=K.gens();q=1+r;P=nor['root_ring'];b,c=P.gens();AP=amp['amplitude_ring'];AA=amp['field'];aa,ar,ay=AA.gens();ab,ac=AP.gens();U=amp['direction_ring'];uy=U.gen()
def const(v):return v.numerator().total_degree()==0 and v.denominator().total_degree()==0
def red(f):
 out={}
 for ex,co in P(f).dict().items():
  e=tuple(map(int,ex));key=(e[0]%2,e[1]%2);out[key]=out.get(key,K.zero())+co*(-q)**(e[0]//2)*(-q/r)**(e[1]//2)
 return P(out)
for ix,(si,ai,ni,ci) in enumerate(zip(src['items'],amp['items'],nor['items'],cl['items'])):
 label=si['label'];M=si['matrix'];rs=si['generic_pivot_rows'][:3];C=M.matrix_from_rows_and_columns(rs,(0,1,2));det=C.det();adj=C.adjugate();ck('critical original block adjugate '+label,C*adj==det*identity_matrix(old,3));rd=vector(old,[M[i,3] for i in rs]);ro=vector(old,[M[i,4] for i in rs]);qb=vector(old,M[6,:3].list());qt=det*M[6,3]-qb*adj*rd;qz=det*M[6,4]-qb*adj*ro
 # Rebuild the actual five-source affine minors and their full coefficient rows.
 matrows=[]
 for rec in ai['sources']:
  j=rec['original_rows'][-1];hb=vector(old,M[j,:3].list());ht=det*M[j,3]-hb*adj*rd;hz=det*M[j,4]-hb*adj*ro;f=(qt*hz-qz*ht)/det;ck('original five-source minor '+str((label,j)),f==rec['original_determinant'])
  mp=rec['mapped_determinant'];MK=mp.parent();ma,mr,mH,mZ,my,mb,mc=MK.gens();xm=[my/mb,my,my*mb][ix];ph=old.ring().hom([ma,mr,mH,mZ,xm,mb,mc],MK);ck('entire amplitude/direction transport',ph(f.numerator())/ph(f.denominator())==mp)
  bh=AA.ring().hom([ma,mr,my],MK)
  def coeff(v):v=AA(v);return bh(v.numerator())/bh(v.denominator())
  def polpull(v):return sum((coeff(co)*mb**int(ex[0])*mc**int(ex[1]) for ex,co in v.dict().items()),MK.zero())
  rebuilt=MK.zero()
  for key,g in rec['base_coefficients']:
   rebuilt+=mH**key[0]*mZ**key[1]*polpull(g*rec['removed_amplitude_unit']);matrows.append([g.monomial_coefficient(v) for v in (AP.one(),ab**2,ac**2)])
  ck('full base coefficient reconstruction',rebuilt==mp.numerator())
 ck('all critical amplitude coefficients',matrix(AA,matrows)==ai['amplitude_matrix']);NN=ai['amplitude_matrix'];ck('canonical amplitude kernel',NN*ai['amplitude_kernel']==0)
 # Only source pairs actually used in the Bezout identity are required.
 terms=[]
 for rec,co in zip(ai['direction_minor_records'],ai['direction_bezout']):
  if not co:continue
  i,j=rec['pair'];f=NN[i,1]*NN[j,2]-NN[i,2]*NN[j,1];ck('used genuine two-minor',f==rec['minor']);uh=U.hom([ay],AA);ck('rank certificate unit replay',f==uh(rec['direction_part'])*rec['removed_material_unit'])
  for fac,e in rec['removed_material_unit'].factor():ck('only protected material units in amplitude rank',const(AA(fac)) or any(const(AA(fac)/v) for v in (aa,ar,ar+1)))
  terms.append(co*rec['direction_part'])
 ck('complete direction Bezout',sum(terms,U.zero())==ai['direction_gcd']);ck('only support direction zero',ai['direction_gcd']==uy**(6 if ix<2 else 4))
 # Independently differentiate before imposing either root relation.
 xf=[-y*b/q,P(y),y*b][ix];ph=old.ring().hom([P(a),P(r),P(H),P(Z),xf,b,c],P)
 def cv(v):
  v=old(v);nn=red(ph(v.numerator()));dd=red(ph(v.denominator()));ck('protected finite-cover denominator',len(dd.dict())==1);ex,co=next(iter(dd.dict().items()));i,j=map(int,ex);return red(nn*b**i*c**j/(co*(-q)**i*(-q/r)**j))
 for i in range(25):
  f=M[i,3];ck('entire original delta column vanishes canonically',cv(f)==0)
  raw=[cv(ob*f.derivative(ob)+[-1,0,1][ix]*ox*f.derivative(ox)),cv(oc*f.derivative(oc))]+[cv(M[i,j]) for j in (0,1,2,4)]
  ri,rj=ni['row_phases'][i];iv=b**ri*c**rj/((-q)**ri*(-q/r)**rj)
  for j,v in enumerate(raw):
   pi,pj=ni['column_phases'][j];ck('literal complete source normal entry',red(v*b**pi*c**pj*iv)==P(ni['matrix'][i,j]))
 N=ni['matrix'];rows=ci['critical_source_rows'];cr=N.matrix_from_rows_and_columns(rows,(2,3,4,5));ck('literal critical four-row determinant',cr.det()==ci['critical']['determinant'])
 ir=ci['super_transverse_rows'];T=N.matrix_from_rows_and_columns(ir,(0,1));td=T.det();V=ci['super_reduced_matrix'];ck('literal supercritical Schur identity',td*V==td*N[:,2:]-N[:,:2]*T.adjugate()*N.matrix_from_rows_and_columns(ir,(2,3,4,5)))
 for rec in ci['super_minors']:
  rr=list(rec['source_rows']);br=rr[2:5];j=rr[-1];CC=V.matrix_from_rows_and_columns(br,(0,1,2));co=vector(K,[V[i,3] for i in br]);ff=td*(CC.det()*V[j,3]-vector(K,V[j,:3].list())*CC.adjugate()*co);ck('literal six-row source determinant',ff==rec['determinant'])
 for rec in [ci['critical']]+ci['super_minors']:
  ck('no discarded determinant factor',K(rec['core'])*rec['removed_unit']==rec['determinant']);phu=cl['direction_ring'].hom([y],K);z0=sum((H**key[0]*Z**key[1]*phu(v) for key,v in rec['base_coefficients'].items()),K.zero());ck('all original base coefficients retained',z0==K(rec['core']))
 CP=cl['direction_ring'];cg=list(ci['critical']['base_coefficients'].values());ck('whole critical no-direction certificate',sum((c*g for c,g in zip(ci['critical_bezout'],cg)),CP.zero())==1);ck('whole supercritical no-direction certificate',sum((c*g for c,g in zip(ci['super_bezout'],ci['super_coefficients'])),CP.zero())==1)
 print('SUPPORT PASS',label,'checks',len(checks),'seconds',time.monotonic()-st,flush=True)
record={'checks':len(checks),'seconds':time.monotonic()-st,'supports':3,'all_complex_amplitudes_retained':True,'all_materials_retained':True,'scope':'All three four-p full-RR supports in the critical and supercritical strata; source/normal criteria combined with original projective-p descent and finite incidence atlas.'}
save({'record':record,'checks':checks,'scope':__doc__},str(p/'independent_RR_torus_missing_p_reader.sobj'),compress=False)
(p/'independent_RR_torus_missing_p_reader.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2),flush=True)
