"""Independent full ORIGINAL pair1 pullback at the resonant coupled corner.
Four root embeddings and the complete compensated second path are retained.
No creator pair Hessian or conjugation-reduction routine is used.
"""
from sage.all import *
from pathlib import Path
from itertools import combinations_with_replacement,permutations
import time,json,resource,os
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(label,z):assert bool(z),label;checks.append(label)
d=load(str(p/'independent_inactive_RRtorus_resonant_second_module.sobj'));raw=load(str(p/'independent_inactive_RRtorus_wall_variation_raw.sobj'));pairs=load(str(p/'independent_ratio_p_boundary_global_pair_syzygies.sobj'));pc=load(str(p/'independent_inactive_RRtorus_resonant_pair_second_module.sobj'));six=load(str(p/'independent_inactive_RRtorus_resonant_six_generator_certificate.sobj'))
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+700*1024**2,)*2)
K=PolynomialRing(QQ,'ir_r').fraction_field();r=K.gen();q=1+r;a=-q/r;m=-r;lc=-q/r;rc=-q
U=PolynomialRing(K,names=('ir_L','ir_R'));ul,ur=U.gens();E=U.quotient(U.ideal(ul**2-lc,ur**2-rc),names=('L','R'));L,R=E.gens();bas=(E(1),L,R,L*R)
def coord(f):return vector(K,[f.lift().monomial_coefficient(v) for v in (U(1),ul,ur,ul*ur)])
inv_cache={}
def inv(f):
 f=E(f)
 if f not in inv_cache:
  mat=matrix(K,4,4,lambda i,j:coord(f*bas[j])[i]);z=mat.solve_right(vector(K,[1,0,0,0]));g=sum((z[j]*bas[j] for j in range(4)),E(0));ck('root inverse multiplyback '+str(len(inv_cache)),f*g==1);inv_cache[f]=g
 return inv_cache[f]
old=d['ring'];of=old.base_ring();oph=of.ring().hom([r],K)
def oldpoly(f):return sum((oph(co.numerator())/oph(co.denominator())*L**int(e[0])*R**int(e[1]) for e,co in f.dict().items()),E(0))
def oe(f):f=old.fraction_field()(f);return oldpoly(f.numerator())*inv(oldpoly(f.denominator()))
rk=raw['ring'].base_ring();rh=rk.ring().hom([a,r,1,-1],K)
def rawpoly(f):return sum((rh(co.numerator())/rh(co.denominator())*L**int(e[0])*R**int(e[1]) for e,co in f.dict().items()),E(0))
def re(f):f=raw['ring'].fraction_field()(f);return rawpoly(f.numerator())*inv(rawpoly(f.denominator()))
ids=d['source_rows'];Nraw=raw['normal'];N=matrix(E,[[re(z) for z in row] for row in Nraw.rows()]);Nk=matrix(E,[[re(z) for z in row] for row in raw['kappa_derivative'].rows()]);Ns=matrix(E,[[oe(z) for z in row] for row in d['normal'].rows()]);ck('full normal independent root conversion',N==Ns)
ks=vector(E,[oe(z) for z in d['kernels'][0]]);kd=vector(E,[oe(z) for z in d['kernels'][1]]);ck('s kernel truly free',ks[0]==1 and ks[5]==0);ck('delta kernel truly free',kd[0]==0 and kd[5]==1);ck('all source normal s kernel',N*ks==0);ck('all source normal delta kernel',N*kd==0)
rr=[ids.index(j) for j in (0,5,7)];ec=[2,3,4]
def matinv(A):
 n=A.nrows();rows=[[E(A[i,j]) for j in range(n)]+[E(i==j) for j in range(n)] for i in range(n)]
 for j in range(n):
  t=next(i for i in range(j,n) if rows[i][j]);rows[j],rows[t]=rows[t],rows[j];z=inv(rows[j][j]);rows[j]=[f*z for f in rows[j]]
  for i in range(n):
   if i!=j:
    z=rows[i][j];rows[i]=[rows[i][k]-z*rows[j][k] for k in range(2*n)]
 return matrix(E,[row[n:] for row in rows])
Ae=N.matrix_from_rows_and_columns(rr,ec);Aei=matinv(Ae);ck('three e range inverse',Ae*Aei==identity_matrix(E,3));ww=[1,2];erow=next(ids.index(j) for j in ww if (N[ids.index(j),1]-N.matrix_from_columns(ec).row(ids.index(j))*Aei*N.matrix_from_rows_and_columns(rr,[1]).column(0))!=0)
rows4=rr+[erow];cols4=[1,2,3,4];A4=N.matrix_from_rows_and_columns(rows4,cols4);Ai4=matinv(A4);ck('four range inverse original normals',A4*Ai4==identity_matrix(E,4));print('FOUR TRUE RANGE SOURCE IDS',[ids[i] for i in rows4],flush=True)
qs,qd,qcross=d['quad_sources'];V=[]
for i,sid in enumerate(ids):
 dl=[re(Nraw[i,j].derivative(raw['ring'].gen(0))) for j in range(6)];dr=[re(Nraw[i,j].derivative(raw['ring'].gen(1))) for j in range(6)]
 V.append([sum((dl[j]*ks[j] for j in range(6)),E(0))*inv(2*L),sum((dr[j]*ks[j] for j in range(6)),E(0))*inv(2*R),Nk.row(i)*ks,sum((dl[j]*kd[j] for j in range(6)),E(0))*inv(2*L),sum((dr[j]*kd[j] for j in range(6)),E(0))*inv(2*R),Nk.row(i)*kd,oe(qs[sid]),oe(qcross[sid]),oe(qd[sid])])
V=matrix(E,V);C=-Ai4*V.matrix_from_rows(rows4);Vred=V+N.matrix_from_columns(cols4)*C;ck('all four second ranges compensated',Vred.matrix_from_rows(rows4)==0)
ck('second full old source module matches creator',Vred==matrix(E,[[oe(z) for z in row] for row in d['quotient_second'].rows()]))
# Polynomial quotient keeps all four square-root embeddings and all five
# independent normal directions through internal degree3. kappa^-1 shifts
# that degree by at most one; final complete second jet is extracted below.
A=PolynomialRing(K,names=('jr_L','jr_R','Wl','WR','Kap','S','D'));al,ar,Wl,WR,Kap,S,D=A.gens();nv=[Wl,WR,Kap,S,D];rels=[al**2-lc,ar**2-rc]
rels += [prod(nv[j] for j in ee) for ee in combinations_with_replacement(range(5),4)]
B=A.quotient(A.ideal(rels),names=('L','R','Wl','WR','Kap','S','D'));Lb,Rb,Wlb,WRb,Kb,Sb,Db=B.gens();eh=U.hom([Lb,Rb],B)
def eb(f):return eh(E(f).lift())
mon=[Sb*Wlb,Sb*WRb,Sb*Kb,Db*Wlb,Db*WRb,Db*Kb,Sb**2,Sb*Db,Db**2]
normal=[eb(ks[j])*Sb+eb(kd[j])*Db for j in range(6)]
for j,cj in enumerate(cols4):normal[cj]+=sum((eb(C[j,h])*mon[h] for h in range(9)),B(0))
s,w,e0,e1,e2,de=normal
ll=Lb*(1+Wlb/(2*lc)-Wlb**2/(8*lc**2)+Wlb**3/(16*lc**3));linv=Lb/lc*(1-Wlb/(2*lc)+3*Wlb**2/(8*lc**2)-5*Wlb**3/(16*lc**3));rhob=Rb*(1+WRb/(2*rc)-WRb**2/(8*rc**2)+WRb**3/(16*rc**3));rhoinv=Rb/rc*(1-WRb/(2*rc)+3*WRb**2/(8*rc**2)-5*WRb**3/(16*rc**3))
ck('l root path complete second order',ll**2==lc+Wlb);ck('rho root path complete second order',rhob**2==rc+WRb);ck('l inverse complete second order',ll*linv==1);ck('rho inverse complete second order',rhob*rhoinv==1)
imgs=[B(1),(s*rhob+Kb*e1)*(linv-Kb*e0*linv**2),(w+Kb*e2)*rhoinv,s,w,ll,ll+Kb*e0,rhob,B(1)]
orig=pairs['ring'].base_ring();pp=orig.ring().hom([a,a*q/m,1/r,-q/r],K);cc={}
def pcv(co):
 if co not in cc:cc[co]=pp(co.numerator())/pp(co.denominator())
 return cc[co]
it=next(x for x in pairs['items'] if x['pair']==1);val={}
for ex,co in it['full'].dict().items():
 ee=tuple(map(int,ex));pd=sum(ee[:6]);rd=sum(ee[6:]);kp=9-2*pd-3*rd;dp=4-pd-rd;ck('original pair full layer degree '+str(len(checks)),dp>=0 and kp>=-1)
 if kp>2:continue
 ff=B(pcv(co))*de**dp
 for j,nn in enumerate(ee):
  if nn:ff*=imgs[j]**nn
 for ex2,co2 in ff.lift().dict().items():
  z=list(map(int,ex2));z[4]+=kp
  if sum(z[2:])>2:continue
  key=tuple(z);val[key]=val.get(key,K(0))+co2
val={ex:co for ex,co in val.items() if co};ck('complete original kappa poles cancel',all(ex[4]>=0 for ex in val));ck('complete original value and all first terms zero',all(sum(ex[2:])==2 for ex in val))
jet=sum((co*prod(B.gen(j)**ex[j] for j in range(7)) for ex,co in val.items()),B(0));labels=['sWl','sWR','sK','dWl','dWR','dK','s2','sd','d2'];pv=[]
for h,label in enumerate(labels):
 exp=mon[h].lift().exponents()[0];coef=E(0)
 for ex,co in val.items():
  if tuple(ex[2:])==tuple(exp[2:]):coef+=co*L**ex[0]*R**ex[1]
 pv.append(coef)
pv=vector(E,pv);ck('complete original compensated pair quotient matches creator',pv==vector(E,[oe(z) for z in pc['items'][1]['quotient']]))
rem=vector(E,pv);redrows=[vector(E,[oe(z) for z in row]) for row in d['row_reduced']]
for j,row in zip(d['pivots'],redrows):rem-=rem[j]*row
ck('sixth coefficient only s2 survives',all(rem[j]==0 for j in range(9) if j!=6));ck('true reduced original pair coefficient',rem[6]==L/(2*q));print('ORIGINAL PAIR REDUCED',rem,flush=True)
sel=six['old_row_indices'];cols=six['old_columns'];oldM=Vred.matrix_from_rows_and_columns(sel,cols)
def det(M):
 n=M.nrows();return sum(((-1)**sum(per[i]>per[j] for i in range(n) for j in range(i+1,n))*prod(M[i,per[i]] for i in range(n)) for per in permutations(range(n))),E(0))
five=det(oldM);ck('old five selected source identities',[ids[i] for i in sel]==[1,17,11,19,13]);ck('old five coefficient minor matches literal creator',five==oe(six['old_minor']));combined=five*rem[6];ck('all four embeddings six determinant protected',combined==E(-5*q/(1152*r**5))*R);ck('six determinant matches creator',combined==oe(six['combined_six_minor']))
rec={'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st,'original_pair1_compensated_multivariate_path':True,'four_roots_retained':True,'no_pair_hessian_used':True,'scope':__doc__}
save({'record':rec,'checks':checks,'field':K,'root_ring':E,'range_source_rows':[ids[i] for i in rows4],'range_second':C,'old_second':Vred,'pair_second':pv,'pair_reduced':rem,'six_minor':combined},str(p/'independent_inactive_fullRR_subcritical_independent_pair_reader.sobj'),compress=False);(p/'independent_inactive_fullRR_subcritical_independent_pair_reader.json').write_text(json.dumps(rec,indent=2)+'\n');print('PASS',rec,flush=True)
