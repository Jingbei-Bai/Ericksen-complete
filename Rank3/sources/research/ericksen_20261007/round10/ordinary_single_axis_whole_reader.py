"""Independent exact minor/Sylvester audit for every ordinary p12-axis point.
The radical Schur lemma, rather than full J9 rank, closes u=0.
"""
from sage.all import *
from pathlib import Path
import time,resource,os,json
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
ax=load(str(p/'ordinary_single_axis_J9.sobj'));ref=load(str(p/'ordinary_axis_v_minus_relative_factor.sobj'));nm=load(str(p/'ordinary_axis_nilpotent_schur_two_minors.sobj'));bb=load(str(p/'ordinary_single_axis_J9_blocks.sobj'));dm=load(str(p/'ordinary_single_axis_dimensionless_blocks.sobj'));bw=load(str(p/'ordinary_single_axis_Bwall.sobj'));kt=load(str(p/'independent_ratio_axis_normal_kernel_certificate.sobj'));nr=load(str(p/'independent_ratio_incoming12_generic_normal.sobj'))
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+500*1024**2,)*2)
U=ax['ring'];K=U.base_ring();uu,vv=U.gens();R=PolynomialRing(K,'v');v=R.gen();OK=nr['ring'];ll,ou,ov=OK.gens()
def ck(n,x):assert bool(x),n;checks.append(n)
def toO(f):return OK({(0,int(e[0]),int(e[1])):c for e,c in U(f).dict().items()})
beta=toO(ax['beta'])
ck('beta -2 is -1',beta.subs({ov:-2})==-1)
for i in range(37):
 f=beta*nr['matrix'][i,0]+nr['matrix'][i,3]
 ck('original whole linear divide '+str(i),f==ll*(ov+2)*ref['whole_divided_coefficients'][i])
for i in range(40):
 for j in range(9):
  ck('refined matrix '+str((i,j)),ax['matrix'][i,j]==ref['matrix'][i,j]*((vv+2) if j==5 else 1))
def u0(f):return R({int(e[1]):c for e,c in U(f).dict().items() if not e[0]})
M=ref['matrix'].apply_map(u0)
for i in (7,11,25,37):ck('pure A block '+str(i),all(not M[i,j] for j in (0,1,2,4,6,8)))
for i in (1,2,3,34,35,36):ck('pure quadratic block '+str(i),all(not M[i,j] for j in (0,2,3,4,5,7)))
BK=PolynomialRing(QQ,names=('a','r'));KK=BK.fraction_field();a,r=KK.gens();q=1+r;k=a*r;m=1+k;cy=K.hom([a,a*q/m,1/r,-q/r],KK);S=PolynomialRing(KK,'v');sv=S.gen()
def cyc(f):return S({int(e):cy(c) for e,c in R(f).dict().items()})
def dets(A):
 if not A.nrows():return A.base_ring().one()
 row=min(range(A.nrows()),key=lambda i:sum(bool(c) for c in A.row(i)))
 return sum(((-1)**(row+j)*A[row,j]*dets(A.delete_rows([row]).delete_columns([j])) for j in range(A.ncols()) if A[row,j]),A.base_ring().zero())
def syl(f,g):
 n=int(f.degree());d=int(g.degree());Kf=f.base_ring();cols=[]
 for z,jmax in ((f,d),(g,n)):
  for j in range(jmax):cols.append([z[t-j] if t>=j else Kf.zero() for t in range(n+d)])
 return matrix(Kf,cols).transpose()
# The generic Sylvester matrices are retained without expanding large determinants.
cert=[]
def certify(name,f,g,ef,eg):
 ck(name+' first degree preserved',f.degree()==ef.degree());ck(name+' second degree preserved',g.degree()==eg.degree())
 ck(name+' first literal specialization',cyc(f)==ef);ck(name+' second literal specialization',cyc(g)==eg)
 A=syl(f,g);Ac=A.apply_map(cy);rr=dets(Ac)
 ck(name+' direct cyclic determinant nonzero',rr!=0)
 cert.append({'name':name,'polynomials':(f,g),'sylvester':A,'cyclic_matrix':Ac,'cyclic_determinant':rr})
 print('SYLVESTER',name,'size',A.nrows(),'factor',rr.factor(),flush=True)
 return rr
# All-v tau/lwk two-column rank, no eta*k column.
f=M[7,3]*M[25,5]-M[7,5]*M[25,3];g=M[7,3]*M[37,5]-M[7,5]*M[37,3]
f0=q/(9216*r)*(sv-6)*(sv**3-200*sv**2-1684*sv-2240);g0=q**2/(384*r**2)*(sv+10)*(3*sv**2+2*sv-152)
rr=certify('u0 tau-lwk',f,g,f0,g0)
ck('known scalar determinant',rr==-QQ(4548425)/31701690482688*(q/r)**11)
# All-v rank of the three quadratic columns, from two true three-minors.
f=dets(M.matrix_from_rows_and_columns((1,2,3),(1,6,8)));g=dets(M.matrix_from_rows_and_columns((1,2,35),(1,6,8)))
certify('u0 quadratic',f,g,-k/(8*m**2)*(sv+2)**2,k/(32*m)*(31*sv**2+220*sv+700))
for i,j,want in [(13,0,-5*k/(6*m)),(19,2,-5/(3*m)),(0,4,KK(3)/56)]:
 ck('single column base only '+str((i,j)),M[i,j].degree()==0)
 ck('single column unit witness '+str((i,j)),cyc(M[i,j])==want)
# Independently check the specific full nine-by-nine endpoint determinant.
rows=(0,1,2,7,13,19,25,35,37);D=matrix(KK,[[cy(K(M[i,j](-2))) for j in range(9)] for i in rows]);dd=dets(D)
ck('v=-2 nine source determinant',dd==-QQ(325)/252*q**2*a**2/m**3)
# u!=0: exact seven-generator matrix and source singleton on lwk.
P=dm['ring'];F=P.base_ring();aa,rr,HH,ZZ=F.gens();qq=1+rr;kk=aa*rr;mm=1+kk;pu,pv=P.gens();B=matrix(P,dm['matrix7'])
for i in range(40):B[i,5]=P(B[i,5]/(pv+2))
Bv=HH*(pv+2)-4*ZZ*(HH+1)
ck('pair S0 only lwk',all(not B[37,j] for j in range(7) if j!=5))
ck('pair S0 exact unit times B',B[37,5]==-qq*HH/(2*rr*ZZ*(HH+1)**2)*Bv)
for i in (1,2,3,34,35,36):ck('u nonzero even pure block '+str(i),all(not B[i,j] for j in range(7) if j not in (1,6)))
f=B[1,1]*B[3,6]-B[3,1]*B[1,6];g=B[1,1]*B[2,6]-B[2,1]*B[1,6];cf=-4*kk*HH**2/(mm**2*ZZ**2*(HH+1)**2);cg=2*kk/(mm*ZZ**2*(HH+1)**2)
w=pv+2;Q2=3*HH**2*w**2-4*HH*ZZ*(HH+1)*w+8*ZZ**2*(HH+1)**2
ck('first even 2minor',f==cf*w**2);ck('second even 2minor',g==cg*Q2)
R0=8*ZZ**2*(HH+1)**2;d0=4*HH*ZZ*(HH+1);c0=3*HH**2;coQ=1/R0+d0/R0**2*w;coW=-(c0/R0-d0**2/R0**2+c0*d0/R0**2*w)
ck('uniform fullbase quadratic Bezout',coQ*Q2+coW*w**2==1)
# On B=0, true block minors in (n23,tau,lwk).
vb=4*ZZ*(HH+1)/HH-2;V=PolynomialRing(F,'u');vu=V.gen();hp=P.hom([vu,V(vb)],V);C=B.apply_map(hp)
f=dets(C.matrix_from_rows_and_columns((7,19,28),(2,3,5)));g=dets(C.matrix_from_rows_and_columns((7,25,28),(2,3,5)))
cy2=F.hom([a,r,KK.one(),-KK.one()],KK)
def cyu(z):return S({int(e):cy2(c) for e,c in V(z).dict().items()})
ef=-sv*(59*sv**2-4864*q)/(576*r*m);eg=-(53*sv**2-1600*q)/(72*r*m)
ck('Bwall first source cyclic',cyu(f)==ef);ck('Bwall second source cyclic',cyu(g)==eg)
ck('Bwall two degree bounds',f.degree()==3 and g.degree()==2)
BA=syl(f,g);Bcy=BA.apply_map(cy2);Bdet=dets(Bcy);ck('Bwall resultant nonzero',Bdet!=0)
print('SYLVESTER Bwall',Bdet.factor(),flush=True)
cert.append({'name':'u nonzero Bwall','polynomials':(f,g),'sylvester':BA,'cyclic_matrix':Bcy,'cyclic_determinant':Bdet})
# The pre-existing exact all-u/v first-normal rank certificate supplies
# rank two for (n23,tau); audit its ALL-material cyclic factors explicitly.
normal_cert=[]
for it in kt['base_nonzero_certificates']:
 val=it['cyclic_value'];par=val.parent();ka,kr=par.gens();ph=par.hom([a,r],KK);z=ph(val);ck('normal base unit nonzero '+it['name'],z!=0);normal_cert.append((it['name'],z));print('NORMAL UNIT',it['name'],z.factor(),flush=True)
for it in kt['resultant_sources']:
 val=it['cyclic_determinant'];par=val.parent();ph=par.hom([a,r],KK);z=ph(val);ck('normal resultant nonzero '+it['chart'],z!=0);normal_cert.append((it['chart'],z));print('NORMAL RESULTANT',it['chart'],z.factor(),flush=True)
save({'checks':checks,'source_files':('ordinary_single_axis_J9.sobj','ordinary_axis_v_minus_relative_factor.sobj','ordinary_single_axis_dimensionless_blocks.sobj','independent_ratio_axis_normal_kernel_certificate.sobj'),'sylvester_certificates':cert,'endpoint_true9matrix':D,'endpoint_true9det':dd,'seven_Bwall_matrix':C,'quadratic_bezout':(coQ,coW,Q2),'normal_rank_certificates':normal_cert,'scope':'Whole relative source ideal plus all finite u/v rank covers; nilpotent Schur required on u=0.'},str(p/'ordinary_single_axis_whole_reader.sobj'),compress=False)
rep={'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st};(p/'ordinary_single_axis_whole_reader.json').write_text(json.dumps(rep,indent=2)+'\n');print(rep,flush=True)
