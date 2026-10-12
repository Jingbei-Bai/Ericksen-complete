"""Independent cyclic v=-2 u-nonzero source proof, retaining all materials."""
from sage.all import *
from pathlib import Path
import time,resource,os,json
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(n,b):
 assert bool(b),n
 checks.append(n)
normal=load(str(p/'independent_ratio_incoming12_generic_normal.sobj'));pair=load(str(p/'ordinary_pair_incoming_full_twojet.sobj'));g=load(str(p/'independent_ratio_p_boundary_second_gauge.sobj'));hs=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'));qs=load(str(p.parent/'round8/rank3_joint_fiber_leading.sobj'));saved=load(str(p/'ordinary_cyclic_positive_wall_rank_probe.sobj'))
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+550*1024**2,)*2)
F=PolynomialRing(QQ,names=('a','r')).fraction_field();a,r=F.gens();k=a*r;m=1+k;q=1+r;DD=k*r+m*q
P=PolynomialRing(F,names=('l','u'));l,u=P.gens();OK=g['coefficient_field'];ph=OK.ring().hom([a,a*q/m,1/r,-q/r],F);cache={}
def cf(c):
 c=OK(c)
 if c not in cache:cache[c]=ph(c.numerator())/ph(c.denominator())
 return cache[c]
def raw(f):
 ans={}
 for ee,c in (f.items() if isinstance(f,dict) else f.dict().items()):
  e=tuple(map(int,ee));kd=e[1]+e[4];dg=sum(e)
  if e[2] or e[3] or kd>3 or dg>5:continue
  ans[(dg,kd)]=ans.get((dg,kd),P.zero())+cf(c)*(-1)**e[1]*(-2)**(e[6]+e[8])*l**(e[5]+e[6])*u**e[7]
 return ans
memo={}
def pool(i):
 if i not in memo:memo[i]=raw(hs['polynomial_pool'][i])
 return memo[i]
labels=[('old',i) for i in (0,1,2,5,7,11,13,17)]+[('pair',2)]
pure={}
PK=pair['coefficient_field'];ka,kr,kH,kZ=PK.gens();pm=PK.ring().hom([a,r,1,-1],F)
def pcf(c):return pm(c.numerator())/pm(c.denominator())
for lab in labels:
 typ,i=lab
 if typ=='pair':
  vals=[P.zero() for j in range(4)]
  for ee,c in pair['items'][i]['full_homogeneous_pullback'].dict().items():
   e=tuple(map(int,ee));kd=e[3]+e[6]
   if e[4] or e[5] or e[7] or e[8] or kd>3:continue
   vals[kd]+=pcf(c)*(-1)**e[3]*(-2)**e[2]*l**e[0]*u**e[1]
 else:
  if i<4:
   ff=qs['quartic_source_rows'][(0,1,5,6)[i]];D=int(ff.total_degree());dd=raw(ff);vals=[dd.get((D,j),P.zero()) for j in range(4)]
  elif i<10:
   dd=raw(g['six_full_sources'][i-4]);vals=[dd.get((4,j),P.zero()) for j in range(4)]
  else:
   vals=[P.zero() for j in range(4)]
   for sign,ii,jj in hs['records'][i-10]['full_signed_product_DAG']:
    aa=pool(ii);bb=pool(jj)
    for (da,ea),ca in aa.items():
     for z in range(ea,4):vals[z]+=sign*ca*bb.get((5-da,z-ea),P.zero())
 pure[lab]=vals
 ck('full raw zero boundary '+str(lab),vals[0]==0)
 ck('full raw zero kernel linear '+str(lab),vals[1]==0)
 ck('full all-u quadratic '+str(lab),vals[2]==(8*(m*l*l-k) if lab==('old',2) else 0))
 print('FULL_SOURCE',lab,'seconds',time.monotonic()-st,flush=True)
OR=normal['ring'];ol,ou,ov=OR.gens()
def np(f):
 return sum((cf(c)*l**int(e[0])*u**int(e[1])*(-2)**int(e[2]) for e,c in OR(f).dict().items()),P.zero())
M=[]
for lab in labels:
 typ,i=lab
 if typ=='old':
  row=list(normal['matrix'].row(i));ck('whole kernel v-plus2 '+str(lab),not (row[3]-row[0]).subs({ov:-2}))
  rr=[np(row[j]) for j in (0,1,2,4,5)];wk=np((row[3]-row[0]).derivative(ov))
 else:
  ff=pair['items'][i]['full_homogeneous_pullback'];UR=pair['jet_ring'];pl,pu,pv=UR.gens();row=[]
  for j in range(6):
   ee=tuple(1 if h==j else 0 for h in range(6));row.append(UR({tuple(map(int,e[:3])):c for e,c in ff.dict().items() if tuple(map(int,e[3:]))==ee}))
  def pp(f):return sum((pcf(c)*l**int(e[0])*u**int(e[1])*(-2)**int(e[2]) for e,c in f.dict().items()),P.zero())
  ck('pair whole kernel v-plus2',not (row[3]-row[0]).subs({pv:-2}));rr=[pp(row[j]) for j in (0,1,2,4,5)];wk=pp((row[3]-row[0]).derivative(pv))
 M.append(rr+[wk,P(8) if lab==('old',2) else P.zero(),pure[lab][3]])
M=matrix(P,M)
SP=saved['ring'];SQ=SP.base_ring();SF=SQ.base_ring();sf=F.hom(SF.gens(),SF);sl=SQ.gen();su=SP.gen()
def wall(f):return sum((sf(c)*sl**int(e[0])*su**int(e[1]) for e,c in f.dict().items()),SP.zero())
W=M.apply_map(wall)
for i,lab in enumerate(labels):
 for j in range(8):ck('original source vs cyclic saved '+str((lab,j)),W[i,j]==saved['matrix'][saved['labels'].index(lab),j])
def determinant(A):
 R=A.base_ring()
 if not A.nrows():return R.one()
 i=min(range(A.nrows()),key=lambda i:sum(bool(x) for x in A.row(i)))
 return sum(((-1)**(i+j)*A[i,j]*determinant(A.delete_rows([i]).delete_columns([j])) for j in range(A.ncols()) if A[i,j]),R.zero())
sets=[(1,5,7,11,17),(1,7,11,13,17)];cert=[]
for idx,rs in enumerate(sets):
 labs=[('old',0),('old',2)]+[('old',j) for j in rs]+[('pair',2)];inds=[labels.index(x) for x in labs];eight=determinant(W.matrix_from_rows(inds));five=determinant(W.matrix_from_rows_and_columns(inds[2:7],(0,1,2,4,5)))
 exp=8*a**3*r**2/(27*m**6)*(m*u*u-40*DD) if idx==0 else a**3*r**2/(54*m**6)*(19*m*u*u-96*DD)
 ck('new true five factor '+str(idx),five==wall(P(exp)))
 ratio=wall(P(8*(3*u/56)*(4*k)))
 ck('full eight Laplace with all range columns '+str(idx),eight==ratio*five or eight==-ratio*five)
 cert.append({'rows':labs,'matrix':W.matrix_from_rows(inds),'eight':eight,'five':five,'five_factor':exp,'laplace_sign':1 if eight==ratio*five else -1})
ck('two five primitive Bezout',19*(m*u*u-40*DD)-(19*m*u*u-96*DD)==-664*DD)
# Independent reconstruction of the old true-minor H coefficients at Z=-1.
d=load(str(p/'independent_ordinary_minus_two_coupled_minor.sobj'));OP=d['ring'];oa,orr,oH,oZ,ol,ou=OP.gens();LV=PolynomialRing(F,names=('L','V'));LL,VV=LV.gens();cfs=[];source_factors=[]
for it in d['items']:
 ck('exact stored factorization product '+str(it['rows']),it['factorization'].value()==it['primitive_determinant'])
 fs=[f for f,e in it['factorization'] if (f.degree(ol)>1 or f.degree(ou)>0) and len(f.dict())>1];ck('one state factor '+str(it['rows']),len(fs)==1);ff=fs[0];byH={}
 for ee,c in ff.dict().items():
  e=tuple(map(int,ee));ck('even exponents '+str((it['rows'],e)),e[4]%2==0 and e[5]%2==0)
  h=e[2]+e[5]//2;byH[h]=byH.get(h,LV.zero())+c*a**e[0]*r**e[1]*(-1)**e[3]*LL**(e[4]//2)*VV**(e[5]//2)
 ck('all H degrees 0..4 '+str(it['rows']),sorted(byH)==list(range(5)))
 source_factors.append(byH);cfs+=list(byH.values())
 muP=m*LL-k;muM=m*LL+k
 ck('outer source factors old units on Z=-1 '+str(it['rows']),all(f in (ol,oH,-oH+oZ,orr,oH+1,oZ-1,oH*oZ-oH+oZ,oa*orr+1,oa) for f,e in it['factorization'] if f!=ff))
ck('off-plus highest H coefficient',source_factors[0][4]==-384*m*VV*LL*muP*muM**2)
ck('negative-wall next coefficient old unit',source_factors[0][1].subs({LL:-k/m})==-512*k**4*VV)
# D positivity is the same exact expression in every positive spectral ordering.
AB=PolynomialRing(QQ,names=('A','B')).fraction_field();aa,bb=AB.gens();rr=(aa-bb)/(aa*(bb-1));mp=F.hom([aa,rr],AB)
ck('D positive original spectral expression',mp(DD)==((aa-bb)**2+bb*(aa-1)**2)/(aa*(bb-1)**2))
# The old u=0 positive-wall determinant is already independently source-audited59.
# Its three remaining base factors at the cyclic point are exact constants below.
HZ=PolynomialRing(F,names=('H','Z'));HH,ZZ=HZ.gens();bf=ZZ*(3*HH+2)+3*HH+1;af=-11*HH*ZZ**2+ZZ**3-12*HH*ZZ-4*ZZ**2-HH+2*ZZ+1
ck('u0 positive B cyclic remainder',bf(HH=1,ZZ=-1) if False else bf.subs({HH:1,ZZ:-1})==-1)
ck('u0 positive A cyclic remainder',af.subs({HH:1,ZZ:-1})==-6)
save({'coefficient_field':F,'ring':P,'labels':labels,'original_matrix':M,'cyclic_wall_matrix':W,'full_pure_coefficients':pure,'new_eight_certificates':cert,'old_H_coefficients':source_factors,'positive_D':DD,'checks':checks,'scope':__doc__},str(p/'ordinary_cyclic_incoming_u_nonzero_reader.sobj'),compress=False)
report={'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st,'scope':'Cyclic incoming v=-2 all u nonzero; exact full sources and all materials. Also positive-wall u0 cyclic factors. E0,u0 remains separate.'};(p/'ordinary_cyclic_incoming_u_nonzero_reader.json').write_text(json.dumps(report,indent=2)+chr(10));print(report,flush=True)
