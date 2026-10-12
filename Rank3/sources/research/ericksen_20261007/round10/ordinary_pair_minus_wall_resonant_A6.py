"""All-material-resonant alternate-base true A6 coefficient determinant.
No state is frozen before whole J11 membership; only its boundary coefficient
matrix is evaluated to witness generic-base invertibility for each material.
"""
from sage.all import *
from pathlib import Path
import time,resource,os,json
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
quad=load(str(p/'ordinary_pair_double_kernel_coupled.sobj'))
pure=load(str(p/'ordinary_pair_plus_wall_cubic_reader.sobj'))
qq=load(str(p/'independent_ratio_incoming12_generic_normal_two.sobj'))
pair=load(str(p/'ordinary_pair_incoming_full_twojet.sobj'))
us=load(str(p/'ordinary_pair_minus_whole_u_square.sobj'))
B4=load(str(p/'ordinary_pair_minus_adapted_B.sobj'))
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+650*1024**2,)*2)
P=quad['ring'];K4=P.base_ring();aa,rr,HH,ZZ=K4.gens();ol=P.gen()
K=FunctionField(QQ,'r');r=K.gen();q=1+r;a=-q/r;k=-q;m=-r
U=PolynomialRing(K,'l');E=U.quotient(U.gen()**2+q/r,'l');l=E.gen()
ph=K4.ring().hom([a,r,K(1),K(-2)],K)
def cc(c):return ph(c.numerator())/ph(c.denominator())
def cv(f):return sum((cc(c)*l**int(e) for e,c in P(f).dict().items()),E.zero())
def multcv(f,baseph):
 return sum((baseph(c)*l**int(e[0])*(-2)**int(e[2]) for e,c in f.dict().items() if not e[1]),E.zero())
def ck(n,z):assert bool(z),n;checks.append(n)
OR=qq['ring'];OK=OR.base_ring()
oc=OK.hom([a,a*q/m,K(1)/r,-2*q/r],K)
exps=[tuple(map(int,e)) for e in qq['normal_exponents']]
def ex(i,j):
 z=[0]*6;z[i]+=1;z[j]+=1;return tuple(z)
PU=pure['ring'];pl,pu=PU.gens()
def purecv(f):
 return sum((cc(c)*l**int(e[0]) for e,c in PU(f).dict().items() if not e[1]),E.zero())
UP=us['ring'];ul,uu=UP.gens()
def ucv(f):
 return sum((cc(c)*l**int(e[0]) for e,c in UP(f).dict().items() if not e[1]),E.zero())
labs=quad['labels'];arows=[5,7,11,13,37,39];rows=[];third={};details={}
for i in arows:
 if i<37:
  two=list(qq['coefficient_matrix'].row(i))
  cross1=two[exps.index(ex(1,3))]-two[exps.index(ex(0,1))]
  cross2=two[exps.index(ex(2,3))]-two[exps.index(ex(0,2))]
  c1=multcv(cross1,oc);c2=multcv(cross2,oc)
  unshift=purecv(pure['full_pure_kernel_restrictions'][('old',i)][3])
 else:
  it=pair['items'][i-37];vals=it['jet_coefficients'];pes=[tuple(map(int,e)) for e in pair['normal_exponents']]
  def paircv(f):
   return sum((cc(c)*l**int(e[0])*(-2)**int(e[2]) for e,c in f.dict().items() if not e[1]),E.zero())
  c1=paircv(vals[pes.index(ex(1,3))]-vals[pes.index(ex(0,1))])
  c2=paircv(vals[pes.index(ex(2,3))]-vals[pes.index(ex(0,2))])
  if i==39:unshift=purecv(pure['full_pure_kernel_restrictions'][('pair',2)][3])
  else:
   unshift=E.zero()
   for ee,c in it['full_homogeneous_pullback'].dict().items():
    e=tuple(map(int,ee))
    if e[4] or e[5] or e[7] or e[8] or e[3]+e[6]!=3 or e[1]:continue
    unshift+=cc(c)*l**e[0]*(-2)**e[2]*(-1)**e[3]
 tt=unshift-m*c1-m*l*c2;third[i]=tt
 oldrow=quad['matrix'].row(i);gf=oldrow[8]-(1+aa*rr)*oldrow[1]-(1+aa*rr)*ol*oldrow[2]
 mu=(1+aa*rr)*ol*ol+aa*rr;muco,rem=gf.quo_rem(mu);ck('generic mu quotient '+str(i),not rem)
 uc=ucv(us['items'][i]['wall_remainder'].derivative(uu).subs({uu:0}))
 full=list(map(cv,oldrow[:8]))+[cv(muco),uc,tt]
 ck('A rows no B-block columns '+str(i),all(not full[j] for j in (1,2,6,8)))
 ck('A rows no u eta '+str(i),not full[4])
 rows.append([full[j] for j in (0,3,5,7,9,10)])
 details[i]={'unshifted':unshift,'cross1':c1,'cross2':c2,'third':tt,'full11':full}
A=matrix(E,rows)
# Determinant through exact square-root parity, no numerical material values.
par={}
for i in range(6):
 for j in range(6):
  nz=[h for h,c in enumerate(A[i,j].list()) if c]
  if nz:ck('pure parity '+str((i,j)),len(nz)==1);par[i,j]=nz[0]
rp=[None]*6;cp=[None]*6;rp[0]=0
for step in range(12):
 for (i,j),pp in par.items():
  if rp[i] is not None:cp[j]=(pp-rp[i])%2
  if cp[j] is not None:rp[i]=(pp-cp[j])%2
ck('all parity assigned',all(z is not None for z in rp+cp))
def scalar(f):
 ff=E(f).lift();ck('scalar coefficient',ff.degree()<=0);return K(ff[0])
M=matrix(K,6,6,lambda i,j:scalar(A[i,j]/l**(rp[i]+cp[j])))
dens=[lcm([z.denominator() for z in row]) for row in M.rows()]
CM=matrix(K.maximal_order(),[[dens[i]*M[i,j] for j in range(6)] for i in range(6)])
dd=K(CM.det(algorithm='df'))/prod(dens);phase=sum(rp)+sum(cp);norm=dd**2*(-1)**phase*(-q/r)**phase
print('ALT BASEDET',dd.factor(),flush=True);print('ALT NORM',norm.factor(),flush=True)
bd=cv(B4['final_B_det']);print('ALT B4',bd,flush=True)
up=cv(quad['matrix'][0,4]);print('ALT UETA',up,flush=True)
ck('alternate A6 nonzero generic resonant material',dd!=0)
ck('alternate B4 nonzero generic resonant material',bd!=0)
ck('alternate u eta nonzero',up!=0)
save({'ring':E,'A_matrix':A,'details':details,'A_base_det':dd,'A_norm':norm,'A_phase':phase,'B4_det':bd,'u_eta_unit':up,'rows':arows,'columns':(0,3,5,7,9,10),'checks':checks},str(p/'ordinary_pair_minus_wall_resonant_A6.sobj'),compress=False)
(p/'ordinary_pair_minus_wall_resonant_A6.json').write_text(json.dumps({'status':'PASS ALTERNATE MATERIAL-RESONANT A6','checks':len(checks),'seconds':time.monotonic()-st},indent=2))
print('DONE',len(checks),time.monotonic()-st,flush=True)
