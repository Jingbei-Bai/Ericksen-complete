from sage.all import *
from pathlib import Path
import time,resource,os
p=Path(__file__).resolve().parent;st=time.monotonic()
d=load(str(p/'ordinary_pair_minus_wall_cyclic_coupled_probe.sobj'));E=d['ring'];lb=E.gen();K=E.base_ring();a,r=K.gens();q=1+r;k=a*r;m=1+k;R=PolynomialRing(K,names=('L','T','D'));xl,xt,xd=R.gens();A=R.quotient(R.ideal(xl**2+k/m,xt**4,xd**6),names=('L','T','D'));l,t,dm=A.gens();fromE=E.hom([l],A)
src=load(str(p/'independent_ratio_p_boundary_second_gauge.sobj'));hs=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'));qs=load(str(p.parent/'round8/rank3_joint_fiber_leading.sobj'));pair=load(str(p/'ordinary_pair_incoming_full_twojet.sobj'));resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+650*1024**2,)*2);OK=src['coefficient_field'];ph=OK.hom([a,a*q/m,1/r,-q/r],K);memo={};cvm={}
def cv(c):
 c=OK(c)
 if c not in cvm:cvm[c]=ph(c)
 return cvm[c]
imgs=[dm,-dm*t,-dm*m*t*t,-dm*m*l*t*t,dm*t,dm*l,-2*dm*l,A(0),-2*dm]
def pull(f):
 out=A(0)
 for ex,c in (f.items() if isinstance(f,dict) else f.dict().items()):
  ee=tuple(map(int,ex))
  if ee[7] or sum(ee)>5 or ee[1]+2*ee[2]+2*ee[3]+ee[4]>3:continue
  v=A(cv(c))
  for j,n in enumerate(ee):
   if n:v*=imgs[j]**n
  out+=v
 return out
def pool(i):
 if i not in memo:memo[i]=pull(hs['polynomial_pool'][i])
 return memo[i]
def coef(f,D,n):return sum((K(c)*lb**int(ex[0]) for ex,c in A(f).lift().dict().items() if int(ex[1])==n and int(ex[2])==D),E(0))
labels=[0,2,1,5,7,11,13,17];third={};raws={}
for ii in labels:
 if ii<4:
  ff=qs['quartic_source_rows'][(0,1,5,6)[ii]];deg=int(ff.total_degree());f=pull(ff)
 elif ii<10:deg=4;f=pull(src['six_full_sources'][ii-4])
 else:deg=5;f=sum((sgn*pool(j)*pool(h) for sgn,j,h in hs['records'][ii-10]['full_signed_product_DAG']),A(0))
 co=[coef(f,deg,n) for n in range(4)];assert all(z==0 for z in co[:3]),(ii,co)
 third[ii]=co[3];raws[ii]=f;print('ROW',ii,'THIRD',co[3],flush=True)
# Full ordinary homogeneous global-pair sources include their true range coordinates.
PF=pair['coefficient_field'];pp=PF.ring().hom([a,r,1,-1],K)
for it in pair['items']:
 f=it['full_homogeneous_pullback'];out=E(0);vals=[l,A(0),A(-2),-t,-m*t*t,-m*l*t*t,t,A(0),A(0)];vv=A(0)
 for ex,c in f.dict().items():
  z=A(pp(c.numerator())/pp(c.denominator()))
  for j,n in enumerate(ex):
   if n:z*=vals[j]**int(n)
  vv+=z
 co=[sum((K(cc)*lb**int(ex[0]) for ex,cc in vv.lift().dict().items() if int(ex[1])==n),E(0)) for n in range(4)]
 assert all(z==0 for z in co[:3]),('pair',it['pair'],co)
 third[37+it['pair']]=co[3];print('PAIR',it['pair'],'RAW_THIRD',co[3],flush=True)
M=d['matrix'];N=M[:,:4];rr=(0,2,5,17);Nc=N.matrix_from_rows(rr);rhs=vector(E,[third[j] for j in rr]);sol=Nc.solve_right(-rhs);print('TRUE_RANGE_THIRD',list(sol),flush=True)
res={j:co+vector(E,list(N.row(j)))*sol for j,co in third.items()};print('TRUE_COMPENSATED',res,flush=True)
save({'ring':E,'M9':M,'third':third,'range_third':sol,'compensated':res,'rows':rr},str(p/'ordinary_pair_minus_wall_true_third_probe.sobj'),compress=False);print('DONE',time.monotonic()-st,flush=True)
