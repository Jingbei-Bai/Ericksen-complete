from sage.all import *
from pathlib import Path
import time,resource,os
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'ordinary_pair_minus_wall_true_eleven_probe.sobj'));E=d['ring'];l=E.gen();K=E.base_ring();a,r=K.gens();resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+650*1024**2,)*2)
ids=d['source_rows'];rr=[ids.index(j) for j in [5,7,11,13,37,39]];cols=[0,3,5,7,9,10];A=d['matrix'].matrix_from_rows_and_columns(rr,cols);par={}
for i in range(6):
 for j in range(6):
  f=A[i,j];vs=f.list();nz=[h for h,c in enumerate(vs) if c]
  if f:assert len(nz)==1,('mixed parity',i,j,vs);par[i,j]=nz[0]
rp=[None]*6;cp=[None]*6;rp[0]=0
for step in range(12):
 for (i,j),pp in par.items():
  if rp[i] is not None:cp[j]=(pp-rp[i])%2
  if cp[j] is not None:rp[i]=(pp-cp[j])%2
assert all(z is not None for z in rp+cp)
def scalar(f):
 ff=E(f).lift();assert ff.degree()<=0;return K(ff[0])
B=matrix(K,6,6,lambda i,j:scalar(A[i,j]/l**(rp[i]+cp[j])));print('parity',rp,cp,flush=True)
# Clear row denominators before Bareiss; only a,r remain polynomial variables.
R=K.ring();dens=[lcm([z.denominator() for z in row]) for row in B.rows()];C=matrix(R,[[R(dens[i]*B[i,j]) for j in range(6)] for i in range(6)]);f=C.det(algorithm='df');dd=K(f)/prod(dens);print('BASEDET',dd.factor(),flush=True);phase=sum(rp)+sum(cp);norm=dd*dd*(-1)**phase*(-a*r/(1+a*r))**phase;print('NORM',norm.factor(),flush=True)
save({'root_ring':E,'matrix':A,'unphased':B,'phase_power':phase,'unit_phase_rows':rp,'unit_phase_cols':cp,'base_det':dd,'norm':norm,'source_rows':[5,7,11,13,37,39],'columns':cols},str(p/'ordinary_pair_minus_wall_A6_determinant.sobj'),compress=False);print('DONE',time.monotonic()-st,flush=True)
