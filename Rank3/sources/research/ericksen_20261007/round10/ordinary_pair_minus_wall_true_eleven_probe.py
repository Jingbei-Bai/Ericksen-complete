from sage.all import *
from pathlib import Path
import resource,os
p=Path(__file__).resolve().parent;resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+550*1024**2,)*2)
quad=load(str(p/'ordinary_pair_double_kernel_coupled.sobj'));third=load(str(p/'ordinary_pair_minus_wall_true_third_with_E29.sobj'));E=third['ring'];lb=E.gen();K=E.base_ring();a,r=K.gens();m=1+a*r;k=a*r;P=quad['ring'];l=P.gen();OK=P.base_ring();ph=OK.ring().hom([a,r,1,-1],K)
def cv(f):return sum((ph(co.numerator())/ph(co.denominator())*lb**int(ex) for ex,co in P(f).dict().items()),E(0))
us=load(str(p/'ordinary_pair_minus_whole_u_square.sobj'));UP=us['ring'];ul,uu=UP.gens();def_not_used=0
def cvu(f):return sum((ph(co.numerator())/ph(co.denominator())*lb**int(ex[0]) for ex,co in UP(f).dict().items() if not ex[1]),E(0))
M=third['M9'];ids=list(third['third']);rows=[]
for j in ids:
 f=quad['matrix'][j,8]-(1+OK.gen(0)*OK.gen(1))*quad['matrix'][j,1]-(1+OK.gen(0)*OK.gen(1))*l*quad['matrix'][j,2];qq,rem=f.quo_rem((1+OK.gen(0)*OK.gen(1))*l*l+OK.gen(0)*OK.gen(1));assert rem==0,('mu-div',j)
 urow=us['items'][j]['wall_remainder'].derivative(uu).subs({uu:0});row=list(M.row(j)[:8])+[cv(qq),cvu(urow),third['third'][j]];rows.append(row);print('SUP',j,[h for h,z in enumerate(row) if z],flush=True)
A=matrix(E,rows);print('sourceRows',ids,flush=True);print('rank11',A.rank(),flush=True);print('pivotrows',A.transpose().pivots(),flush=True);print('kernel',A.right_kernel().basis(),flush=True)
save({'ring':E,'matrix':A,'source_rows':ids,'columns':('range0_shift','n1_shift','n2_shift','tau','u_eta','w_k','eta_square','eta_k','mu_minus_k_square','u_k_square','k_cube'),'rank':A.rank(),'pivots':A.transpose().pivots()},str(p/'ordinary_pair_minus_wall_true_eleven_probe.sobj'),compress=False)
