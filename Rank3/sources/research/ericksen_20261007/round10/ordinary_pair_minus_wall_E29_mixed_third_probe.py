from sage.all import *
from pathlib import Path
p=Path(__file__).resolve().parent
third=load(str(p/'ordinary_pair_minus_wall_true_third_probe.sobj'));E=third['ring'];l=E.gen();K=E.base_ring();a,r=K.gens();q=1+r;k=a*r;m=1+k
qq=load(str(p/'independent_ratio_incoming12_generic_normal_two.sobj'));O=qq['ring'];OK=O.base_ring();ph=OK.hom([a,a*q/m,1/r,-q/r],K);exps=[tuple(map(int,e)) for e in qq['normal_exponents']];row=list(qq['coefficient_matrix'].row(29));ol,ou,ov=O.gens()
def cv(f):return sum((ph(co)*l**int(ee[0])*(-2)**int(ee[2]) for ee,co in O(f).dict().items() if not ee[1]),E(0))
def exp(i,j):
 ee=[0]*6;ee[i]+=1;ee[j]+=1;return tuple(ee)
x1=cv(row[exps.index(exp(1,3))]-row[exps.index(exp(0,1))]);x2=cv(row[exps.index(exp(2,3))]-row[exps.index(exp(0,2))]);tt=-m*x1-m*l*x2;print('E29 shiftedthird',tt,flush=True)
third['third'][29]=tt;N=third['M9'][:,:4];third['compensated'][29]=tt+vector(E,list(N.row(29)))*third['range_third'];print('E29compensated',third['compensated'][29],flush=True);save(third,str(p/'ordinary_pair_minus_wall_true_third_with_E29.sobj'),compress=False)
