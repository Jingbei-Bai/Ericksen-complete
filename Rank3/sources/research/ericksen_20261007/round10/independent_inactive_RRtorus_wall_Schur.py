"""Range-eliminated true inactive normal rows, retaining both spectral-amplitude walls."""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_inactive_RRtorus_normal_block_reader.sobj'));Q=d['ring'];l,rho=Q.gens();K=Q.base_ring();a,r,H,Z=K.gens();m=1+a*r;F=Q.fraction_field();N=d['matrix'];rr=[0,3,4];cc=[2,3,4];oo=[0,1,5];A=N.matrix_from_rows_and_columns(rr,cc);dd=A.det();adj=A.adjugate();T=N.matrix_from_rows_and_columns(rr,oo);ids=[1,2,5,6];S=matrix(F,[dd*vector(F,[N[i,j] for j in oo])-vector(F,[N[i,j] for j in cc])*adj*T for i in ids]);print('formed',time.monotonic()-st,flush=True)
# Reduction in the exact quadratic algebra l^2=wall. Denominators remain powers of l,rho.
items=[]
for sign in [1,-1]:
 c=sign*a*r/m
 def redp(f):
  out=Q.zero()
  for ex,co in f.dict().items():out+=co*c**(int(ex[0])//2)*l**(int(ex[0])%2)*rho**int(ex[1])
  return out
 def red(f):
  f=F(f);num=redp(f.numerator());den=redp(f.denominator());conj=den.subs({l:-l});dd0=redp(den*conj);assert dd0!=0;return F(redp(num*conj))/dd0
 R=matrix(F,4,3,[red(v) for v in S.list()]);ei=next(i for i in [0,1] if R[i,1]);wcoef=-R[ei,0]/R[ei,1];ck=red(wcoef)
 pairs=[]
 for i in [2,3]:
  ps=red(R[i,0]+R[i,1]*ck);pt=red(R[i,2]);pairs.append(vector(F,[ps,pt]));print('WALL',sign,'H',i,'s_nonzero',bool(ps),'delta_nonzero',bool(pt),'s_terms',len(ps.numerator().dict()),'delta_terms',len(pt.numerator().dict()),flush=True)
 cross=red(pairs[0][0]*pairs[1][1]-pairs[0][1]*pairs[1][0]);print('CROSS',sign,bool(cross),'seconds',time.monotonic()-st,flush=True)
 items.append({'sign':sign,'l_square':c,'matrix':R,'even_row':ei,'w_over_s':ck,'H_pairs':pairs,'cross':cross})
save({'source_file':'independent_inactive_RRtorus_normal_block_reader.sobj','ring':Q,'range_rows':rr,'range_columns':cc,'range_determinant':dd,'range_adjugate':adj,'projected_rows':S,'selected_row_positions':ids,'items':items,'scope':__doc__},str(p/'independent_inactive_RRtorus_wall_Schur.sobj'),compress=False)
print('DONE',time.monotonic()-st,flush=True)
