from sage.all import *
from pathlib import Path
from itertools import combinations
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();s=load(str(p/'independent_critical_two_RR_R_scaled.sobj'));QP=s['ring'];T,W,S=QP.gens();lifts=[]
for it in s['items']:
 I=QP.ideal(it['QQ_sources']);items=[]
 for g in it['saturated_basis']:
  target=T**2*g;cf=target.lift(I);assert sum(c*f for c,f in zip(cf,it['QQ_sources']))==target;items.append({'target':target,'cofactors':cf});print('LIFT',it['case'],str(target),'terms',sum(len(c.dict()) for c in cf),flush=True)
 lifts.append({'case':it['case'],'sources':it['QQ_sources'],'items':items})
save({'items':lifts},str(p/'independent_critical_two_RR_direction_certificate.sobj'),compress=False)
outs=[]
for ci in range(3):
 d=load(str(p/('independent_subcritical_two_RR_full_normal_case'+str(ci)+'.sobj' if ci<2 else 'independent_subcritical_two_RR_full_normal.sobj')));P=d['ring'];K=P.base_ring();a,r,H,Z=K.gens();q=1+r;BB=[-q,-q/r,1/r][ci];U=PolynomialRing(K,names=('ct',));tt=U.gen();M=d['matrix'];rows=[]
 def red(f):
  dd={}
  for ex,c in f.dict().items():
   assert ex[3]==0
   pow=int(ex[0])+int(ex[2]);deg=pow%2;dd[deg]=dd.get(deg,K.zero())+c*(-QQ(1)/2)**(int(ex[1])+int(ex[2]))*BB**(pow//2)
  return U(dd)
 N=matrix(U,[[red(M[j,col]) for col in [0,1,3]] for j in range(25)]);found=[]
 for ij in combinations(range(7),2):
  for k in range(7,25):
   rr=list(ij)+[k];ff=N.matrix_from_rows(rr).det();f=ff% (tt**2-BB)
   if f:
    aa=K(f[0]);bb=K(f[1]);norm=aa**2-BB*bb**2
    if norm:
     print('CASE',ci,'ROWS',rr,'DET',str(f),'NORMFACTOR',str(norm.factor()),flush=True);found.append({'rows':rr,'determinant':ff,'remainder':f,'norm':norm});break
  if found:break
 outs.append({'case':ci,'square':BB,'ring':U,'matrix':N,'first_nonzero':found})
save({'items':outs},str(p/'independent_critical_two_RR_direction_H_probe.sobj'),compress=False);print('DONE',time.monotonic()-st,flush=True)
