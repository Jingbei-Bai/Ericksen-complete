from sage.all import *
from pathlib import Path
from itertools import combinations
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();items=[];checks=[]
def ck(n,v):
 if not bool(v):raise AssertionError(n)
 checks.append(n)
B=PolynomialRing(QQ,names=('cr_r','cr_t','cr_w','cr_s'));rr,tt,ww,ss=B.gens();BF=B.fraction_field()
for ci in range(3):
 d=load(str(p/('independent_subcritical_two_RR_full_normal_case'+str(ci)+'.sobj' if ci<2 else 'independent_subcritical_two_RR_full_normal.sobj')));P=d['ring'];K=P.base_ring();a,r,H,Z=K.gens();t,w,s,e=P.gens();M=d['matrix'];coefs=[];records=[]
 for rows in combinations(range(7),3):
  f=M.matrix_from_rows_and_columns(rows,[0,1,3]).det()
  if not f:continue
  # Full denominator clearing, then retain every H/Z coefficient. No state-factor stripping.
  den=lcm([c.denominator() for c in f.coefficients()]);by={}
  for ex,c in f.dict().items():
   ck('no e in critical matrix minor',ex[3]==0);cc=K(c)*den;ck('exact polynomial coefficient',cc.denominator()==1)
   for bex,b in cc.numerator().dict().items():
    ck('no a in original critical R/Q minors',bex[0]==0);key=tuple(map(int,bex[2:]));by[key]=by.get(key,B.zero())+QQ(b)*rr**int(bex[1])*tt**int(ex[0])*ww**int(ex[1])*ss**int(ex[2])
  # remove only an integer content; other material and direction factors retained
  rowco=[]
  for z in by.values():
   if z:rowco.append(z/z.content())
  coefs+=rowco;records.append({'rows':rows,'minor':f,'denominator':den,'base_coefficients':by,'normalized_sources':rowco})
 unique=list(dict.fromkeys(B(f) for f in coefs));print('CASE',ci,'MINORS',len(records),'COEFS',len(unique),'TERMS',[len(z.dict()) for z in unique[:5]],flush=True);items.append({'case':ci,'records':records,'sources':unique})
save({'ring':B,'items':items,'checks':checks},str(p/'independent_critical_two_RR_R_coefficients.sobj'),compress=False);print('PASS',len(checks),'SECONDS',time.monotonic()-st,flush=True)
