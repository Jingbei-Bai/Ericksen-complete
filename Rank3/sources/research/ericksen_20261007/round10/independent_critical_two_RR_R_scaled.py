from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_critical_two_RR_R_coefficients.sobj'));BR=d['ring'];K=PolynomialRing(QQ,names=('kr',)).fraction_field();r=K.gen();P=PolynomialRing(K,names=('dt','dw','ds'));t,w,s=P.gens();Q=PolynomialRing(QQ,names=('T','W','S'));T,W,S=Q.gens();items=[];checks=[]
def ck(n,v):
 if not bool(v):raise AssertionError(n)
 checks.append(n)
def fact(c):
 ee=vector(QQ,[0,0])
 for f,pow in c.factor():
  if f in QQ:continue
  inds=[j for j,g in enumerate([r,1+r]) if f/g in QQ and f/g]
  if not inds:raise ValueError(f)
  ee[inds[0]]+=pow
 con=c/(r**ZZ(ee[0])*(1+r)**ZZ(ee[1]));ck('exact protected material factors',con in QQ and con);return QQ(con),ee
for it in d['items']:
 records=[];AA=[];BB=[];orig=[]
 for f in it['sources']:
  g=P({tuple(map(int,ex[1:])):K(c)*r**int(ex[0]) for ex,c in f.dict().items()}) if False else sum(K(c)*r**int(ex[0])*t**int(ex[1])*w**int(ex[2])*s**int(ex[3]) for ex,c in f.dict().items());orig.append(g);rec=[]
  for ex,c in g.dict().items():
   cc,ee=fact(c);rec.append((tuple(map(int,ex)),cc,ee))
  x0,c0,e0=rec[0]
  for ex,c,ee in rec[1:]:AA.append(vector(QQ,ex)-vector(QQ,x0));BB.append(-(ee-e0))
  records.append(rec)
 A=matrix(QQ,AA);B=matrix(QQ,BB);V=A.solve_right(B);ck('exact scaling system',A*V==B);polys=[];scales=[]
 for g,rec in zip(orig,records):
  x0,c0,e0=rec[0];common=e0+vector(QQ,x0)*V;pol=Q({ex:cc for ex,cc,ee in rec})
  for ex,c,ee in rec:ck('whole coefficient scaling',ee+vector(QQ,ex)*V==common)
  polys.append(pol);scales.append(common)
 print('CASE',it['case'],'SCALES',V,flush=True);I=Q.ideal(polys);JJ,exp=I.saturation(Q.ideal(T));gb=JJ.groebner_basis();print('SAT',exp,list(map(str,gb)),flush=True)
 items.append({'case':it['case'],'original_sources':orig,'QQ_sources':polys,'direction_scales':V,'source_scales':scales,'saturated_basis':gb,'saturation_exponent':exp})
save({'ring':Q,'items':items,'checks':checks},str(p/'independent_critical_two_RR_R_scaled.sobj'),compress=False);print('PASS',len(checks),'SECONDS',time.monotonic()-st,flush=True)
