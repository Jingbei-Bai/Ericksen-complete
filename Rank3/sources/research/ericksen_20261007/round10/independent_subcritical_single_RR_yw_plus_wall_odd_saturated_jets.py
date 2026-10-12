"""Genuine odd kappa-saturated sources: each zero-kappa remainder is an H4 multiple."""
from sage.all import *
from pathlib import Path
import resource,os,time
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_subcritical_single_RR_yw_plus_wall_saturated_cubic_probe.sobj'));resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+350*1024**2,)*2);print('PID',os.getpid(),flush=True);T=d['ring'];K=T.base_ring();r,W=K.gens();x,z,t,kap,u=T.gens();res=d['residuals4'];zero=T.hom([x,z,t,T.zero(),u],T);hom=lambda f,n:T({tuple(map(int,ee)):c for ee,c in f.dict().items() if sum(ee)==n});nt=lambda f:T({tuple(map(int,ee)):c for ee,c in f.dict().items() if sum(ee)<=4});H0=zero(res[11]);q=hom(H0,2);assert q==-2*W**2/5*t*(x+z);out=[]
for ri in (0,1,24,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22):
 f0=zero(res[ri]);cc=T.zero()
 for dd in (2,3,4):
  rem=hom(f0-cc*H0,dd)
  if rem:
   ff=T.fraction_field()(rem)/q;assert ff.denominator().is_constant(),(ri,dd,rem);cc+=T(ff)
 num=nt(res[ri]-cc*res[11]);assert zero(num)==0;assert all(int(ee[3])>=1 for ee in num.dict())
 pn=T({(int(ee[0]),int(ee[1]),int(ee[2]),int(ee[3])-1,int(ee[4])):c for ee,c in num.dict().items()});out.append({'row':ri,'cofactor':cc,'jet3':pn});print('ROW',ri,'COF',cc,'P2',hom(pn,2),'P1',hom(pn,1),'P3 TERMS',len(hom(pn,3).dict()),flush=True)
save({'ring':T,'items':out,'scope':__doc__},str(p/'independent_subcritical_single_RR_yw_plus_wall_odd_saturated_jets.sobj'),compress=False);print('SAVED',time.monotonic()-st,flush=True)
