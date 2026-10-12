"""Exact relative zero-kappa range and odd sheet; retain all H,Z,k,r,w,t."""
from sage.all import *
from pathlib import Path
import time,resource,os
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_subcritical_single_RR_yw_kappa_zero_generic_sources.sobj'));resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+350*1024**2,)*2);print('PID',os.getpid(),flush=True);P=d['ring'];K=P.base_ring();x,z,w,e0,e1,e2,de,t=P.gens();U=PolynomialRing(K,names=('gr_x','gr_z','gr_w','gr_t'));xx,zz,ww,tt=U.gens();hm=P.hom([xx,zz,ww,U.zero(),U.zero(),U.zero(),U.zero(),tt],U);rows=d['range_rows'];A=matrix(U,4,4);B=vector(U,4)
for i,ri in enumerate(rows):
 f=d['sources'][ri];assert all(sum(int(ee[j]) for j in (3,4,5,6))<=1 for ee in f.dict()),ri;B[i]=hm(f)
 for j,ev in enumerate((e0,e1,e2,de)):A[i,j]=hm(f.derivative(ev))
det=A.det();N=-A.adjugate()*B
print('DET',len(det.dict()),str(det.factor())[:1800],'sec',time.monotonic()-st,flush=True)
for i,v in enumerate(N):print('NUM',i,len(v.dict()),str(v.factor())[:4000],'sec',time.monotonic()-st,flush=True)
save({'ring':U,'range_rows':rows,'matrix':A,'constant':B,'det':det,'numerators':N,'sources':d['sources'],'scope':__doc__},str(p/'independent_subcritical_single_RR_yw_kappa_zero_generic_range.sobj'),compress=False);print('SAVED',time.monotonic()-st,flush=True)
