"""Finite exact coefficient-module reader: whole-source membership precedes these unit checks."""
from sage.all import *
from pathlib import Path
import resource,os,time
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_subcritical_single_RR_yw_minus_odd_saturated_cubic_probe.sobj'));full=load(str(p/'independent_subcritical_single_RR_full_regular_chart_H2.sobj'));resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+450*1024**2,)*2);print('PID',os.getpid(),flush=True);T=d['ring'];K=T.base_ring();r,W=K.gens();k=-W**2;m=1+k;x,z,t,kap,u=T.gens();s=x+z
def nt(f):return T({tuple(map(int,ee)):c for ee,c in f.dict().items() if sum(ee)<=4})
P=full['ring'];OK=P.base_ring();ph=OK.ring().hom([k,r],K);ims=(x,z,W+u,(7*x+9*z)*(W+u),-(W+u)/2,16*x+14*z,27*r/(1+r)*(W+u)*s,T.zero(),T.zero());C=T.zero()
for ee,c in full['sources'][27].dict().items():
 val=T(ph(c.numerator())/ph(c.denominator()))
 for im,nn in zip(ims,ee):
  if nn:val=nt(val*im**int(nn))
 C+=val
C=nt(C);g=(kap*s,t*s,s*s,(t+2*kap)*x,C,kap*x**3);print('TRUE C JET4',C,flush=True)
def exponent_tuples(total,n):
 if n==1:yield (total,);return
 for j in range(total+1):
  for rest in exponent_tuples(total-j,n-1):yield (j,)+rest
mg=[]
for gg in g:
 order=min(sum(ee) for ee in gg.dict())
 for degree in range(1,5-order):
  for ee in exponent_tuples(degree,5):mg.append(nt(gg*T.monomial(*ee)))
allpoly=mg+list(g)+d['residuals4'];mons=sorted(set(tuple(map(int,ee)) for ff in allpoly for ee in ff.dict()));vec=lambda ff:vector(K,[ff.monomial_coefficient(T.monomial(*ee)) for ee in mons]);E=matrix(K,[vec(ff) for ff in mg]).echelon_form();EE=[row for row in E if any(row)];piv=[next(i for i,c in enumerate(row) if c) for row in EE];print('mJ finite RANK',len(EE),'COLS',len(mons),'sec',time.monotonic()-st,flush=True)
def red(v):
 v=vector(K,v)
 for row,pi in zip(EE,piv):
  if v[pi]:v-=v[pi]*row
 return v
B=matrix(K,[red(vec(gg)) for gg in g]);print('J/mJ RANK',B.rank(),'sec',time.monotonic()-st,flush=True);sel=B.pivots();assert len(sel)==6;BB=B.matrix_from_columns(sel);M=[]
for j,ff in enumerate(d['residuals4']):
 v=red(vec(ff));co=BB.transpose().solve_right(vector(K,[v[ii] for ii in sel]));assert co*B==v,('source coefficient membership failed',j);M.append(co)
M=matrix(K,M);rank=M.rank();rows=M.transpose().pivots();print('TRUE MODULE RANK',rank,'ROWS',rows,'sec',time.monotonic()-st,flush=True)
if rank==6:print('TRUE DET6',M.matrix_from_rows(rows).det().factor(),flush=True)
else:print('TRUE KERNEL',M.right_kernel().basis_matrix(),flush=True)
save({'ring':T,'generators':g,'mJ_jet_echelon':E,'monomials':mons,'basis':B,'selector':sel,'source_coefficient_matrix':M,'rank':rank,'rows':rows,'scope':__doc__},str(p/'independent_subcritical_single_RR_yw_minus_ten_module_reader.sobj'),compress=False);print('SAVED',time.monotonic()-st,flush=True)
