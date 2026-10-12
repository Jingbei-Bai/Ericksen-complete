"""COMPLETE source l^2 divisibility via exact reduction modulo l^2.
All s,w,e2 powers are retained; this is not a total normal two-jet.
"""
from sage.all import *
from pathlib import Path
import os,resource,time,json,sys,gc
p=Path(__file__).resolve().parent;st=time.time();d=load(str(p/'independent_subcritical_two_RR_inactive_axis_kappa_zero_adapted.sobj'))
vm=int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE');resource.setrlimit(resource.RLIMIT_AS,(vm+500*1024**2,vm+500*1024**2))
print('LinuxPID',os.getpid(),flush=True)
P=d['ring'];F=P.base_ring();l,s,w,E0,e1,e2,de=P.gens();rows=d['rows'];nv=[E0,e1,de];checks=[];items=[]
def tr(f):return P({tuple(map(int,ex)):co for ex,co in P(f).dict().items() if int(ex[0])<2})
def mm(f,g):return tr(tr(f)*tr(g))
def aa(f,g):return tr(f+g)
def det3(A):return aa(aa(mm(A[0][0],aa(mm(A[1][1],A[2][2]),-mm(A[1][2],A[2][1]))),-mm(A[0][1],aa(mm(A[1][0],A[2][2]),-mm(A[1][2],A[2][0])))),mm(A[0][2],aa(mm(A[1][0],A[2][1]),-mm(A[1][1],A[2][0]))))
for hn in ([sys.argv[1]] if len(sys.argv)>1 else ['H9','H15']):
 names=['Q0','R3',hn];A=[[tr(rows[n].derivative(v)) for v in nv] for n in names];rhs=[tr(P(rows[n].subs({v:0 for v in nv}))) for n in names];D=det3(A);print('D modulo l2',hn,len(D.dict()),time.time()-st,flush=True);nums=[]
 for n in names+['R1','H7','H3']:assert all(rows[n].derivative(v).derivative(vv)==0 for v in nv for vv in nv),('complete affine row',n)
 for j in range(3):
  AA=[rr[:] for rr in A]
  for i in range(3):AA[i][j]=-rhs[i]
  nums.append(det3(AA));print('CRAMER modulo l2',j,len(nums[-1].dict()),time.time()-st,flush=True)
 res={}
 for n in ['R1','H7','H3']:
  f=rows[n];ff=mm(tr(P(f.subs({v:0 for v in nv}))),D)
  for j in range(3):ff=aa(ff,mm(f.derivative(nv[j]),nums[j]))
  res[n]=ff;assert ff==0,('FULL SOURCE l2 divisibility failed',hn,n);checks.append('COMPLETE '+hn+' '+n+' numerator divisible l2');print('FULL LOW-L ZERO',hn,n,time.time()-st,flush=True)
  save({'ring':P,'source':d,'chart':hn,'range_low_l_matrix':A,'range_low_l_numerators':nums,'range_low_l_determinant':D,'complete_low_l_residual':res.copy(),'checks':checks.copy(),'scope':'COMPLETE source l2 divisibility checkpoint; no normal truncation'},str(p/('independent_subcritical_two_RR_inactive_axis_generic_l2_'+hn+'_checkpoint.sobj')),compress=False);gc.collect()
 items.append({'range_names':names,'full_source_chart':d,'range_low_l_matrix':A,'range_low_l_numerators':nums,'range_low_l_determinant':D,'complete_low_l_residual':res})
 save({'ring':P,'items':items,'checks':checks,'scope':'Exact full-source l2 divisibility at BOTH original range charts; no s/w/z truncation'},str(p/'independent_subcritical_two_RR_inactive_axis_generic_l2_divisibility.sobj'),compress=False)
rep={'PASS':True,'checks':len(checks),'seconds':time.time()-st,'scope':'Full generic source l2 divisibility, arbitrary s/w/z, both charts'};(p/'independent_subcritical_two_RR_inactive_axis_generic_l2_divisibility.json').write_text(json.dumps(rep,indent=2));print(rep,flush=True)
