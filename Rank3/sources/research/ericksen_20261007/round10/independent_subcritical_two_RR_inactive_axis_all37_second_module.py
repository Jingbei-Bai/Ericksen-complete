from sage.all import *
from pathlib import Path
import time,json,resource,os
p=Path(__file__).resolve().parent;st=time.time();print('LinuxPID',os.getpid(),flush=True);d=load(str(p/'independent_subcritical_two_RR_inactive_axis_second_module.sobj'));P=d['ring'];F=P.base_ring();a,r,rho=F.gens();l,kap,s,z,w,E0,e1,de=P.gens();checks=[]
vm=int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE');resource.setrlimit(resource.RLIMIT_AS,(vm+650*1024**2,vm+650*1024**2))
def ck(n,b):
 if not b:raise AssertionError(n)
 checks.append(n)
old=load(str(p/'independent_ratio_p12_axis_normal_fourjets.sobj'));ok=old['coefficient_field'];ph=ok.ring().hom([a,r],F)
LP=LaurentPolynomialRing(F,names=('K','L','S','Z','W','X','Y','D'));K,L,S,Z,W,X,Y,D=LP.gens();images=[LP(rho)*K**(-1),K**(-1),S,(W+K*Z)/rho,(L*S-K*Y+K*X*S)/rho,W,L,L*K**(-1)+X,K*K*D]
rows=[]
for idx,ff in enumerate(old['normal_four_jets']):
 out=LP.zero();fac=K if idx==0 or idx>=10 else LP.one()
 for ee,c in ff.dict().items():
  # Original normal degree five terms cannot enter the weighted second jet:
  # total ordinary degree is at most five, RR degree at most three, hence
  # an omitted five-normal monomial has >=two p normals and final degree>=3.
  ex=tuple(map(int,ee));cc=ph(c.numerator())/ph(c.denominator());mon=fac*LP(cc)
  for j,n in enumerate(ex):
   if n:mon*=images[j]**n
  out+=LP({ee:co for ee,co in mon.dict().items() if sum(ee)<=2})
 ck('no kappa pole after complete source '+str(idx),all(int(ee[0])>=0 for ee in out.dict()))
 pf=P({(int(ee[1]),int(ee[0]))+tuple(map(int,ee[2:])):co for ee,co in out.dict().items()});rows.append(pf)
 if idx<28:ck('full ordinary fourjet vs direct weighted original jet '+str(idx),pf==d['rows'][idx])
 print('transport',idx,'terms',len(pf.dict()),'sec',round(time.time()-st,2),flush=True)
rows+=d['rows'][28:];norms=[w,E0,e1,de];gens=[l*s,l*z,kap*s,kap*z,s*s,s*z,z*z]
N=matrix(F,[[f.derivative(v).constant_coefficient() for v in norms] for f in rows]);items=[]
for rr in [(0,1,7,19),(0,1,7,25)]:
 A=N.matrix_from_rows(rr);qrows=[P({ee:co for ee,co in P(f.subs({v:0 for v in norms})).dict().items() if sum(ee)==2}) for f in rows];comp=-A.inverse()*vector(P,[qrows[j] for j in rr]);res=[qrows[j]+sum((N[j,h]*comp[h] for h in range(4)),P.zero()) for j in range(len(rows))]
 for j,f in enumerate(res):ck('seven-generator full source span '+str(rr)+'/'+str(j),f==sum((f.monomial_coefficient(g)*g for g in gens),P.zero()))
 M=matrix(F,[[f.monomial_coefficient(g) for g in gens] for f in res]);rank=M.rank();pv=M.transpose().pivots();print('MODULE',rr,'rank',rank,'pivots',pv,flush=True)
 # Preserve the exact matrix, even if a state wall remains.
 minors=[]
 if rank>=5:
  MM=M.matrix_from_columns(M.pivots());rv=MM.transpose().pivots();print('MINOR',rr,'cols',M.pivots(),'rows',rv,'factor',MM.matrix_from_rows(rv).det().factor(),flush=True);minors.append({'rows':rv,'cols':M.pivots(),'value':MM.matrix_from_rows(rv).det()})
 items.append({'range_rows':rr,'range_matrix':A,'range_compensation':comp,'residuals':res,'module_matrix':M,'module_rank':rank,'selected_minors':minors})
save({'ring':P,'rows':rows,'normal_order':norms,'generators':gens,'items':items,'checks':checks,'scope':'Cyclic all-material full37+3pair true weighted second source transport, membership and generic-base certificate pending'},str(p/'independent_subcritical_two_RR_inactive_axis_all37_second_module.sobj'),compress=False)
rep={'PASS':True,'checks':len(checks),'seconds':time.time()-st,'module_ranks':[x['module_rank'] for x in items]};(p/'independent_subcritical_two_RR_inactive_axis_all37_second_module.json').write_text(json.dumps(rep,indent=2));print(rep,flush=True)
