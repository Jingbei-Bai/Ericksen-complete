"""Whole generic paired boundary range and true E0^4 source coefficient.

All coefficients are obtained by evaluating the complete source polynomials,
with exact range corrections of orders two and three.  No Hessian routine.
"""
from sage.all import *
from pathlib import Path
import os,resource,time,json
p=Path(__file__).resolve().parent;st=time.time();d=load(str(p/'independent_subcritical_two_RR_inactive_pair_generic_boundary.sobj'))
vm=int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE');resource.setrlimit(resource.RLIMIT_AS,(vm+500*1024**2,vm+500*1024**2))
print('LinuxPID',os.getpid(),flush=True)
P=d['ring'];F=P.base_ring();a,r,H,Z,l=F.gens();q=1+r;m=1+a*r;s,v,E,e1,e2,de=P.gens();rows=d['rows'];labels=d['labels'];checks=[]
def ck(n,t):
 if not t:raise AssertionError(n)
 checks.append(n)
zero={s:0,v:0,E:0,e1:0,e2:0,de:0}
nv=[e1,e2,de];rr=[0,7,28];A=matrix(F,[[F(rows[j].derivative(x).subs(zero)) for x in nv] for j in rr]);adet=A.det()
ck('true three affine range determinant nonzero',bool(adet))
for j in rr:
 ck('complete affine e1/e2/delta '+labels[j],all(rows[j].derivative(x).derivative(y)==0 for x in nv for y in nv))
 ck('complete '+labels[j]+' zero for arbitrary s/v/E on affine graph',P(rows[j].subs({x:0 for x in nv}))==0)
bf=d['zero_range_boundary'];B=matrix(F,[[F(bf[j].derivative(x).subs(zero)) for x in [s,v]] for j in [1,2]])
ck('two genuine complete even range rows unit generically',bool(B.det()))
T=PolynomialRing(F,'ep');ep=T.gen()
def tr(f):return T({int(i):co for i,co in T(f).dict().items() if int(i)<=4})
def ev(f,ss,vv):
 vals=[ss,vv,ep,0,0,0];out=T.zero()
 for ex,co in P(f).dict().items():
  mon=T.one()
  for j,n in enumerate(ex):
   if n:mon=tr(mon*tr(T(vals[j])**int(n)))
   if not mon:break
  if mon:out+=co*mon
 return tr(out)
path2=[ev(f,T.zero(),T.zero()) for f in bf];c2=-B.inverse()*vector(F,[path2[j][2] for j in [1,2]])
ss=c2[0]*ep**2;vv=c2[1]*ep**2;at2=[ev(f,ss,vv) for f in bf]
ck('all COMPLETE boundary source coefficients through order two zero',all(all(not f[h] for h in range(3)) for f in at2))
c3=-B.inverse()*vector(F,[at2[j][3] for j in [1,2]])
ss+=c3[0]*ep**3;vv+=c3[1]*ep**3;at3=[ev(f,ss,vv) for f in bf]
ck('all COMPLETE boundary source coefficients through order three zero',all(all(not f[h] for h in range(4)) for f in at3))
quartic=at3[29][4]
ck('pair1 pure fourth coefficient nonzero generically',bool(quartic))
cyF=PolynomialRing(QQ,names=('a','r','l')).fraction_field();aa,rr0,ll=cyF.gens();ph=F.ring().hom([aa,rr0,1,-1,ll],cyF)
def cyc(c):return ph(c.numerator())/ph(c.denominator())
ck('pair1 compensated true fourth is cyclic unit',cyc(quartic)==1/(16384*ll*((1+aa*rr0)*ll*ll-aa*rr0)))
print('cyclic range3det',cyc(adet).factor(),flush=True)
ck('protected affine source range unit for all l',cyc(adet)==-9*(1+rr0)**2*ll/(14336*rr0**2))
print('cyclic even range2det',cyc(B.det()).factor(),flush=True)
print('cyclic pure fourth',cyc(quartic).factor(),flush=True)
save({'generic_boundary':d,'affine_range_rows':rr,'affine_range_matrix':A,'even_range_rows':[1,2],'even_range_matrix':B,'true_second_range':c2,'true_third_range':c3,'complete_true_path_to_order4':at3,'pair1_fourth':quartic,'cyclic_pair1_fourth':cyc(quartic),'checks':checks,'scope':'COMPLETE generic whole boundary implies pure E0 residual order >=4; pair1 true fourth unit off two amplitude walls; protected affine range covers all amplitudes'},str(p/'independent_subcritical_two_RR_inactive_pair_boundary_fourth_certificate.sobj'),compress=False)
rep={'PASS':True,'checks':len(checks),'seconds':time.time()-st};(p/'independent_subcritical_two_RR_inactive_pair_boundary_fourth_certificate.json').write_text(json.dumps(rep,indent=2));print(rep,flush=True)
