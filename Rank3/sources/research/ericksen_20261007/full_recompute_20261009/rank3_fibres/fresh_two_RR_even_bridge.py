"""Fresh full-quartic linkage for all three two-RR even supports.

The imported A-coefficients are explicit module witnesses, multiplied
against freshly reconstructed L*RR before any support pullback. No archived
source coefficient or stored remainder is accepted as evidence.
"""
from sage.all import *
from pathlib import Path
import json,hashlib,time
out=Path(__file__).absolute().parent;rank=out.parent/'rank3';research=out.parent.parent;st=time.monotonic();checks=[];imports=[]
def ck(n,v):
 assert bool(v),n
 checks.append(n)
def rd(path):
 imports.append(dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest()));return load(str(path))
def old(n):return rd(research/'round10'/(n+'.sobj'))
joint=rd(rank/'fresh_joint_fields.sobj');fg=rd(rank/'fresh_global_polynomials.sobj')
mod=old('independent_delta_even_L_module');probe=old('independent_subcritical_two_RR_even_probe');charts=old('independent_subcritical_missing_support_probe')
K=probe['field'];a,r,H,Z=K.gens();q=1+r;m=1+a*r
R=PolynomialRing(K,names=tuple('pp'+str(i) for i in range(6))+tuple('rr'+str(i) for i in range(3)));xx=R.gens()
OK=joint['ring'].base_ring();phi=OK.ring().hom([a,a*q/m,H/r,q*Z/r],K)
def cv(c):c=OK(c);return phi(c.numerator())/phi(c.denominator())
def raw(f):return R({tuple(map(int,e)):cv(c) for e,c in f.dict().items()})
L=matrix(R,3,3,[raw(f) for f in fg['L'].list()]);ell=L*vector(R,xx[6:])
def copy(f):return R({tuple(map(int,e))+(0,)*(9-len(e)):K(c) for e,c in f.dict().items()})
def part(f,bi):return R({e:c for e,c in f.dict().items() if (sum(e[:6]),sum(e[6:]))==bi})
sources=[]
for j,it in enumerate(mod['items']):
 qi=int(it['Q_index']);f=raw(joint['quartic_source_rows'][qi]);A=[copy(v) for v in it['A_p_quadratic']]
 sq=[2,0,1][j];top=part(f,(2,2));u=next(iter(top.dict().values()))/next(iter((ell[sq]**2).dict().values()))
 ck('fresh complete even square '+str(qi),top==u*ell[sq]**2)
 ck('fresh complete even linear module '+str(qi),part(f,(3,1))==sum(A[k]*ell[k] for k in range(3)))
 sources.append((u,sq,A,part(f,(4,0))))
S=probe['direction_ring'];E=probe['e_ring'];e=E.gens()
allfresh=[]
for case,it in enumerate(charts['items']):
 pp=it['p_images'];free=[2,0,1][case];active=[j for j in range(3) if j!=free]
 def pull(f):return sum((E(c)*prod(E(pp[k])**int(ex[k]) for k in range(6)) for ex,c in f.dict().items()),E.zero())
 ev=[]
 for j,(u,sq,A,c) in enumerate(sources):
  f=u*e[sq]**2+sum(pull(A[k])*e[k] for k in range(3))+pull(c)
  f=f.subs({e[k]:0 for k in active});ev.append(f)
  ck('fresh all support original even '+str((case,j)),f==probe['items'][case]['even_on_active_e_zero'][j])
 name='independent_subcritical_two_RR_full_normal'+('_case'+str(case) if case<2 else '')
 normal=old(name);P=normal['ring'];tt,ww,ss,ee=P.gens();imgs=[ww,ss,tt]
 def pcopy(f):return sum((P(K(c))*prod(imgs[k]**int(ex[k]) for k in range(3)) for ex,c in S(f).dict().items()),P.zero())
 def coeff(c):
  c=S.fraction_field()(c);nn=pcopy(c.numerator());dd=pcopy(c.denominator());ck('only coefficient-field denominator',dd.total_degree()==0 and dd!=0);return nn/K(dd.constant_coefficient())
 for j,p in enumerate(pp):ck('normal source same support chart '+str((case,j)),coeff(p)==normal['p_images'][j])
 mapped=[]
 for j,f in enumerate(ev):
  g=P.zero()
  for ex,c in f.dict().items():
   ck('single free even coordinate '+str((case,j,ex)),all(int(ex[k])==0 for k in active));g+=coeff(c)*ee**int(ex[free])
  mapped.append(g)
 allfresh.append(mapped)
targets=old('independent_subcritical_two_RR_three_even_transport')
for case in [0,1]:
 P=allfresh[case][0].parent()
 for j,g in enumerate(targets['items'][case]['original_even']):
  ck('transport true even linkage '+str((case,j)),allfresh[case][j]==P({tuple(map(int,e)):K(c) for e,c in g.dict().items()}))
 d=old('independent_three_D_F_even_transport');it=d['items'][case]
 for j,g in enumerate(it['original_even']):ck('F even true original linkage '+str((case,j)),allfresh[case][j]==P({tuple(map(int,e)):K(c) for e,c in g.dict().items()}))
 for rec in it['records']:
  g=P({tuple(map(int,e)):K(c) for e,c in rec['linear_even'].dict().items()});ck('F even selected linear is true original',g in allfresh[case])
d=old('independent_subcritical_two_RR_delta_even_elimination');P=allfresh[2][0].parent()
for j,g in enumerate(d['even_sources']):ck('Dcurve true even linkage '+str(j),allfresh[2][j]==P({tuple(map(int,e)):K(c) for e,c in g.dict().items()}))
report=dict(status='PASS',count=len(checks),checks=checks,imports=imports,seconds=time.monotonic()-st,scope='All three complete two-RR even-quadratic pullbacks regenerated from fresh full quartics; every imported linear-module cofactor verified by multiplication; exact support charts and both D/F downstream source lists compared. Independent of pending25x5 R/Q0/H normal matrix reconstruction.')
(out/'fresh_two_RR_even_bridge.json').write_text(json.dumps(report,indent=2));save(dict(sources=allfresh,checks=checks),str(out/'fresh_two_RR_even_bridge.sobj'),compress=False)
print('PASS',len(checks),report['seconds'])
