"""Seven protected equations before last pressure; find a whole real pivot."""
from sage.all import *
from pathlib import Path
import time,signal,json,resource,sys
O=Path(__file__).resolve().parent;st=time.time();signal.alarm(40);resource.setrlimit(resource.RLIMIT_CORE,(0,0))
delta=int(sys.argv[1]) if len(sys.argv)>1 else 1
tag='plus' if delta==1 else 'minus'
d=load(str(O/f'rr_null_physical_hessian_source_{tag}.sobj'));K=d['field'];PR=K.ring();s,t,p,q,c=K.gens();R=d['ring'];vv=R.gens();hv=vv[2:]
M=d['matrix'][:7,:];off=d['offset'][:7];allowed=[0,1,3,4,6,7,8];free=[2,5]
rows=[]
for j in range(7):
 den=lcm(x.denominator() for x in M.row(j));rows.append([PR(den*x) for x in M.row(j)]+[den if j==i else PR(0) for i in range(7)])
act=list(range(7));piv=[]
for h in range(7):
 _,i,j=min([((rows[i][j].total_degree(),len(rows[i][j].dict())),i,j) for i in act for j in allowed if rows[i][j]])
 act.remove(i);allowed.remove(j);pp=rows[i][j];piv.append((i,j,pp))
 for l in act:
  qq=rows[l][j]
  if not qq:continue
  gg=pp.gcd(qq);nr=[(pp//gg)*rows[l][m]-(qq//gg)*rows[i][m] for m in range(16)]
  gg=gcd([z for z in nr if z]);rows[l]=[z//gg for z in nr] if gg and not gg.is_constant() else nr
T=PolynomialRing(K,names=('alpha','beta','x','y'));aa,bb,x,y=T.gens();sol=[None]*9;sol[2]=x;sol[5]=y
ph=R.hom([aa,bb]+[T(0)]*9,T);oo=list(map(ph,off))
for i,j,pp in reversed(piv):
 rhs=-sum(K(rows[i][9+l])*oo[l] for l in range(7))-sum(K(rows[i][m])*sol[m] for m in range(9) if m!=j and sol[m] is not None)
 sol[j]=T(rhs/K(pp))
phi=R.hom([aa,bb]+sol,T);last=phi(d['first_equations'][-1])
assert all(phi(f)==0 for f in d['first_equations'][:7])
rep={'pivots':[str(z.factor()) for _,_,z in piv], 'last_pressure_coefficients':{str(v):str(last.monomial_coefficient(v).factor()) for v in T.gens()},'seconds':time.time()-st}
lx=last.monomial_coefficient(x)
assert lx==-2*(t*t*p*p+s*s*q*q)/(s*s*t*t*c*(s*s-t*t))
xsol=-last.subs({x:0})/lx
TT=PolynomialRing(K,names=('alpha','beta','chi'));aaa,bbb,chi=TT.gens()
lift=T.hom([aaa,bbb,TT(0),chi],TT)
imgs=[lift(f.subs({x:xsol})) for f in sol]
hom=R.hom([aaa,bbb]+imgs,TT)
checks=dict(d['checks']);checks['all_eight_protected_linear_sources']=all(hom(f)==0 for f in d['first_equations'])
checks['protected_pressure_pivot']=lx==-2*(t*t*p*p+s*s*q*q)/(s*s*t*t*c*(s*s-t*t))
assert all(checks.values())
out={'field':K,'ring':TT,'F':d['F'].change_ring(K).change_ring(TT),'H':{key:vector(TT,map(hom,value)) for key,value in d['H'].items()},
 'axes':vector(TT,map(hom,d['axes'])),'axet':vector(TT,map(hom,d['axet'])),'graph':imgs,'checks':checks,'delta':delta,'pressure_pivot':lx,'pivots':piv}
save(out,str(O/f'rr_null_physical_protected_graph_{tag}.sobj'))
rep={'status':'WHOLE_REAL_PROTECTED_NULL_INFINITY_GRAPH','checks':checks,'delta':delta,'free_variable':'chi=H12_target3',
 'last_pivot':str(lx.factor()),'source_terms':[len(f.dict()) for f in imgs],'seconds':time.time()-st}
(O/f'rr_null_physical_protected_graph_{tag}.json').write_text(json.dumps(rep,indent=2),encoding='utf8');print(json.dumps(rep,indent=2),flush=True)
