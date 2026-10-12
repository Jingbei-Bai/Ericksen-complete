"""Whole main-direction physical chart for projective infinity, null axes.

F=[[s,0,0],[0,t,0],[p,q,c]], s,t,c>0, p*q*(s^2-t^2)!=0.
The two source-null signs are conjugate; delta=+/-1 covers target signs.
All pressures are HIGHEST homogeneous zero, not an actual real finite state.
"""
from sage.all import *
from pathlib import Path
import time,signal,json,resource,sys
O=Path(__file__).resolve().parent;start=time.time();signal.alarm(80)
resource.setrlimit(resource.RLIMIT_CORE,(0,0))
vm=int(Path('/proc/self/statm').read_text().split()[0])*4096
resource.setrlimit(resource.RLIMIT_AS,(vm+1100*1024**2,vm+1100*1024**2))
delta=int(sys.argv[1]) if len(sys.argv)>1 else 1
assert delta in (-1,1)
tag='plus' if delta==1 else 'minus'
QI=QuadraticField(-1,'ii');ii=QI.gen();PR=PolynomialRing(QI,names=('s','t','p','q','c'));K=PR.fraction_field();s,t,p,q,c=K.gens()
R=PolynomialRing(K,names=('alpha','beta')+tuple('h'+str(i) for i in range(9)))
alpha,beta=R.gens()[:2];hv=R.gens()[2:]
F=matrix(R,[[s,0,0],[0,t,0],[p,q,c]]);Fi=matrix(R,[[1/s,0,0],[0,1/t,0],[-p/(s*c),-q/(t*c),1/c]])
assert F*Fi==1
axes=vector(R,[alpha,ii*alpha,0]);axet=vector(R,[beta,delta*ii*beta,0])
e=[vector(R,[int(i==j) for i in range(3)]) for j in range(3)]
HH={}
for j,ids in [(0,(0,1,2)),(1,(3,4,5))]:pass
HH[0,0]=vector(R,hv[:3]);HH[0,1]=HH[1,0]=vector(R,hv[3:6]);HH[1,1]=vector(R,hv[6:9])
for a in range(3):HH[a,2]=HH[2,a]=axet.cross_product(F.column(a))-F*axes.cross_product(e[a])
DF=[matrix(R,3,3,lambda i,j:HH[j,a][i]) for a in range(3)]
C=F.transpose()*F;B=Fi*Fi.transpose();DC=[D.transpose()*F+F.transpose()*D for D in DF]
DB=[-B*D*B for D in DC]
first=[]
for a in range(2):first += [DC[a].trace(),(C*DC[a]).trace(),(Fi*DF[a]).trace()]
KC=sum(DC[a][a,2] for a in range(3));KB=sum(DB[a][a,2] for a in range(3));first += [KC,KB]
checks={'source_null':axes*axes==0,'target_null':axet*axet==0,
 'fixed_direction_spectrum_auto':DC[2].trace()==0 and (C*DC[2]).trace()==0 and (Fi*DF[2]).trace()==0,
 'pressure_C_trace':KC==c*sum(HH[a,a][2] for a in range(3))}
M=matrix(K,[[f.monomial_coefficient(v) for v in hv] for f in first]);off=vector(R,[f.subs(dict(zip(hv,[0]*9))) for f in first])
save({'field':K,'ring':R,'F':F,'Fi':Fi,'H':HH,'first_equations':first,'matrix':M,'offset':off,'axes':axes,'axet':axet,'delta':delta,'checks':checks},str(O/f'rr_null_physical_hessian_source_{tag}.sobj'))
print('raw source',time.time()-start,checks,flush=True)
# Primitive row elimination, with exact transformation to all eight original rows.
rows=[]
for j in range(8):
 den=lcm(x.denominator() for x in M.row(j));rows.append([PR(den*x) for x in M.row(j)]+[den if j==i else PR(0) for i in range(8)])
act=list(range(8));cols=list(range(9));piv=[]
for step in range(8):
 options=[((rows[i][j].total_degree(),len(rows[i][j].dict())),i,j) for i in act for j in cols if rows[i][j]]
 _,i,j=min(options);act.remove(i);cols.remove(j);pp=rows[i][j];piv.append((i,j,pp))
 print('pivot',step,i,j,pp.factor(),time.time()-start,flush=True)
 for l in act:
  qq=rows[l][j]
  if not qq:continue
  gg=pp.gcd(qq);a0=pp//gg;b0=qq//gg
  nr=[a0*rows[l][m]-b0*rows[i][m] for m in range(17)]
  nz=[z for z in nr if z];gg=gcd(nz)
  if gg and not gg.is_constant():nr=[z//gg for z in nr]
  rows[l]=nr
assert len(cols)==1
T=PolynomialRing(K,names=('alpha','beta','chi'));aa,bb,chi=T.gens();sol=[None]*9;sol[cols[0]]=chi
ph=R.hom([aa,bb]+[T(0)]*9,T);offT=vector(T,map(ph,off))
for i,j,pp in reversed(piv):
 rhs=-sum(K(rows[i][9+l])*offT[l] for l in range(8))
 rhs-=sum(K(rows[i][m])*sol[m] for m in range(9) if m!=j and sol[m] is not None)
 sol[j]=T(rhs/K(pp))
hom=R.hom([aa,bb]+sol,T);Hg={key:vector(T,map(hom,val)) for key,val in HH.items()}
checks['all_eight_linear_sources']=all(hom(f)==0 for f in first)
assert all(checks.values()),checks
denoms=sorted(set(str(z.denominator().factor()) for f in sol for z in f.coefficients()))
out={'field':K,'ring':T,'F':F.change_ring(K).change_ring(T),'H':Hg,'graph':sol,'free_index':cols[0],'pivots':piv,
 'axes':vector(T,[aa,ii*aa,0]),'axet':vector(T,[bb,delta*ii*bb,0]),'checks':checks,'delta':delta,'denominators':denoms}
save(out,str(O/f'rr_null_physical_hessian_graph_{tag}.sobj'))
rep={'status':'EXACT_FULL_PARAMETER_NULL_INFINITY_HESSIAN_GRAPH_NOT_FINITE_STATE','checks':checks,'free_index':cols[0],
 'pivots':[str(p.factor()) for _,_,p in piv],'denominators':denoms,'graph_terms':[len(f.dict()) for f in sol],'seconds':time.time()-start}
(O/f'rr_null_physical_hessian_graph_{tag}.json').write_text(json.dumps(rep,indent=2),encoding='utf8');print(json.dumps(rep,indent=2),flush=True)
