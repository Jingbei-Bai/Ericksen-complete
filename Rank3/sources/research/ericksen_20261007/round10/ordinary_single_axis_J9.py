"""Full relative J_axis9 coefficients from the existing genuine two-jets.
No original source is regenerated; all kernel/range couplings are retained.
"""
from sage.all import *
from pathlib import Path
import os,resource,time,json
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
d=load(str(p/'independent_ratio_p12_axis_generic_normal_two.sobj'))
gp=load(str(p/'independent_ratio_p_boundary_global_pair_syzygies.sobj'))
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+500*1024**2,)*2)
U=d['ring'];K=U.base_ring();a,b,h,z=K.gens();u,v=U.gens();ex=[tuple(map(int,e)) for e in d['normal_exponents']];slot={e:i for i,e in enumerate(ex)};coeff=[list(row) for row in d['coefficient_matrix'].rows()]
def ck(n,x):assert bool(x),n;checks.append(n)
ck('all37 rows',len(coeff)==37)
def ee(*ix):
 e=[0]*7
 for i in ix:e[i]+=1
 return tuple(e)
def c(row,*ix):return row[slot[ee(*ix)]]
M=d['linear_matrix'];ce=K(M[0,5].monomial_coefficient(u));lam=U(-M[0,4]/(ce*u));beta=U(-M[13,3]/K(M[13,0]))
for i in range(37):
 ck('axis kernel t '+str(i),M[i,0]*beta+M[i,3]==0)
 ck('incoming amplitude axis kernel '+str(i),M[i,4]+lam*M[i,5]==0)
 ck('eta divided by u '+str(i),U(M[i,5]/u)*u==M[i,5])
print('LOADED',time.monotonic()-st,'beta',beta,'lambda',lam,flush=True)
labels=[('old',i) for i in range(37)]
for j,it in enumerate(gp['items']):
 row=[U.zero() for _ in ex];boundary=U.zero();f=it['full']
 ck('fullpair degree4 '+str(j),f.total_degree()<=4)
 for e0,cc in f.dict().items():
  e=tuple(map(int,e0));ne=tuple(e[k] for k in (1,2,3,4,5,6))+(4-sum(e),)
  if min(ne)<0:raise AssertionError(ne)
  if sum(ne)>2:continue
  vv=K(cc)*u**e[7]*v**e[8]
  if not sum(ne):boundary+=vv
  else:row[slot[ne]]+=vv
 ck('pair axis boundary '+str(j),not boundary)
 ck('pair kernel t '+str(j),c(row,0)*beta+c(row,3)==0)
 ck('pair kernel l '+str(j),c(row,4)+lam*c(row,5)==0)
 ck('pair eta/u '+str(j),U(c(row,5)/u)*u==c(row,5))
 coeff.append(row);labels.append(('pair',j))
rows=[]
for row in coeff:
 ranges=[c(row,j) for j in (0,1,2,6)]
 ueta=U(c(row,5)/u)
 lk=beta*c(row,0,4)+c(row,3,4)+lam*(beta*c(row,0,5)+c(row,3,5))
 et2=c(row,5,5);etk=beta*c(row,0,5)+c(row,3,5)
 k2=beta**2*c(row,0,0)+beta*c(row,0,3)+c(row,3,3)
 rows.append(ranges+[ueta,lk,et2,etk,k2])
A=matrix(U,rows);A0=A.apply_map(lambda f:U(f.subs({u:0})))
A7=matrix(U,[[row[j] for j in range(4)]+[u*row[4],row[5],row[8]] for row in rows])
save({'ring':U,'matrix':A,'matrix_u_zero':A0,'matrix_u_nonzero_seven':A7,'beta':beta,'lambda':lam,'labels':labels,'columns':('range0','n21','n23','tau','u_eta','l_k','eta2','eta_k','k2'),'source_axis_coefficients':coeff,'normal_exponents':ex,'checks':checks,'whole_ideal_proof':'ordinary_single_axis_relative_ideal_reduction.txt','source_files':('independent_ratio_p12_axis_generic_normal_two.sobj','independent_ratio_p_boundary_global_pair_syzygies.sobj')},str(p/'ordinary_single_axis_J9.sobj'),compress=False)
rep={'status':'PASS EXACT J9 COEFFICIENTS; NO RANK CLAIM','checks':len(checks),'seconds':time.monotonic()-st,'shape':list(A.dimensions()),'supports_u_zero':[[j for j,x in enumerate(row) if x] for row in A0.rows()]}
(p/'ordinary_single_axis_J9.json').write_text(json.dumps(rep,indent=2)+'\n');print(rep,flush=True)
