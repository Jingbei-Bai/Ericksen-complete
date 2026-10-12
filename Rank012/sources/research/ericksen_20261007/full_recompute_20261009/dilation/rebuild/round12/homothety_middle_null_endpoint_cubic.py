"""Complete original-source cubic face at p=1,q=0, no division by q.

Only a next-face certificate.  The two real null directions are encoded by
real r!=0 and b=a*(1+r^2)/(a+r^2), with 0<a<1.
"""
from sage.all import *
from pathlib import Path
import pickle,json,time,resource,sys
OUT=Path(__file__).absolute().parent;START=time.time()
PZERO='--pzero' in sys.argv
STEM='homothety_middle_null_endpoint_cubic'+('_pzero' if PZERO else '')
vm=int(Path('/proc/self/statm').read_text().split()[0])*4096
resource.setrlimit(resource.RLIMIT_AS,(vm+650*1024**2,vm+650*1024**2))
d=load(str(OUT/'homothety_axis_plane_forcing_schur.sobj'));dd=d['source']
qd=load(str(OUT/'homothety_axis_plane_raw_quadratic_schur.sobj'))
K=FractionField(PolynomialRing(QQ,names=('a','r','A2','B2','k1' if PZERO else 'k2')))
a,r,A2,B2,kfree=K.gens();b=a*(1+r*r)/(a+r*r)
pval,qval=(0,1) if PZERO else (1,0)
k1,k2=(kfree,K(-1)) if PZERO else (K(-1),kfree)
S=PolynomialRing(K,names=('E','t2','h2','z2'));E,t2,h2,z2=S.gens()
R=PolynomialRing(S,'delta');delta=R.gen()
def tr(f):return R([f[i] for i in range(min(3,f.degree())+1)]) if f else R(0)
def cf4(c):return K(c.numerator()(a,b,pval,qval))/K(c.denominator()(a,b,pval,qval))
def cforig(c):return K(c.numerator()(a,b))/K(c.denominator()(a,b))
def substitute(poly,vals,cfmap):
 ans=R(0)
 for ex,c in poly.dict().items():
  term=R(cfmap(c))
  for i,ej in enumerate(ex):
   if ej:term=tr(term*vals[i]**ej)
  ans+=term
 return tr(ans)
obs=list(map(R,[0,0,A2,0,0,B2,k1,k2]));sv=delta;tv=r*delta+t2*delta**2;ev=E*delta**2
frvals=obs+[h2*delta**2,z2*delta**2,sv,tv,ev]
yfirst=[substitute(z,frvals,cf4) for z in d['state_linear_graph']]
vals=[R(cf4(c))+v for c,v in zip(dd['base_state'],yfirst)]
cols=dd['combined_pivot_columns'];rows=dd['combined_pivot_rows']
for i,z in zip(cols,qd['second_state_pivot_graph']):
 vals[i]+=substitute(z,[delta,r*delta,R(0),R(0)],cf4)
src=pickle.loads((OUT.parent/'round3/generic_degree5_sparse_prepared.sobj').read_bytes())
powers=[[R(1),v] for v in vals]
for j in range(17):
 for k in range(2,6):powers[j].append(tr(powers[j][-1]*vals[j]))
memo={}
def monomial(ex):
 if ex in memo:return memo[ex]
 ans=R(1)
 for i,ej in enumerate(ex):
  if ej:ans=tr(ans*powers[i][ej])
 memo[ex]=ans;return ans
def pull(comp):return sum(cforig(cf)*monomial(ex) for ex,cf in comp.items())
F=[[pull(c) for c in row] for row in src['Fraw']]
weights=[1]*9+[2]*6+[3]*2
H=[tr(F[2][j]+sv*F[0][j]+tv*F[1][j]+weights[j]*ev*vals[j]) for j in range(17)]
pi=matrix(R,3,3,vals[:9]);lam=diagonal_matrix(K,[1,1/a,1/b]);I=matrix.identity(R,3);n=vector(R,[sv,tv,1])
pm=matrix(R,[[pi[1,2]-pi[2,1],pi[0,2],-pi[0,1]],[-pi[1,2],pi[2,0]-pi[0,2],pi[1,0]],[pi[2,1],-pi[2,0],pi[0,1]-pi[1,0]]])
obsrc=list((pi.transpose()-pi.trace()*I/2)*n)+list((pi.transpose()*lam-(lam*pi).trace()*I/2)*n)+list((pm*n)[:2])
raw=[tr(v-ev*w) for v,w in zip(obsrc,obs)]+H
for name,seq in [('C3',src['C3']),('Q4',src['Q4']),('G5',src['basis'][:27])]:
 raw.extend(map(pull,seq));print('pulled endpoint',name,'seconds',time.time()-START,flush=True)
checks={f'full_source_order_{j}_{k}':f[k]==0 for j,f in enumerate(raw) for k in range(3)}
assert all(checks.values()),[k for k,v in checks.items() if not v]
M=dd['observation_plus_normal_matrix'].apply_map(cf4);block=M.matrix_from_rows_and_columns(rows,cols)
cubic=vector(S,[f[3] for f in raw]);state3=-block.inverse()*vector(S,[cubic[i] for i in rows]);cubic+=M.matrix_from_columns(cols)*state3
forcing=d['forcing_matrix'].apply_map(cf4)
frows,fcols=([3,4,11,21,27],[0,1,3,4,7]) if PZERO else ([3,4,9,22,27],[0,1,3,4,6])
forcing_pivot=forcing.matrix_from_rows_and_columns(frows,fcols)
comb=forcing.matrix_from_columns(fcols)*forcing_pivot.inverse()
checks['whole_endpoint_forcing_projection']=comb*forcing.matrix_from_rows(frows)==forcing
assert checks['whole_endpoint_forcing_projection']
projected=cubic-comb*vector(S,[cubic[i] for i in frows])
CM=matrix(K,[[f.monomial_coefficient(v) for v in S.gens()] for f in projected])
rhs=vector(K,[-f.constant_coefficient() for f in projected])
assert all(f==sum(CM[i,j]*S.gen(j) for j in range(4))-rhs[i] for i,f in enumerate(projected))
cr=list(CM.transpose().pivots());cc=list(CM.pivots());AA=CM.matrix_from_rows_and_columns(cr,cc)
solution_piv=AA.inverse()*vector(K,[rhs[i] for i in cr]);res=CM.matrix_from_columns(cc)*solution_piv-rhs
checks['all_cubic_rows_reconstructed']=all(z==0 for z in res)
assert checks['all_cubic_rows_reconstructed']
summary={'status':'FULL_ENDPOINT_NEXT_FACE_ONLY','scope':('p=0,q=1; 0<a<1; real r!=0; arbitrary A2,B2,k1; k2=-1' if PZERO else 'p=1,q=0; 0<a<1; real r!=0; arbitrary A2,B2,k2; k1=-1'),
 'source_checks':sum(checks.values()),'combined_pivot':str(block.det().factor()),
 'forcing_pivot':str(forcing_pivot.det().factor()),'rank':len(cc),'rows':cr,'cols':cc,
 'cubic_minor':str(AA.det().factor()),
 'solution_with_free_columns_zero':{str(S.gen(i)):str(z.factor()) if z else '0' for i,z in zip(cc,solution_piv)},
 'kernel':[[str(z.factor()) if z else '0' for z in vv] for vv in CM.right_kernel().basis()],
 'nonzero_forcing_free_cubic':[(i,str(z.factor())) for i,z in enumerate(projected) if z],
 'seconds':float(time.time()-START)}
save({'ring':S,'parameter_field':K,'material_b':b,'pressure_k1':k1,'pressure_k2':k2,'raw_cubic':cubic,
 'forcing_free_cubic':projected,'base_state':vals,'third_state_pivot_graph':state3,
 'source_checks':checks,'matrix':CM,'rhs':rhs,'summary':summary},str(OUT/(STEM+'.sobj')))
(OUT/(STEM+'.json')).write_text(json.dumps(summary,indent=2))
print(json.dumps(summary),flush=True)
