"""Nonresonant middle null next face with epsilon=E delta^3.

Complete original sources are evaluated before and after the 13-range
Schur elimination. No finite face is described as whole-arc exclusion.
"""
from sage.all import *
from pathlib import Path
import pickle,json,time,resource
OUT=Path(__file__).absolute().parent;START=time.time()
vm=int(Path('/proc/self/statm').read_text().split()[0])*4096
resource.setrlimit(resource.RLIMIT_AS,(vm+1400*1024**2,vm+1400*1024**2))
old=load(str(OUT/'homothety_middle_null_cubic_solution.sobj'));cubic=old['source'];K=cubic['parameter_field'];a,p,q,r,A2,B2,k1=K.gens()
b=a*(1+r*r)/(a+r*r);k2=cubic['pressure_k2'];Zc=old['solution'][3]
Os=cubic['ring'];ov=Os.gens();osub={ov[0]:0,ov[1]:0,ov[2]:0,ov[3]:Zc}
S=PolynomialRing(K,names=('E','t3','h3','z3'));E,t3,h3,z3=S.gens();nv=S.gens()
R=PolynomialRing(S,'delta');delta=R.gen()
def tr(f):return R([f[i] for i in range(min(4,f.degree())+1)]) if f else R(0)
vals=[]
for f in cubic['base_state']:
 ff=R(0)
 for n,co in enumerate(f.list()):
  ff+=K(co.subs(osub))*delta**n
  ff+=sum(K(co.derivative(v).subs(osub))*w for v,w in zip(ov,nv))*delta**(n+1)
 vals.append(tr(ff))
sd=load(str(OUT/'homothety_axis_plane_forcing_schur.sobj'));dd=sd['source']
def cf4(c):return K(c.numerator()(a,b,p,q))/K(c.denominator()(a,b,p,q))
def cforig(c):return K(c.numerator()(a,b))/K(c.denominator()(a,b))
M=dd['observation_plus_normal_matrix'].apply_map(cf4);rows=dd['combined_pivot_rows'];cols=dd['combined_pivot_columns'];BI=M.matrix_from_rows_and_columns(rows,cols).inverse()
src=pickle.loads((OUT.parent/'round3/generic_degree5_sparse_prepared.sobj').read_bytes());weights=[1]*9+[2]*6+[3]*2
sv=delta;tv=r*delta+t3*delta**3;ev=E*delta**3;obs=list(map(R,[0,0,A2,0,0,B2,k1,k2]))
def pull_factory(vs):
 powers=[[R(1),v] for v in vs]
 for j in range(17):
  for k in range(2,6):powers[j].append(tr(powers[j][-1]*vs[j]))
 memo={}
 def mon(ex):
  if ex in memo:return memo[ex]
  ans=R(1)
  for j,ej in enumerate(ex):
   if ej:ans=tr(ans*powers[j][ej])
  memo[ex]=ans;return ans
 def pull(comp):return sum(cforig(c)*mon(ex) for ex,c in comp.items())
 return pull
def low_sources(vs,pull):
 F=[[pull(c) for c in row] for row in src['Fraw']]
 H=[tr(F[2][j]+sv*F[0][j]+tv*F[1][j]+weights[j]*ev*vs[j]) for j in range(17)]
 pi=matrix(R,3,3,vs[:9]);lam=diagonal_matrix(K,[1,1/a,1/b]);I=matrix.identity(R,3);n=vector(R,[sv,tv,1])
 pm=matrix(R,[[pi[1,2]-pi[2,1],pi[0,2],-pi[0,1]],[-pi[1,2],pi[2,0]-pi[0,2],pi[1,0]],[pi[2,1],-pi[2,0],pi[0,1]-pi[1,0]]])
 oo=list((pi.transpose()-pi.trace()*I/2)*n)+list((pi.transpose()*lam-(lam*pi).trace()*I/2)*n)+list((pm*n)[:2])
 return [tr(v-ev*w) for v,w in zip(oo,obs)]+H
pull=pull_factory(vals);low=low_sources(vals,pull)
state3=-BI*vector(S,[low[i][3] for i in rows])
for j,v in zip(cols,state3):vals[j]+=v*delta**3
pull=pull_factory(vals);raw=low_sources(vals,pull)
for name,seq in [('C3',src['C3']),('Q4',src['Q4']),('G5',src['basis'][:27])]:
 raw.extend(map(pull,seq));print('quartic face',name,'seconds',time.time()-START,flush=True)
checks={f'raw_order_{i}_{k}':f[k]==0 for i,f in enumerate(raw) for k in range(4)}
assert all(checks.values()),[k for k,v in checks.items() if not v]
Q=vector(S,[f[4] for f in raw]);state4=-BI*vector(S,[Q[i] for i in rows]);Q+=M.matrix_from_columns(cols)*state4
forcing=sd['forcing_matrix'].apply_map(cf4);fr=[3,4,9,11,27];fc=[0,1,2,3,4]
comb=forcing.matrix_from_columns(fc)*forcing.matrix_from_rows_and_columns(fr,fc).inverse()
projected=Q-comb*vector(S,[Q[i] for i in fr])
save({'ring':S,'field':K,'quartic_raw':Q,'quartic_forcing_free':projected,'checks':checks,
      'state_jet':vals,'state4':state4,'material_b':b,'pressure_k2':k2},str(OUT/'homothety_middle_null_quartic_face.sobj'))
MM=matrix(K,[[f.monomial_coefficient(v) for v in nv] for f in projected]);rhs=vector(K,[-f.constant_coefficient() for f in projected])
assert all(f==sum(MM[i,j]*nv[j] for j in range(4))-rhs[i] for i,f in enumerate(projected))
def sample(c):return QQ(c.numerator()(QQ(1)/2,1,3,2,3,4,5))/QQ(c.denominator()(QQ(1)/2,1,3,2,3,4,5))
sm=MM.apply_map(sample);rr=[];rk=0
order=sorted(range(MM.nrows()),key=lambda i:sum(len(z.numerator().dict())+len(z.denominator().dict()) for z in MM.row(i) if z))
for i in order:
 k=sm.matrix_from_rows(rr+[i]).rank()
 if k>rk:rr.append(i);rk=k
cc=list(sm.matrix_from_rows(rr).pivots());assert len(cc)==4
A=MM.matrix_from_rows_and_columns(rr,cc);sol=A.inverse()*vector(K,[rhs[i] for i in rr]);res=MM*sol-rhs
summary={'status':'EXACT_NONRESONANT_QUARTIC_FACE_ONLY','checks':sum(checks.values()),'rows':rr,
 'minor':str(A.det().factor()),'solution':{str(v):str(z.factor()) if z else '0' for v,z in zip(nv,sol)},
 'nonzero_residuals':[(i,str(z.factor())) for i,z in enumerate(res) if z], 'seconds':float(time.time()-START)}
save({'field':K,'matrix':MM,'rhs':rhs,'solution':sol,'residuals':res,'summary':summary},str(OUT/'homothety_middle_null_quartic_solution.sobj'))
(OUT/'homothety_middle_null_quartic_solution.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
