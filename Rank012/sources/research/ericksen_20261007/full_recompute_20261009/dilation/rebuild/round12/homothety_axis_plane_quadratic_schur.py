"""Full second normal Schur source with epsilon forcing projected out.

All p,q parameters remain free.  The result describes a tangent ideal only;
positivity, higher order estimates, and no-incoming-arc statements require
separate arguments.
"""
from sage.all import *
from pathlib import Path
import pickle,json,time,resource
OUT=Path(__file__).absolute().parent;START=time.time()
vm=int(Path('/proc/self/statm').read_text().split()[0])*4096
resource.setrlimit(resource.RLIMIT_AS,(vm+1400*1024**2,vm+1400*1024**2))
d=load(str(OUT/'homothety_axis_plane_forcing_schur.sobj'));dd=d['source'];K=dd['parameter_field'];a,b,p,q=K.gens()
R=PolynomialRing(K,names=('s','t','h','z'));s,t,h,z=R.gens();gens=R.gens()
FR=d['ring'];A0,A1,A2,B0,B1,B2,k1,k2,hh,zz,ss,tt,ee=FR.gens()
def toR(f):return R(f(0,0,0,0,0,0,0,0,h,z,s,t,0))
y1=list(map(toR,d['state_linear_graph']));base=dd['base_state'];val=[R(base[i])+y1[i] for i in range(17)]
def tr(f):return R({m:c for m,c in R(f).dict().items() if sum(m)<=2})
def hp(f,k):return R({m:c for m,c in R(f).dict().items() if sum(m)==k})
powers=[[R(1),v] for v in val]
for i in range(17):
 for k in range(2,6):powers[i].append(tr(powers[i][-1]*val[i]))
src=pickle.loads((OUT.parent/'round3/generic_degree5_sparse_prepared.sobj').read_bytes());memo={}
def monomial(ex):
 if ex in memo:return memo[ex]
 ans=R(1)
 for j,ej in enumerate(ex):
  if ej:ans=tr(ans*powers[j][ej])
 memo[ex]=ans;return ans
def pull(comp):return sum(K(cf)*monomial(ex) for ex,cf in comp.items())
FF=[[pull(comp) for comp in row] for row in src['Fraw']]
H=[tr(FF[2][i]+s*FF[0][i]+t*FF[1][i]) for i in range(17)]
pi=matrix(R,3,3,val[:9]);lam=diagonal_matrix(K,[1,1/a,1/b]);I=matrix.identity(R,3);n=vector(R,[s,t,1])
pm=matrix(R,[[pi[1,2]-pi[2,1],pi[0,2],-pi[0,1]],[-pi[1,2],pi[2,0]-pi[0,2],pi[1,0]],[pi[2,1],-pi[2,0],pi[0,1]-pi[1,0]]])
obs=list((pi.transpose()-pi.trace()*I/2)*n)+list((pi.transpose()*lam-(lam*pi).trace()*I/2)*n)+list((pm*n)[:2])
raw=list(map(tr,obs))+H
for name,seq in [('C3',src['C3']),('Q4',src['Q4']),('G5',src['basis'][:27])]:
 raw.extend(map(pull,seq));print('pulled',name,'seconds',time.time()-START,flush=True)
assert len(raw)==76
assert all(hp(f,0)==0 and hp(f,1)==0 for f in raw)
M=dd['observation_plus_normal_matrix'];rows=dd['combined_pivot_rows'];cols=dd['combined_pivot_columns'];AA=M.matrix_from_rows_and_columns(rows,cols)
quad=vector(R,[hp(f,2) for f in raw]);y2=AA.inverse()*vector(R,[-quad[i] for i in rows]);Q=list(quad+M.matrix_from_columns(cols)*y2)
assert all(Q[i]==0 for i in rows)
mons=[s*s,s*t,s*h,s*z,t*t,t*h,t*z,h*h,h*z,z*z]
QM=matrix(K,[[f.monomial_coefficient(m) for m in mons] for f in Q])
forcing=d['forcing_matrix']
save({'ring':R,'first_state_graph':y1,'second_state_pivot_graph':y2,'quadratic_residuals':Q,
      'monomials':mons,'quadratic_matrix':QM,'forcing_matrix':forcing},str(OUT/'homothety_axis_plane_raw_quadratic_schur.sobj'))
print('saved raw second Schur',time.time()-START,flush=True)
forcing_rows=[3,4,9,11,27];forcing_cols=[0,1,2,3,4]
FM=forcing.matrix_from_rows_and_columns(forcing_rows,forcing_cols)
comb=forcing.matrix_from_columns(forcing_cols)*FM.inverse()
assert comb*forcing.matrix_from_rows(forcing_rows)==forcing
projected=QM-comb*QM.matrix_from_rows(forcing_rows)
left=None
save({'ring':R,'projected_matrix':projected,'forcing_combinations':comb,'raw_matrix':QM},str(OUT/'homothety_axis_plane_projected_quadratics.sobj'))
print('saved projected quadratic matrix',time.time()-START,flush=True)
def exact_small_basis(Mat):
 # A rational sample ONLY proposes pivot rows. Every symbolic row is then
 # reconstructed exactly from the selected block; no sampled conclusion.
 def sample(c):return QQ(c.numerator()(2,4,1,3))/QQ(c.denominator()(2,4,1,3))
 sampled=Mat.apply_map(sample); rows=[]; rk=0
 order=sorted(range(Mat.nrows()),key=lambda i:sum(len(z.numerator().dict())+len(z.denominator().dict()) for z in Mat.row(i) if z))
 for i in order:
  trial=sampled.matrix_from_rows(rows+[i]);rrr=trial.rank()
  if rrr>rk:rows.append(i);rk=rrr
 cols=list(sampled.matrix_from_rows(rows).pivots());SS=Mat.matrix_from_rows(rows);block=SS.matrix_from_columns(cols)
 BI=block.inverse();BB=BI*SS;CC=Mat.matrix_from_columns(cols)
 assert CC*BB==Mat
 return rows,cols,BB,block.det()
rr,cc,basis,minor=exact_small_basis(projected);polys=[sum(row[i]*mons[i] for i in range(10)) for row in basis.rows()]
rawrank=len(exact_small_basis(QM)[0]);rank=len(rr)
strata={}
for name,sub in [('p_nonzero_q_zero',{q:0}),('p_zero_q_nonzero',{p:0})]:
 FF=forcing.apply_map(lambda c:c.subs(sub)); Qspec=QM.apply_map(lambda c:c.subs(sub))
 frows=list(FF.transpose().pivots());fcols=list(FF.pivots());fm=FF.matrix_from_rows_and_columns(frows,fcols)
 cb=FF.matrix_from_columns(fcols)*fm.inverse();assert cb*FF.matrix_from_rows(frows)==FF
 PP=Qspec-cb*Qspec.matrix_from_rows(frows)
 rrows,rcols,BB,bd=exact_small_basis(PP)
 strata[name]={'rank':len(rrows),'pivot_minor':str(bd.factor()),
  'basis':[str(sum(row[i]*mons[i] for i in range(10)).factor()) for row in BB.rows()]}
summary={'status':'EXACT_QUADRATIC_NORMAL_IDEAL_ONLY','unprojected_rank':rawrank,'forcing_free_rank':rank,
 'generic_pivot_minor':str(minor.factor()),'forcing_free_basis':[str(f.factor()) for f in polys],
 'strata':strata,'seconds':float(time.time()-START)}
save({'ring':R,'first_state_graph':y1,'second_state_pivot_graph':y2,'quadratic_residuals':Q,'monomials':mons,
      'quadratic_matrix':QM,'forcing_row_combinations':comb,'projected_quadratic_matrix':projected,
      'generic_basis':polys,'source':d,'summary':summary},str(OUT/'homothety_axis_plane_quadratic_schur.sobj'))
(OUT/'homothety_axis_plane_quadratic_schur.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary),flush=True)
