"""Middle-axis null cone: full cubic source including epsilon=E delta^2.

This is a next-face calculation, not a no-lift theorem. Both real null
directions are retained by real r!=0; b=a(1+r^2)/(a+r^2). The initial chart
p*q!=0 is explicit and its endpoint limits are not declared covered.
"""
from sage.all import *
from pathlib import Path
import pickle,json,time,resource
OUT=Path(__file__).absolute().parent;START=time.time()
vm=int(Path('/proc/self/statm').read_text().split()[0])*4096
resource.setrlimit(resource.RLIMIT_AS,(vm+1400*1024**2,vm+1400*1024**2))
d=load(str(OUT/'homothety_axis_plane_forcing_schur.sobj'));dd=d['source']
qd=load(str(OUT/'homothety_axis_plane_raw_quadratic_schur.sobj'))
K=FractionField(PolynomialRing(QQ,names=('a','p','q','r','A2','B2','k1')));a,p,q,r,A2,B2,k1=K.gens();b=a*(1+r*r)/(a+r*r)
k2=-1+(2*a*p*q*(b*B2-A2)-a*(b-1)*p*p*(k1+1))/((a-b)*q*q)
S=PolynomialRing(K,names=('E','t2','h2','z2'));E,t2,h2,z2=S.gens()
R=PolynomialRing(S,'delta');delta=R.gen()
def tr(f):return R([f[i] for i in range(min(3,f.degree())+1)]) if f else R(0)
def cf4(c):return K(c.numerator()(a,b,p,q))/K(c.denominator()(a,b,p,q))
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
 raw.extend(map(pull,seq));print('pulled face',name,'seconds',time.time()-START,flush=True)
checks={f'full_source_order_{j}_{k}':f[k]==0 for j,f in enumerate(raw) for k in range(3)}
assert all(checks.values()),[k for k,v in checks.items() if not v]
M=dd['observation_plus_normal_matrix'].apply_map(cf4);block=M.matrix_from_rows_and_columns(rows,cols)
cubic=vector(S,[f[3] for f in raw]);state3=-block.inverse()*vector(S,[cubic[i] for i in rows]);cubic+=M.matrix_from_columns(cols)*state3
forcing=d['forcing_matrix'].apply_map(cf4);frows=[3,4,9,11,27];fcols=[0,1,2,3,4]
comb=forcing.matrix_from_columns(fcols)*forcing.matrix_from_rows_and_columns(frows,fcols).inverse()
assert comb*forcing.matrix_from_rows(frows)==forcing
projected=cubic-comb*vector(S,[cubic[i] for i in frows])
save({'ring':S,'parameter_field':K,'material_b':b,'pressure_k2':k2,'raw_cubic':cubic,
      'forcing_free_cubic':projected,'base_state':vals,'source_checks':checks},str(OUT/'homothety_middle_null_cubic_face.sobj'))
terms=[(i,str(z.factor()) if z else '0') for i,z in enumerate(projected) if z]
summary={'status':'NEXT_FACE_ONLY','scope':'p*q*r!=0; full real null signs; arbitrary compatible A2,B2,k1,E',
 'source_checks':sum(checks.values()),'nonzero_forcing_free_cubic':terms,'seconds':float(time.time()-START)}
(OUT/'homothety_middle_null_cubic_face.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary),flush=True)
