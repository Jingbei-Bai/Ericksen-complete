"""Full endpoint-neighborhood quartic with the plane parameter kept free.

The moving cubic center c(p,q)*delta^2 is exact in the endpoint parameter.
This archive supplies protected units for a separate whole-rate argument.
"""
from sage.all import *
from pathlib import Path
import pickle,json,time,resource,sys
OUT=Path(__file__).absolute().parent;START=time.time();PZERO='--pzero' in sys.argv
STEM='homothety_middle_endpoint_uniform_quartic'+('_pzero' if PZERO else '')
vm=int(Path('/proc/self/statm').read_text().split()[0])*4096
resource.setrlimit(resource.RLIMIT_AS,(vm+900*1024**2,vm+900*1024**2))
sd=load(str(OUT/'homothety_axis_plane_forcing_schur.sobj'));dd=sd['source']
qd=load(str(OUT/'homothety_axis_plane_raw_quadratic_schur.sobj'))
K=FractionField(PolynomialRing(QQ,names=('a','r','A2','B2','k1' if PZERO else 'k2','v')))
a,r,A2,B2,kfree,v=K.gens();b=a*(1+r*r)/(a+r*r)
p,q=(v,K(1)) if PZERO else (K(1),v)
if PZERO:
 k1=kfree;k2=-1+(2*a*p*(b*B2-A2)-a*(b-1)*p*p*(k1+1))/(a-b)
else:
 k2=kfree;k1=-1+(2*a*q*(b*B2-A2)-(a-b)*q*q*(k2+1))/(a*(b-1))
Zc=-2*r*p*p*q*(1+r*r)/(p*r+q)
S=PolynomialRing(K,names=('E','t3','h3','z3'));E,t3,h3,z3=S.gens();nv=S.gens()
R=PolynomialRing(S,'delta');delta=R.gen()
def tr(f):return R([f[i] for i in range(min(4,f.degree())+1)]) if f else R(0)
def cf4(c):return K(c.numerator()(a,b,p,q))/K(c.denominator()(a,b,p,q))
def cforig(c):return K(c.numerator()(a,b))/K(c.denominator()(a,b))
def substitute(poly,vals):
 ans=R(0)
 for ex,c in poly.dict().items():
  term=R(cf4(c))
  for i,ej in enumerate(ex):
   if ej:term=tr(term*vals[i]**ej)
  ans+=term
 return tr(ans)
obs=list(map(R,[0,0,A2,0,0,B2,k1,k2]));sv=delta;tv=r*delta+t3*delta**3;ev=E*delta**3
frvals=obs+[h3*delta**3,Zc*delta**2+z3*delta**3,sv,tv,ev]
vals=[R(cf4(c))+substitute(w,frvals) for c,w in zip(dd['base_state'],sd['state_linear_graph'])]
cols=dd['combined_pivot_columns'];rows=dd['combined_pivot_rows']
for i,z in zip(cols,qd['second_state_pivot_graph']):vals[i]+=substitute(z,[delta,r*delta,R(0),R(0)])
M=dd['observation_plus_normal_matrix'].apply_map(cf4);BI=M.matrix_from_rows_and_columns(rows,cols).inverse()
src=pickle.loads((OUT.parent/'round3/generic_degree5_sparse_prepared.sobj').read_bytes());weights=[1]*9+[2]*6+[3]*2
def pull_factory(vs):
 powers=[[R(1),ww] for ww in vs]
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
 return [tr(ww-ev*oo0) for ww,oo0 in zip(oo,obs)]+H
pull=pull_factory(vals);low=low_sources(vals,pull)
state3=-BI*vector(S,[low[i][3] for i in rows])
for j,ww in zip(cols,state3):vals[j]=tr(vals[j]+ww*delta**3)
pull=pull_factory(vals);raw=low_sources(vals,pull)
for name,seq in [('C3',src['C3']),('Q4',src['Q4']),('G5',src['basis'][:27])]:
 raw.extend(map(pull,seq));print('uniform endpoint quartic',name,'seconds',time.time()-START,flush=True)
checks={f'raw_order_{i}_{kk}':f[kk]==0 for i,f in enumerate(raw) for kk in range(4)}
assert all(checks.values()),[kk for kk,vv in checks.items() if not vv]
quartic=vector(S,[f[4] for f in raw]);state4=-BI*vector(S,[quartic[i] for i in rows]);quartic+=M.matrix_from_columns(cols)*state4
forcing=sd['forcing_matrix'].apply_map(cf4)
fr,fc=([3,4,11,21,27],[0,1,3,4,7]) if PZERO else ([3,4,9,22,27],[0,1,3,4,6])
FP=forcing.matrix_from_rows_and_columns(fr,fc);comb=forcing.matrix_from_columns(fc)*FP.inverse()
checks['whole_forcing_elimination']=comb*forcing.matrix_from_rows(fr)==forcing
projected=quartic-comb*vector(S,[quartic[i] for i in fr])
MM=matrix(K,[[f.monomial_coefficient(ww) for ww in nv] for f in projected]);rhs=vector(K,[-f.constant_coefficient() for f in projected])
assert all(ff==sum(MM[i,j]*nv[j] for j in range(4))-rhs[i] for i,ff in enumerate(projected))
rr=([9,12,16,46] if PZERO else [11,12,16,47]);AA=MM.matrix_from_rows(rr);sol=AA.inverse()*vector(K,[rhs[i] for i in rr]);res=MM*sol-rhs
def at0(c):return c.subs({v:0})
summary={'status':'ENDPOINT_NEIGHBORHOOD_QUARTIC_SOURCE','scope':('p=v,q=1' if PZERO else 'p=1,q=v')+'; independent v retained before all source calculations',
 'checks':sum(checks.values()),'rows':rr,'forcing_minor':str(FP.det().factor()),'minor':str(AA.det().factor()),
 'endpoint_minor':str(at0(AA.det()).factor()),'cubic_center':str(Zc.factor()),
 'solution':{str(ww):str(zz.factor()) if zz else '0' for ww,zz in zip(nv,sol)},
 'nonzero_residuals':[(i,str(zz.factor())) for i,zz in enumerate(res) if zz],
 'endpoint_obstruction':str(at0(res[22 if PZERO else 21]).factor()),'seconds':float(time.time()-START)}
save({'ring':S,'field':K,'quartic_raw':quartic,'quartic_forcing_free':projected,'checks':checks,
 'state_jet':vals,'state3':state3,'state4':state4,'material_b':b,'matrix':MM,'rhs':rhs,
 'forcing':forcing,'forcing_projection':comb,'forcing_rows':fr,'forcing_columns':fc,
 'cubic_center':Zc,'solution':sol,'residuals':res,'summary':summary},str(OUT/(STEM+'.sobj')))
(OUT/(STEM+'.json')).write_text(json.dumps(summary,indent=2))
print(json.dumps({kk:summary[kk] for kk in ('status','scope','checks','endpoint_minor','endpoint_obstruction','seconds')}),flush=True)
