"""Export reconstructed state fields and complete ratio-curve equations."""
from sage.all import *
from pathlib import Path
import json,pickle,time
out=Path(__file__).absolute().parent; root=out.parent.parent; mirror=out/'rebuild'; st=time.monotonic();checks=[]
for r in ('round8','round9','round10'):(mirror/r).mkdir(parents=True,exist_ok=True)
def ck(n,v):
 if not v:raise AssertionError(n)
 checks.append(n)
fresh=load(str(out.parent/'rank3/fresh_joint_fields.sobj'))
R=fresh['ring'];K=R.base_ring();a,b,h,z=K.gens()
v=dict(fresh);v['normalized_pi_fields']=matrix(R,3,8,lambda i,j:0)
for i in range(3):
 for j,k in enumerate([0,1,2,4,5,6]):v['normalized_pi_fields'][i,k]=fresh['fields'][i,j]
 v['normalized_pi_fields'][i,3]=fresh['base_fields'][i,0];v['normalized_pi_fields'][i,7]=fresh['base_fields'][i,1]
ck('complete state field dimensions',fresh['fields'].dimensions()==(3,9) and fresh['base_fields'].dimensions()==(3,2))
save(v,str(mirror/'round10/joint_fields_direct.sobj'),compress=False)
save(fresh,str(mirror/'round8/rank3_joint_fiber_leading.sobj'),compress=False)
curve=load(str(out.parent/'rank3/fresh_curve_sources.sobj'));R=curve['ring'];K=R.base_ring();C=K.ring()
scaling=curve['row_scalings']
ck('seven nonzero row scalings',len(scaling)==7 and all(scaling))
rows=[f/c for f,c in zip(curve['rows'],scaling)]
primitive=[];scales=[]
for f in rows:
 den=lcm([c.denominator() for c in f.coefficients()]);nn={m:C(c*den) for m,c in f.dict().items()};con=gcd(list(nn.values()))
 primitive.append(R({m:K(n/con) for m,n in nn.items()}));scales.append(K(con/den))
 ck('full coefficient content identity '+str(len(primitive)-1),primitive[-1]*scales[-1]==f)
save(dict(ring=R,parameters=K.gens(),coefficient_ring=C,rows=rows,primitive_rows=primitive,scalings=scales),str(mirror/'round10/general_curve_dimensionless_full.sobj'),compress=False)
# Primitive bilinear/trilinear rows derived directly from their degree-2/3 terms.
blocks={};bsc={};cubics=[];csc=[]
def prim(cs):
 den=lcm([v.denominator() for v in cs]);nn=[C(v*den) for v in cs];co=gcd(nn);return [n//co for n in nn],co,den
for f,(i,j) in zip(rows[:3],[(0,1),(0,2),(1,2)]):
 cs=[f.monomial_coefficient(R.gen(2*i+ii)*R.gen(2*j+jj)) for ii in range(2) for jj in range(2)]
 nn,co,den=prim(cs);blocks[i,j]=matrix(C,2,2,nn);bsc[i,j]=(co,den)
for f in rows[3:]:
 cs=[f.monomial_coefficient(R.gen(i)*R.gen(2+j)*R.gen(4+k)) for i in range(2) for j in range(2) for k in range(2)]
 nn,co,den=prim(cs);cubics.append(nn);csc.append((co,den))
for f,ij in zip(rows[:3],[(0,1),(0,2),(1,2)]):
 co,den=bsc[ij]
 represented=sum((K(co/den)*blocks[ij][ii,jj]*R.gen(2*ij[0]+ii)*R.gen(2*ij[1]+jj) for ii in range(2) for jj in range(2)),R.zero())
 leading=R({tuple(e):c for e,c in f.dict().items() if sum(e)==2})
 ck('complete quadratic bilinear face identity '+str(ij),represented==leading)
 residual=f-leading
 ck('complete lower linear remainder '+str(ij),all(sum(e)==1 for e in residual.dict()))
for i,(f,coeffs) in enumerate(zip(rows[3:],cubics)):
 co,den=csc[i]
 represented=sum((K(co/den)*coeffs[4*ii+2*jj+kk]*R.gen(ii)*R.gen(2+jj)*R.gen(4+kk) for ii in range(2) for jj in range(2) for kk in range(2)),R.zero())
 leading=R({tuple(e):c for e,c in f.dict().items() if sum(e)==3})
 ck('complete trilinear coefficient identity '+str(i),represented==leading)
save(dict(coefficient_ring=C,parameters=K.gens(),blocks=blocks,block_scales=bsc,cubic_rows=cubics,cubic_scales=csc),str(mirror/'round10/general_curve_dimensionless_leading.sobj'),compress=False)
# The characteristic equations use the complete seventeen-state source.
f17=load(str(out.parent/'foundations/state17_fresh.sobj'));P=f17['ring']
d=dict(ring=P,F=f17['F'],C3=f17['C3'],Q4=f17['Q4'])
ck('complete characteristic source sizes',P.ngens()==17 and len(d['F'])==3 and all(len(row)==17 for row in d['F']) and len(d['C3'])==3 and len(d['Q4'])==21)
(mirror/'round8/rank3_source_minor_candidates.sobj').write_bytes(pickle.dumps(d))
report=dict(status='PASS',checks=checks,count=len(checks),seconds=time.monotonic()-st)
(out/'fresh_seed_report.json').write_text(json.dumps(report,indent=2));print(report,flush=True)
