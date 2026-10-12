"""All-parameter unit-gradient projective chart and six amplitude rows.

This chart requires all three direction m_i != 0. Spectrum independence
and N elimination are exact, never inferred from a specialization.
"""
from pathlib import Path
import json,time
from sage.all import *
out=Path(__file__).resolve().parent;start=time.monotonic();new=load(str(out/'general_curve_extra_quartic_cubics.sobj'));full=load(str(out/'general_curve_dimensionless_full.sobj'));K=new['ring'].base_ring();a,r,H,Z,S,N=K.gens();q=1+r
R=PolynomialRing(K,names=('t1','t2','t3'));ts=R.gens();images=[ts[0],R.one(),ts[1],R.one(),ts[2],R.one()]
oldguard=a*r*(1+r)*(1+a*r)*H*Z*(H-Z)*(Z-1)*(H+1)*(S-1)*(H*S-Z+1)*(H*S-Z+S)*(1+H-Z)
oldfactors=[f for f,e in K.ring()(oldguard).factor()]
def ev(f):return sum((c*prod(images[j]**int(e) for j,e in enumerate(mon) if e) for mon,c in f.dict().items()),R.zero())
mixed=[]
for f in full['rows'][:3]:
 top=new['ring']({tuple(m):c for m,c in f.dict().items() if sum(m)==2});p=ev(top)*r/q
 assert all(c.derivative(a)==0 and c.derivative(r)==0 for c in p.coefficients())
 dn=R({tuple(m):c.derivative(N) for m,c in p.dict().items()});assert dn==-1
 mixed.append(p)
nvalue=mixed[0]+N;differences=[p-mixed[0] for p in mixed[1:]];assert all(c.derivative(N)==0 for p in [nvalue]+differences for c in p.coefficients())
oldcubics=[];oldcubic_scales=[]
for f in full['rows'][3:]:
 top=new['ring']({tuple(m):c for m,c in f.dict().items() if sum(m)==3});p=ev(top);cs=p.coefficients();den=lcm([c.denominator() for c in cs]);nums=[K.ring()(c*den) for c in cs];con=gcd(nums);p=p*den/con
 for pol in (den,con):
  for fac,e in pol.factor():assert any(not fac.quo_rem(v)[1] for v in oldfactors),str(fac)
 oldcubic_scales.append(K(con/den))
 assert all(c.derivative(a)==0 and c.derivative(r)==0 and c.derivative(N)==0 for c in p.coefficients());oldcubics.append(p)
columns=[1/(a*r),K.one(),1/(1+a*r)];row_scales=[a*r,a*r,K.one(),K.one(),q*(1+a*r)/r,q*(1+a*r)/r];odd=[2,1,2,0,1,0];matrix_rows=[];checks={};records=[]
for i,f in enumerate(new['leading_rows']):
 cols=[R.zero() for j in range(3)]
 for mon,c in f.dict().items():
  deg=[sum(mon[2*j:2*j+2]) for j in range(3)];deg[odd[i]]-=1;j=deg.index(2)
  co=c/(row_scales[i]*columns[j]);assert co.derivative(a)==0 and co.derivative(r)==0
  cols[j]+=co*prod(images[k]**int(e) for k,e in enumerate(mon) if e)
 assert all(c.derivative(N,2)==0 for pp in cols for c in pp.coefficients())
 cols=[sum((c.subs({N:0})*R.monomial(*mon)+c.derivative(N)*nvalue*R.monomial(*mon) for mon,c in pp.dict().items()),R.zero()) for pp in cols]
 assert all(c.derivative(N)==0 for pp in cols for c in pp.coefficients())
 matrix_rows.append(cols);checks['all_parameter_spectral_scaling_'+str(i)]=True
 rec={'Q_index':new['row_specs'][i][0],'amplitude_entry_terms':[pp.number_of_terms() for pp in cols],'amplitude_entry_degrees':[int(pp.total_degree()) for pp in cols]};records.append(rec);print(rec,flush=True)
save({'ring':R,'parameters':K.gens(),'mixed_rows':mixed,'mixed_differences':differences,'curvature_value':nvalue,'old_cubic_rows':oldcubics,'old_cubic_scales':oldcubic_scales,'amplitude_matrix':matrix(R,matrix_rows),'spectral_column_weights':columns,'spectral_row_weights':row_scales,'checks':checks,'records':records,'source':'general_curve_extra_quartic_cubics.sobj and general_curve_dimensionless_full.sobj','scope':'Complete unit-gradient direction chart (all m_i nonzero), N solved without any new divisor; exact all-parameter spectral row/column scaling. No rank or dimension conclusion.'},str(out/'general_curve_unit_gradient_directions.sobj'))
report={'status':'EXACT_SMALL_DIRECTION_SYSTEM','mixed_difference_terms':[p.number_of_terms() for p in differences],'old_cubic_terms':[p.number_of_terms() for p in oldcubics],'checks':checks,'records':records,'elapsed_seconds':time.monotonic()-start};(out/'general_curve_unit_gradient_directions.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
