"""Six actual quartic sources reduced with full diagonal curve tangencies.

No parameter specialization. The diagonal pivots and k_i are old guards.
Both source provenance and all degree cancellations are replayed.
"""
from pathlib import Path
from functools import lru_cache
import json,time
from sage.all import *
out=Path(__file__).resolve().parent;start=time.monotonic()
d=load(str(out/'rank3_quartic_diagonal_velocity_form.sobj'))
g=load(str(out/'rank1_general_slope_quartics.sobj'))
full=load(str(out/'general_curve_dimensionless_full.sobj'))
R=full['ring'];K=R.base_ring();CR=K.ring();a,r,H,Z,S,N=K.gens();q=1+r;b=a*q/(1+a*r)
w1,m1,w2,m2,w3,m3=R.gens()
pi=[w3,w2,((1+H)*S-Z)*m3/(2*Z)+H*w3,w1,-(1+H*S-Z)*m2/(2*H)-Z*w2,-(1-S)*m1/(2*H)-Z*w1/H]
rr=[(Z-H*S)*m1/(H*Z),S*m2/Z,-m3/H];images=pi+rr
args=(a,b,H/r,q*Z/r);args5=args+(q*S,)
@lru_cache(None)
def cc(c):return K(c.numerator()(*args))/K(c.denominator()(*args))
def cc5(c):return K(c.numerator()(*args5))/K(c.denominator()(*args5))
def ev(poly,img=images,conv=cc):
 ans=R.zero()
 for mon,c in poly.dict().items():ans+=conv(c)*prod(img[j]**int(e) for j,e in enumerate(mon) if e)
 return ans
checks={}
for i,f in enumerate(g['ratio_graph']):checks['ratio_graph_'+str(i)]=ev(f,pi,cc5)==rr[i]
assert all(checks.values())
piv=[cc5(v) for v in g['pivots']];ks=[cc(v) for v in d['protected_coefficients']]
Wi=[];Ti=[]
for i in range(3):
 vel=ev(d['diagonal_velocities'][i]);Wi.append(vel)
 ti=full['rows'][4+i]-piv[i]*vel;Ti.append(ti)
 checks['diagonal_remainder_degree_'+str(i)]=ti.total_degree()<=2
 assert checks['diagonal_remainder_degree_'+str(i)]
 print('diagonal',i,'degree',ti.total_degree(),'terms',ti.number_of_terms(),'seconds',time.monotonic()-start,flush=True)
rows=[];leading=[];records=[]
for pos,(qi,i,xi,sgn) in enumerate(d['row_specs']):
 rem=ev(d['remainders'][pos]);f=rem-sgn*pi[xi]*Ti[i]/(ks[i]*piv[i]);rows.append(f)
 top=R({tuple(mon):c for mon,c in f.dict().items() if sum(mon)==3});leading.append(top)
 checks['degree_'+str(qi)]=f.total_degree()<=3
 checks['source_identity_'+str(qi)]=f+sgn*pi[xi]*full['rows'][4+i]/(ks[i]*piv[i])==rem+sgn*pi[xi]*Wi[i]/ks[i]
 assert checks['degree_'+str(qi)] and checks['source_identity_'+str(qi)]
 hist={}
 for mon in f.dict():hist[sum(mon)]=hist.get(sum(mon),0)+1
 rec={'Q_index':qi,'histogram':hist,'top_monomials':[list(map(int,mm)) for mm in top.dict()],'seconds':time.monotonic()-start};records.append(rec)
 print('Q',qi,'hist',hist,'seconds',time.monotonic()-start,flush=True)
 save({'ring':R,'parameters':K.gens(),'images':images,'diagonal_pivots':piv,'k_coefficients':ks,'diagonal_velocities':Wi,'diagonal_lower_remainders':Ti,'row_specs':d['row_specs'],'rows':rows,'leading_rows':leading,'checks':checks,'records':records,'source':'rank3_quartic_diagonal_velocity_form.sobj + general_curve_dimensionless_full.sobj','scope':'Six actual necessary cubic sources. Exact linear combinations of full quartics and diagonal curve tangencies; only old protected k_i/pivot_i inverted.'},str(out/'general_curve_extra_quartic_cubics.sobj'))
report={'status':'EXACT_SIX_EXTRA_CUBIC_SOURCES','checks':checks,'records':records,'elapsed_seconds':time.monotonic()-start};(out/'general_curve_extra_quartic_cubics.json').write_text(json.dumps(report,indent=2)+'\n')
