"""Exact all-37 second-jet compatibility on the true extra-tau wall.

This is an obstruction to a second-order no-lift proof, NOT an actual arc.
Setting eta and tangent first corrections to zero is a chosen jet witness,
not omission of their equations. All 37 source rows are replayed.
"""
from sage.all import *
from pathlib import Path
import time,os,resource,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_ratio_incoming12_generic_nullwall_all.sobj'));c=load(str(p/'independent_ratio_incoming12_generic_nullwall_all_cokernel.sobj'));U=d['ring'];K=U.base_ring();l=U.gen();mod=d['relation'];dd=-mod[0];a,b,h,z=K.gens();items={int(it['source_row']):it for it in d['items']};ci={int(it['source_row']):it for it in c['items']}
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+550*1024**2,)*2)
def mul(a,b):return (a*b)%mod
def norm(a):
 a=U(a)%mod;return K(a[0]*a[0]-dd*a[1]*a[1])
# Set k=l*y, tau=1; remaining eta,dl,du,dv all zero.
f=[]
for row in range(37):
 rr=ci[row]['rhs'];f.append([mul(rr[0],dd),mul(rr[2],l),rr[5]])
ref=f[29];checks=0;cross=[]
for row in range(37):
 res=[(mul(ref[0],f[row][j])-mul(f[row][0],ref[j]))%mod for j in range(3)]
 assert not any(res);checks+=3;cross.append(res)
# The reference quadratic has a common nonzero l parity factor only.
parity=0 if ref[0][0] else 1;assert all(x[1-parity]==0 for x in ref)
Y=PolynomialRing(K,'jet_y');y=Y.gen();refpoly=sum(ref[j][parity]*y**(2-j) for j in range(3))
den=lcm([cc.denominator() for cc in refpoly]);num=refpoly*den;content=gcd([K(cc).numerator() for cc in num]);primitive=Y([K(cc)/content for cc in num]);assert primitive.degree()==2
# Exact range solution numerator, without division by the range determinant.
A=matrix(U,[[items[r]['linear_row'][j] for j in (1,2)] for r in (1,17)]);det=c['two_range_determinant'];adj=A.adjugate();c5=c['odd_range_pivot'];delta=mul(det,c5)
V=PolynomialRing(U,'comp_y');yy=V.gen()
def red(f):return V({i:U(cc)%mod for i,cc in enumerate(f) if U(cc)%mod})
def qrow(row):
 rr=items[row]['rhs'];return V([rr[5],mul(l,rr[2]),mul(dd,rr[0])])
q1,q17,q5=qrow(1),qrow(17),qrow(5)
# M rows 1,17 have no col0; row5 has no cols1,2.
assert items[1]['linear_row'][0]==items[17]['linear_row'][0]==0
assert items[5]['linear_row'][1]==items[5]['linear_row'][2]==0
n2num=[red(-det*q5),red(-c5*(adj[0,0]*q1+adj[0,1]*q17)),red(-c5*(adj[1,0]*q1+adj[1,1]*q17)),V.zero(),V.zero(),V.zero()]
remainders=[]
for row in range(37):
 rem=red(delta*qrow(row)+sum(items[row]['linear_row'][j]*n2num[j] for j in range(6)))
 target=V([f[row][2],f[row][1],f[row][0]])
 assert rem==target;checks+=1;remainders.append(rem)
save({'ring':U,'relation':mod,'kernel_k_scaling':'k=l*y, eta=0, tau=1, dl=du=dv=0','reference_row':29,'reference_coefficients_descending':ref,'reference_polynomial':refpoly,'primitive_quadratic':primitive,'primitive_scaling':K(den)/content,'range_denominator':delta,'range_denominator_norm':norm(delta),'range_second_numerators':n2num,'range_second_remainders':remainders,'cross_identities':cross,'reference_leading_norm':norm(ref[0]),'scope':__doc__},str(p/'independent_ratio_incoming12_nullwall_second_compatibility.sobj'),compress=False)
(p/'independent_ratio_incoming12_nullwall_second_compatibility.json').write_text(json.dumps({'checks':checks,'all_source_rows':37,'quadratic_degree':2,'primitive_coefficient_termcounts':[len(K(cc).numerator().dict()) for cc in primitive],'seconds':time.monotonic()-st,'scope':__doc__,'module_diagnostic':'Prior echelon diagnostic reached memory cap; no result from its unfinished odd block is used.'},indent=2)+'\n');print('DONE',checks,time.monotonic()-st,flush=True)
