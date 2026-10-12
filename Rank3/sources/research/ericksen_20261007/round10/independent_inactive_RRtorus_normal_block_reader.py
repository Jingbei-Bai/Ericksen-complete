"""Independent matching of complete-source transported normals to the short inactive block certificates."""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();n=load(str(p/'independent_inactive_RRtorus_infinity_normal.sobj'));b=load(str(p/'independent_inactive_RRtorus_normal_blocks.sobj'));d=load(str(p/'independent_inactive_RRtorus_incoming_boundary.sobj'));Q=b['ring'];l,rho=Q.gens();K=Q.base_ring();a,r,H,Z=K.gens();QF=Q.fraction_field();old=n['coefficient_field'];phi=old.ring().hom([a,a*(1+r)/(1+a*r),H/r,(1+r)*Z/r],K);cache={};checks=[]
def ck(s,z):assert z,s;checks.append(s)
def coeff(c):
 c=old(c)
 if c not in cache:cache[c]=phi(c.numerator())/phi(c.denominator())
 return cache[c]
def pol(f):return Q({tuple(map(int,e)):coeff(co) for e,co in f.dict().items()})
def cv(f):return QF(pol(f.numerator()))/pol(f.denominator())
Nold=n['boundary_matrix'];rows=[0,1,2,5,7,11,23];N=matrix(QF,7,6,[[cv(Nold[i,j]) for j in range(6)] for i in rows]);items=[]
# Every selected original row is independently matched on the block needed in the actual 6-minor.
for pos,(ri,cols,expected) in enumerate([(1,[0,1],b['even_transverse_matrix'].row(0)),(2,[0,1],b['even_transverse_matrix'].row(1)),(0,[2],[b['short_certificates']['Q0_e0']]),(5,[3,4,5],b['range_delta_block'].row(1)),(7,[3,4,5],b['range_delta_block'].row(3)),(11,[3,4,5],b['range_delta_block'].row(8)),(23,[3,4,5],b['range_delta_block'].row(20))]):
 got=[N[rows.index(ri),j] for j in cols];expected=list(map(QF,expected));j=next(j for j,f in enumerate(expected) if f);unit=got[j]/expected[j]
 ck('row_unit_nonzero',unit!=0);ck('row_unit_only_base',unit.numerator().total_degree()==0 and unit.denominator().total_degree()==0)
 for x,y in zip(got,expected):ck('original_block_reconstruction',x==unit*y)
 items.append({'original_row':ri,'columns':cols,'unit':unit})
for ri in [1,2]:
 for j in [2,3,4,5]:ck('even_other_normal_columns_zero',N[rows.index(ri),j]==0)
for ri in [5,7,11,23]:ck('odd_E0_zero',N[rows.index(ri),2]==0)
q0=N[0,2];ck('Q0_E0_oldunit',q0!=0)
minors=[]
for hi in [11,23]:
 sel=[1,2,0,5,7,hi];sm=matrix(QF,[N.row(rows.index(i)) for i in sel]);evdet=sm.matrix_from_rows_and_columns([0,1],[0,1]).det();hdet=sm.matrix_from_rows_and_columns([3,4,5],[3,4,5]).det();product=evdet*sm[2,2]*hdet
 # Block-triangular identity is checked by a sparse permutation expansion, avoiding large generic elimination.
 import itertools
 exact=QF.zero()
 for perm in itertools.permutations(range(6)):
  if any(sm[i,perm[i]]==0 for i in range(6)):continue
  inv=sum(perm[i]>perm[j] for i in range(6) for j in range(i+1,6));exact+=(-1)**inv*prod(sm[i,perm[i]] for i in range(6))
 ck('actual_six_minor_factorization',exact==product);minors.append({'source_rows':sel,'minor':exact,'block_product':product})
save({'source_files':['independent_inactive_RRtorus_infinity_normal.sobj','independent_inactive_RRtorus_normal_blocks.sobj'],'ring':Q,'selected_original_rows':rows,'matrix':N,'row_units':items,'six_minors':minors,'checks':checks,'scope':__doc__},str(p/'independent_inactive_RRtorus_normal_block_reader.sobj'),compress=False)
(p/'independent_inactive_RRtorus_normal_block_reader.json').write_text(json.dumps({'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n');print('PASS',len(checks),'units',[(i['original_row'],str(K(i['unit']).factor())) for i in items],time.monotonic()-st,flush=True)
