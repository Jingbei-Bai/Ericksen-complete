from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;start=time.monotonic();checks={}
def ck(n,v):
 if not v:raise AssertionError(n)
 checks[n]=True
sec=load(str(p/'independent_ratio_RR_section_plane_blocks.sobj'));cy=load(str(p/'independent_ratio_RR_cyclic_full_matrix.sobj'));pt=load(str(p/'independent_ratio_RR_cyclic_residual_point.sobj'));K=cy['coefficient_field'];a,r=cy['parameters'];q=1+r;guard=cy['guard'];T=cy['RR_ring'];rr=T.gens();ab=(K(a),K(a*q/(1+a*r)))
def cv(f):return T({tuple(map(int,e)):K(c.numerator()(*ab))/K(c.denominator()(*ab)) for e,c in f.dict().items()})
M=cy['matrix']
for i in range(28):
 s=cy['row_scales'][i];ck(f'rowunit_{i}',s!=0 and all(guard.quo_rem(f)[1]==0 for v in [s.numerator(),s.denominator()] for f,e in v.factor()))
 for j in range(7):ck(f'original_source_{i}_{j}',cv(sec['normal_matrix'][i,j])==s*M[i,j])
scales=[rr[2],rr[1],rr[2],rr[0],rr[1],rr[0],rr[0]*rr[1]*rr[2]];xvals=(-q/r,-q,K.one());ck('canonical_rr_squares',tuple(pt['RR_squares'])==xvals);ck('column_scales',list(pt['column_scalings'])==scales)
rows=[]
for i,row in enumerate(M.rows()):
 vals=[f*s for f,s in zip(row,scales)];allmons=[tuple(map(int,e)) for f in vals for e in f.dict()];parity=tuple(v%2 for v in allmons[0]);ck(f'parity_{i}',parity==tuple(pt['row_parities'][i]) and all(tuple(v%2 for v in e)==parity for e in allmons))
 rownew=[]
 for f in vals:
  out=K.zero()
  for e,c in f.dict().items():
   e=tuple(map(int,e));ex=tuple((ei-pi)//2 for ei,pi in zip(e,parity));out+=c*prod(xv**ee for xv,ee in zip(xvals,ex))
  rownew.append(out)
 rows.append(rownew)
N=matrix(K,rows);ck('whole_point_matrix',N==pt['matrix']);expected=[vector(K,[1,0,0,0,0,1,0]),vector(K,[0,1,0,1,0,0,0]),vector(K,[0,0,1,0,1,0,0])]
for i,v in enumerate(expected):ck(f'kernel_{i}',N*v==0)
ck('kernel_independent',matrix(K,expected).rank()==3)
prs=N.transpose().pivots();pcs=N.pivots();ck('rank4',len(prs)==len(pcs)==4);minor=N.matrix_from_rows_and_columns(prs,pcs).det();bad=[f for v in [minor.numerator(),minor.denominator()] for f,e in v.factor() if guard.quo_rem(f)[1]]
ck('minor_protected_all_spectra',not bad)
report={'status':'PASS','count':len(checks),'checks':checks,'pivot_rows':list(prs),'pivot_columns':list(pcs),'protected_minor':str(minor.factor()),'seconds':time.monotonic()-start,'scope':'Original28 sources at cyclic base, complete rational spectral map and protected row units, all parity/column transformations replayed. Every one of four complex RR-square preimages has rank EXACTLY4 for every allowed spectrum, with three-dimensional tau-zero kernel. This is a first-normal obstruction only, not an affine lift or actual solution.'}
(p/'independent_RR_canonical_rank4_reader.json').write_text(json.dumps(report,indent=2)+'\n');save({'checks':checks,'minor':minor,'pivot_rows':prs,'pivot_columns':pcs,'kernel':expected},str(p/'independent_RR_canonical_rank4_reader.sobj'));print({k:v for k,v in report.items() if k!='checks'},flush=True)
