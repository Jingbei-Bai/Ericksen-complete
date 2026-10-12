"""Benchmark exact rational evaluations for a degree-bound identity proof."""
from sage.all import *
from pathlib import Path
from scipy.optimize import linear_sum_assignment
import numpy as np,json,time
out=Path(__file__).resolve().parent;start=time.monotonic();d=load(str(out/'general_curve_ranktwo_base_macaulay_input.sobj'));p=load(str(out/'general_curve_ranktwo_base_flint_preconditioned.sobj'));P=p['integer_matrix_rows'][0][0].parent();A=d['matrix'];b=d['target_vector'];F=A.base_ring();Hp=F.ring();selected=d['selected_rows'];B=matrix(P,p['integer_matrix_rows']);c=vector(P,p['integer_rhs']);n=B.ncols();cols=p['column_contents'];assert all(f.degree()==0 for f in cols)
allrows=[];allrhs=[];contents=[];qscales=[]
for i in range(A.nrows()):
    vals=[Hp(A[i,j]/F(cols[j])) for j in range(n)]+[Hp(b[i])];q=lcm([v.denominator() for f in vals for v in f.list()] or [1]);iv=[P(q*f) for f in vals];g=gcd(iv) or P.one();iv=[v//g for v in iv];allrows.append(iv[:-1]);allrhs.append(iv[-1]);contents.append(g);qscales.append(q)
aug=matrix(P,[list(B.row(i))+[c[i]] for i in range(n)]);weights=np.array([[int(f.degree()) if f else -100000 for f in row] for row in aug.rows()],dtype=np.int64);bounds=[]
for r,t in zip(allrows,allrhs):
    w=np.vstack([weights,np.array([int(f.degree()) if f else -100000 for f in r+[t]],dtype=np.int64)]);ri,ci=linear_sum_assignment(w,maximize=True);v=int(w[ri,ci].sum());bounds.append(v)
D=max(bounds);print('DEGREE_BOUND',D,'row_sum',sum(max(int(f.degree()) if f else 0 for f in row) for row in aug.rows())+max(max(int(f.degree()) if f else 0 for f in row+[t]) for row,t in zip(allrows,allrhs)),flush=True)
payload={'integer_ring':P,'square_matrix':B,'square_rhs':c,'all_rows':matrix(P,allrows),'all_rhs':vector(P,allrhs),'all_row_contents':contents,'all_row_rational_scales':qscales,'column_contents':cols,'degree_bound':D,'individual_degree_bounds':bounds,'input_file':'general_curve_ranktwo_base_macaulay_input.sobj','scope':'A full proof requires D+1 distinct exact-Q evaluations with nonsingular square matrix and all rows replayed. This benchmark is not that proof.'};save(payload,str(out/'general_curve_ranktwo_interpolation_input.sobj'));records=[]
for h in [2,3,5]:
    t0=time.monotonic();M=matrix(ZZ,[[f(h) for f in row] for row in B.rows()]);rhs=vector(ZZ,[f(h) for f in c]);x=M.solve_right(rhs);t1=time.monotonic();AM=matrix(ZZ,[[f(h) for f in row] for row in allrows]);ab=vector(ZZ,[f(h) for f in allrhs]);assert AM*x==ab;rec={'H':h,'solve_seconds':t1-t0,'all_replay_seconds':time.monotonic()-t1,'max_denominator_bits':max(int(v.denominator().nbits()) for v in x)};records.append(rec);print(rec,flush=True)
(out/'general_curve_ranktwo_interpolation_benchmark.json').write_text(json.dumps({'degree_bound':D,'records':records,'seconds':time.monotonic()-start,'scope':payload['scope']},indent=2)+'\n')
