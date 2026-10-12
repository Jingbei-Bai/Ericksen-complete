"""Checkable integer row/column degree potentials, independent of optimizer."""
from sage.all import *
from pathlib import Path
from scipy.optimize import linear_sum_assignment
import numpy as np,json,time
out=Path(__file__).resolve().parent;start=time.monotonic();d=load(str(out/'general_curve_ranktwo_interpolation_input.sobj'));B=d['square_matrix'];c=d['square_rhs'];A=d['all_rows'];b=d['all_rhs'];n=B.ncols();low=-100000
w=np.array([[int(f.degree()) if f else low for f in list(B.row(i))+[c[i]]] for i in range(n)],dtype=np.int64)
last=np.array([max([int(f.degree()) for f in A.column(j) if f] or [low]) for j in range(n)]+[max([int(f.degree()) for f in b if f] or [low])],dtype=np.int64);W=np.vstack([w,last]);ri,ci=linear_sum_assignment(W,maximize=True);assert list(ri)==list(range(n+1));v=np.zeros(n+1,dtype=np.int64)
for step in range(n+1):
    nv=np.maximum(v,np.max(v[ci,None]+W-W[np.arange(n+1),ci,None],axis=0))
    if np.array_equal(nv,v):break
    v=nv
u=np.max(W-v[None,:],axis=1);D=int(u.sum()+v.sum());u=[int(q) for q in u];v=[int(q) for q in v]
checks=[]
for i in range(n):
    for j,f in enumerate(list(B.row(i))+[c[i]]):
        if f:assert int(f.degree())<=u[i]+v[j]
for i in range(A.nrows()):
    for j,f in enumerate(list(A.row(i))+[b[i]]):
        if f:assert int(f.degree())<=u[-1]+v[j]
assert D==sum(u)+sum(v)
report={'selected_row_potentials':u[:-1],'remaining_row_potential':u[-1],'column_potentials':v,'degree_bound':D,'all_integer_degree_inequalities_pass':True,'seconds':time.monotonic()-start,'scope':'For each augmented determinant, every permutation product has degree <= sum(row potentials)+sum(column potentials). The inequality checks, not optimality of a numerical assignment, certify the bound.'}
(out/'general_curve_ranktwo_interpolation_degree_certificate.json').write_text(json.dumps(report,indent=2)+'\n');print('CERTIFIED_DEGREE_BOUND',D,time.monotonic()-start,flush=True)
