"""Exact integer matrix/degree preparation for the complete zero-D layer."""
from sage.all import *
from pathlib import Path
from scipy.optimize import linear_sum_assignment
import numpy as np,json,time
out=Path(__file__).resolve().parent;start=time.monotonic();d=load(str(out/'general_curve_zero_det_B_macaulay_input.sobj'));A=d['matrix'];b=d['target_vector'];selected=d['selected_rows'];F=A.base_ring();Hp=F.ring();P=PolynomialRing(ZZ,'H');n=A.ncols();allrows=[];allrhs=[];contents=[];scales=[]
for i in range(A.nrows()):
    vals=[Hp(v) for v in A.row(i)]+[Hp(b[i])];q=lcm([c.denominator() for v in vals for c in v.list()] or [1]);iv=[P(q*v) for v in vals];g=gcd(iv) or P.one();iv=[v//g for v in iv];allrows.append(iv[:-1]);allrhs.append(iv[-1]);contents.append(g);scales.append(q)
AM=matrix(P,allrows);ab=vector(P,allrhs);B=AM.matrix_from_rows(selected);c=vector(P,[ab[i] for i in selected]);assert B.is_square()
assert all(scales[i]*A[i,j]==F(contents[i])*F(AM[i,j]) for i in range(A.nrows()) for j in range(n));assert all(scales[i]*b[i]==F(contents[i])*F(ab[i]) for i in range(A.nrows()))
low=-100000;w=np.array([[int(v.degree()) if v else low for v in list(B.row(i))+[c[i]]] for i in range(n)],dtype=np.int64);last=np.array([max([int(v.degree()) for v in AM.column(j) if v] or [low]) for j in range(n)]+[max([int(v.degree()) for v in ab if v] or [low])],dtype=np.int64);W=np.vstack([w,last]);ri,ci=linear_sum_assignment(W,maximize=True);v=np.zeros(n+1,dtype=np.int64)
for step in range(n+1):
    nv=np.maximum(v,np.max(v[ci,None]+W-W[np.arange(n+1),ci,None],axis=0))
    if np.array_equal(nv,v):break
    v=nv
u=np.max(W-v[None,:],axis=1);D=int(u.sum()+v.sum());u=[int(a) for a in u];v=[int(a) for a in v]
for i in range(n):
    for j,f in enumerate(list(B.row(i))+[c[i]]):assert not f or f.degree()<=u[i]+v[j]
for i in range(AM.nrows()):
    for j,f in enumerate(list(AM.row(i))+[ab[i]]):assert not f or f.degree()<=u[-1]+v[j]
payload={'integer_ring':P,'square_matrix':B,'square_rhs':c,'all_rows':AM,'all_rhs':ab,'all_row_contents':contents,'all_row_rational_scales':scales,'column_contents':[P.one()]*n,'input_file':'general_curve_zero_det_B_macaulay_input.sobj','preconditioning_identity':True,'scope':'Exact source-preserving integer normalization. A D+1 distinct exact-node certificate is still required.'};save(payload,str(out/'general_curve_zero_det_B_interpolation_input.sobj'))
degree={'selected_row_potentials':u[:-1],'remaining_row_potential':u[-1],'column_potentials':v,'degree_bound':D,'all_integer_degree_inequalities_pass':True,'seconds':time.monotonic()-start};(out/'general_curve_zero_det_B_interpolation_degree_certificate.json').write_text(json.dumps(degree,indent=2)+'\n')
# Reuse the already successful deterministic engine with disjoint artifacts.
source=(out/'general_curve_ranktwo_integer_interpolation.py').read_text().replace('general_curve_ranktwo','general_curve_zero_det_B').replace('B.nrows()==B.ncols()==185 and A.nrows()==276','B.nrows()==B.ncols()==163 and A.nrows()==210').replace("'all_rows':276,'square_size':185","'all_rows':210,'square_size':163")
(out/'general_curve_zero_det_B_integer_interpolation.py').write_text(source)
print({'square_size':n,'all_rows':AM.nrows(),'degree_bound':D,'required_nodes':D+1,'seconds':time.monotonic()-start},flush=True)
