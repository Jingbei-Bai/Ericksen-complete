"""Fraction-free sparse elimination in the p23,tau normal block."""
from sage.all import *
from pathlib import Path
import time
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_ratio_p12_axis_generic_normal.sobj'));M=d['matrix'];U=d['ring'];K=d['coefficient_field'];u,v=U.gens()
R=vector(U,[M[7,2],M[7,6]]);cu=R[0].monomial_coefficient(u);assert R[0]==cu*u and cu
cs=M[19,2].monomial_coefficient(u**2);A=vector(U,[M[19,j] for j in (2,6)])-(cs/cu)*u*R
ac=K(A[0]);assert ac
rels=[]
for i in (28,29,30):
 ci=M[i,2].monomial_coefficient(u**2);B=vector(U,[M[i,j] for j in (2,6)])-(ci/cu)*u*R;bc=K(B[0]);F=ac*B[1]-bc*A[1]
 assert not (ac*B[0]-bc*A[0]);assert all(e[0]>=1 for e in F.dict())
 Q=U({(int(e[0])-1,int(e[1])):c for e,c in F.dict().items()});assert F==u*Q
 rels.append({'source_row':i,'B':B,'B_coefficient':bc,'tau_polynomial':F,'after_u':Q,'R_multiplier':ci/cu})
 print('E',i,'after_u support',list(Q.dict()),flush=True)
G=R[0]*A[1]-A[0]*R[1]
cuH=U({(int(e[0])-1,int(e[1])):c for e,c in M[25,2].dict().items()});assert M[25,2]==cuH*u
J=vector(U,[M[25,j] for j in (2,6)])-(cuH/cu)*R;assert J[0]==0
print('R/A puretau support',list(G.dict()),'H2 puretau',list(J[1].dict()),flush=True)
save({'ring':U,'coefficient_field':K,'source_file':'independent_ratio_p12_axis_generic_normal.sobj','R':R,'cu':cu,'A':A,'A_coefficient':ac,'A_R_multiplier':cs/cu,'E_relations':rels,'RA_tau':G,'H2_tau':J[1],'H2_R_multiplier':cuH/cu,'scope':__doc__,'qualification':'Only cu and ac are base factors; their nonvanishing for each material must be audited before use. No u or v has been inverted.'},str(p/'independent_ratio_p12_axis_tau_elimination.sobj'),compress=False)
print('DONE',time.monotonic()-st,flush=True)
