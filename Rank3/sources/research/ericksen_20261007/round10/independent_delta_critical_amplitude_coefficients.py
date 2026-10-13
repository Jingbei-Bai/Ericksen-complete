"""Exact H/Z coefficient consequences of all critical affine minors off lambda=0.

Only the explicitly selected lambda branch and nonzero B,C amplitudes are
divided; source factors in the remaining directions/materials are kept.
"""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_delta_critical_consistency_minors.sobj'));old=d['field'];a0,r0,H0,Z0,p21,p23,p31,p32=old.gens();lam=p21*p23+p21*p32+p23
KK=PolynomialRing(QQ,names=('xa','xr','xs','xz')).fraction_field();a,r,s,z=KK.gens();P=PolynomialRing(KK,names=('xB','xC'));B,C=P.gens();ph=old.ring().hom([P(a),P(r),P.one(),P.one(),P(s),C*z,B,C],P)
rows=[];labels=[];items=[];checks=[]
for it in d['items']:
 f=it['determinant']/lam;n=f.numerator();groups={}
 for ee,c in n.dict().items():
  e=tuple(map(int,ee));key=e[2:4];g=list(e);g[2]=g[3]=0;groups.setdefault(key,{})[tuple(g)]=c
 gg=[(key,ph(n.parent()(dd))) for key,dd in sorted(groups.items())];mins=tuple(min(int(e[i]) for key,g in gg for e in g.dict()) for i in range(2));unit=B**mins[0]*C**mins[1]
 coeffs=[(key,g//unit) for key,g in gg];mons=sorted(set(tuple(map(int,e)) for key,g in coeffs for e in g.dict()));assert set(mons)<=set(((0,0),(0,2),(2,0))),mons
 for key,g in coeffs:
  rr=[g.monomial_coefficient(w) for w in (P.one(),C*C,B*B)];assert g==rr[0]+rr[1]*C*C+rr[2]*B*B;rows.append(rr);labels.append((it['H_row'],key));checks.append('source coefficient '+str(labels[-1]))
 items.append({'H_row':it['H_row'],'divided_source':f,'removed_amplitude_unit':unit,'base_coefficients':coeffs})
M=matrix(KK,rows);print('MATRIX',M.nrows(),3,'seconds',time.monotonic()-st,flush=True)
# A small exact row basis is structural only. All unchosen rows remain saved.
basis=[];piv=[]
for i in range(M.nrows()):
 vv=M.row(i)
 for pi,rr in zip(piv,basis):vv=vv-vv[pi]*rr
 if vv:
  pi=next(j for j in range(3) if vv[j]);piv.append(pi);basis.append(vv/vv[pi]);print('PIVOT',i,'column',pi,flush=True)
  if len(basis)==3:break
save({'field':KK,'amplitude_ring':P,'matrix':M,'labels':labels,'items':items,'monomial_order':['one','C^2','B^2'],'row_basis':basis,'pivot_columns':piv,'checks':checks,'scope':__doc__},str(p/'independent_delta_critical_amplitude_coefficients.sobj'),compress=False)
(p/'independent_delta_critical_amplitude_coefficients.json').write_text(json.dumps({'rows':M.nrows(),'rank_over_function_field':len(basis),'checks':len(checks),'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n');print('DONE','rank',len(basis),'seconds',time.monotonic()-st,flush=True)
