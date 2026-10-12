"""Genuine H normal rows after old R/e and Q0/line elimination.

Every displayed pair is obtained by exact small-matrix adjugate operations.
New determinant factors remain recorded; no no-lift conclusion is assumed.
"""
from sage.all import *
from pathlib import Path
import time,json,os
p=Path(__file__).resolve().parent;st=time.monotonic()
h=load(str(p/'independent_delta_subcritical_H_two_parameter_jets.sobj'));j=load(str(p/'independent_delta_subcritical_RQ_two_parameter_jets.sobj'));q=load(str(p/'independent_delta_subcritical_Q0_direction.sobj'))
P=j['ring'];old=j['field'];a0,r0,H0,Z0,p21,p23,p31,p32=old.gens()
B=PolynomialRing(QQ,names=('na','nr','nH','nZ','nt','nB','nC'));K=B.fraction_field();a,r,H,Z,t,Bb,Cc=K.gens();AA=-t/(t+1);x=1/AA
ph=old.ring().hom([a,r,H,Z,AA,Cc*t,Bb,Cc],K)
qh=q['field'].ring().hom([a,r,H,Z,x,t],K)
def qc(c):c=q['field'](c);return qh(c.numerator())/qh(c.denominator())
detdir=qc(q['C_determinant']);en=vector(K,[qc(z)*s for z,s in zip(q['e_numerator'],(AA*Cc,Bb*Cc,AA*Bb))]);eg=en/detdir
cache={}
def cv(c):
 c=old(c)
 if c not in cache:cache[c]=ph(c.numerator())/ph(c.denominator())
 return cache[c]
def pull(f):return sum((cv(c)*prod(eg[k]**int(ex[k]) for k in range(3)) for ex,c in f.dict().items()),K.zero())
def dp(f):return P({tuple(map(int,e)):c.derivative(p23) for e,c in f.dict().items()})
def de(f):return vector(K,[pull(f.derivative(g)) for g in P.gens()])
def vv(v):return vector(K,[pull(f) for f in v])
C=matrix(K,3,3,[pull(f) for f in j['C'].list()]);dc=C.det();adj=C.adjugate();Rp=vv([dp(f) for f in j['R_boundary']]);Rk=vv(j['R_kappa']);Rd=vv(j['R_delta']);Qe=de(j['Q0_boundary']);Qp=pull(dp(j['Q0_boundary']));Qk=pull(j['Q0_kappa']);Qd=K.zero()
checks=[]
def ck(n,v):assert v,n;checks.append(n)
ck('three full R boundary vanish',vv(j['R_boundary'])==0)
ck('full Q0 boundary vanish',pull(j['Q0_boundary'])==0)
ck('C determinant amplitude transport',dc==AA*Bb*Cc*detdir)
QL=dc*Qp-Qe*adj*Rp;QK=dc*Qk-Qe*adj*Rk;QD=-Qe*adj*Rd
print('FORMED old Schur','QL terms',len(QL.numerator().dict()),'seconds',time.monotonic()-st,flush=True)
items=[]
for it in h['items']:
 ix=it['row'];F=it['boundary_value'];Fe=de(F);Fp=pull(dp(F));Fk=pull(it['kappa_coefficient']);Fd=pull(it['delta_coefficient'])
 ck('full H boundary vanish '+str(ix),pull(F)==0)
 FL=dc*Fp-Fe*adj*Rp;FK=dc*Fk-Fe*adj*Rk;FD=dc*Fd-Fe*adj*Rd
 pk=QL*FK-FL*QK;pd=QL*FD-FL*QD
 items.append({'row':ix,'H_e_row':Fe,'H_p23':Fp,'H_kappa':Fk,'H_delta':Fd,'H_line_schur':FL,'H_kappa_schur':FK,'H_delta_schur':FD,'projected_kappa':pk,'projected_delta':pd})
 print('PROJECTED',ix,'terms',len(pk.numerator().dict()),len(pd.numerator().dict()),'zero',pk==0,pd==0,'seconds',time.monotonic()-st,flush=True)
 tmp=p/'independent_delta_subcritical_H_projected_normals.tmp.sobj';save({'field':K,'parameters':(a,r,H,Z,t,Bb,Cc),'p21_value':AA,'e_graph':eg,'C':C,'C_determinant':dc,'C_adjugate':adj,'R_p23':Rp,'R_kappa':Rk,'R_delta':Rd,'Q0_e_row':Qe,'Q0_p23':Qp,'Q0_kappa':Qk,'Q0_delta':Qd,'Q0_line_schur':QL,'Q0_kappa_schur':QK,'Q0_delta_schur':QD,'items':items,'checks':checks,'scope':__doc__},str(tmp),compress=False);os.replace(tmp,p/'independent_delta_subcritical_H_projected_normals.sobj')
(p/'independent_delta_subcritical_H_projected_normals.json').write_text(json.dumps({'checks':len(checks),'rows':len(items),'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n')
