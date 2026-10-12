from sage.all import *
from pathlib import Path
import time
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_inactive_RRtorus_resonant_second_module.sobj'));sp=load(str(p/'independent_inactive_RRtorus_pair_regular_sources.sobj'));Q=d['ring'];K=Q.base_ring();F=Q.fraction_field();r=K.gen();q=1+r;a=-q/r;m=1+a*r;l,rho=Q.gens();lc=-q/r;rc=-q;GK=sp['ring'].base_ring();ph=GK.ring().hom((a,r),K)
def cf(z):z=GK(z);return ph(z.numerator())/ph(z.denominator())
def rp(ff):return sum((co*lc**(int(e[0])//2)*rc**(int(e[1])//2)*l**(int(e[0])%2)*rho**(int(e[1])%2) for e,co in ff.dict().items()),Q.zero())
def red(ff):
 ff=F(ff);nn=rp(ff.numerator());de=rp(ff.denominator())
 for v in (l,rho):
  cj=de.subs({v:-v});nn=rp(nn*cj);de=rp(de*cj)
 assert de;return F(nn)/de
Nr=d['normal'];ks,kd=d['kernels'];ids=d['source_rows'];rr=[ids.index(j) for j in (0,5,7)];cols=[2,3,4];outs=[0,1,5];A=Nr.matrix_from_rows_and_columns(rr,cols);Ai=matrix(F,[[red(v/red(A.det())) for v in row] for row in A.adjugate().rows()]);T=Nr.matrix_from_rows_and_columns(rr,outs);E=matrix(F,[[red(v) for v in row] for row in (Nr.matrix_from_columns(outs)-Nr.matrix_from_columns(cols)*Ai*T).rows()]);ei=next(ids.index(j) for j in (1,2) if red(E[ids.index(j),1]));ew=E[ei,1]
# Rebuild the original pre-Schur nine columns from saved first-normal derivative source and exact quadratics.
raw=load(str(p/'independent_inactive_RRtorus_wall_variation_raw.sobj'));RK=raw['ring'].base_ring();rph=RK.ring().hom((a,r,1,-1),K)
def rcv(ff):
 ff=raw['ring'].fraction_field()(ff)
 def pp(g):return sum((rph(co.numerator())/rph(co.denominator())*l**int(e[0])*rho**int(e[1]) for e,co in g.dict().items()),Q.zero())
 return F(pp(ff.numerator()))/pp(ff.denominator())
N=matrix(F,[[rcv(v) for v in row] for row in raw['normal'].rows()]);Nk=matrix(F,[[rcv(v) for v in row] for row in raw['kappa_derivative'].rows()]);qs,qd,qcross=d['quad_sources'];V=matrix(F,[[red(sum(N[i,j].derivative(l)*kv[j] for j in range(6))/(2*l)) for kv in (ks,)] for i in range(len(ids))])
V=matrix(F,[[red(sum(N[i,j].derivative(v)*kv[j] for j in range(6))/den) for kv in (ks,kd) for v,den in ((l,2*l),(rho,2*rho))] for i in range(len(ids))])
# explicitly correct order Wl*s,WR*s,k*s,Wl*d,WR*d,k*d,s2,sd,d2
V=matrix(F,[[red(sum(N[i,j].derivative(l)*ks[j] for j in range(6))/(2*l)),red(sum(N[i,j].derivative(rho)*ks[j] for j in range(6))/(2*rho)),red(sum(Nk[i,j]*ks[j] for j in range(6))),red(sum(N[i,j].derivative(l)*kd[j] for j in range(6))/(2*l)),red(sum(N[i,j].derivative(rho)*kd[j] for j in range(6))/(2*rho)),red(sum(Nk[i,j]*kd[j] for j in range(6))),qs[sid],qcross[sid],qd[sid]] for i,sid in enumerate(ids)])
Vrange=-Ai*V.matrix_from_rows(rr);V1=matrix(F,[[red(v) for v in row] for row in (V-Nr.matrix_from_columns(cols)*Ai*V.matrix_from_rows(rr)).rows()]);wv=-V1.row(ei)/ew;erange=Vrange-Ai*T.column(1).column()*wv.row();items=[]
for it in sp['items']:
 ff=it['full_source'];gvars=sp['ring'].gens();lk,rhok,k,s,w,e0,e1,e2,de=gvars;ns=[s,w,e0,e1,e2,de];zero={v:0 for v in ns};zero[k]=0
 def small(zz):
  zz=sp['ring'].fraction_field()(zz)
  def pp(g):
   assert all(not any(e[2:]) for e in g.dict())
   return sum((cf(co)*l**int(e[0])*rho**int(e[1]) for e,co in g.dict().items()),Q.zero())
  return F(pp(zz.numerator()))/pp(zz.denominator())
 Gn=vector(F,[small(ff.derivative(v).subs(zero)) for v in ns]);Gk=vector(F,[small(ff.derivative(v).derivative(k).subs(zero)) for v in ns]);Hess=matrix(F,[[small(ff.derivative(v).derivative(z).subs(zero)) for z in ns] for v in ns]);assert red(Gn*ks)==red(Gn*kd)==0
 pv=vector(F,[red(sum(Gn[j].derivative(l)*ks[j] for j in range(6))/(2*l)),red(sum(Gn[j].derivative(rho)*ks[j] for j in range(6))/(2*rho)),red(Gk*ks),red(sum(Gn[j].derivative(l)*kd[j] for j in range(6))/(2*l)),red(sum(Gn[j].derivative(rho)*kd[j] for j in range(6))/(2*rho)),red(Gk*kd),red(ks*Hess*ks/2),red(ks*Hess*kd),red(kd*Hess*kd/2)])
 pv=vector(F,[red(v) for v in pv+vector(F,[Gn[j] for j in cols])*erange+Gn[1]*wv]);rem=vector(F,pv);comb=[]
 for j,row in zip(d['pivots'],d['row_reduced']):
  z=red(rem[j]);comb.append(z);rem=vector(F,[red(rem[c]-z*row[c]) for c in range(9)])
 print('PAIR',it['pair'],'REMAINDER',[str(v.factor()) if v else '0' for v in rem],flush=True);items.append(dict(pair=it['pair'],normal=Gn,quotient=pv,old_row_coefficients=comb,remainder=rem))
save(dict(ring=Q,items=items,scope='Exact new pair-source second cokernel at resonance; all old four-range compensations included. No finite-jet no-lift inference.'),str(p/'independent_inactive_RRtorus_resonant_pair_second_module.sobj'),compress=False);print('DONE',time.monotonic()-st,flush=True)
