"""Cyclic-base diagnostic of genuine second variations; all materials retained.
No generic-base no-lift conclusion follows from this diagnostic alone.
"""
from sage.all import *
from pathlib import Path
import time
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_inactive_RRtorus_wall_variation_raw.sobj'));K0=d['ring'].base_ring();K=PolynomialRing(QQ,names=('rs_r',)).fraction_field();r=K.gen();q=1+r;a=-q/r;m=1+a*r;ph=K0.ring().hom((a,r,1,-1),K);Q=PolynomialRing(K,names=('cy_l','cy_rho'));l,rho=Q.gens();F=Q.fraction_field();cache={}
def cf(c):
 if c not in cache:cache[c]=ph(c.numerator())/ph(c.denominator())
 return cache[c]
def cv(f):
 return F(Q({tuple(map(int,e)):cf(c) for e,c in f.numerator().dict().items()}))/Q({tuple(map(int,e)):cf(c) for e,c in f.denominator().dict().items()})
ids=d['source_rows'];N=matrix(F,[[cv(v) for v in row] for row in d['normal'].rows()]);Nk=matrix(F,[[cv(v) for v in row] for row in d['kappa_derivative'].rows()]);print('CYCLIC',time.monotonic()-st,flush=True)
qf=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'));hsrc=load(str(p/'independent_delta_subcritical_H_initials.sobj'));ev=load(str(p/'independent_delta_subcritical_even_full_source.sobj'));old=load(str(p/'independent_delta_H_global_L_module.sobj'))['original_field'];oldph=old.ring().hom((a,a*q/m,1/r,-q/r),K);hc={}
def hcf(co):
 co=old(co)
 if co not in hc:hc[co]=oldph(co.numerator())/oldph(co.denominator())
 return hc[co]
DC=qf['D']*qf['C']
def biterms(f,ij):return {tuple(map(int,e)):co for e,co in f.dict().items() if (sum(e[:6]),sum(e[6:]))==ij}
sp=load(str(p/'independent_supercritical_first_boundary_probe.sobj'));c=-q/r;d2=-q
def rp(f):return sum((co*c**(int(ex[0])//2)*d2**(int(ex[1])//2)*l**(int(ex[0])%2)*rho**(int(ex[1])%2) for ex,co in f.dict().items()),Q.zero())
def red(f):
 f=F(f);nm=rp(f.numerator());de=rp(f.denominator())
 for vv in (l,rho):
  cj=de.subs({vv:-vv});nm=rp(nm*cj);de=rp(de*cj)
 assert de
 return F(nm)/de
Nr=matrix(F,[[red(v) for v in row] for row in N.rows()]);rr=[ids.index(i) for i in (0,5,7)];cols=[2,3,4];out=[0,1,5];A=Nr.matrix_from_rows_and_columns(rr,cols);Ai=matrix(F,[[red(v) for v in row] for row in A.inverse().rows()]);T=Nr.matrix_from_rows_and_columns(rr,out);Ne=Nr.matrix_from_columns(cols);E=matrix(F,[[red(v) for v in row] for row in (Nr.matrix_from_columns(out)-Ne*Ai*T).rows()]);ei=next(ids.index(ii) for ii in (1,2) if red(E[ids.index(ii),1]));ew=E[ei,1]
def kernel(sv,dv):
 wv=red(-(E[ei,0]*sv+E[ei,2]*dv)/ew);er=vector(F,[red(v) for v in -Ai*T*vector(F,[sv,wv,dv])]);kv=vector(F,[sv,wv,er[0],er[1],er[2],dv]);assert all(red(z)==0 for z in Nr*kv);return kv
ks=kernel(1,0);kd=kernel(0,1);print('KERNELS',ks,kd,flush=True)
Vs=matrix(F,len(ids),3,[[red(sum(N[i,j].derivative(l)*ks[j] for j in range(6))/(2*l)),red(sum(N[i,j].derivative(rho)*ks[j] for j in range(6))/(2*rho)),red(sum(Nk[i,j]*ks[j] for j in range(6)))] for i in range(len(ids))]);Vd=matrix(F,len(ids),3,[[red(sum(N[i,j].derivative(l)*kd[j] for j in range(6))/(2*l)),red(sum(N[i,j].derivative(rho)*kd[j] for j in range(6))/(2*rho)),red(sum(Nk[i,j]*kd[j] for j in range(6)))] for i in range(len(ids))]);
def quads(kv):
 er=kv[2:5];pp=[F.one(),rho/l*kv[0],kv[1]/rho,kv[0],kv[1],l];Rv=[l,rho,F.one()]
 def ps(dd,order,mapco=cf):
  ans=F.zero()
  for ex,co in dd.items():
   ex=tuple(map(int,ex))
   if sum(ex[1:5])!=order:continue
   ans+=F(mapco(co))*prod(pp[j]**ex[j] for j in range(6))*prod(Rv[j]**ex[6+j] for j in range(3))
  return red(ans)
 quad=[ps(biterms(qf['full_Q0'],(2,1)),2)]
 for jt in ev['items']:
  val=cf(jt['square_unit'])*er[int(jt['square_component'])]**2+sum(ps(ff.dict(),1)*er[j] for j,ff in enumerate(jt['A_p_quadratic']))+ps(jt['pure_p4'].dict(),2);quad.append(red(val))
 for jj in range(6):
  val=sum(ps(DC[jj,k].dict(),1)*er[k] for k in range(3))+ps(biterms(qf['full_R_sources'][jj],(4,0)),2)+kv[5]*ps(sp['F12'][jj].dict(),1);quad.append(red(val))
 for jj,jt in enumerate(hsrc['items']):
  val=sum(ps(f0,1,hcf)*er[j] for j,f0 in enumerate(jt['L_cofactors_original']))+ps(jt['H41_original'],2,hcf)+kv[5]*ps(sp['H13'][jj].dict(),1);quad.append(red(val))
 return vector(F,quad)
qs=quads(ks);qd=quads(kd);qsum=quads(ks+kd);cross=qsum-qs-qd;V=Vs.augment(Vd).augment(vector(F,[qs[i] for i in ids])).augment(vector(F,[cross[i] for i in ids])).augment(vector(F,[qd[i] for i in ids]));V1=matrix(F,[[red(v) for v in row] for row in (V-Ne*Ai*V.matrix_from_rows(rr)).rows()]);V2=matrix(F,[[red(v) for v in row] for row in (V1-E.matrix_from_columns([1])*(V1.row(ei)/ew).row()).rows()]);print('QUOTIENT_READY',time.monotonic()-st,flush=True)
# Exact row reduction in the biquadratic quotient, not the ambient fraction field.
rows=[list(row) for row in V2.rows()];source_ops=identity_matrix(F,len(rows));piv=[];rk=0
for j in range(9):
 ii=next((i for i in range(rk,len(rows)) if red(rows[i][j])),None)
 if ii is None:continue
 rows[rk],rows[ii]=rows[ii],rows[rk];source_ops.swap_rows(rk,ii);z=red(rows[rk][j]);rows[rk]=[red(v/z) for v in rows[rk]];source_ops.rescale_row(rk,1/z)
 for i in range(len(rows)):
  if i==rk:continue
  z=red(rows[i][j])
  if z:
   rows[i]=[red(rows[i][k]-z*rows[rk][k]) for k in range(9)];source_ops.add_multiple_of_row(i,rk,-z)
 piv.append(j);rk+=1
print('RANK',rk,'PIVOTS',piv,'REDUCED',[row for row in rows[:rk]],flush=True)
save(dict(ring=Q,source_rows=ids,normal=Nr,kernels=[ks,kd],quotient_second=V2,row_reduced=rows,pivots=piv,source_ops=source_ops,quad_sources=[qs,qd,cross],range_pivots=[red(A.det()),ew]),str(p/'independent_inactive_RRtorus_resonant_second_module.sobj'),compress=False)
print('DONE',time.monotonic()-st,flush=True)
