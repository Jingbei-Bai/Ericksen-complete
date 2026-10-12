"""Exact first-normal transport of the complete incoming family to its RR-torus infinity.

Rows0..27 are original Q0, three even Q, six R and eighteen H.
No standalone initial equation is differentiated. Source homogeneizations
are related exactly by the indicated powers of kappa.
"""
from sage.all import *
from pathlib import Path
import json,time
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_ratio_incoming12_generic_normal.sobj'));K=d['coefficient_field'];a,b,h,z=K.gens();r=(a-b)/(a*(b-1));q=1+r;H=r*h;Z=r*z/q;alpha=(Z-H)/(Z*(H+1));beta=H*(Z-1)/(Z*(H+1));shift=-2*H*(Z+1)/(Z*(H+1));P=PolynomialRing(K,names=('ix_l','ix_rho','ix_k'));l,rho,k=P.gens();F=P.fraction_field();checks=[]
def ck(s,v):assert v,s;checks.append(s)
J=matrix(F,6,6);J[0,0]=(alpha/beta)*rho/(l*(alpha+k*shift));J[0,3]=k/(l*(alpha+k*shift));J[1,1]=beta/rho;J[1,4]=k/rho;J[2,0]=1;J[3,1]=1;J[4,2]=1;J[5,5]=k*k
Mold=d['matrix'];M=matrix(F,28,6);N=matrix(F,28,6);records=[]
Q=PolynomialRing(K,names=('il','irho'));ll,rr=Q.gens();QF=Q.fraction_field();N0=matrix(QF,28,6)
def ev(f):return sum((F(co)*l**int(ex[0])*(rho/k)**int(ex[1])*(1/k)**int(ex[2]) for ex,co in f.dict().items()),F.zero())
def atzero(f):
 nm=f.numerator().subs({k:0});dm=f.denominator().subs({k:0});ck('regular_at_k_zero',dm!=0)
 def cv(pp):return Q({tuple(map(int,ex[:2])):co for ex,co in P(pp).dict().items()})
 return QF(cv(nm))/cv(dm)
for i in range(28):
 for j in range(6):M[i,j]=ev(Mold[i,j])
 fac=k if i==0 or i>=10 else F.one()
 nv=fac*M.row(i)*J
 for j in range(6):N[i,j]=nv[j];N0[i,j]=atzero(nv[j])
 records.append({'original_source':d['source_records'][i],'homogeneous_row_factor':fac})
 print('ROW',i,'nz0',sum(bool(x) for x in N0.row(i)),'secs',time.monotonic()-st,flush=True)
save({'source_file':'independent_ratio_incoming12_generic_normal.sobj','coefficient_field':K,'ring':P,'old_evaluation':M,'coordinate_Jacobian':J,'row_records':records,'transported_matrix':N,'boundary_ring':Q,'boundary_matrix':N0,'checks':checks,'scope':__doc__},str(p/'independent_inactive_RRtorus_infinity_normal.sobj'),compress=False)
(p/'independent_inactive_RRtorus_infinity_normal.json').write_text(json.dumps({'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n');print('PASS',len(checks),'seconds',time.monotonic()-st,flush=True)
