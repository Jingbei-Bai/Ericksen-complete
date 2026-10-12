"""Exact finite-extension reduction of the critical rank-two amplitude graph.

B=s*b, C=c, b^2=-q and c^2=-q/r. No complex amplitude is discarded.
"""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_delta_critical_linear_sources.sobj'));old=d['field'];M=d['matrix'];K=PolynomialRing(QQ,names=('ka','kr','kH','kZ','ks','kz')).fraction_field();a,r,H,Z,s,z=K.gens();q=1+r;P=PolynomialRing(K,names=('kb','kc'));b,c=P.gens();ph=old.ring().hom([P(a),P(r),P(H),P(Z),P(s),c*z,s*b,c],P)
def red(f):
 out={}
 for ee,v in P(f).dict().items():
  e=tuple(map(int,ee));key=(e[0]%2,e[1]%2);out[key]=out.get(key,K.zero())+v*(-q)**(e[0]//2)*(-q/r)**(e[1]//2)
 return P(out)
def cv(v):
 v=old(v);nn=red(ph(v.numerator()));dd=red(ph(v.denominator()));assert len(dd.dict())==1
 ee,co=next(iter(dd.dict().items()));i,j=map(int,ee);inv=b**i*c**j/(co*(-q)**i*(-q/r)**j)
 return red(nn*inv)
MM=matrix(P,M.nrows(),5,[cv(v) for v in M.list()]);parities=[]
for i in range(MM.nrows()):
 rr=[]
 for j in range(5):
  pp=MM[i,j];assert len(pp.dict())<=1,(i,j,pp);rr.append(next(iter(pp.dict())) if pp else None)
 parities.append(rr)
rowphase=[tuple(map(int,v[4])) if v[4] is not None else None for v in parities]
colphase=[]
for j in range(5):
 vals=set(((int(parities[i][j][0])+rowphase[i][0])%2,(int(parities[i][j][1])+rowphase[i][1])%2) for i in range(MM.nrows()) if rowphase[i] is not None and parities[i][j] is not None)
 assert len(vals)<=1,(j,vals);colphase.append(next(iter(vals)) if vals else (0,0))
NN=matrix(K,MM.nrows(),5)
for i in range(MM.nrows()):
 if rowphase[i] is None:
  jj=next(j for j in range(5) if parities[i][j] is not None);rowphase[i]=tuple((int(parities[i][jj][k])+colphase[jj][k])%2 for k in range(2))
 ri,rj=rowphase[i];inv=b**ri*c**rj/((-q)**ri*(-q/r)**rj)
 for j in range(5):
  ci,cj=colphase[j];v=red(MM[i,j]*b**ci*c**cj*inv);assert v==P(v.constant_coefficient());NN[i,j]=v.constant_coefficient()
print('PHASES',colphase,'rowphases',rowphase,'seconds',time.monotonic()-st,flush=True)
assert all(MM[i,3]==0 for i in range(MM.nrows()))
block=NN[:4,:4];det=block.det();print('FOUR DET',det.factor() if det else 0,'seconds',time.monotonic()-st,flush=True)
save({'source_file':'independent_delta_critical_linear_sources.sobj','field':K,'extension_ring':P,'radical_squares':(-q,-q/r),'extension_matrix':MM,'row_phases':rowphase,'column_phases':colphase,'material_matrix':NN,'old4_determinant':det,'scope':__doc__},str(p/'independent_delta_critical_canonical_amplitude.sobj'),compress=False)
(p/'independent_delta_critical_canonical_amplitude.json').write_text(json.dumps({'rows':MM.nrows(),'column_phases':colphase,'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n')
