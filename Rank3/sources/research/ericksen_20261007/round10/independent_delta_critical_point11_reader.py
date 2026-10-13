"""Direct-original-row replay of the complete point-(1,1) certificate.

This reader does not trust saved Schur entries or the creator's row
normalizations, and reduces square-root relations by polynomial division.
"""
from sage.all import *
from pathlib import Path
import json,time
p=Path(__file__).resolve().parent;st=time.monotonic()
src=load(str(p/'independent_delta_critical_linear_sources.sobj'))
cert=load(str(p/'independent_delta_critical_point11_certificate.sobj'))
saved=load(str(p/'independent_delta_critical_point11_schur.sobj'))
K=cert['field'];a,r,H,Z,B,C=K.gens();q=1+r;old=src['field'];ph=old.ring().hom([a,r,H,Z,K.one(),C,B,C],K)
def cv(f):f=old(f);return ph(f.numerator())/ph(f.denominator())
M=src['matrix'].apply_map(cv);RR=M[:3,:];C0=RR[:,:3];adj=C0.adjugate();det=C0.det();checks=[]
def ck(label,yes):assert yes,label;checks.append(label)
ck('original C equals saved',C0==saved['C']);ck('original detC equals saved',det==saved['C_determinant'])
rows=[]
for j in range(18):
 raw=M.row(j+4);co=vector(K,raw[:3])*adj
 row=det*raw-sum((co[i]*RR.row(i) for i in range(3)),vector(K,5))
 ck('three eliminated original entries '+str(j),tuple(row[:3])==(0,0,0))
 ck('original source Schur '+str(j),tuple(row[3:])==(saved['items'][j]['delta'],saved['items'][j]['constant']))
 rows.append(row)
for it in cert['witnesses']:
 i,j=it['source_rows'];f=rows[i][3]*rows[j][4]-rows[j][3]*rows[i][4]
 ck('original minor '+str((i,j)),f==it['minor'])
 ck('all-factor replay '+str((i,j)),f==it['target']*it['factorization'].value())
# Independent normal form by a fixed monic Groebner division.
F=PolynomialRing(QQ,names=('va','vr','vH','vZ')).fraction_field();aa,rr,hh,zz=F.gens();P=PolynomialRing(F,names=('vb','vc'));bb,cc=P.gens();I=P.ideal(bb*bb+1+rr,cc*cc+(1+rr)/rr);gb=list(I.groebner_basis());ep=K.ring().hom([P(aa),P(rr),P(hh),P(zz),bb,cc],P)
def nf(f):return P(f).reduce(gb)
def image(f):
 f=K(f);n=nf(ep(f.numerator()));d=nf(ep(f.denominator()));ck('nonzero reduced denominator',bool(d))
 ck('monomial reduced denominator',len(d.dict())==1)
 ee,co=next(iter(d.dict().items()));i,j=map(int,ee)
 inv=bb**i*cc**j/(co*(-1-rr)**i*(-(1+rr)/rr)**j)
 ck('denominator inverse',nf(d*inv)==1)
 return nf(n*inv)
ck('canonical delta vanishes original',image(rows[0][3])==0)
ck('canonical constant is saved original',image(rows[0][4])==image(cert['canonical_H0_constant']))
Fbase=H*Z-H+Z;Pbase=H*H-3*H*Z+Z*Z+3*H-3*Z+1;Gbase=-2*H*H+H*Z-4*H-Z
expected=-81*B*C*q*q*H**3*Z**3*Gbase/((Z-H)**2*(Z-1)**2*(H+1)**2*Fbase**2)
ck('explicit complete constant',cert['canonical_H0_constant']==expected)
ck('Pbase nonzero all amplitude/material',Pbase.numerator().coefficient({K.ring().gens()[2]:2})==1)
ck('Gbase nonzero all amplitude/material',Gbase.numerator().coefficient({K.ring().gens()[2]:2})==-2)
save({'source_files':['independent_delta_critical_linear_sources.sobj','independent_delta_critical_point11_certificate.sobj'],
      'reconstructed_schur_rows':rows,'monic_square_relation_basis':gb,'checks':checks,'scope':__doc__},str(p/'independent_delta_critical_point11_reader.sobj'),compress=False)
(p/'independent_delta_critical_point11_reader.json').write_text(json.dumps({'checks':len(checks),'status':'PASS','seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n')
print('PASS',len(checks),'seconds',time.monotonic()-st,flush=True)
