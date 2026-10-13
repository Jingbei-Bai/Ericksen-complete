"""All seven constant direction walls: inspect the complete old e/delta block."""
from sage.all import *
from pathlib import Path
import time
p=Path(__file__).resolve().parent;st=time.monotonic()
src=load(str(p/'independent_delta_critical_linear_sources.sobj'));old=src['field']
W=PolynomialRing(QQ,'wall_omega');wo=W.gen();L=NumberField(wo**2-QQ(13)/2*wo+1,'wall_omega');omega=L.gen()
R=PolynomialRing(L,names=('cb_a','cb_r','cb_H','cb_Z','cb_B','cb_C'),order='degrevlex');K=R.fraction_field();a,r,H,Z,B,C=K.gens();q=1+r
records=[]
for label,s,z,amp_relation in [('point11',L.one(),L.one(),None),('one_omega',L.one(),omega,'B2=rC2'),('omega_one',omega,L.one(),'B2=-qs2'),('omega_omega',omega,omega,'C2=-q/r')]:
    ph=old.ring().hom([a,r,H,Z,K(s),C*K(z),B,C],K)
    def cv(f):f=old(f);return ph(f.numerator())/ph(f.denominator())
    rows=[[cv(src['matrix'][i,j]) for j in range(5)] for i in range(4)]
    dens=[lcm([f.denominator() for f in row]) for row in rows]
    poly=matrix(R,[[R(f*dd) for f in row[:4]] for row,dd in zip(rows,dens)])
    det=poly.det();print('BLOCK',label,'terms',len(det.dict()),'seconds',time.monotonic()-st,flush=True)
    fact=det.factor() if det else None;print('BLOCK_FACTORS',label,[] if fact is None else [(int(e),len(f.dict()),int(f.total_degree()),f if len(f.dict())<=75 else '') for f,e in fact],flush=True)
    red=K.zero()
    for ee,cc in det.dict().items():
        ee=tuple(map(int,ee));base=K(cc)*a**ee[0]*r**ee[1]*H**ee[2]*Z**ee[3];bexp,cexp=ee[4:]
        if amp_relation=='B2=rC2':base*=r**(bexp//2);cexp+=2*(bexp//2);bexp%=2
        elif amp_relation=='B2=-qs2':base*=(-q*K(s)**2)**(bexp//2);bexp%=2
        elif amp_relation=='C2=-q/r':base*=(-q/r)**(cexp//2);cexp%=2
        red+=base*B**bexp*C**cexp
    print('BLOCK_AFTER_RETAINED_AMPLITUDE_RELATION',label,'zero',red==0,'numterms',len(red.numerator().dict()),flush=True)
    if red:print('REDUCED_BLOCK_FACTORS',label,[(int(e),len(f.dict()),f if len(f.dict())<=75 else '') for f,e in red.numerator().factor()],flush=True)
    records.append({'direction_type':label,'p21':s,'z':z,'amplitude_relation':amp_relation,'rows':rows,'row_denominators':dens,'polynomial_block':poly,'block_determinant':det,'block_factors':fact,'restricted_block':red})
    save({'source_file':'independent_delta_critical_linear_sources.sobj','number_field':L,'ring':R,'field':K,'records':records,'scope':__doc__},str(p/'independent_critical_exceptional_old_block.sobj'),compress=False)
print('SECONDS',time.monotonic()-st,flush=True)
