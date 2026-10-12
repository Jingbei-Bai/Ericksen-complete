"""Full-parameter algebraic rank bound on the necessary L0=0 branch.

Exploration only until an explicit source identity is saved and audited.
No spectrum specialization is used.
"""
from pathlib import Path
from sage.all import load, save, PolynomialRing, QQ, matrix
import json

out=Path(__file__).resolve().parent
d=load(str(out.parent/'round6/balanced_line_generic_laurent_sources.sobj'))
P=d['P']; K=d['K']; a,b,h=P.gens()
report=[]
clean=[]
for i,g in enumerate(d['initial']):
    report.append('SOURCE '+str(i))
    gd={}
    for m,c in sorted(g.dict().items()):
        n,den=c.numerator(),c.denominator()
        gcd=n.gcd(den)
        cc=K(n//gcd)/K(den//gcd)
        gd[tuple(m)]=cc
        report.append(str(tuple(m))+' : '+str(cc.factor()))
    clean.append(gd)
save({'coefficients':clean,'P':P,'K':K},str(out/'rank3_balanced_line_clean.sobj'))
(out/'rank3_balanced_line_clean.txt').write_text('\n'.join(report)+'\n')
R=PolynomialRing(K,names=('u','v','s')); u,v,s=R.gens()
polys=[R({tuple(e+1 for e in m):c for m,c in gd.items()}) for gd in clean]
leading=[R({m:c for m,c in f.dict().items() if sum(m)==f.total_degree()}) for f in polys]
E=K(b-1); F=K(a-1); D=K(a-b); A=K(a*(b-1))
mult=K(b*h)*(A+D*h)/(2*a**2*D**3*E**3)
X=D*E*u*u; Y=D*F*v*v; Z=E*F*s*s
q0=(Y-Z)**2-X*(Y+Z)
q1=(X+Z)**2-Y*(X-Z)
checks={
 'source_0_leading_identity':leading[0]==mult*u*q0/F,
 'source_1_leading_identity':leading[1]==-mult*v*q1/E,
 'source_0_total_degree_5':polys[0].total_degree()==5,
 'source_1_total_degree_5':polys[1].total_degree()==5,
 'u_not_divisor_q1':q1.subs({u:0})!=0,
 'v_not_divisor_q0':q0.subs({v:0})!=0,
}
S=PolynomialRing(QQ,names=('X','Y','Z')); xx,yy,zz=S.gens()
p0=(yy-zz)**2-xx*(yy+zz);p1=(xx+zz)**2-yy*(xx-zz)
quadratic_matrix=matrix(QQ,[[0,-QQ(1)/2,-QQ(1)/2],[-QQ(1)/2,1,-1],[-QQ(1)/2,-1,1]])
checks.update({
 'q0_matrix_identity':sum(quadratic_matrix[i,j]*S.gens()[i]*S.gens()[j] for i in range(3) for j in range(3))==p0,
 'q0_matrix_determinant_minus_one':quadratic_matrix.det()==-1,
 'q1_not_multiple_q0':p0.monomial_coefficient(xx**2)==0 and p1.monomial_coefficient(xx**2)==1,
 'sum_positive_form':p0+p1==(xx-yy+zz/2)**2+QQ(7)/4*zz**2,
 'boundary_u_identity':q1.subs({u:0})==Z*(Z+Y),
 'boundary_v_identity':q0.subs({v:0})==Z*(Z-X),
})
assert all(checks.values()),checks
save({'polys':polys,'leading':leading,'q0':q0,'q1':q1,'multiplier':mult,'checks':checks},str(out/'rank3_balanced_line_polys.sobj'))
result={'status':'PASS_EXACT_ALL_PARAMETER_LEADING_PAIR', 'checks':checks,
 'parameters_retained':['a','b','h'], 'spectral_specializations':0,
 'leading_multiplier':'b*h*(a*(b-1)+(a-b)*h)/(2*a^2*(a-b)^3*(b-1)^3)',
 'scope':'Together with the finite-map and fiber-dimension argument in the proof, excludes rank 3 on every actual L0-identically-zero open patch. Does not exclude L0!=0 rank 3.'}
(out/'rank3_balanced_line_rank_bound.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2),flush=True)
