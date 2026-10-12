"""Exact source and coefficient checks for the complete det=0 lemma."""
from pathlib import Path
import time,json
from sage.all import load,save,QQ,PolynomialRing
out=Path(__file__).resolve().parent;start=time.monotonic();checks={}
d=load(str(out/'rank1_ratio_characteristic_linear_pair.sobj'));R=d['ring'];u,v,x,s,y=R.gens()
cd=load(str(out/'rank1_ratio_characteristic_coefficients.sobj'));co=cd['old_coefficients']
lam,A,B,C,D,E,F,G,J,N=(co[k] for k in ['lam','A','B','C','D','E','F','G','Hcoef','C0']);alpha=J+lam*D
f,g=d['linear_source_rows'];T0=u*x+lam*v*y+N
checks['first_affine_formula']=f==v*T0+s*s*(A*v+B*y)+s*(C*u+D*x)
checks['second_affine_formula']=g==u*T0+s*s*(E*u+F*x)+s*(G*v+J*y)
delta=-B*s*u*u-lam*F*s*v*v+alpha*u*v+s*D*J-B*F*s**3
checks['det_formula']=d['determinant']==s*delta
X=(lam*(E-A)*u*v*v*s+B*E*u*s**3+(B*N-J*C)*u*s
   +lam*G*v**3-lam*C*u*u*v+(B*G-J*A)*v*s*s-J*N*v)
checks['x_cramer_formula']=d['cramer_numerators'][0]==s*X
j0=d['remaining_mixed_source']
checks['mixed_no_x_squared']=not j0.derivative(x,2)
checks['mixed_no_y_squared']=not j0.derivative(y,2)
checks['mixed_xy_coefficient']=j0.derivative(x).derivative(y)==co['Jxy']*s

# Abstract zero-column consistency top coefficients, with no division by
# the individual actual state variables u,v,x,y.
P=PolynomialRing(QQ,names=('la','aa','bb','cc','dd','ee','ff','gg','jj','nn','ep','ss','vzero','uzero'))
la,aa,bb,cc,dd,ee,ff,gg,jj,nn,ep,ss,vzero,uzero=P.gens()
def determinant2(a,b,c,d):return a*d-b*c
uu=ep*ss;vv=vzero
b0=vv*(nn+aa*ss**2)+cc*uu*ss;b1=uu*(nn+ee*ss**2)+gg*vv*ss
compat1=determinant2(la*vv**2+bb*ss**2,la*uu*vv+jj*ss,b0,b1)
def coeff_ss(p,n):
 return P({tuple(e if i!=11 else 0 for i,e in enumerate(m)):c for m,c in p.dict().items() if m[11]==n})
checks['first_zero_column_s5']=coeff_ss(compat1,5)==bb*ep*ee
uu=uzero;vv=ep*ss
b0=vv*(nn+aa*ss**2)+cc*uu*ss;b1=uu*(nn+ee*ss**2)+gg*vv*ss
compat2=determinant2(uu*vv+dd*ss,uu**2+ff*ss**2,b0,b1)
checks['second_zero_column_s5']=coeff_ss(compat2,5)==-ff*ep*aa

CR=cd['coefficient_ring'];K=cd['field'];a,q,T,H=cd['parameters'];r=q-1;m=1+a*r;ccf=cd['coefficients']
guard=CR(a*q*r*m*T*(T-1)*H*(H+T)*(H+1)*(1-H-T))
polys={
 'A':-T*T*H+T*H*H-2*T*T-T*H+H*H+2*T+4*H,
 'B':T*H*H+T*T+3*T*H-2*H*H-T-4*H,
 'D':T**3+2*T*T*H+T*H*H-T*T-2*T*H-2*H*H-T-2*H,
 'E':2*T**3+3*T*T*H+T*H*H+3*T*H+H*H-2*T-4*H,
 'F':T**3+2*T*T*H+T*H*H-4*T*T-7*T*H-2*H*H+3*T+4*H,
 'G':T**3*H+2*T*T*H*H+T**3+T*T*H-T*H*H+H**3-T*T-2*T*H,
 'Hcoef':T*H*H-T*T-2*T*H-2*H*H+T+2*H}
prefactors={}
for name,p in polys.items():
 rat=ccf[name]/p;prefactors[name]=rat;extra=[]
 for pol in [rat.numerator(),rat.denominator()]:
  for fac,power in pol.factor():
   if guard.quo_rem(fac)[1]:extra.append(str(fac))
 checks[name+'_only_protected_prefactor']=not extra
def hc(p,n):
 return K(sum((v*a**mo[0]*q**mo[1]*T**mo[2] for mo,v in CR(p).dict().items() if mo[3]==n),K.zero()))
checks['A_constant']=hc(polys['A'],0)==2*T*(1-T)
checks['B_constant']=hc(polys['B'],0)==T*(T-1)
checks['J_constant']=hc(polys['Hcoef'],0)==T*(1-T)
checks['D_H2']=hc(polys['D'],2)==T-2
checks['D_at_T2_constant']=hc(polys['D'],0).subs({T:2})==2
checks['F_H2']=hc(polys['F'],2)==T-2
checks['F_at_T2_constant']=hc(polys['F'],0).subs({T:2})==-2
checks['E_H2']=hc(polys['E'],2)==T+1
checks['E_at_Tminus1_H']=hc(polys['E'],1).subs({T:-1})==-4
checks['G_H3']=hc(polys['G'],3)==1
checks['lambda_protected']=ccf['lam']==m*H/(T+H)
checks['Jxy_protected']=ccf['Jxy']==-4*q*T/((T-1)*r)
checks['alpha_nonzero_pole']=(H*H*ccf['alpha']).subs({H:0})==2*T*T*m/(a*r*r)
assert all(checks.values()),{k:v for k,v in checks.items() if not v}
save({'checks':checks,'bad_coefficient_polynomials':polys,'protected_prefactors':prefactors,
 'zero_column_compatibilities':[compat1,compat2],'source_file':'rank1_ratio_characteristic_linear_pair.sobj'},
 str(out/'rank1_ratio_characteristic_singular_fiber_audit.sobj'))
report={'status':'PASS_ALL_PARAMETER_SINGULAR_FIBER_IDENTITIES','checks':checks,'check_count':len(checks),
 'elapsed_seconds':time.monotonic()-start,
 'scope':'Source identities and uniform-in-spectrum/c finite-bad-h coefficient witnesses for the complete det=0 proof; det!=0 and general ratio curves remain separate.'}
(out/'rank1_ratio_characteristic_singular_fiber_audit.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2),flush=True)
