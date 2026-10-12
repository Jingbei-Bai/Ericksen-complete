"""Check the two resonant fat quotients using all 28 boundary rows and R7.

The complete relation kap^9 delta^4 R7 is built directly from the bundled
full pair source. No lifted-field or regular-pair cache is required.
"""
from sage.all import *
from pathlib import Path
import time,json
# Locate the extracted bundle (or the source workspace) from this file.
root = next(parent for parent in Path(__file__).resolve().parents
            if (parent / 'research/ericksen_20261007/round10').is_dir())
p = root / 'research/ericksen_20261007/round10'
st=time.monotonic()
K=PolynomialRing(QQ,names=('r','delta')).fraction_field(); r,delta=K.gens();q=1+r;a=-q/r;m=-r
P=PolynomialRing(K,names=('l','rho'));l,rho=P.gens();F=P.fraction_field()
pv=[F(1),F(0),F(0),F(0),F(0),F(l)];rv=[F(l),F(rho),F(1)];vals=pv+rv
flow_s={1:rho/l,3:F(1)};flow_w={2:1/rho,4:F(1)}
src=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'));old=src['ring'].base_ring();ph=old.ring().hom([a,r,1,-1],K)
def cf(co): co=old(co);return ph(co.numerator())/ph(co.denominator())
def pull(dd,flow=None,cv=cf):
 ans=F(0)
 for ex,co in dd.items():
  ee=list(map(int,ex));co=F(cv(co))
  if flow is None:
   ans+=co*prod(vals[j]**ee[j] for j in range(9))
  else:
   for j,vel in flow.items():
    if not ee[j]:continue
    ee1=ee.copy();ee1[j]-=1;ans+=co*ee[j]*vel*prod(vals[k]**ee1[k] for k in range(9))
 return ans
def layer(f,bi):return {tuple(map(int,e)):c for e,c in f.dict().items() if (sum(e[:6]),sum(e[6:]))==bi}
def row(ee,low,delpart=None,cv=cf):
 ds=pull(low,flow_s,cv);dw=pull(low,flow_w,cv)
 if delpart is not None:ds+=delta*pull(delpart,flow_s,cv);dw+=delta*pull(delpart,flow_w,cv)
 return [ds,dw]+[pull(z,cv=cv) for z in ee]
DC=src['D']*src['C']; rows=[]; values=[]
rows.append(row([f.dict() for f in src['Q0_L_cofactors']],layer(src['full_Q0'],(2,1))))
values.append(pull(layer(src['full_Q0'],(2,1))))
ev=load(str(p/'independent_delta_subcritical_even_full_source.sobj'))
for it in ev['items']:
 rows.append(row([f.dict() for f in it['A_p_quadratic']],it['pure_p4'].dict()));values.append(pull(it['pure_p4'].dict()))
sp=load(str(p/'independent_supercritical_first_boundary_probe.sobj'))
for j,ff in enumerate(src['full_R_sources']):
 rows.append(row([DC[j,k].dict() for k in range(3)],layer(ff,(4,0)),sp['F12'][j].dict()));values.append(pull(layer(ff,(4,0)))+delta*pull(sp['F12'][j].dict()))
full=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'));orig=full['coefficient_field'];op=orig.ring().hom([a,a*q/m,1/r,-q/r],K)
def oc(co):co=orig(co);return op(co.numerator())/op(co.denominator())
hi=load(str(p/'independent_delta_subcritical_H_initials.sobj'))
for j,it in enumerate(hi['items']):
 dp={e:oc(c) for e,c in it['H41_original'].items()};ep=[{e:oc(c) for e,c in z.items()} for z in it['L_cofactors_original']]
 # H13 has the same renamed coefficient field as Q0.
 hp={tuple(e):cf(c) for e,c in sp['H13'][j].dict().items()}
 ident=lambda z:K(z)
 rows.append(row(ep,dp,hp,ident));values.append(pull(dp,cv=ident)+delta*pull(hp,cv=ident))
print('all28 genuine boundary rows',len(rows),'seconds',round(time.monotonic()-st,2),flush=True)
lc=-q/r;rc=-q
def rp(pp,both):
 ans=P.zero()
 for ex,co in P(pp).dict().items():
  e0,e1=map(int,ex);ans+=co*lc**(e0//2)*l**(e0%2)*(rc**(e1//2)*rho**(e1%2) if both else rho**e1)
 return ans
def red(ff,both):
 ff=F(ff);nn=rp(ff.numerator(),both);dd=rp(ff.denominator(),both)
 for vv in ([l,rho] if both else [l]):
  conj=dd.subs({vv:-vv});nn=rp(nn*conj,both);dd=rp(dd*conj,both)
 assert dd
 return F(nn)/dd
def ri(A,both):return matrix(F,[[red(x,both) for x in rr] for rr in A.rows()])
dz=K.hom([r,K(0)],K)
def dzero(z):
 z=F(z)
 def pp(g):return P({tuple(e):dz(c) for e,c in g.dict().items()})
 return F(pp(z.numerator()))/pp(z.denominator())
# Construct the complete R7 over the general ratio base before specializing.
pa=load(str(p/'independent_ratio_p_boundary_global_pair_syzygies.sobj'))
pair=next(it for it in pa['items'] if int(it['pair'])==1)
pair_K=pa['ring'].base_ring()
KG=PolynomialRing(QQ,names=('a','r','H','Z')).fraction_field()
ga,gr,H,Z=KG.gens();gq=1+gr;gm=1+ga*gr
pair_map=pair_K.ring().hom([ga,ga*gq/gm,H/gr,gq*Z/gr],KG)
pair_coefficients={}
def pair_cf(co):
 co=pair_K(co)
 if co not in pair_coefficients:
  pair_coefficients[co]=pair_map(co.numerator())/pair_map(co.denominator())
 return pair_coefficients[co]
QG=PolynomialRing(KG,names=('l','rho','kap','s','w','E0','e1','e2','delta'))
ql,qr,kap,s,w,E0,e1,e2,gd=QG.gens()
alpha=(Z-H)/(Z*(H+1));beta=H*(Z-1)/(Z*(H+1))
shift=-2*H*(Z+1)/(Z*(H+1))
v=alpha*ql+kap*(E0+shift*ql)
pn=[QG(1),(alpha/beta)*s*qr+kap*e1,beta*w+kap*e2,s,w,ql]
rv=[v,qr,QG(1)]
terms=pair['full'].dict()
degree=max(max(int(ex[1]),int(ex[2])) for ex in terms)
powers=[(9-2*sum(ex[:6])-3*sum(ex[6:]),4-sum(ex[:6])-sum(ex[6:])) for ex in terms]
min_k=min(int(z[0]) for z in powers);min_delta=min(int(z[1]) for z in powers)
assert min_delta>=0
clear_k=max(0,-min_k)
N=QG.zero()
for ex,co in terms.items():
 ex=tuple(map(int,ex));pd,rd=sum(ex[:6]),sum(ex[6:])
 ke,dd=9-2*pd-3*rd+clear_k,4-pd-rd
 N+=pair_cf(co)*prod(pn[j]**ex[j] for j in range(6))*prod(rv[j]**ex[j+6] for j in range(3))*kap**ke*gd**dd*v**(degree-ex[1])*qr**(degree-ex[2])
cleared_N=N
for power in range(clear_k):
 assert N.subs({kap:0})==0,('nonregular kappa coefficient',power)
 quotient=QG({tuple(ex[j]-(1 if j==2 else 0) for j in range(9)):co for ex,co in N.dict().items()})
 assert N==kap*quotient
 N=quotient
assert cleared_N==kap**clear_k*N
# All coefficient denominators stay nonzero at the resonant centre.
center=KG.ring().hom([-gq/gr,gr,1,-1],KG)
assert all(center(co.denominator())!=0 for co in N.coefficients())
assert N.subs({kap:0,s:0,w:0,E0:0,e1:0,e2:0,gd:0})==0
pair_denominator=(v*qr)**degree
# Specialize the complete regular equation to the two fat families.
family_coefficients=KG.ring().hom([-q/r,r,1,-1],K)
def at_family(poly):
 result=F.zero()
 for ex,co in QG(poly).dict().items():
  if any(ex[2:8]):continue
  co=KG(co)
  coefficient=family_coefficients(co.numerator())/family_coefficients(co.denominator())
  result+=coefficient*l**int(ex[0])*rho**int(ex[1])*delta**int(ex[8])
 return result
N0=at_family(N);D0=at_family(pair_denominator)
assert D0
pair_value=N0/D0
# Differentiate the polynomial numerator and denominator separately. This
# quotient rule avoids normalization of a large rational function in Sage.
pair_normal=vector(F,[(at_family(N.derivative(z))*D0-N0*at_family(pair_denominator.derivative(z)))/D0**2 for z in [s,w,E0,e1,e2]])
print('complete general-base R7 regular',len(terms),'source terms;',len(N.dict()),'regular numerator terms',flush=True)
reports={}
for both,label in [(True,'canonical_thickened_family'),(False,'negative_wall_thickened_family')]:
 nr=matrix(F,rows)
 if not both:
  nr=nr.apply_map(dzero)
 nr=ri(nr,both)
 rr=[0,5,7,2]; cc=[1,2,3,4]
 ran=nr.matrix_from_rows_and_columns(rr,cc);inv=ri(ran.inverse(),both)
 assert ri(ran*inv,both)==identity_matrix(F,4)
 kernel=vector(F,[1]+list(-inv*nr.matrix_from_rows_and_columns(rr,[0]).column(0)))
 assert all(red(z,both)==0 for z in nr*kernel)
 vv=[red(z if both else dzero(z),both) for z in values]
 assert not any(vv)
 print(label,'all28 value and first s vanish','kernel',kernel,'seconds',round(time.monotonic()-st,2),flush=True)
 reports[label]={'whole_boundary_value_zero':True,'whole_boundary_first_s_zero':True,'rows':28}
 # Add the complete selected pair relation P = kap^9 delta^4 R7.
 # Its coefficients and all monomials come from the existing global source.
 value, normal = pair_value, pair_normal
 if not both:
  value=dzero(value); normal=vector(F,[dzero(z) for z in normal])
 assert red(value,both)==0,(label,'R7 value')
 assert red(normal*kernel,both)==0,(label,'R7 first s')
 print(label,'complete R7 value and first s vanish',flush=True)
 reports[label].update(complete_R7_value_zero=True,complete_R7_first_s_zero=True)
report={'status':'PASS','selected_pair_index':1,'pair_relation':'kap^9 delta^4 R7',
        'full_pair_terms':len(pair['full'].dict()),'whole_pair_regular':True,
        'full_general_ratio_base':True,'full_delta_polynomial':True,
        'minimum_kappa_exponent_before_summing':min_k,
        'minimum_delta_exponent_before_summing':min_delta,
        'exact_cancelled_kappa_power':clear_k,'cleared_denominator':'(v1 rho)^'+str(degree),
        'unit_coefficient_denominators_at_centre':True,'regular_numerator_terms':len(N.dict()),
        'fat_quotients':reports,'seconds':time.monotonic()-st}
(p/'check_C_resonant_families.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2),flush=True)
