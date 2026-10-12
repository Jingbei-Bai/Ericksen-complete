"""Independent true eight-generator source reader for v=-2, m*l^2=ar, u!=0.

The old first normal archive has a separate full-source audit.  Here pure
kernel orders zero through three are rebuilt from the full original Q/R
polynomials and signed H products using two independent polynomial variables
(ordinary homogenization degree, kernel degree).  No saved cubic is used.
"""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;start=time.monotonic();checks=[]
def ck(n,x):
    assert bool(x),n
    checks.append(n)
normal=load(str(p/'independent_ratio_incoming12_generic_normal.sobj'))
pair=load(str(p/'ordinary_pair_incoming_full_twojet.sobj'))
src=load(str(p/'independent_ratio_p_boundary_second_gauge.sobj'))
hs=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'))
qs=load(str(p.parent/'round8/rank3_joint_fiber_leading.sobj'))
K=pair['coefficient_field'];a,r,H,Z=K.gens();q=1+r;k=a*r;m=1+k;G=H*Z-H+Z;N=k*r
P=PolynomialRing(K,names=('l','u'));l,u=P.gens();mu=m*l*l-k
OK=src['coefficient_field'];O=PolynomialRing(OK,names=('l','u'));ol,ou=O.gens()
B=PolynomialRing(O,names=('degree_marker','kernel_marker'));d,s=B.gens()
images=(d,-d*s,B.zero(),B.zero(),d*s,d*ol,-2*d*ol,d*ou,-2*d)
memo={}
def literal(ff):
    ans=B.zero()
    for exp,c in (ff.items() if isinstance(ff,dict) else ff.dict().items()):
        exp=tuple(map(int,exp))
        if exp[2] or exp[3] or exp[1]+exp[4]>3 or sum(exp)>5:continue
        term=B(c)
        for j,e in enumerate(exp):
            if e:term*=images[j]**e
        ans+=term
    return ans
def pool(i):
    if i not in memo:memo[i]=literal(hs['polynomial_pool'][i])
    return memo[i]
ph=OK.hom([a,a*q/m,H/r,q*Z/r],K);cache={}
def coeffmap(c):
    c=OK(c)
    if c not in cache:cache[c]=ph(c)
    return cache[c]
def mp(f):return P({tuple(map(int,e)):coeffmap(c) for e,c in O(f).dict().items()})
labels=[('old',0),('old',2),('old',1),('old',5),('old',7),('old',11),('old',13),('pair',2)]
pure={}
for lab in labels:
    typ,i=lab
    if typ=='pair':
        f=pair['items'][i]['full_homogeneous_pullback'];vals=[P.zero() for j in range(4)]
        for ee,c in f.dict().items():
            e=tuple(map(int,ee));j=e[3]+e[6]
            if e[4] or e[5] or e[7] or e[8] or j>3:continue
            vals[j]+=K(c)*(-1)**e[3]*(-2)**e[2]*l**e[0]*u**e[1]
    else:
        if i<4:
            f=qs['quartic_source_rows'][(0,1,5,6)[i]];D=int(f.total_degree());f=literal(f)
        elif i<10:D=4;f=literal(src['six_full_sources'][i-4])
        else:
            D=5;f=sum((sgn*pool(j)*pool(h) for sgn,j,h in hs['records'][i-10]['full_signed_product_DAG']),B.zero())
        vals=[mp(f.monomial_coefficient(d**D*s**j)) for j in range(4)]
    pure[lab]=vals
    ck('full pure boundary '+str(lab),vals[0]==0)
    ck('full pure kernel linear on v=-2 '+str(lab),vals[1]==0)
    ck('full quadratic mu coefficient '+str(lab),vals[2]==(8*mu if lab==('old',2) else 0))
    print('SOURCE',lab,'sec',round(time.monotonic()-start,2),flush=True)
OR=normal['ring'];ol0,ou0,ov0=OR.gens()
def pull(f):
    out=P.zero()
    for e,c in OR(f).dict().items():out+=coeffmap(c)*l**int(e[0])*u**int(e[1])*(-2)**int(e[2])
    return out
M=[]
for lab in labels:
    typ,i=lab
    if typ=='old':
        row=list(normal['matrix'].row(i))
        ck('full relative kernel divided by v+2 '+str(lab),not (row[3]-row[0]).subs({ov0:-2}))
        ans=[pull(row[j]) for j in (0,1,2,4,5)]
        wk=pull((row[3]-row[0]).derivative(ov0))
    else:
        f=pair['items'][i]['full_homogeneous_pullback'];pr=pair['jet_ring'];pl,pu,pv=pr.gens();row=[]
        for j in range(6):
            ex=tuple(1 if h==j else 0 for h in range(6))
            row.append(pr({tuple(map(int,e[:3])):K(c) for e,c in f.dict().items() if tuple(map(int,e[3:]))==ex}))
        def pp(ff):return sum((K(c)*l**int(e[0])*u**int(e[1])*(-2)**int(e[2]) for e,c in ff.dict().items()),P.zero())
        ck('full relative pair kernel divided by v+2',not (row[3]-row[0]).subs({pv:-2}))
        ans=[pp(row[j]) for j in (0,1,2,4,5)];wk=pp((row[3]-row[0]).derivative(pv))
    M.append(ans+[wk,P(8) if lab==('old',2) else P.zero(),pure[lab][3]])
M=matrix(P,M)
U=Z**2*G/(7*(Z-1)**2*(Z-H))
ck('eta column exact support',list(M.column(3))==[U*u]+[P.zero()]*7)
ck('mu kernel-square column exact support',list(M.column(6))==[P.zero(),P(8)]+[P.zero()]*6)
def determinant(A):
    if not A.nrows():return P.one()
    i=min(range(A.nrows()),key=lambda j:sum(bool(v) for v in A.row(j)))
    return sum(((-1)**(i+j)*A[i,j]*determinant(A.delete_rows([i]).delete_columns([j])) for j in range(A.ncols()) if A[i,j]),P.zero())
six=determinant(M.matrix_from_rows_and_columns(range(2,8),(0,1,2,4,5,7)))
eight=determinant(M)
ck('complete eight determinant exact Laplace',eight==U*u*8*six)
def rem(ff):return sum((c*(k/m)**(int(e[0])//2)*l**(int(e[0])%2)*u**int(e[1]) for e,c in ff.dict().items()),P.zero())
saved=load(str(p/'ordinary_pair_plus_wall_coupled_cubic.sobj'))
ck('creator six minor independent reconstruction',six==saved['items'][0]['six_minor'])
D=N+m*q*Z**2
C=-16*a**4*r**3*H*Z*(Z+1)/((H+1)**2*(Z-H)*m**5*G**2*(Z-1)**3)
F=m*Z**3*(H*(Z-7)+3)*u*u+2*(Z-1)**3*((Z**2+2*Z+3)-8*H*(2*Z+1))*D
ck('whole wall factor without new material division',rem(six)==C*u*F)
Vring=PolynomialRing(K,'V');V=Vring.gen()
coeffs=[2*(Z-1)**3*(Z**2+2*Z+3)*D,3*m*Z**3*V-16*(Z-1)**3*(2*Z+1)*D,m*Z**3*(Z-7)*V]
ck('complete H coefficient expansion',m*Z**3*(H*(Z-7)+3)*H*V+2*(Z-1)**3*((Z**2+2*Z+3)-8*H*(2*Z+1))*D==sum(coeffs[j]*H**j for j in range(3)))
ck('constant coefficient D leading Z old unit',D==N+m*q*Z**2 and m*q!=0)
ck('D exactly extra tau wall on mu',rem(r*l*l+q*Z**2)==D/m)
save({'ring':P,'labels':labels,'matrix':M,'full_pure_kernel_restrictions':pure,'six_minor':six,'eight_minor':eight,'factor_unit':C,'factor':F,'D':D,'H_coefficient_ring':Vring,'H_coefficients':coeffs,'checks':checks,'scope':__doc__},str(p/'ordinary_pair_plus_wall_cubic_reader.sobj'),compress=False)
report={'status':'PASS TRUE EIGHT GENERATOR WHOLE WALL CERTIFICATE','checks':len(checks),'seconds':time.monotonic()-start,'scope':'v=-2,u!=0,l!=0,m*l^2=ar; all material a*r*(1+r)*(1+ar)!=0; generic base, all ramification via whole-ideal membership.'}
(p/'ordinary_pair_plus_wall_cubic_reader.json').write_text(json.dumps(report,indent=2)+'\n');print(report,flush=True)
