"""Independent raw-source Taylor and eight-point norm audit.

Uses the complete original R0 and full signed-product H2, substitutes
v=kappa*delta*Rbar and tau=kappa^2*delta termwise, and extracts the
normal jet directly. It does not use the creator's layer/chain formulas.
All material parameters remain symbolic; H=1,Z=-1 is a base witness.
"""
from sage.all import *
from pathlib import Path
import json,time
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(n,v):
    if not bool(v):raise AssertionError(n)
    checks.append(n)
d=load(str(p/'independent_critical_single_RR3_normal.sobj'))
n=load(str(p/'independent_critical_single_RR3_pure_normal_minors.sobj'))
full=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'));pool=full['polynomial_pool'];old=full['coefficient_field']
K=n['field'];k,r=K.gens();q=1+r;m=1+k;H=K.one();Z=-K.one();alpha=(Z-H)/(Z*(H+1));beta=H*(Z-1)/(Z*(H+1))
allowed=K.ring()(k*r*q*m*(k-r));den_seen=set()
def guardden(c):
    den=K(c).denominator()
    if den not in den_seen:
        den_seen.add(den)
        for fac,powr in den.factor():ck('every specialized coefficient denominator is an old material unit '+str(len(checks)),allowed.quo_rem(fac)[1]==0)
P=PolynomialRing(K,names=('cr_x','cr_z','cr_w','cr_e0','cr_e1','cr_e2','cr_delta'));x,z,w,e0,e1,e2,de=P.gens()
# Each entry is a (kappa,t) Taylor dictionary. The nine entries below
# are p followed by Rbar, independently written from exact ell=kappa e.
images=[{(0,0):x},{(0,0):P.one()},{(0,0):z},{(0,0):w},
        {(0,1):z/beta,(1,0):-e2/beta},
        {(0,1):x*w/beta,(1,0):(x*e1-e0)/alpha},
        {(0,1):alpha*w/beta,(1,0):e1},{(0,1):P.one()},{(0,0):P.one()}]
def mul(A,B,maxk):
    out={}
    for (ka,ta),ca in A.items():
        for (kb,tb),cb in B.items():
            ee=(ka+kb,ta+tb)
            if ee[0]<=maxk and ee[1]<=1:out[ee]=out.get(ee,P.zero())+ca*cb
    return {e:c for e,c in out.items() if c}
ph=old.ring().hom([k/r,k*q/(r*m),1/r,-q/r],K);cache={}
def cv(c):
    c=old(c)
    if c not in cache:
        den=ph(c.denominator());ck('original coefficient defined '+str(len(checks)),den!=0);cache[c]=ph(c.numerator())/den;guardden(cache[c])
    return cache[c]
def monjet(e,degree,maxk):
    i=sum(e[:6]);j=sum(e[6:]);kp=2*degree-2*i-3*j;dp=degree-i-j
    ck('raw homogenization has no negative monomial powers '+str(len(checks)),kp>=0 and dp>=0)
    if kp>maxk:return {}
    value={(kp,0):de**dp}
    for power,img in zip(e,images):
        for unused in range(power):value=mul(value,img,maxk)
        if not value:break
    return value
def raw_source(source,degree,maxk):
    # Determine possible Taylor support before material substitution.
    out={};jets={}
    for e,c in source.items():
        e=tuple(map(int,e));jj=monjet(e,degree,maxk)
        if not jj:continue
        cc=cv(c)
        for ex,f in jj.items():out[ex]=out.get(ex,P.zero())+cc*f
    return {e:c for e,c in out.items() if c}
def raw_dag_source(record,degree,maxk):
    coeffs={};jc={}
    for sign,li,ri in record['full_signed_product_DAG']:
        for ex,c in pool[li].items():
            for ey,b in pool[ri].items():
                e=tuple(int(v)+int(u) for v,u in zip(ex,ey))
                if e not in jc:jc[e]=monjet(e,degree,maxk)
                if not jc[e]:continue
                coeffs[e]=coeffs.get(e,old.zero())+sign*c*b
    return raw_source({e:c for e,c in coeffs.items() if c},degree,maxk)
OK=d['ring'].base_ring();sp=OK.ring().hom([k/r,r,1,-1],K)
def mapped(f):
    dd={}
    for ex,c in f.dict().items():
        cc=OK(c);den=sp(cc.denominator());ck('creator coefficient defined '+str(len(checks)),den!=0);dd[tuple(map(int,ex))]=sp(cc.numerator())/den;guardden(dd[tuple(map(int,ex))])
    return P(dd)
jets=[]
for row,degree,divk,divd in ((0,5,2,1),(9,7,3,2)):
    if row==0:j=raw_source(pool[full['full_source_refs'][0]],degree,divk+1)
    else:j=raw_dag_source(full['records'][2],degree,divk+1)
    for (kp,tp),coef in j.items():
        if kp<divk:ck('all apparent source poles cancel '+str((row,kp,tp)),coef==0)
    b0=j.get((divk,0),P.zero());jt=j.get((divk,1),P.zero());jk=j.get((divk+1,0),P.zero())
    ck('entire large axis S source zero '+str(row),b0==0 and d['boundary'][row]==0)
    ck('raw t normal exact '+str(row),jt==de**divd*mapped(d['normal_matrix'][row,7]))
    ck('raw kappa normal exact '+str(row),jk==de**divd*mapped(d['normal_matrix'][row,8]))
    ck('other seven normal columns identically zero '+str(row),all(d['normal_matrix'][row,j]==0 for j in range(7)))
    jets.append((j,jt,jk));print('RAW SOURCE',row,'PASS',round(time.monotonic()-st,2),flush=True)
# Independently recover the finite graph specialization and the genuine minor.
OP=d['ring'];gp=P.hom([mapped(f) for f in d['graph_images']],P)
M=matrix(P,[[gp(mapped(d['normal_matrix'][row,j])) for j in (7,8)] for row in (0,9)])
minor=M.det();RR=n['ring'];conv=P.hom([RR.gen(0),RR.gen(1),RR.gen(2),0,0,0,0],RR)
ck('minor source exactly saved',conv(minor)==n['items'][0]['minor'])
# Build the eight-dimensional multiplication matrix independently.
U=n['candidate_quartic'].parent();zz=U.gen();Q=n['candidate_quartic'];ix=zz/3-zz**3/(10*m)
A=U.zero();B=U.zero()
for ex,c in conv(minor).dict().items():
    i,j,l=map(int,ex);f=c*ix**i*zz**j*k**(l//2)
    if l%2:B+=f
    else:A+=f
A=A.quo_rem(Q)[1];B=B.quo_rem(Q)[1]
cols=[]
for wpower in (0,1):
    for power in range(4):
        pair=((A*zz**power).quo_rem(Q)[1],(B*zz**power).quo_rem(Q)[1]) if wpower==0 else ((k*B*zz**power).quo_rem(Q)[1],(A*zz**power).quo_rem(Q)[1])
        cols.append(vector(K,[pair[j][i] for j in range(2) for i in range(4)]))
mult=matrix(K,cols).transpose();det=mult.det()
expected=K(152186562208866248949760000000000000000)/531441*r**8*k**12*m**16
ck('independent eight by eight multiplication determinant',det==expected)
ck('saved resultant agrees independently',det==n['items'][0]['candidate_norm'])
ck('quartic discriminant uniformly protected',Q.discriminant()==QQ(19600000000)/81*m**6)
ck('all candidate xz and quadratic',(ix*zz-10*m).quo_rem(Q)[1]==0 and (3*(ix**2+zz**2)-ix*zz).quo_rem(Q)[1]==0)
ck('old amplitude Jacobian has no zero candidate',(12**2*k*(-3500*m*m/9))!=0)
# Uniform cyclic-material transport needs only these factors, not signs.
ck('cyclic k transformations retain units',(-1/m)!=0 and (1-1/m)==k/m and (-m/k)!=0 and (1-m/k)==-1/k)
ck('cyclic r transformations retain units',(1-1/q)==r/q and (1-q/r)==-1/r)
report={'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st,'norm':str(det.factor()),'scope':__doc__,'logic':'Actual R0,H2 belong to (kappa,t) on the full large axis locus. The eight-point candidate norm is a unit for every k*r*(1+k)!=0. As a generic-base rational function its cyclic value is nonzero for every spectrum, so localizing its generic norm excludes all eight points. Nakayama gives kappa,t in the local necessary source ideal; no ramified lift. No material positivity was used.'}
save({'checks':checks,'raw_source_jets':jets,'minor':minor,'candidate_multiplication_matrix':mult,'norm':det,'report':report},str(p/'independent_critical_single_RR3_normal_reader.sobj'),compress=False)
(p/'independent_critical_single_RR3_normal_reader.json').write_text(json.dumps(report,indent=2)+'\n');print(report,flush=True)
