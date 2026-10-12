"""Independent original-source replay of the subcritical Newton Q0 relation.

This is a source/degree/coefficient-ideal audit, not a whole-growth exclusion.
No Groebner basis is computed and no determinant-C inverse is introduced.
"""
from sage.all import *
from pathlib import Path
import time, json
p=Path(__file__).resolve().parent;start=time.monotonic();checks=[]
def ck(name,truth):
    if not bool(truth):raise AssertionError(name)
    checks.append(name)
def eq(name,f,g):ck(name,not (f-g))
d=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'))
di=load(str(p/'independent_delta_subcritical_Q0_direction.sobj'))
bc=load(str(p/'independent_delta_subcritical_Q0_base_coefficients.sobj'))
fall=load(str(p/'independent_ratio_weighted_degreefall.sobj'))
raw=load(str(p.parent/'round8/rank3_joint_fiber_leading.sobj'))
rank=load(str(p/'independent_ratio_p_boundary_rank_structure.sobj'))
K=rank['coefficient_field'];a,r,H,Z=K.gens();q=1+r;mm=1+a*r
R=PolynomialRing(K,names=tuple('isub_p%d'%j for j in range(6))+tuple('isub_v%d'%j for j in range(3)))
vs=R.gens();rr=vector(R,vs[6:]);old=raw['ring'].base_ring()
ph=old.ring().hom([a,a*q/mm,H/r,q*Z/r],K);cache={}
def cv(c):
    c=old(c)
    if c not in cache:cache[c]=ph(c.numerator())/ph(c.denominator())
    return cache[c]
def orig(f):return R({tuple(map(int,e)):cv(c) for e,c in f.dict().items() if c})
def cp(f):
    return R({tuple(map(int,e))+(0,)*(9-len(e)):K(c) for e,c in f.dict().items() if c})
Q=orig(raw['quartic_source_rows'][0]);eq('ORIGINAL full Q0',Q,cp(d['full_Q0']))
Rs=[]
for j,it in enumerate(fall['items']):
    idx=int(it['source_index']);oldidentity=raw['quartic_source_rows'][idx]-it['cofactor']*raw['quartic_source_rows'][0]
    eq('ORIGINAL complete degreefall identity %d'%idx,oldidentity,it['full'])
    f=orig(oldidentity);Rs.append(f);eq('ORIGINAL full R payload %d'%idx,f,cp(d['full_R_sources'][j]))
print('ORIGINAL full Q/R PASS',round(time.monotonic()-start,2),flush=True)
L=matrix(R,3,3,[cp(f) for f in rank['L_RR'].list()]);C=matrix(R,3,3,[cp(f) for f in rank['linear_matrix_C'].list()]);off=vector(R,[cp(f) for f in rank['pure_cubic_offset']])
for label,M,saved in [('L',L,d['L']),('C',C,d['C'])]:
    for j,(f,g) in enumerate(zip(M.list(),saved.list())):eq('original '+label+str(j),f,cp(g))
D=matrix(R,6,3);D[0,0]=vs[4];D[1,0]=-vs[5];D[2,1]=vs[2];D[5,1]=-vs[3];D[3,2]=vs[0];D[4,2]=vs[1]
m=vs[4]*vs[2]*vs[0];W=matrix(R,3,6);W[0,0]=vs[2]*vs[0];W[1,2]=vs[4]*vs[0];W[2,3]=vs[4]*vs[2]
for j,(f,g) in enumerate(zip((W*D).list(),(m*identity_matrix(R,3)).list())):eq('support inverse '+str(j),f,g)
for label,M,saved in [('D',D,d['D']),('W',W,d['W'])]:
    for j,(f,g) in enumerate(zip(M.list(),saved.list())):eq('source matrix '+label+str(j),f,cp(g))
def bi(f,i,j):return R({tuple(map(int,e)):c for e,c in f.dict().items() if sum(e[:6])==i and sum(e[6:])==j})
def bis(f):return {(sum(e[:6]),sum(e[6:])) for e in f.dict()}
ell=L*rr
# Recover B directly from original Q0, not the saved range equation.
B=[]
for pj,rj in ((0,0),(1,0),(2,1)):
    dd={}
    for ee,c in bi(Q,1,2).dict().items():
        e=list(map(int,ee))
        if e[pj]:
            ck('B extraction monomial '+str((pj,ee)),e[pj]==1 and sum(e[:6])==1 and e[6+rj]>0)
            e[pj]-=1;e[6+rj]-=1;dd[tuple(e)]=c
    B.append(R(dd))
B=vector(R,B)
for j,f in enumerate(B):eq('original B payload '+str(j),f,cp(d['Q0_L_cofactors'][j]))
eq('ORIGINAL Q0 top module',bi(Q,1,2),B*ell)
for j,f in enumerate(Rs):
    eq('ORIGINAL R31 module '+str(j),bi(f,3,1),(D*C*ell)[j])
    eq('ORIGINAL R40 module '+str(j),bi(f,4,0),(D*off)[j])
for label,fs,w in [('C',C.list(),(1,0)),('d',off,(3,0)),('L',L.list(),(1,0)),('B',B,(0,1)),('W',W.list(),(2,0))]:
    for j,f in enumerate(fs):ck('homogeneous '+label+str(j),not f or bis(f)=={w})
det=C.det();adj=C.adjugate();mult=list(B*adj*W)
ck('first source multiplier homogeneous',bis(m*det)=={(6,0)})
for j,f in enumerate(mult):ck('R source multiplier homogeneous '+str(j),not f or bis(f)=={(4,1)})
eq('complete (7,2) cancellation',m*det*bi(Q,1,2),sum(mult[j]*bi(Rs[j],3,1) for j in range(6)))
P81=m*det*bi(Q,2,1)-sum(mult[j]*bi(Rs[j],4,0) for j in range(6))
eq('ORIGINAL next source (8,1)',P81,m*(det*bi(Q,2,1)-B*adj*off))
eq('saved source initial replay',P81,cp(d['subcritical_initial']))
ck('initial nonzero and homogeneous',P81!=0 and bis(P81)=={(8,1)})
slots=set()
for i,j in bis(Q):slots.add((i+6,j))
for f in Rs:
    for i,j in bis(f):slots.add((i+4,j+1))
ck('full support inventory',slots==set(map(tuple,d['all_full_DAG_bidegrees'])))
for ij in sorted(slots-{(7,2),(8,1)}):
    i,j=ij;g1=9-i-j;g15=QQ(19)/2-i-QQ(3)*j/2
    ck('STRICT every 1<gamma<3/2 '+str(ij),min(g1,g15)>=0 and max(g1,g15)>0)
print('Complete multiplier, layer, and strict-degree audit PASS',round(time.monotonic()-start,2),flush=True)
# Do not put amplitudes into the coefficient field; keep all three polynomial.
SK=di['field'];sa,sr,sH,sZ,sx,sz=SK.gens();kh=K.ring().hom([sa,sr,sH,sZ],SK)
def kc(c):c=K(c);return kh(c.numerator())/kh(c.denominator())
T=PolynomialRing(SK,names=('isub_A','isub_B','isub_C','isub_rho'));A,Ba,Ca,rho=T.gens()
pimgs=(A*sx,Ba*sx*sz,A,Ca*sz,Ba,Ca)
alpha=(sZ-sH)/(sZ*(sH+1));beta=sH*(sZ-1)/(sZ*(sH+1))
vimgs=(Ca*alpha,Ba*beta*sx,A*sx)
def pull(f,images):
    return sum((T(kc(c))*prod(images[j]**int(e[j]) for j in range(9)) for e,c in f.dict().items()),T.zero())
def pp(f):return pull(f,pimgs+(T.zero(),)*3)
Cf=matrix(T,3,3,[pp(f) for f in C.list()]);df=vector(T,[pp(f) for f in off]);Lf=matrix(T,3,3,[pp(f) for f in L.list()])
for j,(f,g) in enumerate(zip(Cf.list(),(matrix(T,di['C'])*diagonal_matrix(T,[Ba,A,Ca])).list())):eq('ALL amplitudes C '+str(j),f,g)
for j in range(3):eq('ALL amplitudes d '+str(j),df[j],A*Ba*Ca*T(di['offset'][j]));eq('ALL amplitudes RR kernel '+str(j),(Lf*vector(T,vimgs))[j],0)
actual=pull(P81,pimgs+tuple(rho*vj for vj in vimgs));expected=rho*(A**2*Ba*sx)*(A*Ba*Ca)**2*T(di['condition'])
eq('ALL amplitudes and arbitrary RR scale source pull',actual,expected)
# This equality contains no C determinant inverse; verify by direct adjugate.
directions=(sx,sx*sz,SK.one(),sz,SK.one(),SK.one())
def direct(f):return sum((kc(c)*prod(directions[j]**int(e[j]) for j in range(6)) for e,c in f.dict().items()),SK.zero())
Cd=matrix(SK,3,3,[direct(f) for f in C.list()]);dd=vector(SK,[direct(f) for f in off]);vv=vector(SK,[alpha,beta*sx,sx]);allimgs=directions+tuple(vv)
directQ=SK.zero()
for e,c in bi(Q,2,1).dict().items():directQ+=kc(c)*prod(allimgs[j]**int(e[j]) for j in range(9))
Br=vector(SK,[sum((kc(c)*prod(allimgs[j]**int(e[j]) for j in range(9)) for e,c in f.dict().items()),SK.zero()) for f in B])
condition=Cd.det()*directQ-Br*Cd.adjugate()*dd
eq('fraction-free original direction condition',condition,di['condition'])
eq('factorization full SAVED numerator multiplication',di['factorization'].value(),di['numerator'])
eq('factorization original condition cross product',SK(di['factorization'].value()),condition*SK(di['denominator']))
F=next(f for f,e in di['factorization'] if f.total_degree()==8)
eq('selected F8 is actual source factor',F,bc['original_factor'])
rest=condition/SK(F)/ (sx*sz+sz+1)
baseguard=SK(sa*sr*(1+sr)*(sa-1)*(1+sa*sr)*sH*sZ*(sH-sZ)*(sH+1)*(sZ-1)*(1+sH-sZ)*(sH*sZ-sH+sZ)).numerator()
direction_factors=[]
for what,f in [('source remaining numerator',rest.numerator()),('source denominator',rest.denominator())]:
    for ff,e in f.factor():
        if ff.total_degree()==0:continue
        if ff==ff.parent()(sx):direction_factors.append((str(ff),int(e)));continue
        ck(what+' guard '+str(ff),baseguard.quo_rem(ff)[1]==0)
ck('no material in F8',all(not e[0] and not e[1] for e in F.dict()))
CR=PolynomialRing(QQ,names=('isub_x','isub_z'));xx,zz=CR.gens();by={}
for ee,c in F.dict().items():
    e=tuple(map(int,ee));k=e[2:4];by[k]=by.get(k,CR.zero())+CR({e[4:6]:QQ(c)})
def cc(f):return CR({tuple(map(int,e)):QQ(c) for e,c in f.dict().items()})
ck('ALL HZ coefficient keys',set(by)==set(bc['base_coefficients']))
for k,f in by.items():eq('original F8 HZ coefficient '+str(k),f,cc(bc['base_coefficients'][k]));ck('reverse coefficient ideal '+str(k),f(1,1)==0)
rows=[by[k] for k in bc['base_coefficients']]
for j,(target,cofs) in enumerate(zip((xx-1,zz-1),bc['cofactors'])):eq('EXPLICIT coefficient Bezout '+str(j),sum((cc(c)*f for c,f in zip(cofs,rows)),CR.zero()),target)
record={'checks':len(checks),'seconds':time.monotonic()-start,'source_DAG':d['source_DAG'],'strict_gamma_interval':'1 < gamma < 3/2','bidegrees':[list(map(int,t)) for t in sorted(slots)],'remaining_source_factor':str(rest.factor()),'F8_degree':int(F.total_degree()),'F8_terms':len(F.dict()),'coefficient_ideal':['x-1','z-1'],'scope':'Source relation on the A B C x rho nonzero transitive Delta chart. Keeps xz+z+1=0 and all detC=0. Does not exclude a whole growth branch.'}
save({'checks':checks,'record':record,'direction_condition':condition,'F8':F,'remaining_factor':rest,'coefficient_rows':by},str(p/'independent_subcritical_Q0_source_reader.sobj'),compress=False)
(p/'independent_subcritical_Q0_source_reader.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2),flush=True)
