"""Independent full-even Newton source bridge, point coefficient rank and line identities."""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(label,v):
    if not bool(v):raise AssertionError(label)
    checks.append(label)
def eq(label,f,g):ck(label,not (f-g))
d=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'));mod=load(str(p/'independent_delta_even_L_module.sobj'));rg=load(str(p/'independent_delta_even_range_growth.sobj'));raw=load(str(p.parent/'round8/rank3_joint_fiber_leading.sobj'))
dirs=load(str(p/'independent_delta_subcritical_even_directions.sobj'));ranks=load(str(p/'independent_delta_subcritical_even_coefficient_rank.sobj'))
line=load(str(p/'independent_subcritical_line_material_free.sobj'));short=load(str(p/'independent_subcritical_line_short_elimination.sobj'))
K=mod['field'];a,r,H,Z=K.gens();q=1+r;mat=1+a*r
R=PolynomialRing(K,names=tuple('iev_p%d'%i for i in range(6))+tuple('iev_R%d'%i for i in range(3)));xs=R.gens()
def cp(f):return R({tuple(map(int,e))+(0,)*(9-len(e)):K(c) for e,c in f.dict().items() if c})
L=matrix(R,3,3,[cp(f) for f in d['L'].list()]);C=matrix(R,3,3,[cp(f) for f in d['C'].list()]);W=matrix(R,3,6,[cp(f) for f in d['W'].list()]);D=matrix(R,6,3,[cp(f) for f in d['D'].list()]);m=cp(d['support_denominator']);ell=L*vector(R,xs[6:]);det=C.det();adj=C.adjugate();Rs=[cp(f) for f in d['full_R_sources']]
comp=load(str(p/'independent_delta_Q0_kernel_compatibility.sobj'));off=vector(R,[cp(f) for f in comp['offset']])
def bi(f,i,j):return R({tuple(map(int,e)):c for e,c in f.dict().items() if sum(e[:6])==i and sum(e[6:])==j})
for j in range(6):eq('complete R31 link '+str(j),bi(Rs[j],3,1),(D*C*ell)[j]);eq('complete R40 link '+str(j),bi(Rs[j],4,0),(D*off)[j])
eq('z highest cancellation',m*det*ell,adj*W*vector(R,[bi(f,3,1) for f in Rs]))
eq('z pure-p8 coefficient',adj*W*vector(R,[bi(f,4,0) for f in Rs]),m*adj*off)
for j,f in enumerate(Rs):
    for ee in f.dict():
        i,k=sum(ee[:6]),sum(ee[6:])
        if (i,k) in ((3,1),(4,0)):continue
        g1=4-i-k;g15=4-i-QQ(3)*k/2
        ck('z lower slot '+str((j,tuple(map(int,ee)))),min(g1,g15)>=0 and max(g1,g15)>0)
# Universal polynomial identity, independent of the source implementation.
Uring=PolynomialRing(QQ,names=('ev_M','ev_U')+tuple('ev_l%d'%j for j in range(3))+tuple('ev_T%d'%j for j in range(3))+tuple('ev_A%d'%j for j in range(3))+('ev_Q',))
zz=Uring.gens();M,U=zz[:2];ll=vector(Uring,zz[2:5]);tv=vector(Uring,zz[5:8]);av=vector(Uring,zz[8:11]);Qsym=zz[11];wv=M*ll-tv
for j in range(3):eq('full necessary source algebra '+str(j),M*M*Qsym-2*U*M*ll[j]*tv[j]+U*tv[j]**2-M*(av*tv),M*M*(Qsym-U*ll[j]**2-av*ll)+U*wv[j]**2+M*(av*wv))
old=raw['ring'].base_ring();ph=old.ring().hom([a,a*q/mat,H/r,q*Z/r],K);cache={}
def cv(c):
    c=old(c)
    if c not in cache:cache[c]=ph(c.numerator())/ph(c.denominator())
    return cache[c]
qdata=[]
for it,rec in zip(mod['items'],rg['records']):
    qi,j=int(rec['Q_index']),int(rec['square_component']);full=R({tuple(map(int,e)):cv(c) for e,c in raw['quartic_source_rows'][qi].dict().items()});AA=vector(R,[cp(f) for f in it['A_p_quadratic']]);unit=K(rec['square_unit'])
    eq('ORIGINAL even square Q'+str(qi),bi(full,2,2),unit*ell[j]**2)
    eq('ORIGINAL even p3RR module Q'+str(qi),bi(full,3,1),AA*ell)
    rem=full-unit*ell[j]**2-AA*ell;q40=bi(full,4,0)
    for ee in (rem-q40).dict():
        i,k=sum(ee[:6]),sum(ee[6:]);g1=4-i-k;g15=4-i-QQ(3)*k/2
        ck('ORIGINAL even strict lower '+str((qi,tuple(map(int,ee)))),min(g1,g15)>=0 and max(g1,g15)>0)
    qdata.append((qi,j,unit,AA,q40))
print('FULL even original-source algebra and strict degrees PASS',round(time.monotonic()-st,2),flush=True)
for ditem,ritem in zip(dirs['items'],ranks['items']):
    branch=int(ditem['branch']);SK=ditem['coefficient_field'];sa,sr,sH,sZ,tt=SK.gens();kh=K.ring().hom([sa,sr,sH,sZ],SK)
    def kc(c):c=K(c);return kh(c.numerator())/kh(c.denominator())
    T=PolynomialRing(SK,names=('iev_A','iev_B','iev_C'));A,B,Ca=T.gens();dx,dz=(SK.one(),SK.one()) if branch==0 else (-(tt+1)/tt,tt);pim=(A*dx,B*dx*dz,A,Ca*dz,B,Ca)
    def pull(f):return sum((T(kc(c))*prod(pim[k]**int(e[k]) for k in range(6)) for e,c in f.dict().items()),T.zero())
    def tc(f):return T({tuple(map(int,e)):SK(c) for e,c in f.dict().items()})
    Cf=matrix(T,3,3,[pull(f) for f in C.list()]);df=vector(T,[pull(f) for f in off]);dc=Cf.det();nc=-Cf.adjugate()*df
    eq('full determinant amplitude factor '+str(branch),dc,A*B*Ca*T(ditem['detC_direction']))
    for j in range(3):eq('full adjugate amplitude factor '+str((branch,j)),nc[j],A*B*Ca*tc(ditem['e_numerators'][j]))
    Es=[]
    for n,(qi,j,unit,AA,q40) in enumerate(qdata):
        ap=vector(T,[pull(f) for f in AA]);fullbracket=kc(unit)*nc[j]**2+dc*(ap*nc)+dc**2*pull(q40)
        E=tc(ditem['even_initial_rows'][n]);eq('GENUINE pure-p16 amplitude link '+str((branch,qi)),fullbracket,(A*B*Ca)**2*E)
        ck('stored E has quartic even amplitudes '+str((branch,qi)),all(sum(e)==4 and all(int(z)%2==0 for z in e) for e in E.dict()))
        Es.append(E)
    # Reconstruct all H,Z coefficient rows in the same recorded order.
    RK=ritem['field'];PR=ritem['polynomial_ring'];mons=list(map(tuple,ritem['amplitude_monomials']));rows=[];labels=[]
    for n,(E,cl) in enumerate(zip(Es,ritem['cleared_sources'])):
        cdict={tuple(map(int,e)):SK(c) for e,c in E.dict().items()};den=SK(cl['denominator']);polys=[]
        for jj,ex in enumerate(mons):
            rawcoef=SK(cdict.get(ex,0))*den;eq('cleared ORIGINAL coefficient '+str((branch,n,jj)),rawcoef,SK(cl['amplitude_coefficients'][jj]));polys.append(SK.ring()(rawcoef))
        by={}
        for jj,f in enumerate(polys):
            for ee,c in f.dict().items():
                e=tuple(map(int,ee));key=e[2:4];by.setdefault(key,[PR.zero() for z in range(6)])[jj]+=PR({(e[0],e[1],e[4]):QQ(c)})
        for key,row in sorted(by.items()):
            if any(row):rows.append(row);labels.append((n,key))
    Mcoef=matrix(RK,rows);ck('ALL coefficient row labels '+str(branch),labels==ritem['labels'])
    for ij,(f,g) in enumerate(zip(Mcoef.list(),ritem['matrix'].list())):eq('ALL coefficient entry '+str((branch,ij)),f,g)
    if branch==0:
        minor=Mcoef.matrix_from_rows(list(map(int,ritem['pivot_rows']))).det();eq('point actual 6minor',minor,ritem['selected_determinant'])
        aa,rr,tx=RK.gens();expected=175693766054100480*aa**8*rr**10*(1+aa*rr)**8;eq('point all-material unit factor',minor,expected)
    print('EVEN branch',branch,'ALL source/coefficient PASS',round(time.monotonic()-st,2),flush=True)
# Reconstruct the line material transport from the audited coefficient rows.
li=ranks['items'][1];LK=li['field'];la,lr,lt=LK.gens();LP=PolynomialRing(LK,names=('im_A','im_B','im_C'));A,B,Ca=LP.gens();m1=1+la*lr
material_map=line['material_free_ring'].hom([LP(lt),la*lr*A**2,la*lr*m1*lt**2*B**2,m1*lt**2*Ca**2],LP)
for j,k in enumerate(line['source_rows']):
    source=sum((LP(c)*A**e[0]*B**e[1]*Ca**e[2] for e,c in zip(li['amplitude_monomials'],li['matrix'].row(int(k)))),LP.zero())
    eq('line complete material row '+str(k),source,LK(line['source_material_units'][j])*material_map(line['material_free_forms'][j]))
    c=LK(line['source_material_units'][j]);allowed=LK.ring()(la*lr*(1+la*lr)*lt*(lt+1))
    for f in (c.numerator(),c.denominator()):
        for fac,e in f.factor():
            if fac.total_degree()>0:ck('line row unit '+str((j,fac)),allowed.quo_rem(fac)[1]==0)
# Only use short exact source identities. Do not rerun or trust a GB display.
SR=PolynomialRing(QQ,names=('ish_t','ish_u','ish_v'));t,u,v=SR.gens();h=line['material_free_ring'].hom([t,u,v,SR.one()],SR);fs=[h(f) for f in line['material_free_forms']]
def sc(f):return SR({tuple(map(int,e)):QQ(c) for e,c in f.dict().items()})
for j,f in enumerate(fs):eq('short source linked to conic '+str(j),f,sc(short['original_conics'][j]))
for j,it in enumerate(short['original_conic_lifts']):
    target=sc(it['target']);mul=sc(it['support_multiplier']);ck('short multiplier exactly v '+str(j),mul==v)
    eq('EXPLICIT conic source identity '+str(j),sum((sc(c)*f for c,f in zip(it['cofactors'],fs)),SR.zero()),mul*target)
    eq('short consequence target '+str(j),target,sc(short['short_consequences'][j]))
Sym=PolynomialRing(QQ,names=('is_s','is_d'));s,dd=Sym.gens();hs=SR.hom([Sym.zero(),(s+dd)/2,(s-dd)/2],Sym)
cubic=(dd-5)*s**2+(3*dd+1)*(2*dd-1);eq('symmetric necessary cubic',hs(sc(short['short_consequences'][3])),cubic)
record={'checks':len(checks),'seconds':time.monotonic()-st,'point_minor':'175693766054100480*a^8*r^10*(1+a*r)^8','correct_full_amplitude_factor':'m^2*(ABC)^2*stored_E_i','scope':'Original necessary even source and unique pure-p16 Newton initial for 1<gamma<3/2. x=z=1 impossible on the A B C nonzero constant-direction chart. Hyperbola xz+z+1=0 retains three material-free conics and their four displayed consequences; no exclusion.'}
save({'checks':checks,'record':record,'universal_source_identity_verified':True},str(p/'independent_subcritical_even_line_reader.sobj'),compress=False)
(p/'independent_subcritical_even_line_reader.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2),flush=True)
