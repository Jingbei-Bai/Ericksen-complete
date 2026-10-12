"""Independent generic-base inactive-row support/source certificate replay."""
from sage.all import *
from pathlib import Path
from itertools import permutations
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(n,v):
    if not bool(v):raise AssertionError(n)
    checks.append(n)
gc=load(str(p/'independent_ratio_inactive_generic_compatibility.sobj'));d=load(str(p/'independent_ratio_inactive_row3_linear_strata.sobj'));qd=load(str(p/'independent_ratio_inactive_row3_RR_derivatives.sobj'));ad=load(str(p/'independent_ratio_inactive_fullsupport_amplitudes.sobj'));g=load(str(p/'independent_ratio_inactive_fullsupport_graph.sobj'));mat=load(str(p/'independent_ratio_inactive_amplitude_material.sobj'));mins=load(str(p/'independent_ratio_inactive_material_minors.sobj'));cert=load(str(p/'independent_ratio_inactive_fullsupport_base_certificate.sobj'));small=load(str(p/'independent_ratio_inactive_small_support_certificate.sobj'));three=load(str(p/'independent_ratio_inactive_three_p_RRzero_certificate.sobj'))
K=d['coefficient_field'];a,r,H,Z=K.gens();q=1+r;m=1+a*r
R=PolynomialRing(K,names=tuple('is_p%d'%i for i in range(6))+tuple('is_R%d'%i for i in range(3)));v=R.gens()
def cp(f):return R({tuple(map(int,e)):K(c) for e,c in f.dict().items() if c})
def sd(f,j):
    ans={}
    for ee,c in f.dict().items():
        e=tuple(map(int,ee))
        if e[j]:ne=list(e);ne[j]-=1;ans[tuple(ne)]=c*e[j]
    return f.parent()(ans)
zs=[cp(f) for f in d['Z']];tt=cp(d['T']);Q={int(it['source_index']):cp(it['restricted_initial']) for it in qd['items']};Ss={it['source'][1]:cp(it['original_initial']) for it in ad['items'] if it['source'][0]=='S'}
orig=load(str(p.parent/'round8/rank3_joint_fiber_leading.sobj'));pairs=load(str(p/'independent_ratio_p_boundary_global_pair_syzygies.sobj'));OK=orig['ring'].base_ring();coeffmap=OK.ring().hom(qd['original_parameter_map'],K);cache={}
def cv(c):
    c=OK(c)
    if c not in cache:cache[c]=coeffmap(c.numerator())/coeffmap(c.denominator())
    return cache[c]
def origpull(f,degree):return R({tuple(map(int,e)):cv(c) for e,c in f.dict().items() if sum(map(int,e))==degree and not e[4] and not e[5]})
for j in (0,1,5,6):ck('ORIGINAL ordinary Q%d'%j,Q[j]==origpull(orig['quartic_source_rows'][j],3 if j==0 else 4))
for it in pairs['items']:ck('ORIGINAL full pair initial %s'%it['pair'],Ss[it['pair']]==origpull(it['ordinary_initial'],4))
def aug(f):
    c0,c1=sd(f,6),sd(f,7);last=f-c0*v[6]-c1*v[7]
    ck('augmented_only_two_RR_'+str(len(checks)),all(not any(e[6:]) for z in (c0,c1,last) for e in z.dict()))
    return [c0,c1,last]
M=matrix(R,[aug(f) for f in zs]);MT=matrix(R,[aug(zs[1]),aug(zs[2]),aug(tt)])
ck('allZ original matrix',M==matrix(R,3,3,[cp(f) for f in gc['all_Z_augmented_matrix'].list()]));ck('ZT original matrix',MT==matrix(R,3,3,[cp(f) for f in gc['ZT_augmented_matrix'].list()]));ck('allZ determinant',M.det()==cp(gc['all_Z_determinant']));ck('ZT determinant',MT.det()==cp(gc['ZT_determinant']))
T=PolynomialRing(K,'is_t');t=T.gen();toT=R.hom([1,1,t,1,0,0,0,0,0],T)
def ct(f):return T({int(e):K(c) for e,c in enumerate(f.list()) if c})
A=matrix(T,2,2,[toT(M[i,j]) for i in (1,2) for j in (0,1)]);bvec=vector(T,[-toT(M[i,2]) for i in (1,2)]);D=A.det();ns=A.adjugate()*bvec
ck('active graph original coefficient matrix',A==matrix(T,2,2,[ct(f) for f in g['matrix'].list()]));ck('graph determinant',D==ct(g['determinant']));ck('graph Cramer all',ns==vector(T,[ct(f) for f in g['RR_ratio_numerators']]));ck('graph identity no division',A*ns==D*bvec)
for label,row,saved in [('missing',M.row(0),g['missing_Z0_numerator']),('T',MT.row(2),g['T_numerator'])]:ck('graph compatibility '+label,toT(row[0])*ns[0]+toT(row[1])*ns[1]+toT(row[2])*D==ct(saved))
# Reconstruct the original full-amplitude polynomial numerators independently.
AP=ad['polynomial_ring'];aa,ar,AH,AZ,at=AP.gens();KR=K.ring()
def liftcoef(c):return AP({tuple(map(int,e))+(0,):QQ(x) for e,x in KR(c).dict().items()})
DP,N0P,N1P=ad['graph_polynomials']
ap_to_T=AP.hom([T(a),T(r),T(H),T(Z),t],T)
DPt,N0t,N1t=map(ap_to_T,(DP,N0P,N1P));graphscale=D.leading_coefficient()/DPt.leading_coefficient()
ck('graph polynomial D linked to original',D==graphscale*DPt)
ck('graph polynomial N0 linked to original',ns[0]==graphscale*N0t)
ck('graph polynomial N1 linked to original',ns[1]==graphscale*N1t)
cr=cert['ring'];cert_to_T=cr.hom([T(H),T(Z),t],T)
for saved,factor in zip((g['missing_Z0_numerator'],g['T_numerator']),cert['direction_polynomials']):
    rat=ct(saved)/((t+1)*cert_to_T(factor));ck('exact direction factor source '+str(len(checks)),rat.denominator()==1 and rat.numerator().degree()==0)
ck('certificate selected minor is original',cert['selected_minor']*cert['selected_minor_constant_scale']==mins['minors'][0]['polynomial'])
for j,f0 in enumerate((DPt,N0t,N1t)):
    ck('primitive graph certificate '+str(j),cert_to_T(cert['graph_D_N'][j])*cert['graph_constant_scales'][j]==f0)
for rowid in (1,2,5,6):
    it=ad['items'][rowid];f=Q[it['source'][1]] if it['source'][0]=='Q' else Ss[it['source'][1]];den=it['coefficient_denominator'];computed={}
    for ee,c in f.dict().items():
        e=tuple(map(int,ee));key=(e[1]+e[7],e[3]+e[6],e[8]);val=liftcoef(K(c*den))*at**e[2]*N0P**e[6]*N1P**e[7]*DP**(2-e[6]-e[7]);computed[key]=computed.get(key,AP.zero())+val
    computed={e:z for e,z in computed.items() if z};saved=it['cleared_amplitude_coefficients'];ck('original amplitude source '+str(rowid),computed==saved)
    co=it['common_polynomial_content'];keys=[(0,0,0),(2,0,0),(0,2,0)] if rowid<5 else ([(0,1,0),(2,1,0),(0,3,0)] if rowid==5 else [(1,0,0),(3,0,0),(1,2,0)])
    ii=(1,2,5,6).index(rowid)
    for jj,e in enumerate(keys):ck('exact content quotient %d %d'%(ii,jj),computed.get(e,AP.zero())==co*mat['matrix'][ii,jj])
    ck('no hidden amplitude coefficient '+str(rowid),set(computed)<=set(keys))
    print('ORIGINAL AMPLITUDE',rowid,'PASS',round(time.monotonic()-st,2),flush=True)
NR=mins['ring'];nh,nz,nt=NR.gens();N=mins['matrix'];nmap=NR.hom([AH,AZ,at],AP);col=(aa*ar,aa*ar*(1+aa*ar),1+aa*ar)
for i in range(4):
    for j in range(3):ck('coherent material matrix %d %d'%(i,j),mat['matrix'][i,j]==nmap(N[i,j])*col[j])
for it in mins['minors']:ck('actual maximal minor '+str(it['rows']),N.matrix_from_rows(it['rows']).det()==it['polynomial'])
# Exact lazy Sylvester definitions: no expanded large determinant required.
def coeffs_t(f):
    B=PolynomialRing(QQ,names=('syl_H','syl_Z'));rows={}
    for ee,c in f.dict().items():
        e=tuple(map(int,ee));rows[e[2]]=rows.get(e[2],B.zero())+QQ(c)*B.gen(0)**e[0]*B.gen(1)**e[1]
    return B,rows
for item in cert['resultant_certificates']:
    ff,gg=item['polynomials'];BB,fc=coeffs_t(ff);_,ggc=coeffs_t(gg);nf=max(fc);ng=max(ggc);size=nf+ng
    if item['name']=='graph_pivot':ck('pivot resultant is actual graph pair',(ff,gg)==tuple(cert['graph_D_N'][:2]))
    else:ck('source resultant is actual necessary pair',ff==prod(cert['direction_polynomials']) and gg==cert['selected_minor'])
    A=matrix(BB,size,size)
    for colid in range(ng):
        for e,c in fc.items():A[e+colid,colid]=c
    for colid in range(nf):
        for e,c in ggc.items():A[e+colid,ng+colid]=BB(c)
    for i in range(size):
        for j in range(size):ck('Sylvester exact %s %s %s'%(item['name'],i,j),A[i,j].dict()==item['sylvester'][i,j].dict())
    bc=BB.hom([1,-1],QQ);Ac=A.apply_map(bc);ck('cyclic Sylvester determinant '+item['name'],Ac.det()==item['cyclic_determinant'] and Ac.det()!=0)
    U=PolynomialRing(QQ,'tc');u=U.gen();ffu=sum(bc(c)*u**e for e,c in fc.items());ggu=sum(bc(BB(c))*u**e for e,c in ggc.items());ck('independent univariate resultant '+item['name'],abs(ffu.resultant(ggu))==abs(Ac.det()))
# Track every potentially canceled factor.  Any extra base factor must be
# material-free; pure material factors must divide the original spectral guard.
guardmat=K.ring()(a*r*q*(a-1)*m);baseexceptions=[]
def factors_allowed(poly,label,allow_t=False):
    for fac,power in poly.factor():
        if fac in QQ:continue
        # AP ring has ordered a,r,H,Z,t; classify coefficients exactly.
        ed=[tuple(map(int,e)) for e in fac.dict()]
        if all(e[2:]==(0,0,0) for e in ed):
            pf=K.ring()({tuple(map(int,e[:2]))+(0,0):c for e,c in fac.dict().items()});ck(label+' spectral '+str(fac),guardmat.quo_rem(pf)[1]==0)
        elif all(e[0]==e[1]==0 for e in ed):
            if fac.degree(at):ck(label+' only balanced t wall '+str(fac),allow_t and fac.degree(at)==1 and fac.subs({at:-1})==0)
            else:baseexceptions.append(str(fac));ck(label+' nonzero material-free base',bool(fac))
        else:raise AssertionError('unaccounted mixed material wall '+label+str(fac))
for i,co in enumerate(mat['row_contents']):factors_allowed(co,'rowcontent '+str(i),True)
for it in cert['balanced_original_sources'][5:7]:
    for e,poly in it['coefficients'].items():factors_allowed(poly,'balanced original '+str(it['source']))
# The actual t=-1 source is a unit times J^2, with the graph determinant
# excluding J=0.  Check all factors explicitly rather than support alone.
J=cert['balanced_J'];CR=cert['ring'];ch,cz,ctt=CR.gens();cm=CR.hom([AH,AZ,at],AP);jap=cm(J)
for it in cert['balanced_original_sources'][5:7]:
    poly=next(iter(it['coefficients'].values()));quo,rem=poly.quo_rem(jap**2);ck('balanced J square factor '+str(it['source']),rem==0 and bool(quo));factors_allowed(quo,'balanced remaining units')
print('FULLSUPPORT RESULTANT CHAIN PASS',round(time.monotonic()-st,2),flush=True)
def rest(f,mask):return R({tuple(map(int,e)):c for e,c in f.dict().items() if all(mask>>j&1 or not e[j] for j in range(4))})
def base_unit_or_exception(c,label):
    c=K(c)
    for typ,poly in [('num',c.numerator()),('den',c.denominator())]:
        poly=K.ring()(poly)
        for fac,power in poly.factor():
            if fac in QQ:continue
            ed=[tuple(map(int,e)) for e in fac.dict()]
            if all(not e[2] and not e[3] for e in ed):ck(label+typ+' old spectral '+str(fac),guardmat.quo_rem(fac)[1]==0)
            elif all(not e[0] and not e[1] for e in ed):baseexceptions.append(str(fac));ck(label+typ+' base factor nonzero',bool(fac))
            else:raise AssertionError('mixed material exception '+label+str(fac))
def monomial_nonzero(F,mask,label):
    F=R(F);ck(label+' one monomial',len(F.dict())==1);e,c=next(iter(F.dict().items()));ck(label+' supported nonzero coordinates',all(not e[j] or (j<4 and (mask>>j)&1) for j in range(9)));base_unit_or_exception(c,label)
base_unit_or_exception(K(graphscale),'graph polynomial scaling')
for j,c in enumerate(cert['direction_prefactors']):base_unit_or_exception(K(c),'direction source prefactor '+str(j))
for it in small['pure_monomial_exclusions']:
    source=rest(Q[it['source'][1]],it['mask']);saved,mon,co=it['identity'];ck('pure support '+str(it['mask']),source==cp(saved)==cp(mon)*K(co));monomial_nonzero(source,it['mask'],'pure support '+str(it['mask']))
for it in small['chain_pair_exclusions']:
    mask=it['mask'];rp=it['RR_position'];f=cp(it['active_Z']);ck('chain active original '+str(mask),f in [rest(z,mask) for z in zs[1:]]);ck('chain RR forcing '+str(mask),f==v[rp]*cp(it['active_coefficient']))
    c0,c1=map(cp,it['RR_coefficient_pair']);f0,f1=map(cp,it['affine_factors']);ck('chain exact cancellation '+str(mask),c0*f1-c1*f0==cp(it['difference']))
    monomial_nonzero(it['active_coefficient'],mask,'chain active coefficient '+str(mask));monomial_nonzero(it['difference'],mask,'chain contradiction '+str(mask))
    im=list(v);im[rp]=R.zero();ph=R.hom(im,R)
    for kk,(qq,saved) in enumerate(zip((5,6) if mask==6 else (1,6),it['original_Q_after_RR0'])):
        actual=ph(rest(Q[qq],mask));ck('chain original Q '+str((mask,qq)),actual==cp(saved));aff=(f0,f1)[kk];seen=False
        for fac,power in actual.factor():
            if fac==aff*(fac.lc()/aff.lc()):seen=True
            else:monomial_nonzero(fac,mask,'chain harmless factor '+str((mask,qq)))
        ck('chain affine factor actually occurs '+str((mask,qq)),seen)
pair=small['opposite_pair'];ck('opposite S1',rest(Ss[1],5)==cp(pair['S1']));ck('opposite S2',rest(Ss[2],5)==cp(pair['S2']))
forced=[]
for key in ('S1','S2'):
    F=cp(pair[key]);ck('opposite forcing single monomial '+key,len(F.dict())==1);ee,cc=next(iter(F.dict().items()));active=[j for j in (6,7,8) if ee[j]];ck('opposite forcing one RR '+key,len(active)==1 and ee[active[0]]==1);forced+=active;monomial_nonzero(F/v[active[0]],5,'opposite RR coefficient '+key)
ck('both opposite RR components forced',set(forced)=={6,7})
ph=R.hom(list(v[:6])+[0,0,v[8]],R)
for j,factor in ((1,v[0]*v[2]*(v[0]+v[2])**2),(6,v[0]*v[2]*(v[0]**2+v[2]**2))):
    ff=ph(rest(Q[j],5));ck('opposite original Q '+str(j),ff==cp(pair['Q%d_reduced'%j]));ratio=ff/factor;ck('opposite final quadratic '+str(j),ratio in K and bool(ratio));base_unit_or_exception(K(ratio),'opposite Q factor '+str(j))
PU=PolynomialRing(QQ,'op_t');ot=PU.gen();ck('opposite quadratic resultant',((ot+1)**2).resultant(ot**2+1)==4)
for it in small['endpoint_three_supports']:
    for k,key in enumerate(('all_Z_determinant','ZT_determinant')):
        ck('three support determinant '+str((it['mask'],k)),rest(cp(gc[key]),it['mask'])==cp(it['determinants'][k]));monomial_nonzero(it['determinants'][k],it['mask'],'endpoint determinant '+str((it['mask'],k)))
incoming=small['shifted_incoming'];U=incoming['ring'];ll,uu,vv=U.gens();ph=R.hom([0,1,0,ll,0,0,incoming['RR1'],uu,vv],U)
for label,source in [('Z1',zs[1]),('Z2',zs[2]),('T',tt)]+[('Q'+str(k),z) for k,z in Q.items()]+[('S'+str(k),z) for k,z in Ss.items()]:ck('shifted incoming original '+label,ph(source)==0)
F10=rest(Q[6],10);T10=rest(tt,10);ratio=F10/T10**2;ck('incoming Q6 equivalent square',ratio in K and bool(ratio));base_unit_or_exception(K(ratio),'incoming square coefficient')
for mask,qi,rp,pi in ((1,5,6,0),(2,6,6,1),(4,1,7,2),(8,6,7,3)):
    pure=rest(Q[qi],mask);ratio=pure/(v[pi]**2*v[rp]**2);ck('singleton necessary RR square '+str(mask),ratio in K and bool(ratio));base_unit_or_exception(K(ratio),'singleton coefficient '+str(mask))
    im=list(v);im[rp]=0;ph=R.hom(im,R)
    for label,source in list(Q.items())+[(10+k,z) for k,z in Ss.items()]:ck('singleton source residual '+str((mask,label)),ph(rest(source,mask))==0)
for it in three['items']:
    t3=it['linear_polynomials'][0].parent().gen();imgs=[1,1,t3,0,0,0,1,0,0] if it['mask']==7 else [1,0,t3,1,0,0,0,1,0];ph=R.hom(imgs,t3.parent());f0,f1=[ph(z) for z in zs[1:]]
    ck('three p coefficient sources '+str(it['mask']),[f0,f1]==it['linear_polynomials']);ck('three p resultant '+str(it['mask']),f0.resultant(f1)==it['resultant'] and bool(it['resultant']));base_unit_or_exception(K(it['resultant']),'three p resultant '+str(it['mask']))
inventory={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};ck('complete support inventory',inventory==set(range(1,16)))
save({'checks':checks,'base_exception_factors':baseexceptions,'scope':'Fullsupport ordinary-initial exclusion on a material-uniform generic base; small-support source inventory and shifted incoming residual graph. Masks7/13 unit exclusion belongs to root, not this audit.'},str(p/'independent_inactive_generic_support_reader.sobj'),compress=False)
report={'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st,'base_factors':sorted(set(baseexceptions)),'scope':'No entire inactive-D or finite-fiber conclusion; shifted mask10 and singleton axes retained.'};(p/'independent_inactive_generic_support_reader.json').write_text(json.dumps(report,indent=2)+'\n');print(report,flush=True)
