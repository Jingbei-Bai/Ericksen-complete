"""Both endpoint orientations, original-source and exact projection audit."""
from pathlib import Path
from itertools import combinations
import time,json
from sage.all import load,save,QQ,PolynomialRing,matrix,vector,prod

out=Path(__file__).resolve().parent;start=time.monotonic();checks={};cert=[]
norm=load(str(out/'general_curve_homogeneous_normalized.sobj'))
lead=load(str(out/'general_curve_dimensionless_leading.sobj'))
R=norm['ring'];K=R.base_ring();C=K.ring();H,Z,S,N=K.gens();J=matrix(K,[[0,-1],[1,0]])
guards=[C(H),C(Z),C(H-Z),C(Z-1),C(H+1),C(H-Z+1),C(S-1),C(H*S-Z+1),C((H+1)*S-Z)]
A=PolynomialRing(K,names=('l','k'));ll,kk=A.gens()
P=PolynomialRing(QQ,names=('H','Z','S'));pH,pZ,pS=P.gens();PK=P.fraction_field()
B=PolynomialRing(QQ,names=('H','Z'));bH,bZ=B.gens()
pg=[bH,bZ,bH-bZ,bZ-1,bH+1,bH-bZ+1]
def isunit(c):
    if not c:return False
    for pol in [c.numerator(),c.denominator()]:
        for f,e in C(pol).factor():
            f=C(f)
            if f.is_constant():continue
            if not any(all((K(f)/g).derivative(v)==0 for v in C.gens()) for g in guards):return False
    return True
def to_plane(p):return B(p(bH,bZ,0))
def sourcecoeffs(p):
    cc=[B.zero()]*(p.degree(pS)+1)
    for mon,c in p.dict().items():cc[mon[2]]+=c*bH**mon[0]*bZ**mon[1]
    return cc
def rem(f,g):return f.quo_rem(g)[1]
def contact_eval(cs,num,den,g):
    np=[B.one()];dp=[B.one()]
    for _ in range(len(cs)-1):np.append(rem(np[-1]*num,g));dp.append(rem(dp[-1]*den,g))
    ans=B.zero()
    for j,c in enumerate(cs):ans=rem(ans+rem(c*np[j],g)*dp[len(cs)-1-j],g)
    return ans

for reverse in [False,True]:
    prefix='general_curve_two_pair_reverse' if reverse else 'general_curve_two_pair_endpoint'
    ep=load(str(out/('general_curve_two_pair_reverse_endpoint.sobj' if reverse else 'general_curve_two_pair_binary_minors.sobj')))
    st=load(str(out/(prefix+'_structure.sobj')))
    big=load(str(out/(prefix+'_big_certificate.sobj')))
    for missing in range(3):
        tag=('rev' if reverse else 'fwd')+str(missing)
        active=[i for i in range(3) if i!=missing]
        BB=matrix(K,2,2,[K(cc(0,0,H,Z,S,N)) for cc in lead['blocks'][tuple(active)].list()])
        if reverse:active.reverse();BB=BB.transpose()
        vv=J*BB.transpose()*vector(K,[1,0]);images=[A.zero()]*6
        images[2*active[0]]=ll
        images[2*active[1]]=kk*vv[0];images[2*active[1]+1]=kk*vv[1]
        ph=R.hom(images,A);rows=[];positions=[]
        for pos,f in enumerate(norm['cubics']):
            if norm['odd_pairs'][pos] not in active:continue
            pol=ph(f);ei,ej=(1,0) if norm['odd_pairs'][pos]==active[0] else (0,1)
            rows.append([K(pol.monomial_coefficient(ll**(ei+2)*kk**ej)),K(pol.monomial_coefficient(ll**ei*kk**(ej+2)))])
            positions.append(norm['cubic_sources'][pos])
        MM=matrix(K,rows)
        saved=ep['items'][missing] if reverse else ep['layers'][missing]
        originalM=saved['matrix'] if reverse else saved['amplitude_matrix'].apply_map(lambda x:K(x(1,0)))
        checks[tag+'_original_four_cubic_rows']=len(rows)==4 and MM==originalM
        checks[tag+'_original_cubic_indices']=positions==(saved['source_indices'] if reverse else saved['sources'])
        raw=[MM[i,0]*MM[j,1]-MM[i,1]*MM[j,0] for i,j in combinations(range(4),2)]
        if reverse:
            prim=saved['primitives']
            checks[tag+'_endpoint_minors']=raw==saved['raw']
            pairs=list(zip(raw,prim))
        else:
            checks[tag+'_endpoint_minors']=raw==[K(f(1,0)) for f in saved['minors']]
            prim=saved['endpoint_primitives'];pairs=list(zip([v for v in raw if v],prim))
            checks[tag+'_nonzero_endpoint_count']=len(pairs)==len(prim)==sum(bool(v) for v in raw)
        for j,(val,p) in enumerate(pairs):
            checks[tag+'_protected_endpoint_strip_'+str(j)]=(not val and not p) or (bool(p) and isunit(val/K(p)))
        structure=st['items'][missing]['rows']
        checks[tag+'_all_Nfree_rows']=[v['source_index'] for v in structure]==[j for j,p in enumerate(prim) if p and p.degree(C.gen(3))==0]
        for j,row in enumerate(structure):
            ff=prod(f**e for f,e in row['factors']);rat=K(row['source'])/K(ff)
            checks[tag+'_Nfree_original_'+str(j)]=row['source']==prim[row['source_index']]
            checks[tag+'_source_factorization_'+str(j)]=bool(rat) and all(rat.derivative(v)==0 for v in C.gens())
        split=next(row for row in structure if len(row['factors'])>1)
        others=[P(row['source'](pH,pZ,pS,0)) for row in structure if row is not split]
        projection=load(str(out/(prefix+'_projection_'+str(missing)+'.sobj')))
        contact=load(str(out/(prefix+'_contact_'+str(missing)+'.sobj')))
        checks[tag+'_complete_branch_count']=projection['complete'] and len(projection['items'])==len(split['factors'])==3
        expected_contact=[]
        for branch,((of,e),item) in enumerate(zip(split['factors'],projection['items'])):
            f=P(of(pH,pZ,pS,0));btag=tag+'_branch'+str(branch)
            checks[btag+'_original_sources']=f==item['source_factor'] and others==item['other_sources']
            rp=f.resultant(others[0],pS);rq=f.resultant(others[1],pS);common=rp.gcd(rq)
            checks[btag+'_two_resultants']=bool(rp) and bool(rq) and [rp,rq]==item['resultants']
            checks[btag+'_gcd']=common==item['gcd']
            aa,ar=rp.quo_rem(common);bb,br=rq.quo_rem(common)
            checks[btag+'_finite_residual_intersection']=not ar and not br and aa.gcd(bb).is_constant()
            product=prod(fac**ee for fac,ee in item['factors']);ratio=PK(common/product)
            checks[btag+'_complete_factorization']=bool(ratio) and all(ratio.derivative(v)==0 for v in P.gens())
            for index,(fac,ee) in enumerate(item['factors']):
                g=to_plane(fac)
                if any(not g.quo_rem(guard)[1] for guard in pg):continue
                if len(g.dict())>100:
                    bi=big['items'][missing]
                    checks[btag+'_big_source_link']=branch==2 and bi['F']*bi['F_unit']==g
                    continue
                expected_contact.append((branch,index))
                pos=next(j for j,rec in enumerate(contact['records']) if rec['branch']==branch and rec['factor_index']==index)
                savedc=contact['certificates'][pos];ck=btag+'_factor'+str(index)
                checks[ck+'_same_factor']=savedc['factor']==g
                if not g.degree(bZ):
                    checks[ck+'_vertical']=savedc.get('vertical',False) and g.degree(bH)>0
                    continue
                den=g.derivative(bZ);num=-g.derivative(bH)
                checks[ck+'_actual_slope']=savedc['num_S']==num and savedc['den_S']==den and g.gcd(den).is_constant()
                gg=g
                for jj,si in enumerate(savedc['source_indices']):
                    residual=contact_eval(sourcecoeffs(([f]+others)[si]),num,den,g)
                    checks[ck+'_residue'+str(si)]=residual==savedc['residues'][jj]
                    gg=gg.gcd(residual)
                checks[ck+'_no_actual_interval']=gg.is_constant()
            print(btag,'source/projection/contact PASS',round(time.monotonic()-start,2),flush=True)
        checks[tag+'_short_contact_complete']=contact['complete'] and expected_contact==[(v['branch'],v['factor_index']) for v in contact['records']]
        cert.append({'orientation':'reverse' if reverse else 'forward','missing':missing,'original_endpoint_matrix':MM,'raw_endpoint_minors':raw})
assert all(checks.values()),[k for k,v in checks.items() if not v]
save({'checks':checks,'certificates':cert,'scope':'All six endpoint orientations reconstructed from original normalized cubics, all 18 exact resultant branches and complete short contact factors. Combine with the separately audited integer good-reduction certificates for six large factors. Applies to any nonzero active bilinear block, not only invertible blocks.'},str(out/'independent_two_pair_endpoint_source_reader.sobj'))
report={'status':'EXACT_BOTH_ENDPOINT_ORIENTATIONS_SOURCE_PASS','check_count':len(checks),'checks':checks,'elapsed_seconds':time.monotonic()-start}
(out/'independent_two_pair_endpoint_source_reader.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='checks'}),flush=True)
