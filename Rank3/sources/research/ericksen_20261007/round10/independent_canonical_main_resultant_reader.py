"""All exact integer nodes independently replayed, with a strict determinant bound."""
from sage.all import *
from pathlib import Path
import json,time,hashlib,sys
sys.set_int_max_str_digits(0)
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(n,v):
    if not bool(v):raise AssertionError(n)
    checks.append(n)
inp=load(str(p/'independent_ratio_canonical_second_resultant_inputs.sobj'))
extra=load(str(p/'independent_ratio_canonical_second_extra_resultant_index31_inputs.sobj'))
B=inp['ring'];br,bv=B.gens();polys=list(inp['polynomials'])+[B(extra['polynomials'][1])]
ck('same_first_polynomial',B(extra['polynomials'][0])==B(polys[0]))
Zr=PolynomialRing(ZZ,'audit_r');r=Zr.gen();Zv=PolynomialRing(ZZ,'audit_v')
coefs=[]
for k,f in enumerate(polys):
    cs=[Zr({int(e[0]):ZZ(c) for e,c in f.dict().items() if int(e[1])==j}) for j in range(int(f.degree(bv))+1)]
    coefs.append(cs)
    ck('integer_reconstruction_'+str(k),B({(i,j):c for j,g in enumerate(cs) for i,c in g.dict().items()})==f)
names=['independent_ratio_canonical_second_resultant12','independent_ratio_canonical_second_resultant13','independent_ratio_canonical_second_extra_resultant_index31']
nodefiles=[names[0]+'_interpolation.nodes',names[1]+'_interpolation.nodes',names[2]+'.nodes']
candidates=[];reports=[]
for j,name in enumerate(names,1):
    payload=load(str(p/(name+'.sobj')));candidate=Zr(payload['resultant']);candidates.append(candidate)
    n,m=len(coefs[0])-1,len(coefs[j])-1
    # Every term in the padded Sylvester determinant uses m coefficients
    # of the first polynomial and n of the second. This bound does not
    # depend on leading coefficients surviving after specialization.
    bound=m*max(f.degree() for f in coefs[0])+n*max(f.degree() for f in coefs[j])
    ck('strict_sylvester_degree_bound_'+str(j),bound==2652)
    ck('candidate_under_bound_'+str(j),candidate.degree()<=bound)
    visited=set();det_hash=hashlib.sha256();rowcount=0
    with (p/nodefiles[j-1]).open() as inp:
        for line in inp:
            nn,yy=line.split();node=ZZ(nn);expected=ZZ(yy)
            ck('distinct_node_%d_%d'%(j,rowcount),node not in visited);visited.add(node)
            fv=[f(node) for f in coefs[0]];gv=[f(node) for f in coefs[j]]
            ff=Zv(fv);gg=Zv(gv)
            ck('no_node_degree_drop_%d_%s'%(j,node),ff.degree()==n and gg.degree()==m)
            actual=ff.resultant(gg)
            ck('exact_integer_resultant_%d_%s'%(j,node),actual==expected)
            ck('entire_candidate_evaluation_%d_%s'%(j,node),candidate(node)==actual)
            # Independent determinant convention at fixed exact nodes;
            # all other nodes above are independently rebuilt resultants.
            if node in (2,17):
                S=matrix(ZZ,n+m,n+m)
                for col in range(m):
                    for k,c in enumerate(fv):S[k+col,col]=c
                for col in range(n):
                    for k,c in enumerate(gv):S[k+col,m+col]=c
                ck('independent_Sylvester_integer_det_%d_%s'%(j,node),S.det()==(-1)**(n*m)*actual)
            det_hash.update(str(node).encode()+b':'+str(actual).encode()+b'\n');rowcount+=1
            if rowcount%100==0:
                print('PAIR',j,'EXACT NODES',rowcount,'seconds',round(time.monotonic()-st,2),flush=True)
                (p/'independent_canonical_main_resultant_reader_progress.json').write_text(json.dumps({'pair':j,'nodes':rowcount,'seconds':time.monotonic()-st})+'\n')
    ck('enough_exact_nodes_'+str(j),len(visited)==bound+1)
    ck('node_interval_'+str(j),visited==set(range(2,2655)))
    rv=candidate.valuation();rem=candidate//r**rv;qv=0
    while rem(-1)==0:rem=rem//(r+1);qv+=1
    ck('protected_valuations_'+str(j),(rv,qv)==(188,218))
    ck('saved_stripped_'+str(j),rem==Zr(payload['stripped']))
    reports.append({'pair':[0,j],'v_degrees':[n,m],'degree_bound':int(bound),'resultant_degree':int(candidate.degree()),'distinct_exact_nodes':rowcount,'replayed_hash':det_hash.hexdigest(),'r_valuation':int(rv),'q_valuation':qv})
    print('PAIR COMPLETE',j,round(time.monotonic()-st,2),flush=True)
q=r+1;strip=[f//(r**188*q**218) for f in candidates];g12=strip[0].gcd(strip[1]);gg=g12.gcd(strip[2]);ck('three_stripped_resultants_coprime',gg.degree()==0 and gg!=0)
fullgcd=candidates[0].gcd(candidates[1]).gcd(candidates[2]);quo,rem=fullgcd.quo_rem(r**188*q**218);ck('full_resultant_gcd_only_protected',not rem and quo.degree()==0 and quo!=0)
report={'status':'PASS_EXACT_ALL_PARAMETER_RESULTANTS','checks':len(checks),'seconds':time.monotonic()-st,'pairs':reports,'first_two_stripped_gcd_degree':int(g12.degree()),'final_stripped_gcd_degree':int(gg.degree()),'scope':'Three complete resultant identities proved by 2653 distinct exact ZZ evaluations each and degree bound2652. Exact gcd leaves only r^188(r+1)^218. Together with independently replayed four original-source certificates and protected column/prefactors, the entire cyclic normal chart t2*(t0+t1+t2)!=0 has full rank6 for every ordered positive simple material. No single-spectrum inference or assumed leading-degree protection is used.'}
save({'checks':checks,'integer_coefficient_polynomials':coefs,'reports':reports,'resultants':candidates,'stripped_gcd':gg,'full_gcd':fullgcd},str(p/'independent_canonical_main_resultant_reader.sobj'),compress=False)
(p/'independent_canonical_main_resultant_reader.json').write_text(json.dumps(report,indent=2)+'\n');print(report,flush=True)
