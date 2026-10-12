"""Independent determinant-pivot and full three-singular-pair trace audit."""
from pathlib import Path
import time,json
from sage.all import load,save,QQ,PolynomialRing,prod,matrix

out=Path(__file__).resolve().parent;start=time.monotonic();checks={}
raw=load(str(out/'independent_scalar_lambda_lift.sobj'))
src=load(str(out/'general_curve_scalar_singular_source.sobj'))
pv=load(str(out/'general_curve_singular_determinant_pivots.sobj'))
gr=load(str(out/'general_curve_singular_trace_graph.sobj'))
ct=load(str(out/'general_curve_singular_trace_contact.sobj'))
O=raw['ring'];h,z,s,n=O.gens();J=matrix(O,[[0,-1],[1,0]])
og=[h,z,h-z,1-z,1+h,1+h-z,1-s,1+h*s-z,(1+h)*s-z]
P=raw['blocks'][0,1];Q=raw['blocks'][0,2];RR=raw['blocks'][1,2]
M=J*Q*J*RR.transpose()*J*P.transpose()
checks['original_monodromy']=M==raw['source_monodromy'] and M.trace()==gr['trace_source']

def isunit(poly,allowed):
    R=allowed[0].parent();poly=R(poly);K=R.fraction_field()
    if not poly:return False
    for f,e in poly.factor():
        f=R(f)
        if f.is_constant():continue
        if not any(all((K(f)/a).derivative(x)==0 for x in R.gens()) for a in allowed):return False
    return True

for ij,block in raw['blocks'].items():
    checks['original_det_'+str(ij)]=block.det()==src['determinants'][ij]*src['removed'][ij]
    checks['det_protected_strip_'+str(ij)]=isunit(src['removed'][ij],og)

B=PolynomialRing(QQ,names=('H','Z'));BH,BZ=B.gens();BK=B.fraction_field()
U=PolynomialRing(BK,'S');US=U.gen()
lead=[BH*(BH-1)*US+2*BH+BZ,(BH+2*BZ)*US-BZ*(BZ+1)]
vel=[(BH*(BH-1),-2*BH-BZ),(BH+2*BZ,BZ*(BZ+1))]
plane_guards=[BH,BZ,BH-BZ,1-BZ,BH+1,1+BH-BZ]
def to_U(p):
    result=U.zero()
    for mon,co in p.dict().items():
        assert not mon[3]
        result+=co*BH**mon[0]*BZ**mon[1]*US**mon[2]
    return result
for i,ij in enumerate([(0,1),(0,2)]):
    f=src['determinants'][ij];aN=to_U(f.derivative(n));bN=to_U(f.subs({n:0}))
    checks['linear_N_'+str(i)]=f==f.derivative(n)*n+f.subs({n:0})
    expected=(-2*(1+BH-BZ)*BZ*BH*(BH+1)*lead[0] if i==0 else
              2*(1+BH-BZ)*(BZ-1)*BZ*BH*lead[1])
    checks['det_own_leading_'+str(i)]=aN==expected
    item=pv['items'][i]
    checks['pivot_original_source_'+str(i)]=item['block']==ij and item['leading_factor']==lead[i] and item['constant']==bN
    resultant=lead[i].resultant(bN)
    checks['pivot_resultant_'+str(i)]=resultant==item['resultant'] and bool(resultant)
    product=prod([f**e for f,e,dp,gcd in item['factor_checks']]);ratio=BK(resultant/product)
    checks['pivot_factorization_'+str(i)]=bool(ratio) and all(ratio.derivative(v)==0 for v in B.gens())
    for j,(f,e,dp,gd) in enumerate(item['factor_checks']):
        f=B(f);actual=vel[i][0]*f.derivative(BH)+vel[i][1]*f.derivative(BZ)
        checks['pivot_contact_'+str(i)+'_'+str(j)]=actual==dp and f.gcd(actual)==gd
        checks['pivot_factor_excluded_'+str(i)+'_'+str(j)]=gd.is_constant() or isunit(f,plane_guards)
checks['pivot_Hone_vertical']=B(BH-1).degree(BH)==1
checks['pivot_Hplus2Z_fixed']=lead[1](S=0).subs({BH:-2*BZ})==-BZ*(BZ+1)

R=gr['ring'];H,Z,S=R.gens();K=R.fraction_field()
old=[R(g(H,Z,S,0)) for g in og]
expected_pairs=[((0,1),(0,2)),((0,1),(1,2)),((0,2),(1,2))]
checks['all_pairs']=len(gr['graphs'])==3 and [v['pair'] for v in gr['graphs']]==expected_pairs
for i,item in enumerate(gr['graphs']):
    first,second=item['pair'];det=src['determinants'][first]
    aa=R(det.derivative(n)(H,Z,S,0));bb=R(det(H,Z,S,0))
    Nvalue=-K(bb)/aa;own=H*(H-1)*S+2*H+Z if i<2 else (H+2*Z)*S-Z*(Z+1)
    allowed=old+[own]
    checks['own_pivot_'+str(i)]=item['N_pivot']==own and item['guards']==allowed
    checks['protected_N_graph_'+str(i)]=Nvalue==item['N'] and isunit(aa,allowed) and isunit(Nvalue.denominator(),allowed)
    dvalue=K(src['determinants'][second](H,Z,S,Nvalue));tvalue=K(M.trace()(H,Z,S,Nvalue))
    checks['det_graph_'+str(i)]=dvalue.numerator()==item['determinant_compatibility']*item['determinant_removed']
    checks['trace_graph_'+str(i)]=tvalue.numerator()==item['trace_row']*item['trace_removed']
    checks['det_graph_units_'+str(i)]=isunit(item['determinant_removed'],allowed) and dvalue.denominator()==item['determinant_denominator'] and isunit(dvalue.denominator(),allowed)
    checks['trace_graph_units_'+str(i)]=isunit(item['trace_removed'],allowed) and tvalue.denominator()==item['trace_denominator'] and isunit(tvalue.denominator(),allowed)
    p=item['determinant_compatibility'];t=item['trace_row'];res=p.resultant(t,S)
    checks['resultant_'+str(i)]=bool(res) and res==item['resultant']
    product=prod([f**e for f,e in item['factors']]);ratio=K(res/product)
    checks['factorization_'+str(i)]=bool(ratio) and all(ratio.derivative(v)==0 for v in R.gens())
    checks['degrees_'+str(i)]=p.degree(S)==5 and t.degree(S)==7
    print('source and projection',i,'PASS',round(time.monotonic()-start,2),flush=True)

def base(p):return B(p(BH,BZ,0))
def coeffs(p):
    cs=[B.zero()]*(p.degree(S)+1)
    for mon,c in p.dict().items():cs[mon[2]]+=c*BH**mon[0]*BZ**mon[1]
    return cs
def rem(p,f):return p.quo_rem(f)[1]
def evaluate_independently(cs,num,den,f):
    np=[B.one()];dp=[B.one()]
    for _ in range(len(cs)-1):np.append(rem(np[-1]*num,f));dp.append(rem(dp[-1]*den,f))
    result=B.zero()
    for j,c in enumerate(cs):result=rem(result+rem(c*np[j],f)*dp[len(cs)-1-j],f)
    return result
seen=[];cert=[];vertical=0
for pair,item in enumerate(gr['graphs']):
    for index,(fac,e) in enumerate(item['factors']):
        f=base(fac)
        if isunit(f,plane_guards):continue
        seen.append((pair,index));saved=next(v for v in ct['certificates'] if v['pair']==pair and v['factor_index']==index)
        key=str(pair)+'_'+str(index);checks['factor_'+key]=saved['factor']==f
        if not f.degree(BZ):
            checks['vertical_'+key]=saved.get('vertical',False) and f.degree(BH)>0;vertical+=1
            cert.append({'pair':pair,'index':index,'factor':f,'vertical':True});continue
        num=-f.derivative(BH);den=f.derivative(BZ)
        checks['contact_slope_'+key]=saved['num_S']==num and saved['den_S']==den and f.gcd(den).is_constant()
        result=evaluate_independently(coeffs(item['determinant_compatibility']),num,den,f)
        checks['original_compatibility_residue_'+key]=saved['residues']==[result] and saved['source_labels']==['det_compatibility']
        checks['no_curvature_extra_'+key]=saved['curvature'] is None
        checks['excludes_interval_'+key]=f.gcd(result).is_constant()
        cert.append({'pair':pair,'index':index,'factor':f,'residue':result,'gcd':f.gcd(result)})
        print('contact',pair,index,'PASS',round(time.monotonic()-start,2),flush=True)
checks['complete_contact_list']=seen==[(v['pair'],v['factor_index']) for v in ct['certificates']] and len(seen)==21 and vertical==2
checks['complete_flags']=gr['complete'] and ct['complete']
assert all(checks.values()),[k for k,v in checks.items() if not v]
save({'checks':checks,'certificates':cert,'scope':'All three pair-singular trace-zero actual-open direction-curve cases. Uses separate determinant-leading pivot exclusions, not the narrower scalar A-zero lemma. No amplitude or zero-pair-boundary dimension conclusion.'},str(out/'independent_singular_trace_reader.sobj'))
report={'status':'EXACT_THREE_SINGULAR_TRACE_BRANCHES_PASS','check_count':len(checks),'checks':checks,'elapsed_seconds':time.monotonic()-start}
(out/'independent_singular_trace_reader.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='checks'}),flush=True)
