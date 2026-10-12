"""Exact all-spectrum necessary sources in a native Laurent representation.

No spectrum specialization, Groebner calculation, or solution exclusion.
The six original graph coordinates are reconstructed directly from saved
all-parameter source data; all twelve source identities are replayed.
"""
from pathlib import Path
from functools import lru_cache
import pickle,json,time,hashlib
import sympy as sy
from sage.all import QQ,PolynomialRing,LaurentPolynomialRing,load,save
OUT=Path(__file__).resolve().parent; BASE=OUT.parent; start=time.time()
line=load(str(OUT/'balanced_rank_one_line.sobj'))
source_path=BASE/'generic_frame_complete_api.pkl'
source=pickle.load(source_path.open('rb'))
vs=source['variables'];aa,bb=source['parameters'];F=source['vectorfields']
P=PolynomialRing(QQ,names=('a','b','h'));a,b,h=P.gens();K=P.fraction_field()
L=LaurentPolynomialRing(K,names=('u','v','s'),implementation='generic');u,v,s=L.gens()
A=a*(b-1);C=b*(a-1);D=a-b
zmap=line['Q'].hom([K(a),K(b),K(h),K(0),K(0),K(0)],K)
z=zmap(line['z'].numerator())/zmap(line['z'].denominator())
assert 1+D*h/A-D*z/C==0
scales=[K(-z/C),K(1/D),K(h/A)]
@lru_cache(None)
def scale_power(j,n):return scales[j]**n
def inverse_unit(q):
    terms=q.dict();assert len(terms)==1,('nonmonomial Laurent divisor',q)
    mon,c=next(iter(terms.items()))
    return L({tuple(-n for n in mon):1/c})
def map_polynomial(q):
    terms={}
    for ex,c in q.dict().items():
        mon=(ex[4],ex[3],ex[5]); coefficient=K(c*a**ex[0]*b**ex[1]*h**ex[2])
        for j in range(3):coefficient*=scale_power(j,ex[3+j])
        terms[mon]=terms.get(mon,K.zero())+coefficient
    return L(terms)
def mapped(q):return map_polynomial(q.numerator())*inverse_unit(map_polynomial(q.denominator()))
ts=[]
for j,q in enumerate(line['solution']):
    ts.append(mapped(q));print('t',j,'terms',len(ts[-1].dict()),'seconds',round(time.time()-start,2),flush=True)
values=[L(1),L(1/D)*u,L(1/D)*v,L(h/A)*u,L(h),L(h/A)*s,L(-z/C)*v,L(-z/C)*s,L(z)]
context=dict(zip(vs[:15],values+ts));context[aa]=L(a);context[bb]=L(b)
@lru_cache(None)
def evaluate(q):
    if q in context:return context[q]
    if q.is_Rational:return L(QQ(int(q.p))/int(q.q))
    if q.is_Add:return sum((evaluate(x) for x in q.args),L.zero())
    if q.is_Mul:
        result=L.one()
        for x in q.args:result*=evaluate(x)
        return result
    if q.is_Pow:
        n=int(q.exp); base=evaluate(q.base)
        return base**n if n>=0 else inverse_unit(base)**(-n)
    raise ValueError(q)
checks={}
def check(label,value):assert value,label;checks[label]=True
AA=aa*(bb-1);CC=bb*(aa-1);DD=aa-bb
B=[vs[8]*vs[5]+CC/AA*vs[4]*vs[7],vs[0]*vs[6]+DD/CC*vs[8]*vs[2],vs[4]*vs[1]-AA/DD*vs[0]*vs[3]]
check('Bcal_parameterization',all(evaluate(q)==0 for q in B))
raw=[sum(sy.diff(q,x)*F[i][j] for j,x in enumerate(vs[:9])) for q in B for i in range(3)]
raw+=source['constraints_degree3']
residuals=[]
for j,q in enumerate(raw):
    actual=evaluate(q)
    check('original_residual_%d'%j,actual==mapped(line['residuals'][j]))
    residuals.append(actual)
    print('residual',j,'terms',len(actual.dict()),'seconds',round(time.time()-start,2),flush=True)
for j in [9,10,11,0,4,8]:check('selected_graph_equation_%d'%j,residuals[j]==0)
det=mapped(line['det'])
check('graph_determinant_protected',det==L(-4*C*h**2*(A+D*h)**2/(A**5*D))*u*v*s)
fields=[]
for i in range(3):
    sigma=evaluate(F[i][0]);vel=[evaluate(F[i][j])-values[j]*sigma for j in range(9)]
    check('ratio_line_tangency_%d'%i,vel[8]==L(z.derivative(h))*vel[4])
    fields.append([vel[4],L(D)*vel[1],-L(C/z)*vel[6]+L(C/z**2)*values[6]*vel[8],
                   L(A/h)*vel[5]-L(A/h**2)*values[5]*vel[4]])
    print('field',i,'terms',[len(q.dict()) for q in fields[-1]],'seconds',round(time.time()-start,2),flush=True)
def dh(q):return L({tuple(mon):c.derivative(h) for mon,c in q.dict().items()})
initial=[residuals[j] for j in [1,2,5]]
derivatives=[]
for j,q in enumerate(initial):
    partials=[dh(q)]+[q.derivative(c) for c in L.gens()]
    for i,field in enumerate(fields):
        derivative=sum((p*f for p,f in zip(partials,field)),L.zero())
        derivatives.append(derivative)
        print('derivative',j,i,'terms',len(derivative.dict()),'seconds',round(time.time()-start,2),flush=True)
known=[a,b,h,a-1,b-1,a-b,A+D*h];support=set()
@lru_cache(None)
def protected(den):
    q=P(den)
    for factor in known:
        while q%factor==0:
            q=q//factor;support.add(str(factor))
    assert q.is_constant(),('unprotected divisor',q)
for q in ts+values+sum(fields,[])+residuals+derivatives:
    for c in q.coefficients():protected(c.denominator())
check('all_coefficient_denominators_protected',True)
payload={'P':P,'K':K,'L':L,'z':z,'values':values,'solution_t':ts,'determinant':det,
         'fields':fields,'source_residuals':residuals,'source_indices':[1,2,5],
         'initial':initial,'derivatives':derivatives,'checks':checks,
         'known_parameter_factors':known,'scope':'Actual Bcal=0,L0=0,uvw!=0 open branch only; necessary conditions, not an exclusion.'}
save(payload,str(OUT/'balanced_line_generic_laurent_sources.sobj'))
report={'status':'PASS_ALL_PARAMETER_ORIGINAL_SOURCE_LAURENT',
        'source_sha256':hashlib.sha256(source_path.read_bytes()).hexdigest(),
        'checks':checks,'initial_terms':[len(q.dict()) for q in initial],
        'derivative_terms':[len(q.dict()) for q in derivatives],
        'parameter_denominator_factors':sorted(support),
        'scope':payload['scope'],'seconds':time.time()-start}
(OUT/'balanced_line_generic_laurent_sources.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k!='checks'},indent=2),flush=True)
