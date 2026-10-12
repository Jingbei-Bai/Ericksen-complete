from sage.all import *
from pathlib import Path
import json,time

HERE=Path(__file__).absolute().parent
SOURCE=HERE/'source_inputs'
start=time.monotonic()
d=load(str(SOURCE/'independent_subcritical_single_RR_full_regular_chart_H1.sobj'))
fresh=load(str(HERE/'generic_G_relative_regular.sobj'))
K=FractionField(PolynomialRing(QQ,names=('k','r','w')))
k,r,w=K.gens();m=1+k;q=1+r
P=PolynomialRing(K,names=('x','z','ww','e0','e1','e2','delta','t','kappa'))
def convert(f,images):
    F=f.parent().base_ring();h=F.ring().hom(images,K)
    return P({tuple(map(int,e)):h(c.numerator())/h(c.denominator()) for e,c in f.dict().items()})
sources=[convert(f,[k,r]) for f in d['sources']]
sources.append(convert(fresh['source'],[k,r,1,-1]))
X=PolynomialRing(K,'X');E=X.quotient(X.gen()**2-2*m*(k+w*w)/k,'xb');xb=E.gen()
base=[xb,-xb,E(w),-2*xb*w,E.zero(),2*xb,E.zero(),E.zero(),E.zero()]
cache={}
def ev_monomial(ex):
    ex=tuple(map(int,ex))
    if ex not in cache:
        if any(n and not base[j] for j,n in enumerate(ex)):cache[ex]=E.zero()
        else:cache[ex]=prod(base[j]**n for j,n in enumerate(ex) if n)
    return cache[ex]
def derivative_value(f,indices=()):
    result=E.zero()
    for ex,c in f.dict().items():
        ex=list(map(int,ex));mult=1
        for j in indices:
            if not ex[j]:mult=0;break
            mult*=ex[j];ex[j]-=1
        if mult:result+=c*mult*ev_monomial(ex)
    return result
rows=[2,3,6,28,21,25,26,27]
gradient={i:vector(E,[derivative_value(sources[i],(j,)) for j in range(9)]) for i in rows}
dx=vector(E,[1,-1,0,-2*w,0,2,0,0,0])
cx=gradient[27]*dx
assert cx==-16*xb*(k+w*w)/m
xs=(3*w*w-k)/(2*(k+w*w))
ds=vector(E,[xs,1-xs,0,-2*xs*w+12*w,0,2*xs+10,32*r/q*w,0,0])
def second(i):
    active=[j for j,z in enumerate(ds) if z]
    value=E.zero()
    for a,j in enumerate(active):
        value+=ds[j]*ds[j]*derivative_value(sources[i],(j,j))/2
        for l in active[a+1:]:value+=ds[j]*ds[l]*derivative_value(sources[i],(j,l))
    return value
assert gradient[27]*ds==0
xss=-second(27)/cx
matrix_rows=[]
for i in rows:
    assert derivative_value(sources[i])==0 and gradient[i]*ds==0,i
    matrix_rows.append([gradient[i][j] for j in [3,4,5,6]]+
                       [(gradient[i]*dx)/cx,second(i)+xss*(gradient[i]*dx),gradient[i][7],gradient[i][8]])
M=matrix(E,matrix_rows)
det_cache={}
def sparse_det(row,columns):
    if not columns:return E.one()
    key=(row,columns)
    if key not in det_cache:
        ans=E.zero()
        for place,j in enumerate(columns):
            if M[row,j]:
                ans+=(-1)**place*M[row,j]*sparse_det(row+1,columns[:place]+columns[place+1:])
        det_cache[key]=ans
    return det_cache[key]
det=sparse_det(0,tuple(range(8)))
expected=E(QQ(3)/14/k**5*w*r*q*m**3*(w*w-k)**2*(w*w+k)**5)*xb
assert det==expected or det==-expected
norm=K(det[0]**2-2*m*(k+w*w)/k*det[1]**2)
target=-QQ(9)/98/k**11*w*w*r*r*q*q*m**7*(w*w-k)**4*(w*w+k)**11
assert norm==target
payload={'status':'OPEN_CONIC_RELATIVE_NORMAL_MINOR_CONFIRMED',
         'sources':['R2','R3','Q0','G(p21,Q6)','H14','Q1','Q5','Q6'],
         'normal_generators':['e0-f0','e1-f1','e2-f2','delta-d','C','s^2','t','kappa'],
         'determinant':str(K(det[1]).factor())+' * xb','norm':str(norm.factor()),
         'nonzero_locus':'w*(w^2-k)*(w^2+k) != 0',
         'fresh_sources_used':True,'seconds':time.monotonic()-start}
(HERE/'conic_open_minor.json').write_text(json.dumps(payload,indent=2)+'\n')
save({'matrix':M,'ring':E,'determinant':det,'norm':norm,'report':payload},str(HERE/'conic_open_minor.sobj'),compress=False)
print(json.dumps(payload,indent=2))
