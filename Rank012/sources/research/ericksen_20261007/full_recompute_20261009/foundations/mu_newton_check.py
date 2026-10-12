"""Current mu-wall cubic sum and all three Laurent-root calculations.

The arithmetic is retained from the reviewed mu_newton_check.py. Paths point
into the extracted minimal proof bundle; paper/main.tex contains Appendix B.
--inspect-tex checks only the eight coefficient declarations with ordinary
Python. It does not run Sage or certify polynomial identities.
"""
from pathlib import Path
import argparse, json, re, hashlib

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--bundle-root', type=Path)
parser.add_argument('--main-tex', type=Path)
parser.add_argument('--report', type=Path)
parser.add_argument('--inspect-tex', action='store_true')
args=parser.parse_args()
bundle=args.bundle_root.resolve() if args.bundle_root is not None else next(
    (p for p in Path(__file__).resolve().parents if (p/'paper/main.tex').is_file()), None)
if bundle is None and not args.inspect_tex:
    parser.error('Run inside the extracted proof bundle or pass --bundle-root')
paper=args.main_tex.resolve() if args.main_tex is not None else bundle/'paper/main.tex'
tex=paper.read_text(encoding='utf-8')
blocks=re.findall(r'\\begin\{equation\}(.*?)\\end\{equation\}',tex,re.S)
coefficient_declarations=[]
for block in blocks:
    hit=re.search(r'f_\{([12]),([0-3])\}=(.*)',block,re.S)
    if hit: coefficient_declarations.append(tuple(map(int,hit.group(1,2))))
if set(coefficient_declarations)!={(n,j) for n in (1,2) for j in range(4)} or len(coefficient_declarations)!=8:
    raise AssertionError('Expected all eight current mu-wall cubic coefficients exactly once')
for label in ('fnd:mu-short-cubic','fnd:mu-newton-roots','fnd:mu-newton-values'):
    assert '\\label{'+label+'}' in tex, label
if args.inspect_tex:
    print(json.dumps(dict(status='PARSE_PASS',cubic_coefficients=8,
                         current_formula_labels=3,arithmetic_executed=False),indent=2))
    raise SystemExit(0)

from sage.all import *
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr,standard_transformations,implicit_multiplication_application
out=bundle/'research/ericksen_20261007/full_recompute_20261009/foundations'
report_path=args.report or out/'mu_newton_checks.json'
src=load(str(out/'mu_wall_fresh.sobj'))
F,G=src['printed_cubics'];R=F.parent().gen();A=F.base_ring();a,d,z=A.gens();e=a*d-a-1
printed={1:[None]*4,2:[None]*4}
for block in re.findall(r'\\begin\{equation\}(.*?)\\end\{equation\}',tex,re.S):
    hit=re.search(r'f_\{([12]),([0-3])\}=(.*)',block,re.S)
    if not hit: continue
    n,j=map(int,hit.group(1,2));expr=hit.group(3)
    expr=expr.replace(r'\end{aligned}','').replace(r'\\','').replace('&','')
    expr=re.sub(r'\^\{(\d+)\}',r'**\1',expr).strip()
    syms=sp.symbols('a d z')
    parsed=parse_expr(expr,local_dict=dict(zip(('a','d','z'),syms)),transformations=standard_transformations+(implicit_multiplication_application,))
    printed[n][j]=A({mon:QQ(str(c)) for mon,c in sp.Poly(parsed,*syms).terms()})
assert all(c is not None for row in printed.values() for c in row)
assert F==sum(c*R**j for j,c in enumerate(printed[1]))
assert G==sum(c*R**j for j,c in enumerate(printed[2]))
print('All eight maintained cubic coefficients re-parsed and matched.',flush=True)
K=2*(d*d+3*d*z-z*z)*R**2+(2*d**3+4*d*d*z-9*d*z*z+z**3)*R+3*d*z*z*(z-d)
H=4*R*(R+d)**2+e*K
assert F+G==4*a*d*z*(d+z)*H
# Change coefficient order to inspect z asymptotics over fixed material field.
AD=PolynomialRing(QQ,names=('a','d'));aa,dd=AD.gens();KK=AD.fraction_field();ee=aa*dd-aa-1
Z=PolynomialRing(KK,names=('z',));zz=Z.gen();ZR=PolynomialRing(Z,names=('R',));rr=ZR.gen()
def transform(f):
    return sum(Z(c(aa,dd,zz))*rr**j for j,c in enumerate(f.list()))
f,h=map(transform,(F,H));checks=[]
def check(name,ok,value):
    checks.append(dict(name=name,passed=bool(ok),value=str(value)));print(name,ok,value,flush=True);assert ok
aux=PolynomialRing(KK,names=('r','c','y'));r,c,y=aux.gens()
laur=LaurentPolynomialRing(aux,names=('t',));tt=laur.gen()
def rescaled_h(argument,power):
    return tt**power*sum(laur(cc)*tt**(-k)*argument**j for j,zzpoly in enumerate(h.list()) for k,cc in enumerate(zzpoly.list()))
for n,argument,power,target,var,root,expectedderiv in [
    (0,laur(r),3,ee*(r+3*dd),r,-3*dd,ee),
    (1,1/(2*tt)+c,3,aux(1)/2-ee*c,c,1/(2*ee),-ee),
    (2,y/tt**2,6,2*y*y*(2*y-ee),y,ee/2,ee**2)]:
    scaled=rescaled_h(argument,power)
    check('scaled implicit equation has no negative t powers '+str(n),scaled.valuation()>=0,scaled.valuation())
    limit=aux(scaled[0])
    check('scaled implicit limiting equation '+str(n),limit==target,limit)
    value=limit.subs({var:root})
    deriv=limit.derivative(var).subs({var:root})
    check('limiting root and derivative '+str(n),value==0 and deriv==expectedderiv,deriv)
q0=h(-3*dd)
q1=h(zz/2+1/(2*ee))
q2=h(ee*zz**2/2)
check('small root residual has degree at most 2',q0.degree()<=2,q0.degree())
check('middle root residual has degree at most 2',q1.degree()<=2,q1.degree())
check('large root residual has degree at most 5',q2.degree()<=5,q2.degree())
v0=f(-3*dd);v1=f(zz/2+1/(2*ee));v2=f(ee*zz**2/2)
for n,v,degree,expected in [(0,v0,6,12*aa*dd*ee),(1,v1,6,aa*(dd-1)*(aa*dd-1)),(2,v2,9,aa*ee**3)]:
    check('F at root '+str(n)+' degree',v.degree()==degree,v.degree())
    check('F at root '+str(n)+' leading',v.leading_coefficient()==expected,v.leading_coefficient())
# Error terms O(z^-1), O(z^-1), O(z) have strictly smaller F order.
df=f.derivative()
for n,r0,order,leaddegree in [(0,-3*dd,-1,6),(1,zz/2+1/(2*ee),-1,6),(2,ee*zz**2/2,1,9)]:
    bounds=[]
    for k in (1,2,3):
        val=f.derivative(k)(r0)
        deg=-100 if not val else val.degree()
        bounds.append(deg+k*order)
    check('all Taylor root errors lower order '+str(n),max(bounds)<leaddegree,bounds)
res=F.resultant(H)
resz=Z(res(aa,dd,zz))
expected=-768*aa**3*dd*(dd-1)*(aa*dd-1)*ee**4
check('resultant z degree',resz.degree()==21,resz.degree())
check('resultant leading coefficient',resz.leading_coefficient()==expected,resz.leading_coefficient())
report_path.write_text(json.dumps(dict(status='PASS',checks=checks,new_identity='F1+F2=4adz(z+d)H',sources=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in (paper,)],all_eight_printed_coefficients_independently_reparsed=True,scope='all allowed fixed material parameters; exceptional e=0 positive; generic e!=0 three distinct Newton slopes'),indent=2))
