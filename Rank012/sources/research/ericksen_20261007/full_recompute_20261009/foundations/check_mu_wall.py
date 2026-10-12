"""Fresh mu-wall restriction and Sylvester computation from rebuilt fields.

The arithmetic is retained from the original source. The independent eight
printed cubic coefficients are read from current expanded paper/main.tex.
--inspect-tex reads declarations with ordinary Python; it performs no algebra.
"""
from pathlib import Path
import argparse, json, time, re, hashlib

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--bundle-root', type=Path)
parser.add_argument('--main-tex', type=Path)
parser.add_argument('--output-dir', type=Path)
parser.add_argument('--inspect-tex', action='store_true')
args=parser.parse_args()
bundle=args.bundle_root.resolve() if args.bundle_root is not None else next(
    (p for p in Path(__file__).resolve().parents if (p/'paper/main.tex').is_file()), None)
if bundle is None and (not args.inspect_tex or args.main_tex is None):
    parser.error('Run inside the extracted proof bundle or pass --bundle-root')
paper=args.main_tex.resolve() if args.main_tex is not None else (bundle/'paper/main.tex').resolve()
tex=paper.read_text(encoding='utf-8')
blocks=re.findall(r'\\begin\{equation\}(.*?)\\end\{equation\}',tex,re.S)
coefficient_declarations=[]
for block in blocks:
    hit=re.search(r'f_\{([12]),([0-3])\}=(.*)',block,re.S)
    if hit: coefficient_declarations.append(tuple(map(int,hit.group(1,2))))
if set(coefficient_declarations)!={(n,j) for n in (1,2) for j in range(4)} or len(coefficient_declarations)!=8:
    raise AssertionError('Expected all eight current mu-wall cubic coefficients exactly once')
if args.inspect_tex:
    print(json.dumps(dict(status='PARSE_PASS',cubic_coefficients=8,
                         arithmetic_executed=False),indent=2))
    raise SystemExit(0)

from sage.all import *
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr,standard_transformations,implicit_multiplication_application
SOURCE_OUT=bundle/'research/ericksen_20261007/full_recompute_20261009/foundations'
OUT=args.output_dir.resolve() if args.output_dir is not None else SOURCE_OUT
OUT.mkdir(parents=True,exist_ok=True)
START=time.monotonic();checks=[]
def check(name,ok,**extra):
    checks.append(dict(name=name,passed=bool(ok),**extra)); print(name,bool(ok),flush=True)
    (OUT/'mu_wall_checks.json').write_text(json.dumps(dict(status='RUNNING',checks=checks),indent=2))
    assert ok,name
src=load(str(SOURCE_OUT/'state17_fresh.sobj')); oldP=src['ring'];old=oldP.gens()
AD=PolynomialRing(QQ,names=('a','d'));a,d=AD.gens();KK=AD.fraction_field();bd=KK(a*(1-d)/(1-a*d))
R=PolynomialRing(KK,names=('z','X','Y','U','V','W','T','J','K'),order='lex');z,X,Y,U,V,W,T,J,K=R.gens();q=1-d;g=1-a*d
point=[z,X,Y,U,1,W,V,-W,q,0,-X*V+W*z*(1+2*z/d),2*(U-X),0,T,0,J,K]
def coefficient(c):
    return KK(c.numerator()(a,bd))/KK(c.denominator()(a,bd))
def restrict(f):return sum((coefficient(c)*prod(point[i]**e for i,e in enumerate(mon)) for mon,c in f.dict().items()),R.zero())
F=[[restrict(f) for f in row] for row in src['F']]
C3=[restrict(f) for f in src['C3']]
E2=(X*X+d)*Y+2*z*(1+z/d)*X*W
E3=(d*q-Y*Y)*X+2*z*(1-z/d)*Y*W
pair=a*d/(4*g)*(-(1+a*q)*(z*z+d*d)+2*(1+q*g)*z-4*X*U-4*g*Y*V)
check('three differentiated mu equations',all(F[i][8]-q*F[i][4]==0 for i in range(3)))
check('opposite-pair differentiated equations',[F[i][5]+F[i][7] for i in range(3)]==[pair,0,0])
check('mu cubic sources exactly printed',[C3[0],C3[1]+E2,C3[2]-E3]==[0,0,0])
target=2*X*X*Y*Y+(d-z)*Y*Y-q*(d+z)*X*X
check('undivided biquadratic identity',(1-z/d)*Y*E2-(1+z/d)*X*E3==target)
eq=[F[i][j] for i in range(3) for j in (9,12,14)]
check('D3t0 equals J',eq[6]==J)
check('D3t3 has protected coefficient two',eq[7].derivative(T)==2)
check('D2t5 has protected coefficient q/2',eq[5].derivative(K)==q/2)
tvalue=-(eq[7]-2*T)/2
kvalue=-2*(eq[5]-q*K/2).subs({T:tvalue})/q
check('two protected higher-jet eliminations',eq[7].subs({T:tvalue})==0 and eq[5].subs({T:tvalue,K:kvalue})==0)
remaining=[eq[idx]-eq[idx].derivative(T)*eq[7]/2 for idx in (3,8)]
check('remaining equations contain no eliminated jets',all(all(f.degree(v)==0 for v in (T,J,K)) for f in remaining))
L=R.fraction_field();ux=-g*W**2/(a*d)-(z-d)**2/(4*d);vy=W**2/(a*d)+q*(z+d)**2/(4*d)
check('UX and VY remove opposite-pair quadratic',L(pair)(U=ux/X,V=vy/Y)==0)
wval=-d*(X**2+d)*Y/(2*z*(z+d)*X);y2=q*(z+d)*X**2/(2*X**2+d-z)
reduced=[]
for f in remaining:
    ff=L(f)(U=ux/X,V=vy/Y)(W=wval)
    num=ff.numerator();den=ff.denominator()
    # Exact quotient reduction modulo Y^2-y2; the denominator stays guarded.
    rem=sum((L(c)*prod(R.gen(i)**e for i,e in enumerate(mon) if i!=2)*Y**(mon[2]%2)*y2**(mon[2]//2) for mon,c in num.dict().items()),L.zero())
    reduced.append(rem/den)
print('fresh regular tangency restrictions complete',flush=True)
AZ=PolynomialRing(QQ,names=('a','d','z'));az,dz,zz=AZ.gens();RR=PolynomialRing(AZ,names=('R',));rvar=RR.gen()
printed={1:[None]*4,2:[None]*4}
for block in blocks:
    hit=re.search(r'f_\{([12]),([0-3])\}=(.*)',block,re.S)
    if not hit:continue
    n,j=map(int,hit.group(1,2));expr=hit.group(3)
    expr=expr.replace(r'\end{aligned}','').replace(r'\\','').replace('&','')
    expr=re.sub(r'\^\{(\d+)\}',r'**\1',expr).strip()
    parsed=parse_expr(expr,local_dict=dict(zip(('a','d','z'),sp.symbols('a d z'))),transformations=standard_transformations+(implicit_multiplication_application,))
    printed[n][j]=AZ({mon:QQ(str(c)) for mon,c in sp.Poly(parsed,*sp.symbols('a d z')).terms()})
assert all(c is not None for coeffs in printed.values() for c in coeffs)
Fprinted=[sum(c*rvar**j for j,c in enumerate(printed[n])) for n in (1,2)]
for i,ff in enumerate(reduced):
    # Compare complete rational expressions, allowing only an explicit nonzero multiplier.
    fp=sum(R(c(a,d,z)) * X**(2*j) for j,c in enumerate(printed[i+1]))
    ratio=L(ff)/L(fp)
    protected=R(X*Y*z*(z+d)*(2*X**2+d-z))
    factors=[f for pol in (ratio.numerator(),ratio.denominator()) for f,e in pol.factor()]
    check('fresh source %d equals printed cubic by a protected unit'%(i+1),bool(ratio) and all(protected%R(f)==0 for f in factors),multiplier=str(ratio),factors=list(map(str,factors)))
    explicit=(d-1)*X/(4*a*z**2*Y*(z+d)*(2*X**2+d-z)**2)
    if i==0:explicit/=a*d-1
    check('explicit source-to-cubic multiplier %d'%(i+1),ratio==explicit,formula=str(explicit))
coeffs=[[printed[n][j] for j in (3,2,1,0)] for n in (1,2)]
Syl=matrix(AZ,[[AZ.zero()]*shift+co+[AZ.zero()]*(2-shift) for co in coeffs for shift in range(3)])
res=Syl.det();factor=-1024*(dz-1)*(az*dz-1)*az**3*dz**4*zz**6*(zz+dz)**9
P12,rem=res.quo_rem(factor)
check('full printed six by six resultant factorization',rem==0 and P12.degree(zz)==12,terms=len(P12.dict()))
check('generic leading coefficient',P12.coefficient({zz:12})==48*az**3*(az*dz-az-1)**4)
KAZ=AZ.fraction_field();special=KAZ(P12)(d=1+1/KAZ(az))
check('entire exceptional specialization',special==55296*zz**4*(az+1)**2*(az*zz-az-1)*(az*zz+az+1)/az)
save({'printed_cubics':Fprinted,'resultant':res,'P12':P12,'fresh_remainders':reduced},str(OUT/'mu_wall_fresh.sobj'))
report=dict(status='PASS',checks=checks,elapsed_seconds=time.monotonic()-START,sample_spectra_used=0,
            maintained_tex=str(paper.relative_to(bundle)) if paper.is_relative_to(bundle) else str(paper),
            tex_sha256=hashlib.sha256(paper.read_bytes()).hexdigest())
(OUT/'mu_wall_checks.json').write_text(json.dumps(report,indent=2))
print('COMPLETE',report['elapsed_seconds'],flush=True)
