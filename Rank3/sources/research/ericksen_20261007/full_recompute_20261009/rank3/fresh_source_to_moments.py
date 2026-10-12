"""Derive the printed five-row moment matrix from fresh actual curve sources.
The six nonzero-gradient directions are normalized only after their full source
equations have been constructed. All stripped row factors are checked against
the explicit guard; no coefficient or moment determinant is inverted.
"""
from sage.all import *
from pathlib import Path
import json,time
out=Path(__file__).absolute().parent;research=out.parent.parent;st=time.monotonic();checks=[]
d=load(str(out/'fresh_curve_sources.sobj'));oldR=d['ring'];oldK=oldR.base_ring();a,r,H,Z,S,N=oldK.gens();q=1+r
T=PolynomialRing(oldK,names=('t1','t2','t3'));tt=T.gens()
def ck(name,value):
 if not value:raise AssertionError(name)
 checks.append(name)
guardold=oldK.ring()(a*r*(1+r)*(1+a*r)*H*Z*(H-Z)*(Z-1)*(H+1)*(S-1)*(H*S-Z+1)*((H+1)*S-Z)*(1+H-Z)*S*(H*S-Z))
guardf=[g for g,e in guardold.factor()]
def guardunit(c,label):
 c=oldK(c)
 for part in [c.numerator(),c.denominator()]:
  for factor,e in part.factor():ck('protected row factor '+label,any(factor.quo_rem(g)[1]==0 for g in guardf))
def leading(f,degree):return oldR({tuple(e):c for e,c in f.dict().items() if sum(e)==degree})
ev=oldR.hom([tt[0],T.one(),tt[1],T.one(),tt[2],T.one()],T)
mixed=[]
for i in range(3):
 p=ev(leading(d['rows'][i],2));cn=T({e:c.derivative(N) for e,c in p.dict().items()})
 ck('mixed curvature coefficient constant '+str(i),cn.is_constant() and cn!=0)
 unit=-cn.constant_coefficient();guardunit(unit,'mixed'+str(i));p=p/unit
 ck('mixed curvature coefficient minusone '+str(i),T({e:c.derivative(N) for e,c in p.dict().items()})==-1);mixed.append(p)
dirs=[mixed[1]-mixed[0],mixed[2]-mixed[0]];scales=[]
for i in range(3,7):
 p=ev(leading(d['rows'][i],3));den=lcm([c.denominator() for c in p.coefficients()]);con=gcd([oldK.ring()(c*den) for c in p.coefficients()]);unit=oldK(con/den);guardunit(unit,'cubic'+str(i));p=p/unit
 ck('spectral independence '+str(i),all(c.derivative(a)==0 and c.derivative(r)==0 and c.derivative(N)==0 for c in p.coefficients()));dirs.append(p);scales.append(unit)
for i,p in enumerate(dirs):ck('all six direction source spectral independence '+str(i),all(c.derivative(a)==0 and c.derivative(r)==0 and c.derivative(N)==0 for c in p.coefficients()))
CC=PolynomialRing(QQ,names=('H','Z','S'));HH,ZZ,SS=CC.gens();KK=CC.fraction_field();RR=PolynomialRing(KK,names=('U','V','W'));U,V,W=RR.gens();LL=HH*SS-ZZ
mp=oldK.ring().hom([0,0,HH,ZZ,SS,0],KK)
def reduced(c):
 c=oldK(c);ck('spectral removal denominator',mp(c.denominator())!=0);return mp(c.numerator())/mp(c.denominator())
images=[(V-HH*W)/LL+SS/(2*ZZ),(U+W)/SS-LL/(2*HH*ZZ),W-LL/(2*HH*ZZ)]
full=[sum((reduced(c)*prod(images[j]**int(e) for j,e in enumerate(mon) if e) for mon,c in p.dict().items()),RR.zero()) for p in dirs]
AA=PolynomialRing(KK,'w');w=AA.gen()
def univariate(p):
 ck('only W coefficient',all(e[0]==e[1]==0 for e in p.dict()));return AA({int(e[2]):c for e,c in p.dict().items()})
mom=[]
for i,f in enumerate(full):
 ca=f.coefficient({U:1,V:1});cb=f.coefficient({U:1,V:0});cc=f.coefficient({U:0,V:1});const=f.subs({U:0,V:0})
 ck('constant divisible by W '+str(i),all(e[2]>0 for e in const.dict()))
 dd=RR({(0,0,int(e[2])-1):c for e,c in const.dict().items()})
 ck('full moment decomposition '+str(i),f==U*V*ca+U*cb+V*cc+W*dd);mom.append([univariate(v) for v in [ca,cb,cc,dd]])
M6=matrix(AA,mom);cap=1+HH-ZZ;Qstar=2*ZZ*w*(2*HH*ZZ*cap*w-(LL*LL+HH*SS*SS-ZZ));qdiv=Qstar//w;rows=[]
for i in range(2,6):
 div,rem=M6[i,3].quo_rem(qdiv);ck('Qstar exact division '+str(i),rem==0);rows.append(list(M6.row(i))[:3]+[div])
comb=(LL+1)*M6.row(0)+(SS-1)*M6.row(1);div,rem=comb[3].quo_rem(qdiv);ck('fifth row Qstar division',rem==0);rows.append(list(comb)[:3]+[div])
Ms=matrix(AA,rows)
inputfile=research/'rank3_simplification_20261009/original_moment_definitions.json';defs=json.loads(inputfile.read_text());env=dict(H=AA(HH),Z=AA(ZZ),S=AA(SS),w=w)
target=matrix(AA,[[sage_eval(e,locals=env) for e in row] for row in defs['general_curve_mixed_Q_matrix_plain']]);ratios=[]
guard=CC(HH*ZZ*(HH-ZZ)*(ZZ-1)*(HH+1)*(1+HH-ZZ)*(SS-1)*(HH*SS-ZZ+1)*((HH+1)*SS-ZZ)*SS*(HH*SS-ZZ))
for i in range(5):
 ratio=Ms[i,0]/target[i,0]
 ck('row multiplier independent w '+str(i),ratio.numerator().degree()==0 and ratio.denominator().degree()==0)
 ratio=KK(ratio);ck('all four printed entries are fresh source '+str(i),all(Ms[i,j]==ratio*target[i,j] for j in range(4)))
 for part in [ratio.numerator(),ratio.denominator()]:
  for factor,power in part.factor():ck('printed row old-guard multiplier '+str(i),guard.quo_rem(factor)[1]==0)
 ratios.append(ratio)
save(dict(source_matrix=Ms,printed_matrix=target,row_multipliers=ratios,six_directions=dirs,full_moment_matrix=M6),str(out/'fresh_source_to_moments.sobj'),compress=False)
rep=dict(status='PASS',checks=checks,count=len(checks),seconds=time.monotonic()-st,printed_row_multipliers=[str(x.factor()) for x in ratios],fresh_derivation_chain=['Section2 state17','cubic and u reconstruction','RR coordinate graph','actual joint vector fields','all six second curve tangencies','unit-gradient direction equations','centered moments','five printed source rows'],material_specializations=0,
 scope='All five printed moment rows are exact old-guard combinations of freshly regenerated actual sources. Other projective direction supports and amplitude finiteness are separate obligations.')
(out/'fresh_source_to_moments.json').write_text(json.dumps(rep,indent=2));print(rep,flush=True)
