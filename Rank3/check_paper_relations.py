"""Compare current printed formulas with freshly reconstructed algebra."""
from pathlib import Path
import argparse, hashlib, json, re

def expression(tex):
    tex=tex.strip().rstrip('.,')
    for token in [r'\left',r'\right',r'\quad',r'\qquad',r'\,',r'\!']:
        tex=tex.replace(token,' ')
    tex=tex.replace(r'\tfrac',r'\frac')
    tex=re.sub(r'\\frac([0-9])([0-9])',r'\\frac{\1}{\2}',tex)
    aliases={r'\kappa_0':'r',r'\varsigma':'S',r'\eta':'H',r'\zeta':'Z',r'k_{\mathrm{dir}}':'k'}
    for a,b in aliases.items():tex=tex.replace(a,' '+b+' ')
    tex=re.sub(r'([dfh])_\{?(\d)\}?',r' \1\2 ',tex)
    def group(t,start):
        if t[start]!='{':raise ValueError('Unbraced fraction')
        depth=1;i=start+1
        while depth:
            depth+=(t[i]=='{')-(t[i]=='}');i+=1
        return t[start+1:i-1],i
    while r'\frac' in tex:
        start=tex.index(r'\frac');a,end=group(tex,start+5);b,end=group(tex,end)
        tex=tex[:start]+'('+expression(a)+')/('+expression(b)+')'+tex[end:]
    tex=re.sub(r'\^\{(\d+)\}',r'**\1',tex).replace('^','**')
    if re.search(r'[\\{}]',tex):raise ValueError('Unsupported formula syntax: '+tex)
    tokens=re.findall(r'\*\*|[A-Za-z][A-Za-z0-9_]*|\d+|[()+\-*/]',tex)
    if ''.join(tokens)!=re.sub(r'\s+','',tex):raise ValueError('Incomplete formula tokenization')
    out=[]
    for t in tokens:
        if out and (re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*|\d+',out[-1]) or out[-1]==')') and (re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*|\d+',t) or t=='('):out.append('*')
        out.append(t)
    return ''.join(out)

def value(tex,env):
    from sage.all import sage_eval, ZZ
    code=expression(tex)
    names=set(re.findall(r'[A-Za-z][A-Za-z0-9_]*',code))
    if not names<=set(env):raise ValueError('Unknown formula variables: '+str(names-set(env)))
    code=re.sub(r'\b\d+\b',lambda m:'ZZ('+m[0]+')',code)
    return sage_eval(code,locals=dict(env,ZZ=ZZ),preparse=False)

def verify(package,root):
    from sage.all import PolynomialRing, QQ, load, matrix, vector, gcd, prod
    mapping=json.loads((package/'PROOF_MAP.json').read_text())
    targets=mapping['printed_relations'];checks=[]
    def ck(name,test):
        if not test:raise AssertionError(name)
        checks.append(name)
    p=root/'research/ericksen_20261007/full_recompute_20261009'
    fg=load(str(p/'rank3/fresh_FG_sources.sobj'))
    C=fg['ring'];H,Z,S=C.gens();env={'H':H,'Z':Z,'S':S,'L':H*S-Z}
    guard=C(value(targets['gstar'],env));ck('printed localization equals reconstructed guard',guard==fg['guard'])
    P0=C(value(targets['P0'],env));remaining=gcd(fg['resultants'])//fg['G']
    for f,_ in guard.factor():
        while remaining and remaining.quo_rem(f)[1]==0:remaining=remaining//f
    ck('printed P0 equals the reconstructed residual common factor',remaining/remaining.leading_coefficient()==P0/P0.leading_coefficient())
    # All seven source rows and the full normal matrix were rebuilt by the
    # original raw-source reader, including the compensated quadratic term.
    data=load(str(p/'rank3/source_replays/independent_supercritical_inactive_root_source_reader_cyclic/generated/round10/independent_supercritical_inactive_root_source_reader_cyclic.sobj'))
    K=data['field'];A,R=K.gens()
    F=PolynomialRing(QQ,'r').fraction_field();r=F.gen();q=1+r
    phi=K.ring().hom([-q/r,r],F)
    def cv(v):return phi(v.numerator())/phi(v.denominator())
    ker=vector(F,[1,-1,-4,-2,-2,0]);corr=vector(F,[-6,13,38,38,0,-4*q])
    J=matrix(F,[[cv(c) for c in row] for row in data['J'].rows()])
    ck('complete source kernel on the resonance',J*ker==0)
    quadratic=F.zero()
    ck('normal source has no residual root phase',all(mon[0]==mon[1]==0 for mon in data['sources'][20].dict()))
    for mon,co in data['sources'][20].dict().items():
        if sum(mon[2:])==2:quadratic+=cv(co)*prod(ker[j]**int(mon[j+2]) for j in range(6))
    actual=quadratic+(J*corr)[20]
    ck('printed normal coefficient equals the compensated full-source coefficient',actual==value(targets['normal_coefficient'],{'r':r,'q':q}))
    # The conics are evaluated from the freshly source-derived quartics.
    data=load(str(p/'rank3_fibres/rebuild/round10/independent_simple_direction_reader.sobj'))
    forms=data['quartic_amplitudes'];R=forms[0].parent();K=R.base_ring()
    a,r,H,Z,S,N=K.gens();L=H*S-Z;env={'H':H,'Z':Z,'S':S,'k':H*S/Z}
    scales=[S/(8*H**3*Z),L/(8*H**4*Z),-S*L/(8*H**3*Z**2)]
    forms=[f/scales[i] for i,f in enumerate(forms)]
    x,y,z=R.gens()
    minor=forms[0].monomial_coefficient(y*y)*forms[1].monomial_coefficient(z*z)-forms[0].monomial_coefficient(z*z)*forms[1].monomial_coefficient(y*y)
    ck('printed conic coefficient minor equals the source minor',minor==value(targets['conic_minor'],env))
    for i,f in enumerate(forms):
        M=matrix(K,3,3,lambda j,k:f.monomial_coefficient(R.gen(j)**2) if j==k else f.monomial_coefficient(R.gen(j)*R.gen(k))/2)
        ck('printed conic determinant '+str(i+1),M.det()==value(targets['conic_determinants'][i],env))
    # Evaluate the two actual printed pseudoremainders, then the fixed degree
    # Sylvester determinant. This is a polynomial identity before specialization.
    C=PolynomialRing(QQ,names=('d0','d1','d2','f0','f1','f2','f3','f4','w'))
    d0,d1,d2,f0,f1,f2,f3,f4,w=C.gens();env=dict(zip(C.variable_names(),C.gens()))
    af=C(value(targets['norm_af'],env));bf=C(value(targets['norm_bf'],env))
    d=d2*w*w+d1*w+d0;f=f0+f1*w+f2*w*w+f3*w**3+f4*w**4
    syl=matrix(C,6,6)
    for i in range(4):
        for j,co in enumerate([d2,d1,d0]):syl[i,i+j]=co
    for i in range(2):
        for j,co in enumerate([f4,f3,f2,f1,f0]):syl[4+i,i+j]=co
    ck('printed universal pseudoremainders and fixed degree norm identity',d2*bf*bf-d1*af*bf+d0*af*af==d2**3*syl.det())
    report={'status':'PASS','checks':checks,'count':len(checks),'paper_sha256':mapping['paper']['sha256'],
            'scope':'Current printed formulas compared with fresh exact source objects and universal polynomial arithmetic.'}
    dest=root/'run/current_paper_relations.json';dest.write_text(json.dumps(report,indent=2)+'\n')
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--package',type=Path,required=True);p.add_argument('--work-root',type=Path,required=True)
    a=p.parse_args();print(json.dumps(verify(a.package,a.work_root),indent=2))

if __name__=='__main__':main()
