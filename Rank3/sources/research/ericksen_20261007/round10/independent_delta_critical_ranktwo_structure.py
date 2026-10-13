"""Exact spectral separation and small rank-drop ideal of the critical matrix."""
from sage.all import *
from pathlib import Path
import time,json,itertools
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_delta_critical_amplitude_coefficients.sobj'));K=d['field'];a,r,s,z=K.gens();M=d['matrix'];w=vector(K,[-r,1+r,r*(1+r)*s*s]);assert M*w==0;N=M*diagonal_matrix(K,w);rows=[];items=[];checks=[];T=PolynomialRing(QQ,names=('cs','cz'));ss,zz=T.gens();ph=K.ring().hom([T.zero(),T.zero(),ss,zz],T)
protected=[a,r,1+r,1+a*r]
for i in range(N.nrows()):
 rr=list(N.row(i));den=lcm([v.denominator() for v in rr]);pp=[K.ring()(v*den) for v in rr];gg=gcd(pp);unit=K(gg.lc())
 for f,e in gg.factor():
  if f.degree(K.ring().gens()[2])==0 and f.degree(K.ring().gens()[3])==0:
   assert any(K(f)/v in QQ and K(f)/v!=0 for v in protected),(i,f)
   unit*=K(f/f.lc())**e
 nr=[K(v)/unit for v in pp];assert all(v.derivative(a)==0 and v.derivative(r)==0 for v in nr),(i,nr)
 out=tuple(ph(v.numerator())/QQ(v.denominator()) for v in nr);assert sum(out)==0
 items.append({'row':i,'clearing_denominator':den,'removed_protected_material_unit':unit,'normalized_row':out});checks.append(i)
 if out not in rows:rows.append(out)
print('ROWS',len(rows),'degrees',sorted(set(max(v.total_degree() for v in row) for row in rows)),'seconds',time.monotonic()-st,flush=True)
minors=[];links=[]
for i,j in itertools.combinations(range(len(rows)),2):
 f=rows[i][0]*rows[j][1]-rows[i][1]*rows[j][0]
 if f:minors.append(f);links.append((i,j))
gg=gcd(minors);print('MINOR GCD',gg.factor() if gg else 0,'count',len(minors),flush=True)
I=T.ideal(minors);gb=I.groebner_basis();print('BASIS',list(gb),flush=True)
save({'source_file':'independent_delta_critical_amplitude_coefficients.sobj','field':K,'short_kernel':w,'material_transform_items':items,'ring':T,'distinct_rows':rows,'minors':minors,'minor_links':links,'minor_gcd':gg,'minor_ideal_basis':gb,'checks':checks,'scope':__doc__},str(p/'independent_delta_critical_ranktwo_structure.sobj'),compress=False)
(p/'independent_delta_critical_ranktwo_structure.json').write_text(json.dumps({'checks':len(checks),'rows':len(rows),'minors':len(minors),'basis':[str(v) for v in gb],'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n')
