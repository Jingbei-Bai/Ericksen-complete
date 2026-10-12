"""Link the actual line pivot and equality-wall polynomials to their sources."""
from sage.all import *
from pathlib import Path
import json,hashlib
out=Path(__file__).absolute().parent;research=out.parent.parent;checks=[];imports=[]
def old(n):
 p=research/'round10'/(n+'.sobj');imports.append(dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()));return load(str(p))
def ck(n,v):assert bool(v),n;checks.append(n)
line=old('independent_supercritical_line_even');ex=old('independent_supercritical_exceptional_material');mat=old('independent_supercritical_line_material')
K=line['field'];a,r,H,Z,t=K.gens();V=line['normal_matrix'];det=V[:3,:3].det()
ck('literal original line pivot determinant',det==line['determinant'])
pa,pb=ex['line_pivot_polynomials'];ph=pa.parent().hom([t],K);P=H*ph(pa)+ph(pb)
ck('printed line pivot polynomial',P==H*(t**4+5*t+3)-2*t**4-5*t**3+5*t+2)
u=det/P;allowed=[a,r,1+r,1+a*r,a-1,H,Z,H+1,Z-1,Z-H,H*Z-H+Z,t,t+1]
for f,e in u.factor():
 ck('every actual line pivot unit '+str(f),any((K(f)/g).numerator().total_degree()==0 and (K(f)/g).denominator().total_degree()==0 for g in allowed) or K(f) in QQ)
c1,c2=ex['line_pivot_bezout'];ck('actual base coefficients exact Bezout',c1*pa+c2*pb==1)
U=ex['line_wall_ring'];tt=U.gen();n=U.base_ring().gen();hh=mat['ring'].hom([U(n),U(n),tt],U)
for j,f in enumerate(mat['forms']):
 actual=hh(f);power=actual.valuation();reduced=actual//(tt**power)
 ck('wall removes only support t '+str(j),reduced*tt**power==actual)
 ck('all original wall polynomial source '+str(j),reduced==ex['line_wall_forms'][j])
ck('original wall Sylvester matches',ex['line_wall_forms'][0].sylvester_matrix(ex['line_wall_forms'][1]).det()==ex['line_wall_resultant'])
report=dict(status='PASS',count=len(checks),checks=checks,imports=imports,scope='Actual source-derived line pivot has exactly the printed coprime H coefficients up to protected units; equality wall forms are substitutions of all original three line forms with only t-support powers removed.')
(out/'supercritical_line_guards.json').write_text(json.dumps(report,indent=2));print('PASS',len(checks))
