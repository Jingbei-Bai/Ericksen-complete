"""Necessary linear incidence on the Z0!=0 inactive-row branch; all p support kept."""
from pathlib import Path
from sage.all import load,save,PolynomialRing,matrix
import time,json
out=Path(__file__).resolve().parent;start=time.monotonic()
d=load(str(out/'independent_ratio_inactive_row3_RR_derivatives.sobj'))
gf=load(str(out/'independent_ratio_p_boundary_second_gauge.sobj'))
R=d['ring'];K=d['coefficient_field'];a,r,H,Z=K.gens();v=R.gens();oldK=gf['coefficient_field'];ph=oldK.ring().hom(d['original_parameter_map'],K);cache={}
def cv(c):
 c=oldK(c)
 if c not in cache:cache[c]=ph(c.numerator())/ph(c.denominator())
 return cache[c]
def pull(f):return R({tuple(map(int,e)):cv(c) for e,c in f.dict().items() if sum(map(int,e))==3 and not e[4] and not e[5]})
units=(2*(H-Z)/(H*Z),-2*(Z-1)/Z,r*(H+1)/(H*(1+r)))
zs=[d['missing_Z0'],pull(gf['first_gauged_fields'][1,7])/units[1],pull(gf['first_gauged_fields'][0,6])/units[0]]
der=d['items'][0]['derivatives'][0]['derivative'];factors=list(der.factor());T=next(p for p,n in factors if p.total_degree()>0)
assert der/T in K
def row(f):
 rr=[]
 for j in (6,7):
  c=f.derivative(v[j]);assert all(not e[6] and not e[7] and not e[8] for e in c.dict());rr.append(c)
 c=f-rr[0]*v[6]-rr[1]*v[7];assert all(not e[6] and not e[7] and not e[8] for e in c.dict());return rr+[c]
M=matrix(R,[row(zs[1]),row(zs[2]),row(T)]);det=M.det();fac=det.factor() if det else None
print('DET',fac,'terms',len(det.dict()),'sec',time.monotonic()-start,flush=True)
records=[];support=[]
for mask in range(1,16):
 def restrict(f):return R({tuple(map(int,e)):c for e,c in f.dict().items() if all((mask>>j)&1 or not e[j] for j in range(4))})
 active=(bool(mask&12),bool(mask&3))
 rows=[restrict(T)]+[restrict(zs[i+1]) for i in range(2) if active[i]]
 # Keep missing Z0 itself, never infer it vanishes on this branch.
 z0=restrict(zs[0]);qdeg=[]
 for it in d['items']:
  f=restrict(it['restricted_initial']);qdeg.append({'index':it['source_index'],'terms':len(f.dict()),'factorization':f.factor() if f else None})
 rec={'support_mask':mask,'active_rows':[i+1 for i,b in enumerate(active) if b],'Z0_zero_identically':not z0,'T_zero_identically':not rows[0],'Q_term_counts':[s['terms'] for s in qdeg]};records.append(rec)
 support.append({'record':rec,'linear_sources':rows,'missing_Z0':z0,'quartics':qdeg})
save({'ring':R,'coefficient_field':K,'Z':zs,'T':T,'T_original_derivative':der,'linear_augmented_matrix':M,'compatibility_determinant':det,'determinant_factorization':fac,'support_strata':support,'scope':__doc__,'qualification':'Only on an initial point with Z0!=0 does the genuine principal derivative source force T=0. Z1/Z2 are imposed only when their D columns are active. No actual invariance of a p support stratum is assumed.'},str(out/'independent_ratio_inactive_row3_linear_strata.sobj'),compress=False)
(out/'independent_ratio_inactive_row3_linear_strata.json').write_text(json.dumps({'determinant_factorization':str(fac),'records':records,'seconds':time.monotonic()-start},indent=2)+'\n');print('DONE',time.monotonic()-start,flush=True)
