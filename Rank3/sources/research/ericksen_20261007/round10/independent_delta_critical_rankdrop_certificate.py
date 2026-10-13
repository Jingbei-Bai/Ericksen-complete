"""Small exact source certificate for the full critical rank-drop locus."""
from sage.all import *
from pathlib import Path
import time,json,itertools
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_delta_critical_ranktwo_structure.sobj'));R=d['ring'];s,z=R.gens();rows=d['distinct_rows'];mons=sorted(set((i,tuple(map(int,e))) for row in rows for i in range(2) for e in row[i].dict()));dic=[{(i,tuple(map(int,e))):c for i in range(2) for e,c in row[i].dict().items()} for row in rows];M=matrix(QQ,[[v.get(e,0) for e in mons] for v in dic]);inds=list(M.transpose().pivots());selected=[rows[i] for i in inds];print('CONSTANT_ROW_RANK',len(inds),'rows',inds,flush=True)
src=[];links=[]
for i,j in itertools.combinations(range(len(selected)),2):
 f=selected[i][0]*selected[j][1]-selected[i][1]*selected[j][0]
 if f:
  q,rem=f.quo_rem(s**4);assert rem==0;src.append(q);links.append((inds[i],inds[j]))
I=R.ideal(src);gb=I.groebner_basis();print('SMALLBASIS',list(gb),'sources',len(src),'seconds',time.monotonic()-st,flush=True)
cert=[]
for f in gb:
 cc=list(f.lift(I));assert sum((a*b for a,b in zip(cc,src)),R.zero())==f;cert.append({'target':f,'cofactors':cc});print('LIFT',f,'nonzero',sum(bool(v) for v in cc),'seconds',time.monotonic()-st,flush=True)
save({'source_file':'independent_delta_critical_ranktwo_structure.sobj','ring':R,'constant_coefficient_matrix':M,'constant_basis_indices':inds,'minor_sources_divided_s4':src,'minor_links':links,'basis':gb,'certificates':cert,'scope':__doc__},str(p/'independent_delta_critical_rankdrop_certificate.sobj'),compress=False)
(p/'independent_delta_critical_rankdrop_certificate.json').write_text(json.dumps({'basis_rows':len(inds),'source_minors':len(src),'identities':len(cert),'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n')
