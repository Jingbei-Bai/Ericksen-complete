"""Exact amplitude forms with all six homogeneous direction coordinates retained."""
from pathlib import Path
import json,time
from sage.all import load,save,matrix
out=Path(__file__).resolve().parent;start=time.monotonic()
d=load(str(out/'general_curve_two_pair_offdiagonal.sobj'))
meta=load(str(out/'general_curve_homogeneous_normalized.sobj'));R=d['ring'];xx=R.gens()
A=matrix(R,6,3);Q=matrix(R,3,6);pairs=[(0,0),(1,1),(2,2),(0,1),(0,2),(1,2)]
checks={}
for row,f in enumerate(d['cubics']):
 odd=meta['odd_pairs'][row]
 for m,c in f.dict().items():
  deg=[int(m[2*i]+m[2*i+1]) for i in range(3)]
  deg[odd]-=1;assert min(deg)>=0 and sorted(deg)==[0,0,2]
  j=deg.index(2);A[row,j]+=R({tuple(m):c})
 checks['cubic_reconstruction_'+str(row)]=sum(A.row(row))==f
for row,f in enumerate(d['quartics']):
 for m,c in f.dict().items():
  deg=[int(m[2*i]+m[2*i+1]) for i in range(3)];assert all(v%2==0 for v in deg)
  ids=[i for i in range(3) for _ in range(deg[i]//2)];assert len(ids)==2
  j=pairs.index(tuple(ids));Q[row,j]+=R({tuple(m):c})
 checks['quartic_reconstruction_'+str(row)]=sum(Q.row(row))==f
assert all(checks.values())
records={'pure_Q_factors':[[str(Q[i,j].factor()) for j in range(3)] for i in range(3)],
 'A_terms':[[len(f.dict()) for f in row] for row in A.rows()],
 'Q_terms':[[len(f.dict()) for f in row] for row in Q.rows()],'checks':checks,'elapsed_seconds':time.monotonic()-start}
save({'ring':R,'A':A,'Q':Q,'Q_monomial_pairs':pairs,'odd_pairs':meta['odd_pairs'],
 'source_cubics':d['cubics'],'source_quartics':d['quartics'],'checks':checks,
 'scope':'All homogeneous pair directions; no direction or amplitude coordinate is inverted.'},str(out/'general_curve_full_amplitude_forms.sobj'))
(out/'general_curve_full_amplitude_forms.json').write_text(json.dumps(records,indent=2)+'\n');print(json.dumps(records,indent=2),flush=True)
