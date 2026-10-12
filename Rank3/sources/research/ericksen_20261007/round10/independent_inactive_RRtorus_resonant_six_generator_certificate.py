from sage.all import *
from pathlib import Path
from itertools import permutations
p=Path(__file__).resolve().parent;d=load(str(p/'independent_inactive_RRtorus_resonant_second_module.sobj'));sp=load(str(p/'independent_inactive_RRtorus_resonant_pair_second_module.sobj'));Q=d['ring'];F=Q.fraction_field();r=Q.base_ring().gen();q=1+r;l,rho=Q.gens();lc=-q/r;rc=-q

def rp(ff):return sum((co*lc**(int(e[0])//2)*rc**(int(e[1])//2)*l**(int(e[0])%2)*rho**(int(e[1])%2) for e,co in ff.dict().items()),Q.zero())
def red(ff):
 ff=F(ff);nn=rp(ff.numerator());de=rp(ff.denominator())
 for v in (l,rho):
  cj=de.subs({v:-v});nn=rp(nn*cj);de=rp(de*cj)
 assert de;return F(nn)/de
old=d['quotient_second'];cs=[0,2,3,4,5];rows=[list(old.row(i)[j] for j in cs) for i in range(old.nrows())];ids=list(range(old.nrows()));rank=0;sel=[];pivs=[]
for j in range(5):
 ii=next(i for i in range(rank,len(rows)) if red(rows[i][j]));rows[rank],rows[ii]=rows[ii],rows[rank];ids[rank],ids[ii]=ids[ii],ids[rank];sel.append(ids[rank]);z=red(rows[rank][j]);pivs.append(z);rows[rank]=[red(v/z) for v in rows[rank]]
 for i in range(rank+1,len(rows)):
  z=red(rows[i][j]);rows[i]=[red(rows[i][h]-z*rows[rank][h]) for h in range(5)]
 rank+=1
M=old.matrix_from_rows_and_columns(sel,cs);ans=F(0)
for per in permutations(range(5)):
 inv=sum(per[i]>per[j] for i in range(5) for j in range(i+1,5));ans=red(ans+(-1)**inv*prod(M[i,per[i]] for i in range(5)))
c=sp['items'][1]['remainder'][6];assert c
v=red(ans*c);print('SELECTED',[d['source_rows'][i] for i in sel],'MINOR5',ans.factor(),'NEW6',v.factor(),flush=True)
assert all(sp['items'][1]['remainder'][j]==0 for j in range(9) if j!=6)
save(dict(ring=Q,old_row_indices=sel,old_source_rows=[d['source_rows'][i] for i in sel],old_columns=cs,old_minor=ans,new_pair=1,new_remainder=c,combined_six_minor=v,columns_in_order=['sWl','sK','dWl','dWR','dK','s2']),str(p/'independent_inactive_RRtorus_resonant_six_generator_certificate.sobj'),compress=False)
