"""All-material exclusion of the v=2H ordinary sheet for l nonzero.
Only generic-base polynomial nonvanishing is used; all u and ramification
are retained through true full normals and the coupled eta ideal.
"""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
d=load(str(p/'ordinary_pair_eta_algebraicity_certificate.sobj'));P=d['ring'];K=P.base_ring();a,r,H,Z=K.gens();q=1+r;m=1+a*r;l,v=P.gens()
def ck(n,x):assert bool(x),n;checks.append(n)
f=d['FA'].subs({v:2*H});top=P.zero();maxz=-1
for e,c in f.dict().items():
 c=K(c);ck('polynomial coefficient '+str(e),c.denominator()==1)
 for ce,n in c.numerator().dict().items():
  ce=tuple(map(int,ce));maxz=max(maxz,ce[3])
  if ce[3]==5:top+=n*a**ce[0]*r**ce[1]*H**ce[2]*l**int(e[0])
ck('highest Z degree five',maxz==5)
ck('unit highest Z coefficient',top==-16*m*q*(H+1)**5*l*l)
# Check the complete first-normal source rows at u=0 against the already
# reconstructed coupled matrix. This certifies nonzero constant terms in u
# of genuine five-minors, without expanding the unnecessary higher terms.
full=load(str(p/'ordinary_pair_full_normal_blocks.sobj'));FP=full['ring'];fl,fu,fv=FP.gens();MM=full['matrix'];labs=full['labels'];CC=d['matrix'];clabs=d['labels']
def mp(z):return P({(int(e[0]),int(e[2])):c for e,c in z.dict().items() if not e[1]})
selected=[]
for i,b in enumerate(d['two_rows']):
 rs=b+d['three_rows'][0];actual=[]
 for lab in rs:
  row=MM.row(labs.index(lab));ck('only first source eta '+str((i,lab)),row[4]==0)
  got=[mp(row[j]) for j in (0,1,2,3,5)]
  want=list(CC.row(clabs.index(lab))[:5]);ck('full true-normal constant term '+str((i,lab)),got==want)
  actual.append(got)
 A=matrix(P,actual)
 # Its block pattern makes the constant determinant a signed product.
 ck('upper block support '+str(i),all(A[x,y]==0 for x in (0,1) for y in (0,3,4)))
 ck('lower block support '+str(i),all(A[x,y]==0 for x in (2,3,4) for y in (1,2)))
 z=A.det();prod0=d['two_minors'][i]*d['three_minors'][0]
 sign=1 if z==prod0 else -1
 ck('genuine five minor constant product '+str(i),z==sign*prod0)
 selected.append({'rows':rs,'constant_in_u':z,'sign':sign})
save({'ring':P,'FA_on_sheet':f,'highest_Z_coefficient':top,'selected_full_normal_minors':selected,'checks':checks,'scope':__doc__},str(p/'ordinary_pair_v_plus_sheet_reader.sobj'),compress=False)
rec={'status':'PASS ALL L NONZERO V=2H SHEET EXCLUSION CERTIFICATE','checks':len(checks),'seconds':time.monotonic()-st}
(p/'ordinary_pair_v_plus_sheet_reader.json').write_text(json.dumps(rec,indent=2)+'\n');print(rec,flush=True)
