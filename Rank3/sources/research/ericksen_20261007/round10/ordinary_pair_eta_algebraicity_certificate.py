"""All-material, all nonzero l certificate: surviving u=0 forces algebraic v.

Independent reconstruction from original first normals, original eta-square
coefficients, and the new literal full pair source. No Groebner computation.
True seven-column minors retain range corrections through the coupled ideal.
"""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
old=load(str(p/'independent_ratio_incoming12_generic_normal.sobj'))
sq=load(str(p/'independent_incoming_simple_wall_quadratics.sobj'))
pair=load(str(p/'ordinary_pair_incoming_full_twojet.sobj'))
saved=load(str(p/'ordinary_pair_eta_coupled_matrix.sobj'))
P=saved['ring'];K=P.base_ring();a,r,H,Z=K.gens();l,v=P.gens();q=1+r;k=a*r;m=1+k;G=H*Z-H+Z
O=old['ring'];ol,ou,ov=O.gens();OK=O.base_ring();ph=OK.hom([a,a*q/m,H/r,q*Z/r],K)
def ck(n,x):
 assert bool(x),n
 checks.append(n)
cache={}
def cv(x):
 x=OK(x)
 if x not in cache:cache[x]=ph(x)
 return cache[x]
def omap(f):
 return P({(int(e[0]),int(e[2])):cv(c) for e,c in O(f).dict().items() if not e[1]})
eta=next(z for z in sq['items'] if z['name']=='eta');eta2={int(z['row']):z['second'] for z in eta['values']}
es=pair['normal_exponents'];PU=pair['jet_ring'];pl,pu,pv=PU.gens()
def pmap(f):return P({(int(e[0]),int(e[2])):K(c) for e,c in f.dict().items() if not e[1]})
rows=[]
for typ,i in saved['labels']:
 if typ=='old':
  ck('actual eta kernel '+str(i),not old['matrix'][i,4].subs({ou:0}))
  row=[omap(old['matrix'][i,j]) for j in (0,1,2,3,5)]+[omap(old['matrix'][i,4].derivative(ou)),omap(eta2.get(i,O.zero()))]
 else:
  it=next(z for z in pair['items'] if z['pair']==i);vals=it['jet_coefficients'];f=it['full_homogeneous_pullback']
  ck('whole pair boundary '+str(i),not vals[0])
  # Read first/second coefficients from the whole polynomial again.
  extracted=[]
  for ex in es:
   extracted.append(PU({tuple(map(int,e[:3])):c for e,c in f.dict().items() if tuple(map(int,e[3:]))==tuple(ex)}))
  ck('all literal pair jets '+str(i),extracted==vals)
  ff=extracted[1:7];row=[pmap(ff[j]) for j in (0,1,2,3,5)]+[pmap(ff[4].derivative(pu)),pmap(extracted[es.index((0,0,0,0,2,0))])]
 rows.append(row)
M=matrix(P,rows);ck('complete coupled matrix',M==saved['matrix']);labels=saved['labels']
def rid(x):return labels.index(x)
def det(A):
 if not A.nrows():return P.one()
 i=min(range(A.nrows()),key=lambda i:sum(bool(x) for x in A.row(i)))
 return sum(((-1)**(i+j)*A[i,j]*det(A.delete_rows([i]).delete_columns([j])) for j in range(A.ncols()) if A[i,j]),P.zero())
def minor(rs,cs):return det(M.matrix_from_rows_and_columns([rid(z) for z in rs],cs))
ueta=Z**2*G/(7*(Z-1)**2*(Z-H));eta_square=2*Z**2*H**2/((Z-H)**2*(Z-1)**2)
ck('u eta protected pivot',M[rid(('old',0)),5]==ueta)
ck('eta square protected pivot',M[rid(('old',2)),6]==eta_square)
two_rows=[(('old',1),('old',17)),(('old',17),('old',29))]
three_rows=[(('old',5),('old',7),('old',11)),(('old',5),('old',7),('pair',2))]
two=[minor(rs,(1,2)) for rs in two_rows];three=[minor(rs,(0,3,4)) for rs in three_rows]
expected2=-4/(r**3*a**3*m**2*G)*l*(m*l*l+k)**2*(m*l*l*Z*(3*H+2)+k*(3*H+1))
ck('two block generic short identity',two[0]==expected2)
F6=P(two[1]/(-4/(r**2*a**2*m**2*G)*l*l))
R=PolynomialRing(K,'L');L=R.gen()
def squares(f):
 dd={}
 for e,c in P(f).dict().items():
  ck('even l and no v '+str(e),e[0]%2==0 and e[1]==0)
  dd[int(e[0])//2]=c
 return R(dd)
wall=squares(F6)(-k/m)
ck('two block remaining wall identity',wall==-48*k**3*(4*H*Z-4*H+3*Z-1))
units=[-H/(4*(H+1)**6*(Z-H)*r*G**2*m),H/(16*(H+1)**6*(Z-H)*r*G*m)]
FA=P(three[0]/(units[0]*(v+2)*l**3));FB=P(three[1]/(units[1]*(v+2)*l**2))
def vcoef(f,n):return P({(int(e[0]),0):c for e,c in f.dict().items() if e[1]==n})
ck('first v degree four',FA.degree(v)==4)
ck('second v degree four',FB.degree(v)==4)
ca=m*H**4*(3*H+2)*l*l*(5*r*l*l+4*q)
cb=3*m*H**4*Z*l*l*(4*r*l*l+5*q)
ck('first exact highest coefficient',vcoef(FA,4)==ca)
ck('second exact highest coefficient',vcoef(FB,4)==cb)
ck('protected highest coefficient Bezout',5*(4*r*l*l+5*q)-4*(5*r*l*l+4*q)==9*q)
full=[]
for i,rs2 in enumerate(two_rows):
 for j,rs3 in enumerate(three_rows):
  rs=(('old',0),('old',2))+rs2+rs3
  z=minor(rs,tuple(range(7)));prod0=ueta*eta_square*two[i]*three[j]
  sign=1 if z==prod0 else -1
  ck('true full seven determinant '+str((i,j)),z==sign*prod0)
  full.append({'rows':rs,'sign':sign,'two_index':i,'three_index':j})
save({'ring':P,'matrix':M,'labels':labels,'two_rows':two_rows,'three_rows':three_rows,'two_minors':two,'three_minors':three,'u_eta_pivot':ueta,'eta_square_pivot':eta_square,'FA':FA,'FB':FB,'highest_coefficients':[ca,cb],'two_wall':wall,'full_seven_factorizations':full,'checks':checks,'scope':__doc__},str(p/'ordinary_pair_eta_algebraicity_certificate.sobj'),compress=False)
rec={'status':'EXACT ALL MATERIAL ALL NONZERO L ALGEBRAICITY CERTIFICATE','checks':len(checks),'seconds':time.monotonic()-st,'scope':'Surviving actual u=0 ordinary boundary with l nonzero has algebraic v. Does not exclude sheets v=-2 or v=2H, or l=0.'}
(p/'ordinary_pair_eta_algebraicity_certificate.json').write_text(json.dumps(rec,indent=2)+'\n');print(rec,flush=True)
