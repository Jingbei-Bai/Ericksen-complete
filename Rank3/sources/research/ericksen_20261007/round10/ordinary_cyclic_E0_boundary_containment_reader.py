
"""Independent whole-source matrix linkage and sparse Leibniz replay.
No rank, determinant library algorithm, or Groebner basis is used.
"""
from sage.all import *
from pathlib import Path
import time,json,resource,os
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(s,x):assert bool(x),s;checks.append(s)
def lap(A):
 if A.nrows()==0:return A.base_ring().one()
 i=min(range(A.nrows()),key=lambda j:sum(bool(c) for c in A.row(j)))
 return sum(((-1)**(i+j)*A[i,j]*lap(A.delete_rows([i]).delete_columns([j])) for j in range(A.ncols()) if A[i,j]),A.base_ring().zero())
D=load(str(p/'ordinary_cyclic_E0_boundary_containment_matrix.sobj'));Alt=load(str(p/'ordinary_cyclic_E0_boundary_alternate_minors.sobj'))['items'][0];S=load(str(p/'ordinary_pair_double_kernel_coupled.sobj'));E=D['ring'];l=E.gen();K=E.base_ring();a,r=K.gens();q=1+r;k=a*r;m=1+k;N=k*r;lam=m*q;OR=S['ring'];OK=OR.base_ring();ph=OK.ring().hom([a,r,1,-1],K)
def cv(f):return sum((ph(c.numerator())/ph(c.denominator())*l**int(e) for e,c in OR(f).dict().items()),E.zero())
used=int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE');resource.setrlimit(resource.RLIMIT_AS,(used+500*1024**2,)*2)
ck('all40whole sources',S['matrix'].nrows()==40)
for i in range(40):
 for j in range(9):ck('original40coefficient_'+str((i,j)),cv(S['matrix'][i,j])==D['matrix'][i,j])
ck('literal tau column deleted',tuple(D['boundary_columns'])==(0,1,2,4,5,6,7,8))
formulas=[l**5*(lam-N)**2*(N+lam)*(5*lam-2*N)*(5*lam+4*N)/(42*r**8*a**3*m**3),-l**6*(lam-N)**2*(N+lam)*(11*lam-5*N)*(5*lam+4*N)/(168*r**7*a**3*q*m**3)]
for idx,it in enumerate((D,Alt)):
 rows=list(it['pivot_rows'] if idx==0 else it['rows']);A=D['matrix'].matrix_from_rows_and_columns(rows,D['boundary_columns']);saved=it['minor_matrix'] if idx==0 else it['matrix'];ck('full source rows '+str(idx),A==saved)
 rp=it['phase_rows'] if idx==0 else it['row_phase'];cp=it['phase_cols'] if idx==0 else it['column_phase'];B=it['unphased']
 for i in range(8):
  for j in range(8):ck('actual full l-phase '+str((idx,i,j)),A[i,j]==B[i,j]*l**(rp[i]+cp[j]))
 R=K.ring();den=[lcm([c.denominator() for c in row]) for row in B.rows()];C=matrix(R,[[R(den[i]*B[i,j]) for j in range(8)] for i in range(8)]);f=K(lap(C))/prod(den)*l**(sum(rp)+sum(cp));ck('independent sparse Leibniz exact '+str(idx),f==formulas[idx]);ck('saved determinant exact '+str(idx),f==(it['minor'] if idx==0 else it['scalar_det']*l**it['phase']))
ck('two amplitude factors exact Bezout',5*(5*lam-2*N)-2*(11*lam-5*N)==3*lam)
# All positive label orders, not merely 1<a<b.
V=PolynomialRing(QQ,names=('aa','bb')).fraction_field();aa,bb=V.gens();rv=(aa-bb)/(aa*(bb-1));mv=(aa-1)/(bb-1);qv=bb*(aa-1)/(aa*(bb-1));ck('N all positive spectra formula',aa*rv**2==(aa-bb)**2/(aa*(bb-1)**2));ck('lambda all positive spectra formula',mv*qv==bb*(aa-1)**2/(aa*(bb-1)**2));ck('GP equivalence',aa*rv+qv==(aa**2-bb)/(aa*(bb-1)))
# The GP block is pulled afresh from the already source-audited full A6.
G=load(str(p/'ordinary_cyclic_GP_boundary_five_minors.sobj'));Orig=load(str(p/'ordinary_pair_minus_wall_A6_determinant.sobj'));C=G['matrix'];F=C.base_ring();rr=F.gen();qq=1+rr;OLD=Orig['unphased'].base_ring();oa,orr=OLD.gens();gph=OLD.ring().hom([-qq/rr,rr],F)
for i in range(6):
 for j in range(6):
  z=Orig['unphased'][i,j];ck('GP original A6 '+str((i,j)),C[i,j]==gph(z.numerator())/gph(z.denominator()))
ck('GP complete one-dimensional kernel source',C*vector(F,[1,F(1)/2,0,0,0,0])==vector(F,6))
rs=(0,1,2,3,4);cc=(0,2,3,4,5);B=C.matrix_from_rows_and_columns(rs,cc);R=F.ring();den=[lcm([z.denominator() for z in row]) for row in B.rows()];BC=matrix(R,[[R(den[i]*B[i,j]) for j in range(5)] for i in range(5)]);det5=F(lap(BC))/prod(den);ck('GP true five-column Leibniz',det5==F(160)/9*qq**3/rr)
phase=sum(Orig['unit_phase_rows'][i] for i in rs)+sum(Orig['unit_phase_cols'][j] for j in cc);ck('GP retained phase four',phase==4)
# Verify final B from the actual four coefficient rows, retaining all ranges.
BD=load(str(p/'ordinary_pair_minus_adapted_B.sobj'));BP=BD['ring'];BK=BP.base_ring();ga,gr,gH,gZ=BK.gens();bp=BK.ring().hom([-qq/rr,rr,1,-1],F);PE=PolynomialRing(F,'ell');ell=PE.gen();RE=PE.quotient(ell**2+qq/rr,'e');ee=RE.gen()
def gcv(z):return sum((bp(c.numerator())/bp(c.denominator())*ee**int(ex) for ex,c in BP(z).dict().items()),RE.zero())
BB=matrix(RE,[[gcv(z) for z in row] for row in BD['final_B_matrix'].rows()]);db=lap(BB);ck('GP original four B-row determinant',db==-1536*qq**2/rr**3)
unit=F(3)/56;product=unit*db*det5*ee**4;ck('full ten-generator block product',product==-F(10240)/7*qq**7/rr**6)
# Complete zero blocks are original first/even/third source identities in
# the previous whole-negative-wall certificates; here also check the
# source labels and deleted column, so no alternative A block is used.
ck('A source row identities',Orig['source_rows']==[5,7,11,13,37,39]);ck('A actual columns',Orig['columns']==[0,3,5,7,9,10]);ck('B original row identities',tuple(BD['final_B_rows'])==(2,17,29,35))

# Literal original coefficient zero blocks needed for the ten-column
# Laplace argument (no assumption about B rows' later A entries).
ol=OR.gen();oa,orr,oH,oZ=OK.gens();om=1+oa*orr;onu=om*ol*ol+oa*orr;pp=OK.ring().hom([-qq/rr,rr,1,-1],F)
def ocv(z):return sum((pp(c.numerator())/pp(c.denominator())*ee**int(ex) for ex,c in OR(z).dict().items()),RE.zero())
ck('original u eta pivot',ocv(S['matrix'][0,4])==F(3)/56)
for row in (2,17,29,35,5,7,11,13,37):ck('original lower ueta zero '+str(row),ocv(S['matrix'][row,4])==0)
for row in (5,7,11,13,37):
 ff=S['matrix'][row,8]-om*S['matrix'][row,1]-om*ol*S['matrix'][row,2];quo,rem=ff.quo_rem(onu);ck('original exact nu polynomial quotient '+str(row),rem==0)
 for j in (1,2,6):ck('literal A/B zero '+str((row,j)),ocv(S['matrix'][row,j])==0)
 ck('literal A nu quotient zero '+str(row),ocv(quo)==0)

rec={'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st,'scope':'All-spectrum cyclic E0 boundary locally contained in whole incoming family; no claim that all tau-nonzero arcs through special points are excluded.'};save({'checks':checks,'nonGP_formulas':formulas,'GP_block_product':product,'scope':rec['scope']},str(p/'ordinary_cyclic_E0_boundary_containment_reader.sobj'),compress=False);(p/'ordinary_cyclic_E0_boundary_containment_reader.json').write_text(json.dumps(rec,indent=2)+'\n');print(rec,flush=True)
