"""Actual characteristic-curve tangencies, keeping every parameter c.
No conclusion on the characteristic branch is asserted by this source file.
"""
from pathlib import Path
import pickle,json,time
from sage.all import QQ,PolynomialRing,matrix,vector,save
out=Path(__file__).resolve().parent;base=out.parent;start=time.monotonic()
old=pickle.loads((base/'round8/rank3_source_minor_candidates.sobj').read_bytes())
C=PolynomialRing(QQ,names=('a','b','c','h'));K=C.fraction_field();a,b,c,h=map(K,C.gens())
R=PolynomialRing(K,names=('p12','p13','p21','p23','p31','t0','t1','t2','t3','t4','t5'))
p12,p13,p21,p23,p31,*tt=R.gens();r=(a-b)/(a*(b-1));q=1+r;z=q*h+c;p32=-z*p23/(q*h)
imgs=[R.one(),p12,p13,p21,R(h),p23,p31,p32,R(z),*tt,R.zero(),R.zero()]
def ck(co):
 def pp(p):return sum((K(v)*a**m[0]*b**m[1] for m,v in p.dict().items()),K.zero())
 return pp(co.numerator())/pp(co.denominator())
def cv(p):
 ans=R.zero()
 for m,co in p.dict().items():
  term=R(ck(co))
  for x,e in zip(imgs,m):
   if e:term*=x**e
  ans+=term
 return ans
F=[[cv(old['F'][i][j]) for j in range(9)] for i in range(3)]
DL=[F[i][8]-q*F[i][4]-c*F[i][0] for i in range(3)]
DB=[F[i][8]*p23+z*F[i][5]+q*F[i][4]*p32+q*h*F[i][7] for i in range(3)]
C3=[cv(f) for f in old['C3']]
assert not DL[0]
sources=C3+DL[1:]+DB
zero={t:R.zero() for t in tt}
M=matrix(R,[[p.derivative(t) for t in tt] for p in sources]);constants=vector(R,[p.subs(zero) for p in sources])
assert all(sources[i]==constants[i]+sum(M[i,j]*tt[j] for j in range(6)) for i in range(8))
minor=M.matrix_from_rows(range(6)).det()
payload={'ring':R,'parameters':(a,b,c,h),'variables':R.gens(),'images':imgs,'z':z,'p32':p32,
 'sources':sources,'source_names':['C3_0','C3_1','C3_2','F2L','F3L','F1B1','F2B1','F3B1'],
 't_coefficient_matrix':M,'constant_terms':constants,'selected_minor':minor,
 'Fpi_on_graph':F,'DL':DL,'DB':DB,'scope':'All-spectral source tangencies for L=z-qh-c=0, B1=0. No t-pivot division has been made.'}
save(payload,str(out/'rank1_ratio_characteristic_source.sobj'))
print('saved_sources',round(time.monotonic()-start,2),flush=True)
fac=minor.factor()
report={'status':'PASS_EXACT_SOURCE_GRAPH','first_direction_L_identity':True,'all_eight_sources_affine_in_t':True,
 'selected_minor_factorization':str(fac),'source_term_counts':[len(p.dict()) for p in sources],
 'elapsed_seconds':time.monotonic()-start,'scope':payload['scope']}
(out/'rank1_ratio_characteristic_source.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2),flush=True)
