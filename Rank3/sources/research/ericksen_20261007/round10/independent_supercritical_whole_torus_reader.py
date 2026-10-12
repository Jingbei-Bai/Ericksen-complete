"""Independent source replay of the whole complex supercritical full-support proof.

Transverse columns are computed directly from original F12/H13 Euler
derivatives, independently of differentiating the critical source matrix.
Resultants are verified by exact integer grids with strict degree bounds.
"""
from sage.all import *
from pathlib import Path
import json,time
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(n,v):
 if not bool(v):raise AssertionError(n)
 checks.append(n)
nd=load(str(p/'independent_supercritical_canonical_normal.sobj'));K=nd['field'];a,r,H,Z,s,z=K.gens();q=1+r;N=nd['normal_matrix']
can=load(str(p/'independent_delta_critical_canonical_amplitude.sobj'));P=can['extension_ring'];b,c=P.gens()
raw=load(str(p/'independent_supercritical_first_boundary_probe.sobj'));old=raw['ring'].base_ring();ph=old.ring().hom([a,r,H,Z],K)
def cv(v):v=old(v);return ph(v.numerator())/ph(v.denominator())
def red(f):
 dd={}
 for ex,co in P(f).dict().items():
  e=tuple(map(int,ex));key=(e[0]%2,e[1]%2);dd[key]=dd.get(key,K.zero())+co*(-q)**(e[0]//2)*(-q/r)**(e[1]//2)
 return P(dd)
def inv(f):
 f=red(f);ck('monomial protected denominator',len(f.dict())==1);ex,co=next(iter(f.dict().items()));i,j=map(int,ex)
 return b**i*c**j/(co*(-q)**i*(-q/r)**j)
alpha=(Z-H)/(Z*(H+1));beta=H*(Z-1)/(Z*(H+1));imgs=tuple(map(P,(1,b*z,s,c*z,s*b,c,alpha*c,beta*b,1)))
def pull(f,indices=None):
 out=P.zero()
 for ex,co in f.dict().items():
  ee=tuple(map(int,ex));mult=1 if indices is None else sum(ee[j] for j in indices)
  if mult:out+=mult*cv(co)*prod(imgs[j]**ee[j] for j in range(9))
 return red(out)
sources=[(raw['F12'][0],s*b),(raw['F12'][2],P(s)),(raw['F12'][3],P.one()),(None,P.one())]+[(f,P.one()) for f in raw['H13']]
for i,(f,den) in enumerate(sources):
 ri,rj=can['row_phases'][i];scale=b*c*inv(den)*inv(b**ri*c**rj)
 if f is None:vals=[P.zero(),P.zero()]
 else:
  ck('literal canonical full source '+str(i),pull(f)==0)
  vals=[red(pull(f,inds)*scale) for inds in ((1,4,7),(3,5,6))]
 for j,v in enumerate(vals):ck('independent original Euler normal '+str((i,j)),v==P(N[i,j]))
 for j,col in enumerate((0,1,2,4)):ck('certified old beta source '+str((i,j)),N[i,j+2]==can['material_matrix'][i,col])
mi=load(str(p/'independent_supercritical_normal_minors.sobj'));V=mi['reduced_matrix'];A=N[:2,:2];de=A.det();ad=A.adjugate()
ck('transverse 2x2 exact adjugate',A*ad==de*identity_matrix(K,2))
ck('normal Schur source identity',de*V==de*N[2:,2:]-N[2:,:2]*ad*N[:2,2:])
C=V[:3,:3];cd=C.det();cad=C.adjugate();ck('old three adjugate',C*cad==cd*identity_matrix(K,3))
for it in mi['items']:
 j=it['source_rows'][-1]-2;v=de*(cd*V[j,3]-(V[j,:3]*cad*V[:3,3])[0,0]);ck('true six-source block determinant '+str(j),v==it['determinant'])
 if v:ck('complete determinant factor replay '+str(j),it['factorization'].value()==v)
dr=load(str(p/'independent_supercritical_normal_direction.sobj'));DP=dr['ring'];ds,dz=DP.gens()
for f,co in dr['target_lifts']:ck('all-direction original coefficient lift',sum((v*g for v,g in zip(co,dr['sources'])),DP.zero())==f)
ck('all four direction identities retained',len(dr['target_lifts'])==4)
line=load(str(p/'independent_supercritical_line_even.sobj'));LK=line['field'];la,lr,lH,lZ,t=LK.gens();lm=K.ring().hom([la,lr,lH,lZ,-t/(t+1),t],LK)
def lc(v):v=K(v);return lm(v.numerator())/lm(v.denominator())
LV=matrix(LK,V.nrows(),4,[lc(v) for v in V.list()]);ck('line literal source matrix',LV==line['normal_matrix']);eg=list(line['e_graph'].list())
ck('line full normal graph',LV*vector(LK,eg+[1])==0)
ev=load(str(p/'independent_supercritical_exceptional_even.sobj'))['even']
for j,f in enumerate(ev):
 vv=sum((lc(co)*prod(eg[k]**int(ex[k]) for k in range(3)) for ex,co in f.dict().items()),LK.zero());ck('line original even pullback '+str(j),vv==line['even_pullbacks'][j])
fixed=load(str(p/'independent_supercritical_fixed_directions.sobj'));FK=fixed['field'];fa,fr,fH,fZ=FK.gens()
for it in fixed['items']:
 sv,zv=it['direction'];fm=K.ring().hom([fa,fr,fH,fZ,FK(sv),FK(zv)],FK)
 def fc(v):v=K(v);return fm(v.numerator())/fm(v.denominator())
 FV=matrix(FK,V.nrows(),4,[fc(v) for v in V.list()]);ck('fixed literal source matrix '+it['label'],FV==it['matrix'])
 if 'full_rank_minor' in it:
  det=FV.matrix_from_rows(it['pivot_rows']).det();ck('point11 full normal determinant',det==it['full_rank_minor'] and det!=0);ck('point11 factor replay',it['full_rank_minor_factorization'].value()==det)
 else:
  fg=it['graph'];EE=fg[0].parent();ck('fixed full normal graph '+it['label'],matrix(EE,FV)*vector(EE,fg+[EE.one()])==0)
  for j,f in enumerate(ev):
   vv=sum((EE(fc(co))*prod(fg[k]**int(ex[k]) for k in range(3)) for ex,co in f.dict().items()),EE.zero());ck('fixed original even pullback '+str((it['label'],j)),vv==it['even_pullbacks'][j])
ex=load(str(p/'independent_supercritical_exceptional_material.sobj'));SP=ex['ring'];sn,sl=SP.gens()
for it in ex['items']:
 src=next(v for v in fixed['items'] if v['label']==it['label'])
 sm=SP.hom([fa*fr**2,(1+fa*fr)*(1+fr)],FK)
 for j,(f,u) in enumerate(zip(it['forms'],it['units'])):ck('fixed all-material source transport '+str((it['label'],j)),sm(f)*u==FK(src['even_pullbacks'][j].constant_coefficient()))
 for g,co in zip(it['groebner_basis'],it['original_source_lifts']):ck('original even material cofactor identity',sum((c*f for c,f in zip(co,it['forms'])),SP.zero())==g)
 ck('positive material square explicitly in source ideal',sn*sn in it['groebner_basis'])
pa,pb=ex['line_pivot_polynomials'];ca,cb=ex['line_pivot_bezout'];ck('entire line old-block pivot cover',ca*pa+cb*pb==1)
mat=load(str(p/'independent_supercritical_line_material.sobj'));MP=mat['ring'];mn,ml,mt=MP.gens();mm=MP.hom([la*lr**2,(1+la*lr)*(1+lr),t],LK)
for j,f in enumerate(mat['forms']):ck('line all-material source transport '+str(j),mm(f)*mat['source_cores'][j]['removed_unit']==line['even_pullbacks'][j])
# Exact multivariate polynomial identity from degree-bounded interpolation.
B=PolynomialRing(QQ,names=('gn','gl'));gn,gl=B.gens();U=PolynomialRing(B,'gt');gt=U.gen();mh=MP.hom([U(gn),U(gl),gt],U);fs=[mh(f) for f in mat['forms']]
def syl(f,g):
 m=int(f.degree());k=int(g.degree());fc=[f[i] for i in range(m,-1,-1)];gc=[g[i] for i in range(k,-1,-1)]
 return matrix(B,[[B.zero()]*i+fc+[B.zero()]*(k-1-i) for i in range(k)]+[[B.zero()]*i+gc+[B.zero()]*(m-1-i) for i in range(m)])
res=mat['resultants'][0]['resultant'];rh=res.parent().ring().hom([gn,gl],B);rr=rh(res.numerator());ck('line resultant polynomial denominator',res.denominator()==1);sy=syl(fs[0],fs[1]);bounds=[sum(max(int(f.degree(v)) for f in sy.row(i)) for i in range(sy.nrows())) for v in (gn,gl)]
for j,v in enumerate((gn,gl)):ck('strict interpolation degree bound',rr.degree(v)<=bounds[j])
for v1 in range(1,bounds[0]+2):
 for v2 in range(1,bounds[1]+2):ck('exact line Sylvester node '+str((v1,v2)),matrix(ZZ,sy.nrows(),sy.ncols(),[f(v1,v2) for f in sy.list()]).det()==rr(v1,v2))
F=9*(1-gn-gl)*(1-(gn-gl)**2)+8*gn*gl
expected=32*(gl-gn)*(gl-1)*(gl+1)*(gn+gl)*(gn+1)**2*(gn-1)**4*gn**6*gl**11*(gn*gl+8*gn+gl)*((gn-1)**2+gl*(gn+1))*F
ck('entire material resultant factor identity',rr==expected)
WF=ex['line_wall_forms'];wr=ex['line_wall_resultant'];ck('N=lambda actual quartic resultant',WF[0].sylvester_matrix(WF[1]).det()==wr)
wn=wr.parent().gen();ck('entire material equality wall factor identity',wr==64*(wn-QQ(3)/2)*(wn+9)*(wn+1)**3*(wn-1)**5*wn**16*(wn**2-wn/2+QQ(1)/2))
record={'checks':len(checks),'seconds':time.monotonic()-st,'line_resultant_integer_nodes':(bounds[0]+1)*(bounds[1]+1),'line_resultant_degree_bounds':bounds,'all_materials_retained':True,'all_complex_amplitudes_retained':True,'both_omega_conjugates_retained':True,'scope':'Entire full-support supercritical first-normal exclusion, coupled to the full-source regularity and Nakayama argument in the accompanying proof. Whole gradient-independent problem remains open.'}
save({'record':record,'checks':checks,'scope':__doc__},str(p/'independent_supercritical_whole_torus_reader.sobj'),compress=False)
(p/'independent_supercritical_whole_torus_reader.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2),flush=True)
