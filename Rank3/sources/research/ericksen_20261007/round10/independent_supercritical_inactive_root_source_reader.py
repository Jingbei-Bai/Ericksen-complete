"""Root independent replay from full R/Q and raw signed-product H DAGs.

Evaluate in the canonical quadratic root algebra, without the creator's
phase-pulling routine. Generic first-five, shifted whole line, and cyclic
full two-jet are checked separately. Every lower Q0 layer is inspected.
"""
from sage.all import *
from pathlib import Path
import sys,time,json,os,resource
p=Path(__file__).resolve().parent;st=time.monotonic();mode=sys.argv[1];assert mode in ('generic','cyclic');N=1 if mode=='generic' else 2
checks=[]
def ck(s,b):
 assert b,s
 checks.append(s)
qdata=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'));full=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'));init=load(str(p/'independent_delta_subcritical_H_initials.sobj'));ev=load(str(p/'independent_delta_subcritical_even_full_source.sobj'))
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+850*1024**2,)*2)
K=PolynomialRing(QQ,names=('audit_a','audit_r','audit_H','audit_Z')).fraction_field() if mode=='generic' else PolynomialRing(QQ,names=('audit_a','audit_r')).fraction_field()
a,r=K.gens()[:2];H,Z=K.gens()[2:] if mode=='generic' else (K(1),K(-1));q=1+r;m=1+a*r;al=(Z-H)/(Z*(H+1));bb=H*(Z-1)/(Z*(H+1));d=-2*H*(Z+1)/(Z*(H+1))
P=PolynomialRing(K,names=('root_b','root_c','ns','nw','ne0','ne1','ne2','nsigma'));b,c,s,w,e0,e1,e2,sig=P.gens();normals=P.gens()[2:]
def red(f):
 out={}
 for ee,co in P(f).dict().items():
  e=tuple(map(int,ee))
  if sum(e[2:])>N:continue
  ex=(e[0]%2,e[1]%2)+e[2:];out[ex]=out.get(ex,K(0))+co*(-q)**(e[0]//2)*(-q/r)**(e[1]//2)
 return P(out)
def mul(x,y):return red(x*y)
def invroot(f):
 ff=red(f);ck('invert coherent root unit '+str(len(checks)),len(ff.dict())==1)
 e,co=next(iter(ff.dict().items()));e=tuple(map(int,e));ck('root unit has no normals '+str(len(checks)),sum(e[2:])==0)
 return red(b**e[0]*c**e[1]/(co*(-q)**e[0]*(-q/r)**e[1]))
sigactual=mul(sig,invroot(b*c));es=(c*(e0+d),b*c*e1,b*e2)
imgs=(P(1),b*s,w,c*s,b*w,c,al*c,bb*b,P(1))
nk=qdata['ring'].base_ring();ph=nk.ring().hom([a,r,H,Z],K);old=full['coefficient_field'];oph=old.ring().hom([a,a*q/m,H/r,q*Z/r],K);cc={};oc={};mc={}
def cv(co,orig=False):
 f=(old if orig else nk)(co);cache=oc if orig else cc;hp=oph if orig else ph
 if f not in cache:cache[f]=hp(f.numerator())/hp(f.denominator())
 return cache[f]
def mon(e):
 e=tuple(map(int,e))
 if e not in mc:
  v=P(1)
  for j,ex in enumerate(e):
   if ex:v=mul(v,imgs[j]**ex)
   if not v:break
  mc[e]=v
 return mc[e]
def pull(dd,orig=False,euler=None,derivative=None):
 out=P(0)
 for ee,co in dd.items():
  e=tuple(map(int,ee));vv=P(0)
  if derivative is None:vv=mon(e)
  else:
   for j,df in derivative.items():
    if e[j]:vv+=mul(e[j]*df,mon(tuple(e[k]-(k==j) for k in range(9))))
  if not vv:continue
  if euler is not None:vv*=sum(e[j] for j in euler)
  if vv:out+=cv(co,orig)*vv
 return red(out)
def dd(f):return f.dict() if hasattr(f,'dict') else f
def face(f,bi):return {tuple(map(int,e)):co for e,co in dd(f).items() if (sum(e[:6]),sum(e[6:]))==bi}
def module(cfs,orig=False):return red(sum((mul(pull(dd(f),orig),es[j]) for j,f in enumerate(cfs)),P(0)))
def logmodule(cfs,orig,lg):
 weights=(0,1,1) if lg==0 else (1,1,0)
 return red(sum((mul(pull(dd(f),orig,euler=logs[lg])+weights[j]*pull(dd(f),orig),es[j]) for j,f in enumerate(cfs)),P(0)))
pool=full['polynomial_pool'];pc={};degrees={i:set((sum(e[:6]),sum(e[6:])) for e in f) for i,f in enumerate(pool)} if isinstance(pool,list) else {i:set((sum(e[:6]),sum(e[6:])) for e in f) for i,f in pool.items()}
def poolpart(i,bi,euler=None,derivative=None):
 key=(i,bi,tuple(euler) if euler else None,tuple(derivative) if derivative else None)
 if key not in pc:pc[key]=pull(face(pool[i],bi),True,euler,derivative)
 return pc[key]
def hface(j,bi,euler=None,derivative=None):
 out=P(0)
 for sign,li,ri in full['records'][j]['full_signed_product_DAG']:
  for dg in degrees[li]:
   rg=(bi[0]-dg[0],bi[1]-dg[1])
   if rg not in degrees[ri]:continue
   if euler or derivative:
    out+=sign*(mul(poolpart(li,dg,euler,derivative),poolpart(ri,rg))+mul(poolpart(li,dg),poolpart(ri,rg,euler,derivative)))
   else:out+=sign*mul(poolpart(li,dg),poolpart(ri,rg))
 return red(out)
def const(f):return red(P({e:co for e,co in red(f).dict().items() if sum(e[2:])==0}))
def strip(f):
 f=red(f);phases=set(tuple(e[:2]) for e in f.dict());ck('coherent source '+str(len(checks)),len(phases)<=1)
 if not phases:return P(0)
 i,j=next(iter(phases));return red(mul(f,invroot(b**i*c**j)))
logs=((1,4,7),(3,5,6));DC=qdata['D']*qdata['C'];AJ=[];Rbeta=[]
dk={1:b*(e1-s*e0),2:e2,6:c*e0}
for j in (1,3):
 f12=face(qdata['full_R_sources'][j],(1,2));AJ.append([pull(f12,euler=lg) for lg in logs]);v=mul(b*c,module(DC.row(j))+pull(face(qdata['full_R_sources'][j],(4,0))))
 if mode=='cyclic':v+=mul(sig,pull(face(qdata['full_R_sources'][j],(2,1)))+pull(f12,derivative=dk))
 Rbeta.append(red(v))
det=red(AJ[0][0]*AJ[1][1]-AJ[0][1]*AJ[1][0]);ide=invroot(det);Rsol=[-mul(ide,mul(AJ[1][1],Rbeta[0])-mul(AJ[0][1],Rbeta[1])),-mul(ide,-mul(AJ[1][0],Rbeta[0])+mul(AJ[0][0],Rbeta[1]))]
ck('true R range solution0',red(mul(AJ[0][0],Rsol[0])+mul(AJ[0][1],Rsol[1])+Rbeta[0])==0)
ck('true R range solution1',red(mul(AJ[1][0],Rsol[0])+mul(AJ[1][1],Rsol[1])+Rbeta[1])==0)
Q0=module(qdata['Q0_L_cofactors'])+pull(face(qdata['full_Q0'],(2,1)))
q02=pull(face(qdata['full_Q0'],(0,2)));q30=pull(face(qdata['full_Q0'],(3,0)))
if mode=='cyclic':
 ck('Q0 pure RR2 cyclic exactly zero',q02==0);Q0+=mul(sigactual,q02)
else:ck('Q0 pure RR2 generalbase retained as sigma layer',q02!=0)
fs=[strip(Q0)];labels=['Q0'];linechecks=[]
if mode=='generic':
 ck('Q0 shifted whole line',const(Q0)==0);ck('Q0 pure p3 whole line',const(q30)==0)
 for j,f in enumerate(qdata['full_R_sources']):
  ck('raw R12 line '+str(j),const(pull(face(f,(1,2))))==0);ck('raw R40 line '+str(j),const(pull(face(f,(4,0))))==0);ck('raw DCe shifted line '+str(j),const(module(DC.row(j)))==0)
  for lg in range(2):ck('R beta transverse line derivative %d %d'%(j,lg),const(pull(face(f,(4,0)),euler=logs[lg])+logmodule(DC.row(j),False,lg))==0)
for j,it in enumerate(init['items']):
 h13=hface(j,(1,3));HJ=[hface(j,(1,3),euler=lg) for lg in logs];h41=hface(j,(4,1));hm=module(it['L_cofactors_original'],True)
 ck('raw DAG H41 agrees cached face '+str(j),h41==pull(it['H41_original'],True))
 hh=mul(b*c,hm+h41)
 if mode=='cyclic':hh+=mul(sig,hface(j,(2,2))+hface(j,(1,3),derivative=dk))+mul(sig**2,mul(invroot(b*c),hface(j,(0,3))))
 hh+=sum((mul(HJ[k],Rsol[k]) for k in range(2)),P(0));fs.append(strip(hh));labels.append('H'+str(j))
 if mode=='generic':
  ck('raw H13 shifted line '+str(j),const(h13)==0);ck('raw H41+module shifted line '+str(j),const(h41+hm)==0);ck('raw H50 whole line '+str(j),const(hface(j,(5,0)))==0)
  ck('raw H13 entire canonical amplitudes '+str(j),h13==0)
  for lg in range(2):ck('H beta transverse line derivative %d %d'%(j,lg),const(hface(j,(4,1),euler=logs[lg])+logmodule(it['L_cofactors_original'],True,lg))==0)
 print('H audited',mode,j,'seconds',round(time.monotonic()-st,2),flush=True)
for it in ev['items']:
 j=int(it['square_component']);v=red(cv(it['square_unit'])*es[j]**2+module(it['A_p_quadratic'])+pull(dd(it['pure_p4'])))
 if mode=='cyclic':
  f=it.get('full_Q',it.get('full_source'))
  for ee,co in f.dict().items():
   e=tuple(map(int,ee));rd=sum(e[6:])
   if sum(e[:6])+2*rd==4 and rd:v+=mul(sigactual**rd,pull({e:co}))
 if mode=='generic':ck('raw even shifted line '+str(it['Q_index']),const(v)==0)
 fs.append(strip(v));labels.append('Q'+str(it['Q_index']))
creator=load(str(p/('independent_supercritical_inactive_generic_first5.sobj' if mode=='generic' else 'independent_supercritical_inactive_saturated_twojet.sobj')));cp=creator['ring'];ck('source labels exact',labels==creator['labels']);ckf=creator['field'];cph=ckf.ring().hom(list(K.gens()),K)
def creator_pull(f):
 return P({(0,0)+tuple(map(int,e)):cph(co.numerator())/cph(co.denominator()) for e,co in f.dict().items()})
for j,(x,y) in enumerate(zip(fs,creator['sources'])):ck('independent fullsource jet equals creator '+labels[j],x==creator_pull(y))
if mode=='cyclic':
 J=matrix(K,[[f.monomial_coefficient(x) for x in normals] for f in fs]);rows=[0,2,4,8,19,20];six=J.matrix_from_rows(rows);det6=six.det();ck('sixnormal determinant factors',det6!=0)
 rrK=PolynomialRing(QQ,'audit_res_r').fraction_field();rr=rrK.gen();rp=K.ring().hom([-(1+rr)/rr,rr],rrK)
 def rc(z):return rp(z.numerator())/rp(z.denominator())
 RJ=matrix(rrK,[[rc(z) for z in row] for row in J.rows()]);ker=vector(rrK,[1,-1,-4,-2,-2,0]);corr=vector(rrK,[-6,13,38,38,0,-4*(1+rr)]);ck('all22 true source kernel',RJ*ker==0)
 back=rrK.ring().hom([r],K)
 kback=[back(t.numerator())/back(t.denominator()) for t in ker]
 quadratic=[]
 for f in fs:
  z=K(0)
  for ee,co in f.dict().items():
   e=tuple(map(int,ee));assert e[0]==e[1]==0
   if sum(e[2:])==2:z+=co*prod(kback[j]**e[j+2] for j in range(6))
  quadratic.append(rc(z))
 residue=vector(rrK,quadratic)+RJ*corr;range5=[0,2,4,8,19];cols5=[2,0,1,3,5];ck('range5 all quadratic zero',all(residue[j]==0 for j in range5));ck('range5 protected determinant',RJ.matrix_from_rows_and_columns(range5,cols5).det()==5*(1+rr)**6/(5376*rr**5));ck('true Q5 compensated k2 unit',residue[20]==-16*(1+rr)*(rr**2+rr+QQ(13)/4)/rr)
 print('SIX NORMAL DET',det6.factor(),flush=True);print('TRUE RESONANT Q5',residue[20].factor(),flush=True)
else:
 J=matrix(K,[[f.monomial_coefficient(x) for x in normals[:5]] for f in fs]);rK=PolynomialRing(QQ,names=('root_res_r','root_res_H','root_res_Z')).fraction_field();rr,h,z=rK.gens();rp=K.ring().hom([-(1+rr)/rr,rr,h,z],rK);RJ=matrix(rK,[[rp(t.numerator())/rp(t.denominator()) for t in row] for row in J.rows()]);B=h*(z-1)/(z*(h+1));ker=vector(rK,[1/B,-1/B,-2*(h*z-h+2*z)/(h*(z-1)),-2/B,-2]);ck('generic all22 fullsource firstfive kernel',RJ*ker==0)
record={'status':'PASS','checks':len(checks),'mode':mode,'seconds':time.monotonic()-st,'full_raw_DAG_replay':True,'Q02_cyclic_zero':mode=='cyclic','generic_Q02_sigma_layer_retained':mode=='generic'}
save({'record':record,'checks':checks,'field':K,'ring':P,'sources':fs,'labels':labels,'J':J},str(p/('independent_supercritical_inactive_root_source_reader_'+mode+'.sobj')),compress=False)
(p/('independent_supercritical_inactive_root_source_reader_'+mode+'.json')).write_text(json.dumps(record,indent=2)+'\n');print('PASS',record,flush=True)
