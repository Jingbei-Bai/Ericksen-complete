"""COMPLETE generic H/Z seven-original-source kappa-zero identity.

Build only original source layers contributing kappa^0 or kappa^1; all
higher kappa layers remain in the explicitly specified COMPLETE original
source DAG. The entire kappa-zero identity makes the full division legal.
No state jet or material specialization is used to prove that identity.
"""
from sage.all import *
from pathlib import Path
import os,resource,time,json,gc
p=Path(__file__).resolve().parent;st=time.time();checks=[]
qp=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'))
hs=load(str(p/'independent_delta_subcritical_H_initials.sobj'))
full=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'))
H1=load(str(p/'pure_x_root_cyclic_h_source_H1.sobj'))
H2=load(str(p/'subcritical_single_pure_x_linear_kappa_syzygies_H2.sobj'))
vm=int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE');resource.setrlimit(resource.RLIMIT_AS,(vm+600*1024**2,vm+600*1024**2))
print('LinuxPID',os.getpid(),'VMcap_MB',vm//1024**2+600,flush=True)
K=PolynomialRing(QQ,names=('k','r','H','Z')).fraction_field();k,r,H,Z=K.gens();m=1+k;q=1+r
P=PolynomialRing(K,names=('y','z','w','e0','e1','e2','delta','t','v','kap'));y,z,w,e0,e1,e2,de,t,v,kap=P.gens()
al=(Z-H)/(Z*(H+1));be=H*(Z-1)/(Z*(H+1));ga=al/be
inc0=y*v-ga*w*t;inc=inc0-kap*e1
images0=(P.one(),y,z,w,t*z/be,v/al,v,t,P.one())
es=(e0,e1,e2);QO=qp['ring'];qph=QO.base_ring().ring().hom([k/r,r,H,Z],K)
old=full['coefficient_field'];oph=old.ring().hom([k/r,k*q/(r*m),H/r,q*Z/r],K);converted={}
def ck(n,b):
 if not b:raise AssertionError(n)
 checks.append(n)
def cvOld(c):
 if c not in converted:converted[c]=oph(c.numerator())/oph(c.denominator())
 return converted[c]
def cvQ(c):return qph(c.numerator())/qph(c.denominator())
def nt(f):return f.parent()({tuple(map(int,e)):K(c) for e,c in f.dict().items() if c})
def project(dd,converter,normal_shift=0):
 """Actual FULL p-chart pullback through kappa1, all remaining states kept."""
 out=P.zero()
 for ex,c in dd.items():
  ex=tuple(map(int,ex));cc=converter(c)
  if not cc:continue
  val=P(cc)*prod(images0[j]**ex[j] for j in range(9) if ex[j])
  if normal_shift==0:
   out+=val
   for j,direction in [(4,-e2/be),(5,-e0/al)]:
    if ex[j]:out+=kap*ex[j]*direction*P(cc)*prod(images0[h]**(ex[h]-(1 if h==j else 0)) for h in range(9) if ex[h]-(1 if h==j else 0))
  elif normal_shift==1:out+=kap*val
  else:raise AssertionError(normal_shift)
 return nt(out)
def split(dd):
 out={}
 for ee,c in dd.items():out.setdefault((sum(ee[:6]),sum(ee[6:])),{})[tuple(map(int,ee))]=c
 return out
raws={};firsttwo={}
DC=matrix(QO,qp['D'])*matrix(QO,qp['C'])
R5=sum((project(DC[5,j].dict(),cvQ)*es[j] for j in range(3)),P.zero())
rfaces=split(qp['full_R_sources'][5].dict())
for face,dd in rfaces.items():
 if face==(3,1):continue
 ke=8-2*face[0]-3*face[1];dp=4-face[0]-face[1]
 ck('original R5 source layer regular '+str(face),ke>=0 and dp>=0)
 if ke<=1:R5+=de**dp*project(dd,cvQ,ke)
firsttwo[5]=nt(R5);raws[5]={'full_R5':qp['full_R_sources'][5],'DC_row5':DC.row(5)}
print('ORIGINAL R5 first-two complete layers',len(R5.dict()),round(time.time()-st,2),flush=True)

# All six H layers come from the COMPLETE signed product DAG. Degree
# (3,2) is replaced by its independently source-certified exact L module.
pool=full['polynomial_pool'];blocks={}
def grouped(ref):
 if ref not in blocks:blocks[ref]=split(pool[ref])
 return blocks[ref]
targets={(4,1),(1,3),(5,0),(2,2)}
for j in [4,9,11,13,16,17]:
 dd={x:{} for x in targets}
 for sign,left,right in full['records'][j]['full_signed_product_DAG']:
  lb=grouped(left);rb=grouped(right)
  for lx,lv in lb.items():
   for tx in targets:
    rx=(tx[0]-lx[0],tx[1]-lx[1])
    if rx not in rb:continue
    out=dd[tx]
    for ex,c in lv.items():
     mc=cvOld(c)
     for ey,d in rb[rx].items():
      ee=tuple(int(x)+int(y) for x,y in zip(ex,ey));out[ee]=out.get(ee,K.zero())+sign*mc*cvOld(d)
 dd={face:{e:c for e,c in di.items() if c} for face,di in dd.items()}
 actual=sum((project(hs['items'][j]['L_cofactors_original'][h],cvOld)*es[h] for h in range(3)),P.zero())
 for face,di in dd.items():
  ke=11-2*face[0]-3*face[1];dp=5-face[0]-face[1]
  ck('H needed full layer regular '+str((j,face)),ke in [0,1] and dp>=0)
  actual+=de**dp*project(di,K,ke)
 originalH41=project(hs['items'][j]['H41_original'],cvOld)
 ck('COMPLETE signed DAG H41 original L archive multiplyback '+str(j),project(dd[(4,1)],K)==originalH41)
 idx=7+j;firsttwo[idx]=nt(actual)
 raws[idx]={'original_signed_record':full['records'][j],'actual_layers_first_two':dd,'original_L_cofactors':hs['items'][j]['L_cofactors_original']}
 save({'ring':P,'source_first_two':firsttwo,'raw_sources':raws,'incidence':inc,'checks':checks},str(p/'pure_x_super_reader_generic_h_sources.partial.sobj'),compress=False)
 print('COMPLETE ORIGINAL H',j,'firsttwo terms',len(actual.dict()),'seconds',round(time.time()-st,2),flush=True)
 gc.collect()

zero=P.hom(list(P.gens()[:9])+[P.zero()],P)
support=[(5,t),(11,P.one()),(16,w),(18,P.one()),(20,y),(23,v),(24,t)]
polys=[nt(mult*firsttwo[idx]) for idx,mult in support]
def nfco(f):
 nf=P.zero();co=P.zero()
 for ex,c in zero(f).dict().items():
  ee=list(map(int,ex));n=min(ee[0],ee[8]);gg=ee.copy();gg[0]-=n;gg[8]-=n;gg[2]+=n;gg[7]+=n;nf+=P({tuple(gg):c*ga**n})
  for h in range(n):
   gg=ee.copy();gg[0]-=1+h;gg[8]-=1+h;gg[2]+=h;gg[7]+=h;co+=P({tuple(gg):c*ga**h})
 ck('COMPLETE incidence normal form identity '+str(len(checks)),zero(f)==nf+co*inc0)
 return nt(nf),nt(co)
pairs=[nfco(f) for f in polys]
save({'ring':P,'source_first_two':firsttwo,'raw_sources':raws,'support':support,'seven_full_first_two':polys,'seven_kappa0_NF_CF':pairs,'incidence':inc,'checks':checks},str(p/'pure_x_super_reader_generic_h_sources.sobj'),compress=False)

# Choose a small genuine generic minor using its H2 fiber ONLY as a
# nonzero-minor witness. The final equation is verified as a FULL generic
# polynomial, so no frozen-fiber ideal or saturation is used.
CB=PolynomialRing(QQ,names=('k','r')).fraction_field();ckk,crr=CB.gens();ph2=K.ring().hom([ckk,crr,2,-1],CB)
def cyc2(c):c=K(c);return ph2(c.numerator())/ph2(c.denominator())
mons=sorted(set(tuple(map(int,e)) for nf,cf in pairs for e in nf.dict()))
M2=matrix(CB,[[cyc2(pairs[j][0].monomial_coefficient(P.monomial(*e))) for e in mons] for j in range(1,7)])
piv=list(M2.pivots());ck('six other sources independent H2',len(piv)==6)
A=matrix(K,[[pairs[j][0].monomial_coefficient(P.monomial(*mons[h])) for j in range(1,7)] for h in piv])
rhs=-vector(K,[pairs[0][0].monomial_coefficient(P.monomial(*mons[h])) for h in piv])
print('GENERIC 6x6 sparse minor ready',round(time.time()-st,2),flush=True)
vec=vector(K,[K.one()]+list(A.solve_right(rhs)))
nf=nt(sum((co*it[0] for co,it in zip(vec,pairs)),P.zero()))
ck('ENTIRE generic HZ kappa-zero syzygy NF exactly ZERO',nf==0)
co=nt(sum((c*it[1] for c,it in zip(vec,pairs)),P.zero()))
num=nt(sum((c*f for c,f in zip(vec,polys)),P.zero())-co*inc)
ck('COMPLETE generic source identity mod kappa',zero(num)==0)
quot=P({tuple(int(n)-(1 if h==9 else 0) for h,n in enumerate(e)):c for e,c in num.dict().items() if e[9]==1})
linear=vector(K,[quot.monomial_coefficient(x) for x in P.gens()]);quadratic=P({tuple(map(int,e)):c for e,c in quot.dict().items() if sum(e)==2})
ck('generic h divided source linear ZERO',not any(linear))
expected2=[CB.one(),-10740*crr/10049,-12960*crr/10049,CB(2960)/10049,476*crr/10049,-11692*crr/10049,CB(3489)/10049]
ck('H2 complete kernel coefficient witness',list(map(cyc2,vec))==expected2)
ph1=K.ring().hom([ckk,crr,1,-1],CB)
def cyc1(c):c=K(c);return ph1(c.numerator())/ph1(c.denominator())
expected1=[CB.one(),-2864*crr/2415,-3832*crr/2415,CB(968)/2415,344*crr/7245,-1856*crr/1449,CB(3592)/7245]
ck('H1 complete kernel coefficient witness',list(map(cyc1,vec))==expected1)
save({'ring':P,'source_first_two':firsttwo,'raw_sources':raws,'support':support,'seven_kappa0_NF_CF':pairs,'generic_combination':vec,'incidence_cofactor':co,'true_h_kappa_zero_divided_source':quot,'linear':linear,'quadratic':quadratic,'selected_minor':A,'selected_monomials':[mons[h] for h in piv],'checks':checks,'scope':__doc__},str(p/'pure_x_super_reader_generic_h_kappa_bridge.sobj'),compress=False)
report={'PASS':True,'checks':len(checks),'seconds':time.time()-st,'full_generic_kappa0_identity':True,'H1_H2_witnesses':True,'only_original_R5_H4_H9_H11_H13_H16_H17':True,'higher_kappa_source_layers_retained_as_full_original_DAG':True}
(p/'pure_x_super_reader_generic_h_kappa_bridge.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS',report,flush=True)
