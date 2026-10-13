"""Independent support, material and whole-base scheme-input audit.
Original D masks, all three cyclic labels and every root are retained.
"""
from sage.all import *
from pathlib import Path
from itertools import product
import time,json,resource,os
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(label,z):assert bool(z),label;checks.append(label)
qf=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'));raw=load(str(p/'independent_inactive_RRtorus_wall_variation_raw.sobj'));pairs=load(str(p/'independent_ratio_p_boundary_global_pair_syzygies.sobj'));can=load(str(p/'independent_supercritical_first_boundary_probe.sobj'));negpair=load(str(p/'independent_inactive_RRtorus_pair_negative_generic_kernel.sobj'));dec=load(str(p/'independent_inactive_RRtorus_canonical_delta_decoupling.sobj'))
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+650*1024**2,)*2)
D=qf['D'];patterns=[]
for row in D.rows():ck('original D disjoint source row '+str(len(checks)),sum(bool(f) for f in row)<=1)
for mask in product((0,1),repeat=6):
 if not any(mask):continue
 if not(mask[0]==mask[5] and mask[1]==mask[3] and mask[2]==mask[4]):continue
 active=[]
 for j in range(3):
  nz=any(any(all(not ex[k] or mask[k] for k in range(6)) for ex in D[i,j].dict()) for i in range(6));active.append(nz)
 rank=sum(active);amps=sum((mask[0],mask[1],mask[2]));ck('all original RRtorus paired D ranks '+str(mask),rank==(2 if amps==1 else 3));patterns.append((mask,rank))
ck('exactly three inactive full RR supports',sum(rank<=2 for mask,rank in patterns)==3)
F=PolynomialRing(QQ,names=('cm_a','cm_b','cm_h','cm_z')).fraction_field();a,b,h,z=F.gens();cycles=[];per=(3,2,5,4,0,1);pa,pb,ph,pz=a,b,h,z
for j in range(3):
 r=(pa-pb)/(pa*(pb-1));q=1+r;m=1+pa*r;N=pa*r*r;lam=m*q
 ck('positive-label q identity '+str(j),q==pb*(pa-1)/(pa*(pb-1)));ck('positive-label m identity '+str(j),m==(pa-1)/(pb-1));ck('lambda positive square identity '+str(j),lam==(pb/pa)*m*m);ck('material difference identity '+str(j),N-lam==-(pa*r+r+1))
 cycles.append((pa,pb,r,q,m,N,lam,ph,pz));pa,pb,ph,pz=pb/pa,1/pa,pz/ph,1/ph
ck('all original material/base cycles order three',(pa,pb,ph,pz)==(a,b,h,z));ck('projective permutation order three',all(per[per[per[i]]]==i for i in range(6)))
incoming=({0,5},{2,4},{1,3})
for j,ss in enumerate(incoming):ck('paired orientation transport '+str(j),{i for i in range(6) if per[i] in ss}==incoming[(j+1)%3])
R=PolynomialRing(QQ,'u');u=R.gen();ck('resonance real bound identity',QQ(1)/4+u*(1+u)==(u+QQ(1)/2)**2)
# Full pair layer inspection proves delta-independence at kappa=0 without
# using a derivative routine or assuming the added pair list is invariant.
for it in pairs['items']:
 counts={}
 for ex in it['full'].dict():
  pd=sum(ex[:6]);rd=sum(ex[6:]);kp=9-2*pd-3*rd;dp=4-pd-rd;counts[(pd,rd)]=counts.get((pd,rd),0)+1
  ck('complete pair exponent '+str((it['pair'],len(checks))),dp>=0 and kp>=-1)
  if kp<=0:ck('whole canonical pair delta absent '+str((it['pair'],len(checks))),dp==0)
 ck('sole regularization pole original (2,2) '+str(it['pair']),all((pd,rd)==(2,2) for pd,rd in counts if 9-2*pd-3*rd<0))
# The 24 full canonical RR initials are independently evaluated in a four-root
# algebra, rather than the creator's rational conjugation implementation.
K=PolynomialRing(QQ,names=('dc_a','dc_r','dc_H','dc_Z','dc_s','dc_w')).fraction_field();aa,rr,HH,ZZ,ss,ww=K.gens();qq=1+rr;alpha=(ZZ-HH)/(ZZ*(HH+1));beta=HH*(ZZ-1)/(ZZ*(HH+1));U=PolynomialRing(K,names=('dc_L','dc_B'));ul,ub=U.gens();E=U.quotient(U.ideal(ul**2+qq/rr,ub**2+qq),names=('L','B'));ll,bb=E.gens();pv=(E(1),bb*ss*ll/(-qq/rr),ww*bb/(-qq),E(ss),E(ww),ll);rv=(alpha*ll,beta*bb,E(1));ck('canonical incidence L equations',pv[0]*rv[0]-alpha*pv[5]*rv[2]==0 and pv[1]*rv[0]-alpha/beta*pv[3]*rv[1]==0 and pv[2]*rv[1]-beta*pv[4]*rv[2]==0)
ok=can['ring'].base_ring();cph=ok.ring().hom([aa,rr,HH,ZZ],K)
def co(z):return cph(z.numerator())/cph(z.denominator())
for name in ('F12','H13'):
 for j,f in enumerate(can[name]):
  v=sum((co(coef)*prod(pv[i]**int(ex[i]) for i in range(6))*prod(rv[i]**int(ex[i+6]) for i in range(3)) for ex,coef in f.dict().items()),E(0));ck('original canonical whole-amplitude '+name+str(j),v==0)
ck('creator full canonical source count matches',len(dec['items'])==24)
print('SUPPORT+MATERIAL+CANONICAL PASS',len(checks),'seconds',time.monotonic()-st,flush=True)
# Positive/negative first kernels: independently solve the true three-e range
# in the one-root algebra; no cached conjugation reduction is reused.
KK=PolynomialRing(QQ,names=('ak_a','ak_r','ak_H','ak_Z','ak_Rho')).fraction_field();aa,rr,HH,ZZ,Rho=KK.gens();qq=1+rr;mm=1+aa*rr;be=HH*(ZZ-1)/(ZZ*(HH+1));srcK=raw['ring'].base_ring();mk=srcK.ring().hom([aa,rr,HH,ZZ],KK);ids=raw['source_rows'];kernel_results=[]
for eps,name in ((1,'positive'),(-1,'negative')):
 lc=eps*aa*rr/mm;Q=PolynomialRing(KK,'Al');al=Q.gen();EE=Q.quotient(Q.ideal(al**2-lc),names=('L',));L=EE.gen()
 def iv(f):
  f=EE(f);h=f.lift();c0=h[0];c1=h[1];den=c0*c0-lc*c1*c1;ck('generic root inversion nonzero '+str(len(checks)),den!=0);return (c0-c1*L)/den
 def pp(g):return sum((mk(co.numerator())/mk(co.denominator())*L**int(ex[0])*Rho**int(ex[1]) for ex,co in g.dict().items()),EE(0))
 def ff(z):z=raw['ring'].fraction_field()(z);return pp(z.numerator())*iv(pp(z.denominator()))
 N=matrix(EE,[[ff(z) for z in row] for row in raw['normal'].rows()]);ir=[ids.index(j) for j in (0,5,7)];cols=[2,3,4];M=N.matrix_from_rows_and_columns(ir,cols);A=M.augment(identity_matrix(EE,3));rows=[list(row) for row in A.rows()]
 for j in range(3):
  ix=next(i for i in range(j,3) if rows[i][j]);rows[j],rows[ix]=rows[ix],rows[j];ivv=iv(rows[j][j]);rows[j]=[x*ivv for x in rows[j]]
  for i in range(3):
   if i!=j:
    zz=rows[i][j];rows[i]=[rows[i][k]-zz*rows[j][k] for k in range(6)]
 Ai=matrix(EE,[row[3:] for row in rows]);ck('independent genuine e range inverse '+name,M*Ai==identity_matrix(EE,3));wv=-eps*mm*L*Rho/(aa*rr*be);ev=-Ai*N.matrix_from_rows_and_columns(ir,[0,1,5])*vector(EE,[1,wv,0]);kv=vector(EE,[1,wv,*ev,0]);ck('all12 original wall normals whole-base '+name,N*kv==0)
 oldk=load(str(p/('independent_inactive_RRtorus_'+name+'_generic_kernel.sobj')));ck('creator wall kernel fully matches '+name,kv==vector(EE,[ff(z) for z in oldk['kernel']]))
 if eps==-1:
  for it in negpair['items']:
   ck('all3 added pair whole-base first normals '+str(it['pair']),vector(EE,[ff(z) for z in it['normal']])*kv==0)
 kernel_results.append((name,kv));print('GENERIC '+name+' KERNEL PASS',time.monotonic()-st,flush=True)
# Exact monomial intersection used in the relative Nakayama argument.
T=PolynomialRing(QQ,names=('s','d','Wl','WR','k'));s,d0,wl,wr,k=T.gens();J=T.ideal(s*wl,s*k,s*s,d0*wl,d0*wr,d0*k);I=T.ideal(s,d0).intersection(T.ideal(wl,wr,k,s*s)).intersection(T.ideal(d0,wl,k,s*s));ck('relative three-scheme intersection exact',I==J)
rec={'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st,'all_supports_and_three_material_labels':True,'whole_amplitude_canonical_delta_inputs':True,'whole_base_wall_kernels':True,'scope':__doc__}
save({'record':rec,'checks':checks,'patterns':patterns,'cycles':cycles,'kernel_results':kernel_results,'intersection':J},str(p/'independent_inactive_fullRR_subcritical_independent_coverage_reader.sobj'),compress=False);(p/'independent_inactive_fullRR_subcritical_independent_coverage_reader.json').write_text(json.dumps(rec,indent=2)+'\n');print('PASS',rec,flush=True)
