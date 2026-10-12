"""Replay the complete generic divided source, retaining every positive kappa layer."""
from sage.all import *
from pathlib import Path
import resource,os,time,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'pure_x_root_generic_w_bridge.sobj'));old=load(str(p/'subcritical_single_pure_x_linear_kappa_syzygies_H2.sobj'))
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+500*1024**2,)*2)
P=d['ring'];K=P.base_ring();k,r,H,Z=K.gens();y,z,w,e0,e1,e2,de,t,v,kap=P.gens();al=(Z-H)/(Z*(H+1));be=H*(Z-1)/(Z*(H+1));inc=y*v-al/be*t*w-kap*e1
images=(P.one(),y,z,w,(t*z-kap*e2)/be,(v-kap*e0)/al,v,t,P.one());checks=[]
def ck(n,vv):assert bool(vv),n;checks.append(n)
def nk(c):
 c=K(c)
 if not c:return K.zero()
 n,den=c.numerator(),c.denominator();gg=n.gcd(den);n,den=n//gg,den//gg;lc=den.leading_coefficient();return K(n/lc)/K(den/lc)
def nt(f):return P({tuple(map(int,e)):nk(c) for e,c in f.dict().items() if c})
full=[]
for it in d['items']:
 R=it['original_derivative'].parent();hom=R.hom(images,P);layers={};offset,do=it['offset'];known={(int(a),int(b)):co for a,b,c,co in it['division_certificates']}
 for ee,c in it['original_derivative'].dict().items():layers.setdefault((sum(ee),sum(ee[6:])),{})[tuple(map(int,ee))]=c
 out=P.zero()
 for (deg,rd),dd in layers.items():
  power=offset-2*deg-rd;dp=do-deg;ck('nonnegative delta '+str((it['name'],deg,rd)),dp>=0);ff=nt(hom(R(dd)))
  if power<0:
   cofs=known[(deg,rd)];ck('all negative divisions retained '+str((it['name'],deg,rd)),len(cofs)==-power)
   for j,co in enumerate(cofs):
    ff=nt(ff-co*inc);ck('literal complete incidence division '+str((it['name'],deg,rd,j)),all(e[9]>=1 for e in ff.dict()))
    ff=P({tuple(int(a)-(1 if ii==9 else 0) for ii,a in enumerate(e)):c for e,c in ff.dict().items()})
  else:ff*=kap**power
  out=nt(out+de**dp*ff)
 ck('complete first two layers '+it['name'],P({tuple(map(int,e)):c for e,c in out.dict().items() if e[9]<=1})==it['first_two_kappa']);full.append(out)
 save({'ring':P,'full_G_sources':full,'checks':checks},str(p/'pure_x_root_generic_w_reader.partial.sobj'),compress=False)
 print(it['name'],'complete terms',len(out.dict()),flush=True)
vec=d['combination'];acc={}
for scalar,shift,poly in [(vec[0],None,full[0]),(vec[1],2,full[1]),(vec[2],8,full[1]),(vec[3],0,full[2])]:
 for ee,c in poly.dict().items():
  ex=list(map(int,ee))
  if shift is not None:ex[shift]+=1
  ex=tuple(ex);acc[ex]=nk(acc.get(ex,K.zero())+nk(scalar*c))
for ee,c in d['incidence_cofactor'].dict().items():
 for ie,ic in inc.dict().items():
  ex=tuple(int(a)+int(b) for a,b in zip(ee,ie));acc[ex]=nk(acc.get(ex,K.zero())-nk(c*ic))
numerator=P({e:c for e,c in acc.items() if c})
ck('whole source exact kappa divisibility',all(e[9]>=1 for e in numerator.dict()))
source=P({tuple(int(a)-(1 if ii==9 else 0) for ii,a in enumerate(e)):c for e,c in numerator.dict().items()});zero=P.hom(list(P.gens()[:9])+[P.zero()],P)
ck('whole source division multiplyback',kap*source==numerator)
ck('whole kappa-zero source replay',zero(source)==d['kappa_zero_divided_source'])
ck('whole material-independent-base normal',vector(K,[source.monomial_coefficient(vv) for vv in P.gens()])==vector(K,[0,0,2*k/(1+k),0,0,0,0,0,0,0]))
HP=old['ring'];HK=HP.base_ring();ph=K.ring().hom([HK.gen(0),HK.gen(1),HK(2),HK(-1)],HK)
special=HP({tuple(map(int,e)):ph(c.numerator())/ph(c.denominator()) for e,c in source.dict().items()})
ck('independent full H2 syzygy source equality',special==old['items'][6]['source'])
save({'ring':P,'source':source,'full_G_sources':full,'combination':vec,'incidence_cofactor':d['incidence_cofactor'],'checks':checks,'scope':__doc__},str(p/'pure_x_root_generic_w_reader.sobj'),compress=False)
report={'status':'PASS','checks':len(checks),'source_terms':len(source.dict()),'seconds':time.monotonic()-st,'scope':__doc__};(p/'pure_x_root_generic_w_reader.json').write_text(json.dumps(report,indent=2)+'\n');print(report,flush=True)
