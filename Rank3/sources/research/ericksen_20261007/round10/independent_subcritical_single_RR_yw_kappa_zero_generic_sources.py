"""Whole-base genuine kappa-zero sources for the plus yw saturation bridge.

Only zero-kappa layers are built; H41, H13 and certified L-cofactors
give the entire genuine H at kappa=0, with t left arbitrary.
"""
from sage.all import *
from pathlib import Path
import resource,os,time,json
p=Path(__file__).resolve().parent;st=time.monotonic()
qf=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'))
hs=load(str(p/'independent_delta_subcritical_H_initials.sobj'))
hh=load(str(p/'independent_supercritical_RR_small_support_sources.sobj'))
ev=load(str(p/'independent_delta_subcritical_even_full_source.sobj'))
gau=load(str(p/'independent_partial_D_Newton_gauges_reader.sobj'))
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+550*1024**2,)*2)
print('PID',os.getpid(),'VM_CAP',resource.getrlimit(resource.RLIMIT_AS),flush=True)
OK=qf['ring'].base_ring();old=hh['original_field']
K=PolynomialRing(QQ,names=('kg_k','kg_r','kg_H','kg_Z')).fraction_field();k,r,H,Z=K.gens();m=1+k;q=1+r
ph=OK.ring().hom([k/r,r,H,Z],K);oh=old.ring().hom([k/r,k*q/(r*m),H/r,q*Z/r],K);cache={}
P=PolynomialRing(K,names=('kg_x','kg_z','kg_w','kg_e0','kg_e1','kg_e2','kg_delta','kg_t'))
x,z,w,e0,e1,e2,de,t=P.gens();alpha=(Z-H)/(Z*(H+1));beta=H*(Z-1)/(Z*(H+1));es=(e0,e1,e2)
U=PolynomialRing(K,names=('kg2_x','kg2_z','kg2_w','kg2_e0','kg2_e1','kg2_e2','kg2_delta','kg2_t','kg2_kap'))
ux,uz,uw,ue0,ue1,ue2,ud,ut,uk=U.gens();ims=(ux,U.one(),uz,uw,(ut*uz-uk*ue2)/beta,ux*uw*ut/beta+uk*(ux*ue1-ue0)/alpha,alpha/beta*uw*ut+uk*ue1,ut,U.one());zero=U.hom(list(P.gens())+[P.zero()],P);images=tuple(zero(f) for f in ims)
def cv(c,field=OK):
 c=field(c);key=(str(field),c)
 if key not in cache:
  h=ph if field==OK else oh;cache[key]=h(c.numerator())/h(c.denominator())
 return cache[key]
def norm(f,off,discard=()):
 out=P.zero()
 for ee,c in f.dict().items():
  ee=tuple(map(int,ee));pd=sum(ee[:6]);rd=sum(ee[6:]);ex=off[0]-2*pd-3*rd;dd=off[1]-pd-rd
  if (pd,rd) in discard or ex>0:continue
  assert ex==0 and dd>=0,(pd,rd,ex,dd)
  out+=cv(c)*de**dd*prod(a**n for a,n in zip(images,ee))
 return out
def mp(f):return sum((cv(c)*prod(a**int(n) for a,n in zip(images,ee)) for ee,c in f.dict().items()),P.zero())
def om(dd):return sum((cv(c,old)*prod(a**int(n) for a,n in zip(images,ee)) for ee,c in dd.items()),P.zero())
R=qf['ring'];DC=matrix(R,qf['D'])*matrix(R,qf['C']);sources={}
for j in (0,4,5):sources[j]=norm(qf['full_R_sources'][j],(8,4),((3,1),))+sum((mp(DC[j,i])*es[i] for i in range(3)),P.zero())
sources[6]=norm(qf['full_Q0'],(7,3),((1,2),))+sum((mp(qf['Q0_L_cofactors'][i])*es[i] for i in range(3)),P.zero())
for j in (4,16,17):
 it=hs['items'][j];sources[7+j]=om(it['H41_original'])+de*om(hh['H13_original'][j])+sum((om(it['L_cofactors_original'][i])*es[i] for i in range(3)),P.zero())
 print('BUILT H',j,len(sources[7+j].dict()),'sec',time.monotonic()-st,flush=True)
ell=matrix(R,qf['L'])*vector(R,R.gens()[6:])
for j,it in enumerate(ev['items']):
 rem=R(it['full_source'])-OK(it['square_unit'])*ell[int(it['square_component'])]**2-sum((R(it['A_p_quadratic'][i])*ell[i] for i in range(3)),R.zero())
 sources[25+j]=norm(rem,(8,4))+cv(it['square_unit'])*es[int(it['square_component'])]**2+sum((mp(it['A_p_quadratic'][i])*es[i] for i in range(3)),P.zero())
 print('BUILT EVEN',j,len(sources[25+j].dict()),'sec',time.monotonic()-st,flush=True)
save({'ring':P,'sources':sources,'range_rows':(4,5,6,23),'scope':__doc__},str(p/'independent_subcritical_single_RR_yw_kappa_zero_generic_sources.sobj'),compress=False)
print('SAVED BASE',time.monotonic()-st,flush=True)
# Assemble the original complete p13-active numerator once, including
# both original h/z coefficient derivatives before specialization.
O=ev['ring'];OQ=next(it['full_source'] for it in ev['items'] if int(it['Q_index'])==6);it=gau['items'][1];of=gau['ring'].base_ring();GO=PolynomialRing(K,names=gau['ring'].variable_names());vs=GO.gens()
def mapped(f):return GO({tuple(map(int,ee)):cv(c,f.parent().base_ring()) for ee,c in f.dict().items()})
Q=mapped(OQ);fv=list(map(mapped,it['field_numerators']));bv=list(map(mapped,it['base_numerators']));G=sum((fv[j]*Q.derivative(vs[j]) for j in range(9)),GO.zero())
for j,param in enumerate((OK.gen(2),OK.gen(3))):
 scale=OK.gen(1) if j==0 else OK.gen(1)/(1+OK.gen(1));dd=O({tuple(map(int,ee)):scale*c.derivative(param) for ee,c in OQ.dict().items() if c.derivative(param)});G+=bv[j]*mapped(dd)
print('ORIGINAL G TERMS',len(G.dict()),'sec',time.monotonic()-st,flush=True)
save({'ring':GO,'G':G,'scope':__doc__},str(p/'independent_subcritical_single_RR_G30_original_generic.sobj'),compress=False)
answer=P.zero();layers={}
for ee,c in G.dict().items():
 deg=sum(ee);rd=sum(ee[6:]);power=13-2*deg-rd;dp=6-deg
 if power>0:continue
 layers.setdefault((power,dp),{})[tuple(map(int,ee))]=c
for (power,dp),coeffs in sorted(layers.items()):
 assert dp>=0
 f=U.zero()
 for ee,c in coeffs.items():f+=c*prod(a**n for a,n in zip(ims,ee))
 assert all(int(ee[8])>=-power for ee in f.dict()),(power,dp)
 ff=P({tuple(int(n) for n in ee[:8]):c for ee,c in f.dict().items() if int(ee[8])==-power})
 answer+=de**dp*ff;print('G LAYER',power,dp,len(ff.dict()),'sec',time.monotonic()-st,flush=True)
sources[30]=answer;save({'ring':P,'sources':sources,'range_rows':(4,5,6,23),'scope':__doc__},str(p/'independent_subcritical_single_RR_yw_kappa_zero_generic_sources.sobj'),compress=False)
print('SAVED ALL',len(answer.dict()),time.monotonic()-st,flush=True)
