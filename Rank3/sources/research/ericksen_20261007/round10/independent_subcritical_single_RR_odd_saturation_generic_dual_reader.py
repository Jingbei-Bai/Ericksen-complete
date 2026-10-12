"""Exact ambient R0/H17 saturation bridge and whole-base/all-w J10 dual replay."""
from sage.all import *
from pathlib import Path
import resource,os,time,json
p=Path(__file__).resolve().parent;st=time.monotonic();ra=load(str(p/'independent_subcritical_single_RR_yw_kappa_zero_generic_range.sobj'));d=load(str(p/'independent_subcritical_single_RR_yw_kappa_zero_generic_sources.sobj'));eh=load(str(p/'independent_subcritical_single_RR_yw_kappa_zero_generic_exact.sobj'));normal=load(str(p/'independent_critical_single_RR3_normal.sobj'));resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+550*1024**2,)*2);print('PID',os.getpid(),flush=True);U=ra['ring'];OK=U.base_ring();ok,orr,oH,oZ=OK.gens();x,z,w,t=U.gens();P=d['ring'];px,pz,pw,pe0,pe1,pe2,pde,pt=P.gens();D=eh['D'];N=eh['N'];Sigma=eh['Sigma'];BH=eh['H4_tSigma_quotient'];A=ra['matrix'];det=ra['det'];checks=[]
def ck(n,v):assert bool(v),n;checks.append(n)
def nk(c,field=OK):
 c=field(c)
 if not c:return field.zero()
 nn,dd=c.numerator(),c.denominator();gg=nn.gcd(dd);nn,dd=nn//gg,dd//gg;lc=dd.leading_coefficient();return field(nn/lc)/field(dd/lc)
def nt(f):return U({tuple(map(int,ee)):nk(c) for ee,c in f.dict().items() if c})
def image(f):
 md=max((sum(int(ee[j]) for j in (3,4,5,6)) for ee in f.dict()),default=0);ck('complete zero-kappa source range-affine',md==1);out=U.zero()
 for ee,c in f.dict().items():
  val=U(c)*x**int(ee[0])*z**int(ee[1])*w**int(ee[2])*t**int(ee[7])*D**(1-sum(int(ee[j]) for j in (3,4,5,6)))
  for j,im in zip((3,4,5,6),N):val*=im**int(ee[j])
  out+=val
 return nt(out)
hm=P.hom([x,z,w,U.zero(),U.zero(),U.zero(),U.zero(),t],U);he=vector(U,[hm(d['sources'][11].derivative(v)) for v in (pe0,pe1,pe2,pde)]);items=[]
for ri in (0,24):
 fn=image(d['sources'][ri]);bf=U.fraction_field()(fn)/(t*Sigma);ck('exact t Sigma factor row '+str(ri),bf.denominator().is_constant());BF=nt(U(bf));ck('exact residual multiplyback row '+str(ri),nt(fn-t*Sigma*BF)==0);fe=vector(U,[hm(d['sources'][ri].derivative(v)) for v in (pe0,pe1,pe2,pde)]);items.append({'row':ri,'BF':BF,'source_range_gradient':fe,'cofactor_range_factored_column':BH*fe-BF*he});print('BRIDGE',ri,'BF TERMS',len(BF.dict()),'sec',time.monotonic()-st,flush=True)
K=PolynomialRing(QQ,names=('od_k','od_r','od_W','od_H','od_Z')).fraction_field();k,r,W,H,Z=K.gens();m=1+k;q=1+r;ph=OK.ring().hom([k,r,H,Z],K);T=PolynomialRing(K,names=('od_x','od_s'));xx,ss=T.gens()
def cf0(f):
 return nk(sum((ph(c.numerator())/ph(c.denominator())*W**int(ee[2]) for ee,c in f.dict().items() if int(ee[0])+int(ee[1])+int(ee[3])==0),K.zero()),K)
A0=matrix(K,4,4,[cf0(ai) for ai in A.list()]);BH0=cf0(BH);det0=nk(A0.det(),K);common0=nk(det0*BH0,K);ck('common source denominator is locally nonzero',common0!=0);print('COMMON0 FACT',common0.factor(),flush=True)
NP=normal['ring'];NK=NP.base_ring();nph=NK.ring().hom([k/r,r,H,Z],K);ims=(xx,ss-xx,T(W),-2*xx*W+(4*H*Z-2*H-4*Z**2+2*Z)/(H*Z)*ss*W,T(-2*Z*(H-1)/(H*(Z-1))*W),2*xx+(2*H*Z-2*H-4*Z**2+2*Z)/Z*ss,4*r*(H+1)**3/(q*H**2)*ss*W);cache={}
def pull(f):
 out=T.zero()
 for ee,c in f.dict().items():
  if c not in cache:cache[c]=nk(nph(c.numerator())/nph(c.denominator()),K)
  val=T(cache[c])
  for im,nn in zip(ims,ee):val*=im**int(nn)
  out+=val
 return T({tuple(map(int,ee)):nk(c,K) for ee,c in out.dict().items() if sum(ee)<=1})
rows=ra['range_rows'];nn={ri:pull(normal['normal_matrix'][ri,8]) for ri in (0,24,11)+rows};tests=[]
for it in items:
 ri=it['row'];aN0=A0.adjugate().transpose()*vector(K,[cf0(fa) for fa in it['cofactor_range_factored_column']]);ans=common0*nn[ri]-det0*cf0(it['BF'])*nn[11]-sum((aN0[j]*nn[rj] for j,rj in enumerate(rows)),T.zero());ans=T({tuple(map(int,ee)):nk(c,K) for ee,c in ans.dict().items() if c});ck('whole all-w kappa first saturated source has no linear x/s '+str(ri),ans==0);tests.append(ans);it['cofactor_range_origin_numerators']=aN0;print('DUAL',ri,ans,'sec',time.monotonic()-st,flush=True)
save({'ring':U,'items':items,'range_matrix':A,'determinant':det,'H4_gradient':he,'BH':BH,'dual_ring':T,'kappa_first_dual_remainders':tests,'common_unit_origin':common0,'checks':checks,'scope':__doc__},str(p/'independent_subcritical_single_RR_odd_saturation_generic_dual_reader.sobj'),compress=False);(p/'independent_subcritical_single_RR_odd_saturation_generic_dual_reader.json').write_text(json.dumps({'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st},indent=2)+'\n');print('PASS',len(checks),time.monotonic()-st,flush=True)
