"""Independent original-pair/binomial-source replay and protected ten-unit checks, generic and cyclic."""
from sage.all import *
from pathlib import Path
import resource,os,time,json,sys
p=Path(__file__).resolve().parent;st=time.monotonic();mode=sys.argv[1];assert mode in ('H2','H1');full=load(str(p/('independent_subcritical_single_RR_full_regular_chart_'+mode+'.sobj')));pa=load(str(p/('independent_subcritical_single_RR_full_pair_chart_'+mode+'.sobj')));raw=load(str(p/'independent_ratio_p_boundary_global_pair_syzygies.sobj'));cyc=mode=='H1';jj=load(str(p/('independent_subcritical_single_RR_yw_minus_cyclic_source_jets.sobj' if cyc else 'independent_subcritical_single_RR_yw_minus_odd_saturated_cubic_probe.sobj')));mm=load(str(p/('independent_subcritical_single_RR_yw_minus_cyclic_ten_module.sobj' if cyc else 'independent_subcritical_single_RR_yw_minus_ten_module_reader.sobj')));rel=load(str(p/('independent_subcritical_single_RR_yw_minus_cyclic_relative_ten_ideal.sobj' if cyc else 'independent_subcritical_single_RR_yw_minus_generic_ten_ideal.sobj')));resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+450*1024**2,)*2);print('PID',os.getpid(),'MODE',mode,flush=True);checks=[]
def ck(s,b):assert bool(b),s;checks.append(s)
T=jj['ring'];K=T.base_ring();r,W=K.gens();k=-W**2;m=1+k;q=1+r;H,Z=(K(1),K(-1)) if cyc else (K(2),K(-1));x,z,t,kap,u=T.gens();s=x+z;N=PolynomialRing(K,names=tuple('mr_v'+str(j) for j in range(9)));nv=N.gens();OK=full['ring'].base_ring();ph=OK.ring().hom([k,r],K);cc={}
def nk(c):
 c=K(c)
 if not c:return K.zero()
 nn,dd=c.numerator(),c.denominator();gg=nn.gcd(dd);nn,dd=nn//gg,dd//gg;lc=dd.leading_coefficient();return K(nn/lc)/K(dd/lc)
def cv(c):
 if c not in cc:cc[c]=nk(ph(c.numerator())/ph(c.denominator()))
 return cc[c]
sources=[N({tuple(map(int,ee)):cv(c) for ee,c in f.dict().items()}) for f in full['sources']]
# Reconstruct every normalized pair directly from its original full polynomial.
RK=raw['ring'].base_ring();rp=RK.ring().hom([k/r,k*q/(r*m),H/r,Z*q/r],K);al=(Z-H)/(Z*(H+1));be=H*(Z-1)/(Z*(H+1));nx,nz,nw,ne0,ne1,ne2,nd,nt,nkap=nv;imgs=(nx,N.one(),nz,nw,(nt*nz-nkap*ne2)/be,nx*nw*nt/be+nkap*(nx*ne1-ne0)/al,al/be*nw*nt+nkap*ne1,nt,N.one())
for it,itpa in zip(raw['items'],pa['items']):
 ff=N.zero();negative=N.zero()
 for ee,c in it['full'].dict().items():
  pd,rd=sum(ee[:6]),sum(ee[6:]);ke,dp=9-2*pd-3*rd,4-pd-rd;ck('pair layer regular or sole -1 '+str(len(checks)),ke>=-1 and dp>=0);a=nk(rp(c.numerator())/rp(c.denominator()));val=a*nd**dp*prod(im**int(nn) for im,nn in zip(imgs,ee))
  if ke<0:negative+=val
  else:ff+=nkap**ke*val
 ck('original negative pair layer exactly kappa divisible '+str(it['pair']),all(int(ee[8])>=1 for ee in negative.dict()));ff+=N({tuple(int(nn)-(1 if j==8 else 0) for j,nn in enumerate(ee)):nk(c) for ee,c in negative.dict().items()});ff=N({tuple(map(int,ee)):nk(c) for ee,c in ff.dict().items()});old=N({tuple(map(int,ee)):cv(c) for ee,c in itpa['source'].dict().items()});ck('original full pair agrees archive '+str(it['pair']),ff==old);sources.append(ff)
print('ORIGINAL PAIRS REPLAY',len(checks),'sec',time.monotonic()-st,flush=True)
base_e1=K.zero() if cyc else -W/2
def taylor(ff):
 out={}
 for ee,c in ff.dict().items():
  ext=sum(int(nn) for j,nn in enumerate(ee) if j not in (2,4))
  if ext>4:continue
  for iw in range(min(int(ee[2]),4-ext)+1):
   for ie in range(min(int(ee[4]),4-ext-iw)+1):
    if not base_e1 and ie!=int(ee[4]):continue
    ex=list(map(int,ee));ex[2]=iw;ex[4]=ie;val=c*binomial(int(ee[2]),iw)*W**(int(ee[2])-iw)*binomial(int(ee[4]),ie)*base_e1**(int(ee[4])-ie);ex=tuple(ex);out[ex]=out.get(ex,K.zero())+val
 return N({ee:nk(c) for ee,c in out.items() if c})
jets=[taylor(ff) for ff in sources]
def trunc(ff):return T({tuple(map(int,ee)):nk(c) for ee,c in ff.dict().items() if sum(ee)<=4 and c})
def compose(ff,im):
 pp={}
 for j,v in enumerate(im):
  pp[j,0]=T.one()
  for nn in range(1,5):pp[j,nn]=trunc(pp[j,nn-1]*v)
 out=T.zero()
 for ee,c in ff.dict().items():
  val=T(c)
  for j,nn in enumerate(ee):
   if nn:val=trunc(val*pp[j,int(nn)])
   if not val:break
  out+=val
 return trunc(out)
im=list(jj['range_images4']);im[2]-=W;im[4]-=base_e1;replay=[]
for j,jf in enumerate(jets):
 rr=compose(jf,im);ck('independent binomial true complete jet '+str(j),rr==jj['residuals4'][j]);replay.append(rr)
rows=(4,5,6,23);ck('four true range jets zero',all(replay[j]==0 for j in rows));A=matrix(K,4,4,lambda i,j:jets[rows[i]].monomial_coefficient(nv[(3,4,5,6)[j]]));ck('true range determinant nonzero',A.det()!=0)
C0=(4*H*Z-2*H-4*Z**2+2*Z)/(H*Z);C2=(2*H*Z-2*H-4*Z**2+2*Z)/Z;decoef=4*r*(H+1)**3/(q*H**2);oldim=(x,z,u,-2*x*(W+u)+C0*s*(W+u),base_e1/W*u,2*x+C2*s,decoef*s*(W+u),T.zero(),T.zero());C=compose(jets[27],oldim);gg=(kap*s,t*s,s*s,(t+2*kap)*x,C,kap*x**3);ck('target C is actual old full graph jet',C==mm['generators'][4]);ck('all six target generator jets match',gg==tuple(mm['generators']))
mons=mm['monomials'];vec=lambda ff:vector(K,[ff.monomial_coefficient(T.monomial(*ee)) for ee in mons]);E=mm['mJ_jet_echelon'];EE=[row for row in E if any(row)];piv=[next(i for i,a in enumerate(row) if a) for row in EE]
def red(v):
 v=vector(K,v)
 for row,pi in zip(EE,piv):
  if v[pi]:v-=v[pi]*row
 return v
def tuples(total,nn):
 if nn==1:yield(total,);return
 for i in range(total+1):
  for rest in tuples(total-i,nn-1):yield(i,)+rest
mg=[]
for gf in gg:
 order=min(sum(ee) for ee in gf.dict())
 for deg in range(1,5-order):
  for ee in tuples(deg,5):mg.append(vec(trunc(gf*T.monomial(*ee))))
MG=matrix(K,mg);ck('independent mJ finite rowspace equality',MG.row_space()==matrix(K,EE).row_space());B=matrix(K,[red(vec(gf)) for gf in gg]);ck('J/mJ coefficient class injective into fourjets',B.rank()==6);ck('actual generator basis matches',B==mm['basis']);sel=B.pivots();BB=B.matrix_from_columns(sel);M=[]
for j,ff in enumerate(replay):
 v=red(vec(ff));co=BB.transpose().solve_right(vector(K,[v[ii] for ii in sel]));ck('full relative source class matches coefficient row '+str(j),co*B==v);M.append(co)
M=matrix(K,M);ck('independent actual coefficient matrix matches',M==mm['source_coefficient_matrix']);ck('true six source normals span all target classes',M.rank()==6);ids=M.transpose().pivots();D6=nk(M.matrix_from_rows(ids).det());D4=nk(A.det());D10=nk(D4*D6);expect=QQ(-256)/3*r*W**8/(W**2-1) if cyc else QQ(384)/25*r*W**7/(W**2-1);ck('exact six protected determinant',D6==expect);guard=K.ring()(W*r*(1-W**2)*(1+r))
for label,part in (('numerator',D10.numerator()),('denominator',D10.denominator())):
 for fac,pow in part.factor():ck('all material factors protected '+label+' '+str(fac),guard.quo_rem(fac)[1]==0)
flip=K.ring().hom([r,-W],K);DU=nk(D10*flip(D10.numerator())/flip(D10.denominator()));L=PolynomialRing(QQ,names=('mn_r','mn_k')).fraction_field();lr,lk=L.gens()
def evenmap(pp):
 out=L.zero()
 for ee,c in pp.dict().items():ck('both-root norm even in W '+str(len(checks)),int(ee[1])%2==0);out+=c*lr**int(ee[0])*(-lk)**(int(ee[1])//2)
 return out
norm=evenmap(DU.numerator())/evenmap(DU.denominator());print('DET4',D4.factor(),'DET6',D6.factor(),'DET10',D10.factor(),'BOTH ROOT NORM',norm.factor(),flush=True)
# Whole relative identity checks precede the finite coefficient argument.
RT=rel['ring'];RK=RT.base_ring();rv=RK.gens();rk,rr,rw=(rv[0],rv[1],rv[-1]);rm=1+rk;rx,rs,rt,rkap=RT.gens();RC=rel['ideal_generators'][4];expectedC=-4*rk/rm**2*rx**4+8*rk/rm**2*rx**3*rs+8*(rk+rw**2)/rm*rx**2-8*(3*rk-rw**2)/rm*rx*rs;ck('exact whole soft conic C',RC==expectedC);ck('C-kappa torsion identity without D division',(rk+rw**2)*rkap*rx**2-(rm/8*rkap*RC+rk/(2*rm)*rkap*rx**4-rk/rm*rkap*rx**3*rs+(3*rk-rw**2)*rkap*rx*rs)==0);ck('all31 whole source normal forms',all(ff==0 for i,ff in enumerate(rel['source_normal_forms']) if i!=27) and rel['source_normal_forms'][27]==RC);ck('exact inactive p31 p32 zero in target quotient',rel['chart_p'][4]==0 and rel['chart_p'][5]==0);ck('delta graph divisible by s',all(int(ee[1])>=1 for ee in rel['graph'][6].dict()))
save({'mode':mode,'ring':T,'independent_true_source_jets':replay,'matrix':M,'rows':ids,'range_matrix':A,'det6':D6,'det4':D4,'det10':D10,'norm_ring':L,'both_root_norm':norm,'checks':checks,'scope':__doc__},str(p/('independent_subcritical_single_RR_yw_minus_closed_reader_'+mode+'.sobj')),compress=False);(p/('independent_subcritical_single_RR_yw_minus_closed_reader_'+mode+'.json')).write_text(json.dumps({'status':'PASS','mode':mode,'checks':len(checks),'seconds':time.monotonic()-st},indent=2)+'\n');print('PASS',mode,len(checks),time.monotonic()-st,flush=True)
