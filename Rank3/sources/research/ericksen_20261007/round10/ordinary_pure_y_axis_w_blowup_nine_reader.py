from sage.all import *
from pathlib import Path
import os,resource,time,json,itertools
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(name,v):assert bool(v),name;checks.append(name)
def mul(d,e,cap=2):
 out={}
 for a,c in d.items():
  for b,v in e.items():
   ee=tuple(x+y for x,y in zip(a,b))
   if sum(ee)<=cap:out[ee]=out.get(ee,0)+c*v
 return {ee:v for ee,v in out.items() if v}
def sparse_det(M,rows=None,cols=None):
 if rows is None:rows=tuple(range(M.nrows()));cols=tuple(range(M.ncols()))
 if not rows:return M.base_ring().one()
 i=min(rows,key=lambda i:sum(bool(M[i,j]) for j in cols));ii=rows.index(i);ans=M.base_ring().zero()
 for jj,j in enumerate(cols):
  if M[i,j]:ans+=(-1)**(ii+jj)*M[i,j]*sparse_det(M,tuple(h for h in rows if h!=i),tuple(h for h in cols if h!=j))
 return ans
saved=load(str(p/'ordinary_pure_y_axis_w_blowup_nine_matrix.sobj'));outputs=[]
for it in saved['items']:
 hv=int(it['H']);d=load(str(p/('independent_subcritical_single_RR_full_regular_chart_H%s.sobj'%hv)));pa=load(str(p/('independent_subcritical_single_RR_full_pair_chart_H%s.sobj'%hv)));P=d['ring'];K=P.base_ring();k,r=K.gens();q=1+r;m=1+k
 if hv==1:resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+450*1024**2,)*2)
 sources=d['sources']+[xx['source'] for xx in pa['items']];rows=(4,5,6,23,11,29,25,26,27);wp=(1,2,1,1,0,1,0,0,0);MM=matrix(K,9,9);zero=(0,)*8
 def lin(cs):
  dd={}
  for j,cc in cs.items():
   ee=list(zero);ee[j]=1;dd[tuple(ee)]=K(cc)
  return dd
 xx=lin({4:1});zz=lin({5:1});tt=lin({6:1});kk=lin({7:1})
 baseB=K.zero() if hv==1 else -K.one()/2
 Aimage=lin({0:1,4:10 if hv==1 else 7,5:12 if hv==1 else 9});Bimage=lin({1:1});Bimage[zero]=baseB
 Eimage=lin({2:1,4:12 if hv==1 else 16,5:10 if hv==1 else 14});Dimage=lin({3:1,4:(32 if hv==1 else 27)*r/q,5:(32 if hv==1 else 27)*r/q});im=(xx,zz,None,Aimage,Bimage,Eimage,Dimage,tt,kk)
 reconstructed=[]
 for ii,(row,poww) in enumerate(zip(rows,wp)):
  ans={}
  for ee,c0 in sources[row].dict().items():
   wd=sum(int(ee[j]) for j in (2,3,4,6));ck('all raw w powers H%s row%s'%(hv,row),wd>=poww)
   if wd!=poww:continue
   vv={zero:K(c0)}
   for j,n in enumerate(ee):
    if j==2:continue
    for unused in range(int(n)):vv=mul(vv,im[j])
   for ex,cc in vv.items():ans[ex]=ans.get(ex,K.zero())+cc
  ans={ee:cc for ee,cc in ans.items() if cc};ck('constant zero H%s row%s'%(hv,row),not ans.get(zero,0))
  for j in range(4):
   ex=list(zero);ex[j]=1;MM[ii,j]=ans.get(tuple(ex),K.zero())
  for j in range(4,8):
   ex=list(zero);ex[j]=1;ck('free linear zero H%s row%s col%s'%(hv,row,j),not ans.get(tuple(ex),0))
  def coef(i,j):
   ex=list(zero);ex[i]+=1;ex[j]+=1;return ans.get(tuple(ex),K.zero())
  vk=(coef(4,7),coef(5,7),coef(4,6),coef(5,6),coef(4,4),coef(4,5),coef(5,5))
  c1=vk[1];c4=(vk[0]-c1)/2;c2=vk[3];c3=vk[6];c5=(2*c3-vk[5])/3
  ck('five generator x*t H%s row%s'%(hv,row),vk[2]==c2+c4);ck('five generator x2 H%s row%s'%(hv,row),vk[4]==c3-2*c5)
  for j,cc in enumerate((c1,c2,c3,c4,c5),4):MM[ii,j]=cc
  reconstructed.append(ans)
 ck('all nine matrix entries H%s'%hv,MM==it['matrix']);det=sparse_det(MM);target=(-QQ(25)/4032 if hv==1 else -QQ(544)/10935)*q*k**3;ck('independent sparse Leibniz H%s'%hv,det==target)
 R3=MM.matrix_from_rows_and_columns((0,2,3),(0,1,3));det3=sparse_det(R3);print('H',hv,'DET9',det.factor(),'DET3',det3.factor(),flush=True);outputs.append({'H':hv,'matrix':MM,'range3det':det3,'determinant':det,'target':target,'raw_dict_jets':reconstructed})
# Generic pair source, complete polynomial, original strong J9 normal form.
pa=load(str(p/'independent_subcritical_single_RR_full_pair_chart.sobj'));OK=pa['ring'].base_ring();K=PolynomialRing(QQ,names=('rd_k','rd_r','rd_H','rd_Z','rd_W')).fraction_field();k,r,H,Z,W=K.gens();q=1+r;ph=OK.ring().hom([k,r,H,Z],K);T=PolynomialRing(K,names=('rd_x','rd_s','rd_t','rd_kap'));x,s,t,kap=T.gens();c=(3*k-W*W)/(k+W*W)
def nf(f):
 out={}
 for ee,a in T(f).dict().items():
  i,j,l,n=map(int,ee)
  if j>=2 or i>=3 or (j and l+n) or (i==2 and (j or l+n)):continue
  if i==2:key=(1,1,0,0);a*=c
  elif i==1 and j==0:key=(1,0,0,l+n);a*=(-2)**l
  else:key=(i,j,l,n)
  out[key]=out.get(key,K.zero())+a
 return T(out)
f0=-2*x*W+(4*H*Z-2*H-4*Z*Z+2*Z)/(H*Z)*s*W;f1=-2*Z*(H-1)/(H*(Z-1))*W;f2=2*x+(2*H*Z-2*H-4*Z*Z+2*Z)/Z*s;fd=4*r*(H+1)**3/(q*H*H)*s*W
ims=(x,s-x,T(W),f0,T(f1),f2,fd,t,kap);pairs=[]
for row in range(3):
 ans=T.zero()
 for ee,cc in pa['items'][row]['source'].dict().items():
  vv=T.one()
  for val,n in zip(ims,ee):
   for unused in range(int(n)):vv=nf(vv*val)
  if vv:ans+=ph(cc.numerator())/ph(cc.denominator())*vv
 ans=nf(ans);ck('complete generic pair in strong J9 '+str(row),not ans);pairs.append(ans)
print('GENERIC PAIRS ALL ZERO',time.monotonic()-st,flush=True)
save({'items':outputs,'generic_pair_normal_forms':pairs,'checks':checks,'scope':'Independent original dictionary jets + sparse Leibniz; complete generic pair identities; no source creator rerun'},str(p/'ordinary_pure_y_axis_w_blowup_nine_reader.sobj'),compress=False);(p/'ordinary_pure_y_axis_w_blowup_nine_reader.json').write_text(json.dumps({'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st},indent=2)+'\n');print('PASS',len(checks),'seconds',time.monotonic()-st,flush=True)
