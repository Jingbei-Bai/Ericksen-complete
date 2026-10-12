"""Full ordinary normal four-jets on the entire cyclic p12-axis boundary.

Chart p12=1. Tangent RR2=u, RR3=v are unrestricted. Normal coordinates
are (p13,p21,p23,p31,p32,RR1,tau). No tangent variable is inverted.
Sources: four original Q, six full R, eighteen genuine H, nine genuine E.
H/E degree six cancellation is global before homogenization to degree five.
"""
from pathlib import Path
import time,json,os
from sage.all import load,save,PolynomialRing,QQ,matrix,vector
out=Path(__file__).resolve().parent;start=time.monotonic()
g=load(str(out/'independent_ratio_p_boundary_second_gauge.sobj'))
hdata=load(str(out/'independent_ratio_global_R_derivative_degreefall.sobj'))
edata=load(str(out/'independent_ratio_even_derivative_cyclic_direct.sobj'))
oldq=load(str(out.parent/'round8/rank3_joint_fiber_leading.sobj'))
K0=g['coefficient_field'];oa,ob,oh,oz=K0.gens();O=PolynomialRing(K0,names=tuple('ax_p%d'%i for i in range(6))+tuple('ax_R%d'%i for i in range(3)));ov=O.gens()
K=PolynomialRing(QQ,names=('ax_a','ax_r')).fraction_field();a,r=K.gens();P=K.ring();aa,rr=P.gens();qq=1+rr;mm=1+aa*rr;phcache={};pows={};cfcache={};dcache={}
def pw(i,n):
 key=(int(i),int(n))
 if key not in pows:pows[key]=(aa,rr,qq,mm)[i]**int(n)
 return pows[key]
def ph(f):
 f=K0.ring()(f)
 if f not in phcache:
  dd=f.dict();db=max((int(e[1]) for e in dd),default=0);dr=max((int(e[2])+int(e[3]) for e in dd),default=0)
  n=sum((P(c)*(-1)**int(e[3])*pw(0,int(e[0])+int(e[1]))*pw(2,int(e[1])+int(e[3]))*pw(1,dr-int(e[2])-int(e[3]))*pw(3,db-int(e[1])) for e,c in dd.items()),P.zero())
  phcache[f]=K(n)/K(pw(1,dr)*pw(3,db))
 return phcache[f]
def cf(c):
 c=K0(c)
 if c not in cfcache:cfcache[c]=ph(c.numerator())/ph(c.denominator())
 return cfcache[c]
def dcf(c,p):
 c=K0(c);p=K0.ring()(p);key=(c,str(p))
 if key not in dcache:
  n,de=c.numerator(),c.denominator();nn,dd=ph(n),ph(de)
  dcache[key]=(ph(n.derivative(p))*dd-nn*ph(de.derivative(p)))/(dd*dd)
 return dcache[key]
def cp(f):return O({tuple(map(int,e)):K0(c) for e,c in f.dict().items()})
normal_positions=(1,2,3,4,5,6)
# A factor is indexed by (ordinary_degree,normal_exponents_6,u_degree,v_degree).
# Keeping its normal degree <=4 cannot lose any term of a product normal 4-jet.
def factor(dd,der=None):
 ans={}
 for ee,c in dd.items():
  e=tuple(map(int,ee));ne=tuple(e[j] for j in normal_positions)
  if sum(ne)>4:continue
  cc=cf(c) if der is None else dcf(c,der)
  key=(sum(e),ne,e[7],e[8]);ans[key]=ans.get(key,K.zero())+cc
 return {key:c for key,c in ans.items() if c}
J=PolynomialRing(K,names=('ax_u','ax_v','ax_n13','ax_n21','ax_n23','ax_n31','ax_n32','ax_nR1','ax_tau'));jv=J.gens()
def one_jet(fd,degree):
 terms={}
 for (d,n,up,vp),c in fd.items():
  tp=degree-d
  if tp<0 or sum(n)+tp>4:continue
  exp=(up,vp)+n+(tp,);terms[exp]=terms.get(exp,K.zero())+c
 return J(terms)
def product_jet(left,right,degree):
 terms={}
 for (dl,nl,ul,vl),cl in left.items():
  for (dr,nr,ur,vr),cr in right.items():
   tp=degree-dl-dr
   if tp<0:continue # global degree-six cancellation was proved separately
   n=tuple(nl[i]+nr[i] for i in range(6))
   if sum(n)+tp>4:continue
   exp=(ul+ur,vl+vr)+n+(tp,);terms[exp]=terms.get(exp,K.zero())+cl*cr
 return J({e:c for e,c in terms.items() if c})
records=[];jets=[]
for qi in (0,1,5,6):
 f=cp(oldq['quartic_source_rows'][qi]);deg=max(sum(map(int,e)) for e in f.dict());jet=one_jet(factor(f.dict()),deg)
 records.append({'kind':'Q','index':qi,'ordinary_degree':int(deg)});jets.append(jet)
for i,f in enumerate(g['six_full_sources']):
 jet=one_jet(factor(f.dict()),4);records.append({'kind':'R','position':i,'ordinary_degree':4});jets.append(jet)
print('ORIGINALS',len(jets),'sec',time.monotonic()-start,flush=True)
hc={}
def hf(ref):
 if ref not in hc:hc[ref]=factor(hdata['polynomial_pool'][ref])
 return hc[ref]
for item in hdata['records']:
 jet=sum((sign*product_jet(hf(l),hf(rr),5) for sign,l,rr in item['full_signed_product_DAG']),J.zero())
 records.append({'kind':'H','direction':item['direction'],'R_position':item['source_position'],'ordinary_degree':5});jets.append(jet)
 print('H',item['direction'],item['source_position'],'terms',len(jet.dict()),'sec',round(time.monotonic()-start,1),flush=True)
src=[cp(f) for f in g['six_full_sources']];unit=g['source_unit_coefficients'];efcache={}
def ef(f,der=None):
 key=(id(f),str(der))
 if key not in efcache:efcache[key]=factor(f.dict(),der)
 return efcache[key]
hat=[[cp(g['first_gauged_fields'][i,k]) for k in range(9)] for i in range(3)];base=[[cp(g['base_fields'][i,k]) for k in range(2)] for i in range(3)];Qs={qi:cp(oldq['quartic_source_rows'][qi]) for qi in (1,5,6)}
for item in edata['items']:
 i=item['direction'];qi=item['Q_source_index'];f=Qs[qi];co=[O({tuple(map(int,e)):K0(c) for e,c in dd.items()}) for dd in item['global_support_coefficient_dicts']]
 jet=J.zero();derivatives=[]
 for k in range(9):
  df=f.derivative(ov[k]);derivatives.append(df);jet+=product_jet(ef(hat[i][k]),factor(df.dict()),5)
 for k,p in enumerate((oh,oz)):jet+=product_jet(ef(base[i][k]),ef(f,p),5)
 for k,(_,pos,sgn) in enumerate(edata['supports'][i]):jet-=sgn*product_jet(factor((unit[i]*co[k]).dict()),ef(src[pos]),5)
 records.append({'kind':'E','direction':i,'Q_index':qi,'ordinary_degree':5});jets.append(jet)
 print('E',i,qi,'terms',len(jet.dict()),'sec',round(time.monotonic()-start,1),flush=True)
B=PolynomialRing(K,names=('ax_u','ax_v'));U,V=B.gens();lin=[];quad=[]
for f in jets:
 assert all(sum(map(int,e[2:])) in (1,2,3,4) for e in f.dict()),'nonzero boundary value or excess normal degree'
 lin.append([B({(int(e[0]),int(e[1])):c for e,c in f.dict().items() if sum(map(int,e[2:]))==1 and int(e[2+j])==1}) for j in range(7)])
 quad.append(J({tuple(map(int,e)):c for e,c in f.dict().items() if sum(map(int,e[2:]))==2}))
M=matrix(B,lin)
payload={'ring':J,'coefficient_field':K,'tangent_ring':B,'normal_variables':jv[2:],'tangent_variables':jv[:2],'base_section':(1,-1),'source_records':records,'normal_four_jets':jets,'normal_layers':[[J({tuple(map(int,e)):c for e,c in fq.dict().items() if sum(map(int,e[2:]))==degree}) for fq in jets] for degree in range(1,5)],'normal_jacobian':M,'normal_quadratic_parts':quad,'checks':len(jets),'source_files':['round8/rank3_joint_fiber_leading.sobj','independent_ratio_p_boundary_second_gauge.sobj','independent_ratio_global_R_derivative_degreefall.sobj','independent_ratio_even_derivative_cyclic_direct.sobj'],'scope':__doc__,'coefficient_derivatives':'H includes both original h,z terms in the stored full DAG. E evaluates original quotient-rule h,z derivatives before cyclic specialization.','filtration':'H and E are degree <=5 only after their globally proved degree-six cancellation. Discarded individual degree-six terms sum identically zero before any axis restriction.'}
stem='independent_ratio_p12_axis_normal_fourjets';tmp=out/(stem+'.tmp.sobj');save(payload,str(tmp),compress=False);os.replace(tmp,out/(stem+'.sobj'))
rec={'sources':len(jets),'jet_term_counts':[len(f.dict()) for f in jets],'matrix_shape':list(M.dimensions()),'seconds':time.monotonic()-start,'scope':__doc__};(out/(stem+'.json')).write_text(json.dumps(rec,indent=2)+'\n');print('SAVED',len(jets),time.monotonic()-start,flush=True)
# No pivot inversion is used in the source payload. Rank is diagnostic only;
# any rank-defect locus must be retained in subsequent geometry.
rank=M.change_ring(B.fraction_field()).rank();print('GENERIC_RANK',rank,flush=True)

