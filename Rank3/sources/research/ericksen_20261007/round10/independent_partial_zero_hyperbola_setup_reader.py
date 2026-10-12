from sage.all import *
from pathlib import Path
from itertools import permutations
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(n,v):
 if not bool(v):raise AssertionError(n)
 checks.append(n)
old=load(str(p/'independent_ratio_canonical_amplitude_normal.sobj'));new=load(str(p/'independent_second_prolongation_reader.sobj'));sec=load(str(p/'independent_ratio_RR_cyclic_full_matrix.sobj'));base=load(str(p/'independent_ratio_canonical_fullbase_flat_minor.sobj'));short=load(str(p/'independent_ratio_canonical_fullbase_short_tau.sobj'));charts=load(str(p/'independent_ratio_canonical_partial_zero_hyperbolas.sobj'));K=old['coefficient_field'];a,r=old['parameters'];q=1+r;m=1+a*r;R=PolynomialRing(K,names=('audit_t0','audit_t1','audit_t2'));t=R.gens();GK=base['ring'].base_ring();GA,GR,GH,GZ=base['parameters'];mp=GK.ring().hom([a,r,1,-1],K)
def cv(c):c=GK(c);return mp(c.numerator())/mp(c.denominator())
def cyc(f):return R({tuple(map(int,e)):cv(c) for e,c in f.dict().items()})
def oldp(f):return R({tuple(map(int,e)):K(c) for e,c in f.dict().items()})
cs=[-4*(a*r+r+2)/(a*r*q),4*r*r*(a+1)/(q*q*m),-4*r*r*(2*a*r+a+1)/q];Psi=t[0]*t[1]+t[0]*t[2]+t[1]*t[2]+cs[0]*t[0]**2*t[1]**2+cs[1]*t[0]**2*t[2]**2+cs[2]*t[1]**2*t[2]**2;Lam=sum(t);M=matrix(R,[[cyc(v) for v in row] for row in short['matrix'].rows()]);S=matrix(R,[[cyc(v) for v in row] for row in short['short_rows'].rows()]);v=vector(R,[cyc(v) for v in short['left_syzygy']]);Pcy=cyc(short['pure_tau_coefficient']);G=K.ring();baseguards=[G(a),G(r),G(q),G(m),G(a-1)]
def guard_only(f,guards):
 f=G(f)
 for g in guards:
  while f and f.degree()>0:
   quo,rem=f.quo_rem(g)
   if rem:break
   f=quo
 return bool(f) and f.degree()<=0
for i in range(7):
 u=K(sec['row_scales'][i]);ck('original_old_source_unit_'+str(i),guard_only(u.numerator(),baseguards) and guard_only(u.denominator(),baseguards))
 for j in range(6):ck('fullbase_to_old_source_%d_%d'%(i,j),M[i,j]==u*oldp(old['matrix'][i,j]))
for ii,(i,j,col) in enumerate(short['source_pairs']):
 scale=cv(short['row_scales'][ii]);ck('cyclic_source_combination_unit_'+str(ii),guard_only(scale.numerator(),baseguards) and guard_only(scale.denominator(),baseguards))
 for k in range(6):ck('cyclic_actual_row_combination_%d_%d'%(ii,k),scale*S[ii,k]==M[j,col]*M[i,k]-M[i,col]*M[j,k])
ck('cyclic_true_syzygy',list(v*S)==[0,0,0,Pcy,0,0]);pu=K(Pcy.monomial_coefficient(t[0]*t[1]));lu=K(cyc(base['Lambda']).monomial_coefficient(t[0]));ck('Psi_coefficients_from_true_source',Pcy==pu*Psi);ck('Lambda_from_true_source',cyc(base['Lambda'])==lu*Lam);ck('Psi_normalizing_unit',guard_only(pu.numerator(),baseguards) and guard_only(pu.denominator(),baseguards));ck('Lambda_normalizing_unit',guard_only(lu.numerator(),baseguards) and guard_only(lu.denominator(),baseguards))
minor_records=[]
for j,it in enumerate(base['items']):
 i,k=it['original_direction_rows'];Y=matrix(R,[S.row(0),S.row(1),S.row(2),M.row(0),M.row(i),M.row(k)]);N4=Y[:4,:4];det4=R.zero()
 for perm in permutations(range(4)):
  sign=(-1)**sum(perm[i]>perm[j] for i in range(4) for j in range(i+1,4));val=R(sign)
  for i in range(4):val*=N4[i,perm[i]]
  det4+=val
 det=det4*(Y[4,4]*Y[5,5]-Y[4,5]*Y[5,4]);unit=-pu*lu*cv(it['direction_unit']);idx=int(it['amplitude']);ck('actual_old7_minor_'+str(j),det==unit*t[idx]**2*Lam*Psi);ck('actual_old7_minor_saved_'+str(j),det==cyc(it['row_combination_determinant']));ck('actual_old7_minor_unit_'+str(j),guard_only(unit.numerator(),baseguards) and guard_only(unit.denominator(),baseguards));minor_records.append({'amplitude':idx,'unit':str(unit.factor())})
ck('flat_g2_at_cyclic',cv(base['flat_g_square'])==1);ck('flat_B1_at_cyclic',cv(base['flat_B1'])==1)
ck('positive_c0_numerator_identity',a*r+r+2==m+q);ck('positive_c2_numerator_identity',2*a*r+a+1==2*m+a-1)
print('SOURCE PSI/THREE TRUE MINORS PASS',len(checks),round(time.monotonic()-st,2),flush=True)
raw=list(old['matrix'].rows())+list(new['normal_matrix'].rows());U=PolynomialRing(K,names=('audit_x','audit_unused'));x=U.gen(0);records=[];ck('three_charts_complete',len(charts['items'])==3)
for item in charts['items']:
 case=int(item['case']);i,j=map(int,item['remaining_amplitude_indices']);missing=3-i-j;c=cs[case];ck('hyperbola_coefficient_'+str(case),K(item['coefficient'])==c);ck('coefficient_nonzero_rational_function_'+str(case),c!=0)
 pts=list(t);pts[missing]=R.zero();restricted=R.zero()
 for e,co in Psi.dict().items():
  if not e[missing]:restricted+=co*R.monomial(*tuple(map(int,e)))
 ck('hyperbola_necessary_factorization_'+str(case),restricted==t[i]*t[j]*(1+c*t[i]*t[j]));den_guards=baseguards+[G(a*r+r+2),G(a+1),G(2*a*r+a+1)];den_seen=set()
 for rn,row in enumerate(raw):
  entries=[]
  for f in row:
   dd={}
   for ex,co in f.dict().items():
    ex=tuple(map(int,ex))
    if ex[missing]:continue
    power=ex[i]-ex[j];dd[power]=dd.get(power,K.zero())+K(co)*(-1/c)**ex[j]
   entries.append({pw:val for pw,val in dd.items() if val})
  low=min((pw for dd in entries for pw in dd),default=0);ck('hyperbola_row_shift_%d_%d'%(case,rn),int(item['row_shifts'][rn])==-low)
  for col,dd in enumerate(entries):
   f=U({(pw-low,0):co for pw,co in dd.items()});saved=U({tuple(map(int,e)):K(co) for e,co in item['matrix'][rn,col].dict().items()});ck('hyperbola_fullsource_%d_%d_%d'%(case,rn,col),f==saved)
   for co in dd.values():den_seen.add(co.denominator()/co.denominator().lc())
 for ind,den in enumerate(den_seen):ck('hyperbola_spectral_den_%d_%d'%(case,ind),guard_only(den,den_guards))
 records.append({'case':case,'pair':(i,j),'missing':missing,'coefficient':str(c.factor()),'spectral_extra_factors':['a*r+r+2','a+1','2*a*r+a+1'],'distinct_coefficient_denominators':len(den_seen),'rows_rebuilt':len(raw)})
 print('HYPERBOLA',case,'ALL91x6 SOURCE ENTRIES PASS',round(time.monotonic()-st,2),flush=True)
report={'status':'PASS','count':len(checks),'seconds':time.monotonic()-st,'Psi_coefficients':[str(c.factor()) for c in cs],'actual_old7_minor_units':minor_records,'charts':records,'scope':'Genuine selected old7 six-minors force Lambda*Psi=0 for nonzero surviving amplitudes; after sealed flat-sumzero layer, off-axis partial-zero survivors necessarily lie on these three hyperbolas. Entire91x6 hyperbola matrices reconstructed with recorded row units x^shift and all spectral denominators. Flat coordinates here are normalized to cyclic t via t=g*(xi0,B1xi1,xi2), g2=B1=1 at cyclic; bare geometric flat k-coordinates differ by a material-dependent common constant. No hyperbola rank exclusion asserted.'};save({'checks':checks,'Psi':Psi,'coefficients':cs,'records':records,'minor_records':minor_records},str(p/'independent_partial_zero_hyperbola_setup_reader.sobj'));(p/'independent_partial_zero_hyperbola_setup_reader.json').write_text(json.dumps(report,indent=2)+'\n');print(report,flush=True)
