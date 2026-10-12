from sage.all import *
from pathlib import Path
import time,json,re,hashlib
p=Path(__file__).resolve().parent;st=time.monotonic();base=load(str(p/'independent_sumzero_PID_adjugate_certificate.sobj'));third=load(str(p/'independent_sumzero_PID_resultant13_input.sobj'));cand=load(str(p/'independent_sumzero_PID_resultant_gcd.sobj'));Rr=cand['ring'];rp=Rr.gen();Hcand=cand['resultants'];checks=[]
def ck(n,v):
 if not bool(v):raise AssertionError(n)
 checks.append(n)
f,g=base['univariate_polynomials'];f3,g3=third['polynomials'];ck('shared_first_polynomial',f3==f);P=base['big_factors'][0].parent();a,r=P.gens();source=load(str(p/'independent_ratio_canonical_sumzero_tau_chart.sobj'));C=source['coefficient_matrix'];clear=base['source_matrix_clear_denominator'];N3=vector(P,[P(v) for v in third['third_numerators']]);D3=P(third['original_denominator']);Cm=matrix(P,C.nrows(),C.ncols(),[P(c*clear) for c in C.list()]);ans=Cm*N3
for i,v in enumerate(ans):ck('third_original47_row_'+str(i),v==clear*D3*base['source_target_coefficients'][i])
print('THIRD ORIGINAL47 PASS',round(time.monotonic()-st,2),flush=True)
polys=[f,g,g3];icoefs=[];scales=[]
for j,F in enumerate(polys):
 cs=[Rr(c) for c in F.list()];D=lcm([c.denominator() for cc in cs for c in cc.coefficients()]);scales.append(ZZ(D));icoefs.append([ZZ['r'](c*D) for c in cs]);ck('integer_polynomial_lift_'+str(j),all(Rr(icoefs[-1][i])/D==cs[i] for i in range(len(cs))))
n,m=int(f.degree()),int(g.degree());ck('common_degrees',(n,m)==(24,21) and g3.degree()==21);bound=max(m*max(c.degree() for c in icoefs[0])+n*max(c.degree() for c in icoefs[j]) for j in (1,2));ck('declared_degree_bound',bound==1008)
for j,h in enumerate(Hcand):ck('candidate_degree_bound_'+str(j),h.degree()<=bound)
# Verify one-variable resultant convention by the exact multiplication map;
# det S=(-1)^(n*m) Res(f,g), which is +1 in degrees24,21.
ck('determinant_sign_even',(n*m)%2==0)
results=[]
for node in range(2,2+int(bound)+1):
 vals=[[c(node) for c in cs] for cs in icoefs];row=[]
 for j in (1,2):
  mat=matrix(ZZ,n+m,n+m)
  for k in range(m):
   for i,c in enumerate(vals[0]):mat[i+k,k]=c
  for k in range(n):
   for i,c in enumerate(vals[j]):mat[i+k,m+k]=c
  det=mat.det();expected=scales[0]**m*scales[j]**n*Hcand[j-1](node);ck('exact_ZZ_det_%d_%d'%(node,j),det==expected);row.append(hashlib.sha256(str(det).encode()).hexdigest())
 results.append({'node':node,'sha256_integer_determinants':row})
 if len(results)%50==0:
  print('EXACT NODES',len(results),'/1009',round(time.monotonic()-st,2),flush=True);(p/'independent_sumzero_PID_resultant_reader_progress.json').write_text(json.dumps({'nodes':len(results),'seconds':time.monotonic()-st})+'\n')
ck('distinct_exact_nodes',len(results)==1009 and len({it['node'] for it in results})==1009)
h0,h1=Hcand;gcd=h0.gcd(h1);q=rp+1;rem=gcd//(rp**190*q**192);fA=rp**2-QQ(1193)/27*rp+QQ(292)/27;fB=rp**2+rp-1;unit,rest=rem.quo_rem(fA*fB);ck('final_gcd_factorization',not rest and unit.degree()==0 and unit!=0);ck('positive_interval_first_identity',fA==rp**2+QQ(1193)/27*(-rp)+QQ(292)/27);ck('negative_interval_second_identity',fB==-(1+(-rp)*(rp+1)))
report={'status':'PASS_EXACT_ALL_PARAMETER_RESULTANTS','count':len(checks),'seconds':time.monotonic()-st,'degree_bound':int(bound),'distinct_exact_integer_nodes':len(results),'determinants_per_node':2,'candidate_degrees':[int(h.degree()) for h in Hcand],'gcd_factorization':'nonzero rational * r^190*(r+1)^192*(r^2-1193*r/27+292/27)*(r^2+r-1)','scope':'Each native resultant replayed via1009 distinct exact ZZ Sylvester determinants and degree<=1008, not finite-field testing. Third original47 source rows independently replayed. For ordered positive spectrum1<a<b, -1<r<0 and both extra quadratic factors are nonzero, so the two exact PID resultants never vanish simultaneously there. Unordered material r outside this interval is not covered without relabeling.'};save({'checks':checks,'results':results,'integer_polynomial_coefficients':icoefs,'integer_scales':scales,'resultants':Hcand,'gcd':gcd,'factor_units':unit},str(p/'independent_sumzero_PID_resultant_reader.sobj'),compress=False);(p/'independent_sumzero_PID_resultant_reader.json').write_text(json.dumps(report,indent=2)+'\n');print(report,flush=True)
