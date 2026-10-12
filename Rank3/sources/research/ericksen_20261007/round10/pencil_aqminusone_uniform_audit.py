"""Replay the source specialization, closed formulas, and all real sign claims."""
from pathlib import Path
import json,time
from sage.all import load,save,QQ,PolynomialRing,LaurentPolynomialRing
out=Path(__file__).resolve().parent;start=time.monotonic();checks={}
src=load(str(out/'pencil_vzero_sone_cubics.sobj'));wall=load(str(out/'pencil_aqminusone_newton.sobj'))
L=LaurentPolynomialRing(QQ,names=('q','T','Hb'));qL,TL,HL=L.gens()
source_map=src['base_ring'].hom([-qL**-1,qL,TL,HL],L)
wall_map=wall['base_ring'].hom([qL,TL,HL],L)
for i in range(2):
 e=wall['records'][i]['q_clearing_power']
 for j in range(4):
  checks[f'original_source_specialization_{i}_{j}']=source_map(src['cubics'][i][j])*qL**e==wall_map(wall['cubics'][i][j])
print('source replay',round(time.monotonic()-start,2),flush=True)
up=load(str(out/'pencil_aqminusone_upper_walls.sobj'));C=up['coefficients'][0].parent();q,T=C.gens();r=q-1;A=q*q+6*q+1;NN=A*T-4*r*r;JJ=(T-1)*NN-8*q*T
expected=16384*T**6*q**8*(T-1)**9*(T-2)**13*r**25*NN*JJ**4
checks['upper_coefficient_4_formula']=up['coefficients'][4]==expected
checks['upper_coefficients_0_to_3_zero']=all(not x for x in up['coefficients'][:4])
nf=load(str(out/'pencil_aqminusone_wall_factors.sobj'));K=nf['N']['norm'].parent();qq=K.gen() if K.ngens()==1 else K.ring().gen()
# Use the actual rational-field generator, avoiding display-unit ambiguities.
qq=K(K.ring().gen());rr=qq-1;AA=qq*qq+6*qq+1
en=-QQ(60514543933804184076288)*qq**12*rr**49*(qq*qq-QQ(14)/3*qq+1)**8*(qq*qq-10*qq+1)**14*(qq**4-140*qq**3+150*qq*qq-140*qq+1)/AA**35
ej=-QQ(361680002717660915922878529496292660128776192)*qq**38*rr**90*(qq**4-QQ(4)/27*qq**3+QQ(722)/27*qq*qq-QQ(4)/27*qq+1)**4/AA**39
checks['N_remainder_formula']=nf['N']['norm']==en
checks['J_norm_formula']=nf['J']['norm']==ej
az=load(str(out/'pencil_aqminusone_Azero.sobj'));gg,u,v=az['bezout'];checks['Azero_bezout']=gg==1 and u*az['wall']+v*az['remainder']==1
tt=load(str(out/'pencil_aqminusone_Ttwo_sylvester.sobj'));R=tt['coefficients'][0].parent();x=R.gen()
checks['Ttwo_coefficient_formula']=tt['coefficients'][2]==QQ(274877906944)*x**8*(x-1)**23*(x*x-10*x+1)*(x*x-4*x+1)*(x*x-18*x+1)**4
checks['Ttwo_lower_coefficients_zero']=all(not y for y in tt['coefficients'][:2])
signs=[x*x-QQ(14)/3*x+1,x*x-10*x+1,x**4-140*x**3+150*x*x-140*x+1,x**4-QQ(4)/27*x**3+QQ(722)/27*x*x-QQ(4)/27*x+1,x*x-4*x+1,x*x-18*x+1]
for i,p in enumerate(signs):checks[f'positive_at_negative_q_{i}']=all(c>0 for c in p(-x).coefficients()) and p(0)>0
assert all(checks.values()),{k:v for k,v in checks.items() if not v}
report={'status':'PASS_UNIFORM_AQPLUSONE_WALL','checks':checks,'check_count':len(checks),'elapsed_seconds':time.monotonic()-start,'scope':'Every fixed real q<0 and T!=0,1 has a nonzero Hb coefficient in the original cubic resultant. Covers all temporary A=0 and T=2 branches. Does not claim the entire Ericksen problem is solved.'}
(out/'pencil_aqminusone_uniform_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
