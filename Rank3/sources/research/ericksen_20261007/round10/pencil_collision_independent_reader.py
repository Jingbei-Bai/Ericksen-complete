"""Audit collision polynomial, global/local resultant coefficient, and linear remainder."""
from pathlib import Path
import time,json
from sage.all import load,save
out=Path(__file__).resolve().parent;start=time.monotonic()
d=load(str(out/'pencil_collision_newton.sobj'));Z=d['ring'];K=Z.base_ring();C=K.ring();q,T=C.gens();r=q-1
f,g=d['quadratics'];F0,G0=d['series'][0][0],d['series'][1][0]
checks={}
def check(key,val):checks[key]=bool(val);assert checks[key],key
factor_candidates=[p for p,e in d['resultant'].numerator().factor() if p.degree(q)==6]
assert len(factor_candidates)==1
P=factor_candidates[0]
leading_q=C({(0,mon[1]):cc for mon,cc in P.dict().items() if mon[0]==6})
check('P_q6_coefficient',leading_q==16*T**4*(T-1)**2*(T-2)**2)
at_test=P(-1,3);check('P_minus1_3_negative',at_test<0)
# q=-1,T=3 lies in the physical positive-simple-spectral domain on Psi=0.
av=d['a_on_wall'];bv=av*K(q)/(1+av*K(r))
aa=av(q=-1,T=3);bb=bv(q=-1,T=3)
check('negative_witness_positive_simple_spectrum',aa>0 and bb>0 and aa!=1 and bb!=1 and aa!=bb)
# In the branch T!=2 the noncolliding roots and the quadratic leading
# coefficients are nonzero by the original guards. The rational formula
# for [t^4] is a polynomial identity after denominators are cleared.
prefactor=F0.leading_coefficient()**3*G0(0)/f.leading_coefficient()**2
global_t4=prefactor*d['resultant']
fm=f/f.leading_coefficient();gm=g/g.leading_coefficient();linear=gm-fm
check('normalized_difference_degree_at_most1',linear.degree()<=1)
L=linear[1];M=linear[0]
check('monic_quadratic_resultant_formula',fm.resultant(gm)==M*M-fm[1]*L*M+fm[0]*L*L)
report={'status':'PASS_COLLISION_SOURCE_AND_FACTOR_AUDIT','checks':checks,'P_terms':P.number_of_terms(),
 'P_at_qminus1_T3':str(at_test),'witness_a_b':[str(aa),str(bb)],
 'global_local_t4_prefactor':str(prefactor.factor()),'global_t4_factorization':str(global_t4.factor()),
 'monic_difference_L_factorization':str(L.factor()),'monic_difference_M_factorization':str(M.factor()),
 'elapsed_seconds':time.monotonic()-start,
 'scope':'P is allowed to vanish. The global/local t4 relation uses the standard simple noncolliding root product; an independent full Sylvester t4 replay may strengthen its certificate.'}
save({'P':P,'leading_q':leading_q,'prefactor':prefactor,'global_t4':global_t4,
      'monic_quadratics':(fm,gm),'linear_difference':linear,'L':L,'M':M,'checks':checks},str(out/'pencil_collision_independent_reader.sobj'))
(out/'pencil_collision_independent_reader.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k not in ('global_t4_factorization','monic_difference_M_factorization')},indent=2),flush=True)
